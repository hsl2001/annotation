/**
 * @file ann.c
 * @brief Standalone fixed-shape 1D UNet runtime.
 *
 * Tensors are channel-major.  Every forward pass uses the same fixed window
 * size
 */
#include "unet.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

enum {
  MIN_DIMENSION = 1,
  OUTPUT_KERNEL = 1,
  INITIAL_DEPTH = 4,
  INITIAL_KERNEL = 9,
  INITIAL_POOL = 2,
  INITIAL_CHANNEL_0 = 16,
  INITIAL_CHANNEL_1 = 32,
  INITIAL_CHANNEL_2 = 64,
  INITIAL_CHANNEL_3 = 128,
  UNET_LAYER_COUNT = 15,
  UNET_NORMALIZED_LAYER_COUNT = UNET_LAYER_COUNT - 1,
  WORKSPACE_COUNT = 22,
};

#define NORM_EPSILON 1e-5 /**< variance floor of the channel normalization */

typedef struct {
  int in;
  int out;
  int k;
  float *w;
  float *b;
  float *mw;
  float *vw;
  float *mb;
  float *vb;
} Layer;
struct Unet {
  UnetConfig config;
  int nlayer;
  uint64_t optimizer_step;
  Layer layer[UNET_LAYER_COUNT];
};

/**
 * @brief Validate a UNet configuration.
 * @param c Configuration to validate.
 * @return Nonzero when the configuration is valid, otherwise zero.
 */
static int good_config(const UnetConfig *c) {
  if (!c || c->depth != INITIAL_DEPTH ||
      c->kernel < MIN_DIMENSION || !(c->kernel & MIN_DIMENSION) ||
      c->pool != INITIAL_POOL)
    return 0;
  for (int i = 0; i < c->depth; ++i)
    if (c->channels[i] < 1)
      return 0;
  return 1;
}

/**
 * @brief Return the default UNet configuration.
 * @return Default model configuration.
 */
UnetConfig unet_default_config(void) {
  return (UnetConfig){INITIAL_DEPTH,
                      INITIAL_KERNEL,
                      INITIAL_POOL,
                      {INITIAL_CHANNEL_0, INITIAL_CHANNEL_1, INITIAL_CHANNEL_2,
                       INITIAL_CHANNEL_3}};
}

/**
 * @brief Release all memory owned by a UNet.
 * @param u Model to release; NULL is accepted.
 * @return None.
 */
static void release(Unet *u) {
  if (!u)
    return;
  for (int i = 0; i < u->nlayer; ++i) {
    free(u->layer[i].w);
    free(u->layer[i].b);
    free(u->layer[i].mw);
    free(u->layer[i].vw);
    free(u->layer[i].mb);
    free(u->layer[i].vb);
  }
  free(u);
}

/**
 * @brief Append and initialize a convolution layer.
 * @param u Model receiving the layer.
 * @param in Number of input channels.
 * @param out Number of output channels.
 * @param k Convolution kernel width.
 * @param seed Mutable random-number generator state.
 * @return Nonzero on success, otherwise zero on allocation failure.
 */
static int add_layer(Unet *u, int in, int out, int k, uint64_t *seed) {
  Layer *l = &u->layer[u->nlayer++];
  l->in = in;
  l->out = out;
  l->k = k;
  size_t weights = (size_t)in * out * k;
  l->w = (float *)malloc(weights * sizeof(float));
  l->b = (float *)calloc((size_t)out, sizeof(float));
  l->mw = (float *)calloc(weights, sizeof(float));
  l->vw = (float *)calloc(weights, sizeof(float));
  l->mb = (float *)calloc((size_t)out, sizeof(float));
  l->vb = (float *)calloc((size_t)out, sizeof(float));
  if (!l->w || !l->b || !l->mw || !l->vw || !l->mb || !l->vb)
    return 0;
  /* He-uniform: Var(w) = 2 / fan_in keeps ReLU activations at unit scale;
   * the previous fixed 0.01 bound shrank the signal ~20x per layer. */
  float bound = sqrtf(6.0f / (float)((size_t)in * k));
  for (size_t i = 0; i < weights; ++i)
    l->w[i] = bound * anno_random_symmetric(seed);
  return 1;
}

/**
 * @brief Apply one Adam update to a parameter vector.
 * @param parameter Parameters to update.
 * @param gradient Accumulated gradients.
 * @param first_moment Adam first-moment state.
 * @param second_moment Adam second-moment state.
 * @param count Divisor applied to the accumulated gradient.
 * @param size Number of parameters.
 * @param step_size Bias-corrected Adam step size.
 * @return None.
 */
static void adam_update(float *parameter, const float *gradient,
                        float *first_moment, float *second_moment,
                        size_t count, size_t size, double step_size) {
  for (size_t i = 0; i < size; ++i) {
    float g = gradient[i] / (float)count;
    first_moment[i] =
        ANNO_ADAM_BETA1 * first_moment[i] + (1.0f - ANNO_ADAM_BETA1) * g;
    second_moment[i] =
        ANNO_ADAM_BETA2 * second_moment[i] + (1.0f - ANNO_ADAM_BETA2) * g * g;
    parameter[i] -= (float)(step_size * first_moment[i] /
                            (sqrt((double)second_moment[i]) +
                             ANNO_ADAM_EPSILON));
  }
}

static void free_gradients(float **gradient, float **bias, int layers) {
  if (!gradient || !bias)
    return;
  for (int l = 0; l < layers; ++l) {
    free(gradient[l]);
    free(bias[l]);
  }
  free(gradient);
  free(bias);
}

/**
 * @brief Create and initialize a UNet model.
 * @param c Model configuration.
 * @return Newly allocated model, or NULL for invalid input or allocation failure.
 */
Unet *unet_create(const UnetConfig *c) {
  if (!good_config(c))
    return NULL;
  Unet *u = (Unet *)calloc(1, sizeof(*u));
  if (!u)
    return NULL;
  u->config = *c;
  uint64_t seed = anno_random_seed();
  int prev = ANNO_INPUT_CHANNELS;
  for (int l = 0; l < c->depth; ++l) {
    if (!add_layer(u, prev, c->channels[l], c->kernel, &seed) ||
        !add_layer(u, c->channels[l], c->channels[l], c->kernel, &seed)) {
      release(u);
      return NULL;
    }
    prev = c->channels[l];
  }
  for (int l = c->depth - 2; l >= 0; --l) {
    if (!add_layer(u, c->channels[l + 1] + c->channels[l], c->channels[l],
                   c->kernel, &seed) ||
        !add_layer(u, c->channels[l], c->channels[l], c->kernel, &seed)) {
      release(u);
      return NULL;
    }
  }
  if (!add_layer(u, c->channels[0], ANNO_CLASSES, OUTPUT_KERNEL, &seed)) {
    release(u);
    return NULL;
  }
  return u;
}

/**
 * @brief Destroy a UNet model.
 * @param u Model to destroy; NULL is accepted.
 * @return None.
 */
void unet_destroy(Unet *u) { release(u); }
/**
 * @brief Apply one convolution to the channel concatenation [a ; b].
 *
 * Input channel i < split is read from a, the others from b.  The position
 * range of every tap is clipped once, so the innermost loop is a plain
 * contiguous multiply-add that the compiler vectorizes.
 * @param l Layer parameters.
 * @param a First channel-major input (split channels).
 * @param b Second channel-major input (l->in - split channels), or NULL.
 * @param split Number of channels read from a.
 * @param y Channel-major output activation.
 * @param n Number of positions.
 * @param stat NULL for a linear output; otherwise ReLU followed by
 *        per-window channel normalization, storing (mean, sd) per channel.
 * @return None.
 */
static void convolve(const Layer *l, const float *a, const float *b, int split,
                     float *y, int n, int valid, float *stat) {
  int r = l->k / 2;
#pragma omp parallel for schedule(static)
  for (int o = 0; o < l->out; ++o) {
    float *yo = y + (size_t)o * n;
    for (int p = 0; p < n; ++p)
      yo[p] = l->b[o];
    for (int i = 0; i < l->in; ++i) {
      const float *xi =
          i < split ? a + (size_t)i * n : b + (size_t)(i - split) * n;
      const float *w = l->w + ((size_t)o * l->in + i) * l->k;
      for (int q = 0; q < l->k; ++q) {
        int s = q - r, lo = s < 0 ? -s : 0, hi = s > 0 ? valid - s : valid;
        float wq = w[q];
        for (int p = lo; p < hi; ++p)
          yo[p] += wq * xi[p + s];
      }
    }
    if (stat) {
      double sum = 0, square = 0;
      for (int p = 0; p < valid; ++p) {
        yo[p] = yo[p] < 0 ? 0 : yo[p];
        sum += yo[p];
        square += (double)yo[p] * yo[p];
      }
      double mean = sum / valid, var = square / valid - mean * mean;
      float mu = (float)mean,
            sd = (float)sqrt((var > 0 ? var : 0) + NORM_EPSILON);
      for (int p = 0; p < valid; ++p)
        yo[p] = (yo[p] - mu) / sd;
      stat[2 * o] = mu;
      stat[2 * o + 1] = sd;
    }
    for (int p = valid; p < n; ++p)
      yo[p] = 0.0f;
  }
}

/** @brief Plain convolution layer; see convolve(). */
static void conv(const Layer *l, const float *x, float *y, int n, int valid,
                 float *stat) {
  convolve(l, x, NULL, l->in, y, n, valid, stat);
}

/**
 * @brief Apply average pooling to a channel-major tensor.
 * @param x Channel-major input tensor.
 * @param y Channel-major output tensor.
 * @param c Number of channels.
 * @param n Number of input positions.
 * @param f Pooling factor.
 * @return None.
 */
static void pool(const float *x, float *y, int c, int n, int f) {
#pragma omp parallel for collapse(2)
  for (int ch = 0; ch < c; ++ch)
    for (int p = 0; p < n / f; ++p) {
      float s = 0;
      for (int q = 0; q < f; ++q)
        s += x[(size_t)ch * n + p * f + q];
      y[(size_t)ch * (n / f) + p] = s / f;
    }
}

/**
 * @brief Expand each position by repeating it.
 * @param x Channel-major input tensor.
 * @param y Channel-major output tensor.
 * @param c Number of channels.
 * @param n Number of input positions.
 * @param f Upsampling factor.
 * @return None.
 */
static void upsample(const float *x, float *y, int c, int n, int f) {
#pragma omp parallel for collapse(2)
  for (int ch = 0; ch < c; ++ch)
    for (int p = 0; p < n; ++p)
      for (int q = 0; q < f; ++q)
        y[(size_t)ch * n * f + p * f + q] = x[(size_t)ch * n + p];
}

static void zero_tail(float *x, int channels, int n, int valid) {
  if (valid >= n)
    return;
#pragma omp parallel for schedule(static)
  for (int ch = 0; ch < channels; ++ch)
    memset(x + (size_t)ch * n + valid, 0,
           (size_t)(n - valid) * sizeof(*x));
}

static void valid_levels(int valid, int pool, int level[ANNO_MAX_DEPTH]) {
  level[0] = valid;
  for (int i = 1; i < ANNO_MAX_DEPTH; ++i)
    level[i] = (level[i - 1] + pool - 1) / pool;
}

typedef struct {
  float *a[WORKSPACE_COUNT], *d[WORKSPACE_COUNT];
  float *s[UNET_NORMALIZED_LAYER_COUNT]; /**< per-layer normalization (mean, sd) pairs */
  int n[WORKSPACE_COUNT], c[WORKSPACE_COUNT];
} Workspace;

static const int workspace_lengths[WORKSPACE_COUNT] = {
    ANNO_UNET_WINDOW,     ANNO_UNET_WINDOW,     ANNO_UNET_WINDOW,
    ANNO_UNET_WINDOW / 2, ANNO_UNET_WINDOW / 2, ANNO_UNET_WINDOW / 2,
    ANNO_UNET_WINDOW / 4, ANNO_UNET_WINDOW / 4, ANNO_UNET_WINDOW / 4,
    ANNO_UNET_WINDOW / 8, ANNO_UNET_WINDOW / 8, ANNO_UNET_WINDOW / 8,
    ANNO_UNET_WINDOW / 4, ANNO_UNET_WINDOW / 4, ANNO_UNET_WINDOW / 4,
    ANNO_UNET_WINDOW / 2, ANNO_UNET_WINDOW / 2, ANNO_UNET_WINDOW / 2,
    ANNO_UNET_WINDOW,     ANNO_UNET_WINDOW,     ANNO_UNET_WINDOW,
    ANNO_UNET_WINDOW};

static int workspace_init(const Unet *u, Workspace *a) {
  memset(a, 0, sizeof(*a));
  for (int i = 0; i < WORKSPACE_COUNT; ++i) {
    static const int channel_level[WORKSPACE_COUNT] = {
        -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3,
        3,  2, 2, 2, 1, 1, 1, 0, 0, -2};
    int level = channel_level[i];
    a->c[i] = level == -1 ? ANNO_INPUT_CHANNELS
                           : level == -2 ? ANNO_CLASSES
                                         : u->config.channels[level];
    a->n[i] = workspace_lengths[i];
    size_t size = (size_t)a->c[i] * a->n[i];
    a->a[i] = (float *)calloc(size, sizeof(float));
    if (i)
      a->d[i] = (float *)calloc(size, sizeof(float));
    if (!a->a[i] || (i && !a->d[i]))
      return 0;
  }
  for (int l = 0; l < UNET_NORMALIZED_LAYER_COUNT; ++l)
    if (!(a->s[l] = (float *)calloc(2 * (size_t)u->layer[l].out, sizeof(float))))
      return 0;
  return 1;
}

/**
 * @brief Free all tensors and operation metadata in a workspace.
 * @param a Workspace to clear.
 * @return None.
 */
static void workspace_free(Workspace *a) {
  for (int i = 0; i < WORKSPACE_COUNT; ++i) {
    free(a->a[i]);
    free(a->d[i]);
  }
  for (int l = 0; l < UNET_NORMALIZED_LAYER_COUNT; ++l)
    free(a->s[l]);
  memset(a, 0, sizeof(*a));
}


/** @brief Convolution over the concatenation [up ; skip]; see convolve(). */
static void conv_concat(const Layer *l, const float *up, const float *skip,
                        float *y, int skip_channels, int n, int valid,
                        float *stat) {
  convolve(l, up, skip, l->in - skip_channels, y, n, valid, stat);
}

/**
 * @brief Run the UNet forward pass and softmax output.
 * @param u Model to evaluate.
 * @param x Channel-major input tensor.
 * @param len Number of valid input positions.
 * @param a Workspace used by the pass.
 * @param prob Position-major class probabilities.
 * @return None.
 */
static void run_forward(const Unet *u, const float *x, int64_t len, Workspace *a,
                        float *prob) {
  int valid[ANNO_MAX_DEPTH];
  valid_levels((int)len, u->config.pool, valid);
  memcpy(a->a[0], x, (size_t)ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW *
                         sizeof(float));
  zero_tail(a->a[0], ANNO_INPUT_CHANNELS, a->n[0], valid[0]);
  conv(&u->layer[0], a->a[0], a->a[1], a->n[1], valid[0], a->s[0]);
  conv(&u->layer[1], a->a[1], a->a[2], a->n[2], valid[0], a->s[1]);
  pool(a->a[2], a->a[3], a->c[2], a->n[2], u->config.pool);
  conv(&u->layer[2], a->a[3], a->a[4], a->n[4], valid[1], a->s[2]);
  conv(&u->layer[3], a->a[4], a->a[5], a->n[5], valid[1], a->s[3]);
  pool(a->a[5], a->a[6], a->c[5], a->n[5], u->config.pool);
  conv(&u->layer[4], a->a[6], a->a[7], a->n[7], valid[2], a->s[4]);
  conv(&u->layer[5], a->a[7], a->a[8], a->n[8], valid[2], a->s[5]);
  pool(a->a[8], a->a[9], a->c[8], a->n[8], u->config.pool);
  conv(&u->layer[6], a->a[9], a->a[10], a->n[10], valid[3], a->s[6]);
  conv(&u->layer[7], a->a[10], a->a[11], a->n[11], valid[3], a->s[7]);
  upsample(a->a[11], a->a[12], a->c[11], a->n[11], u->config.pool);
  zero_tail(a->a[12], a->c[12], a->n[12], valid[2]);
  conv_concat(&u->layer[8], a->a[12], a->a[8], a->a[13], a->c[8], a->n[13],
              valid[2], a->s[8]);
  conv(&u->layer[9], a->a[13], a->a[14], a->n[14], valid[2], a->s[9]);
  upsample(a->a[14], a->a[15], a->c[14], a->n[14], u->config.pool);
  zero_tail(a->a[15], a->c[15], a->n[15], valid[1]);
  conv_concat(&u->layer[10], a->a[15], a->a[5], a->a[16], a->c[5], a->n[16],
              valid[1], a->s[10]);
  conv(&u->layer[11], a->a[16], a->a[17], a->n[17], valid[1], a->s[11]);
  upsample(a->a[17], a->a[18], a->c[17], a->n[17], u->config.pool);
  zero_tail(a->a[18], a->c[18], a->n[18], valid[0]);
  conv_concat(&u->layer[12], a->a[18], a->a[2], a->a[19], a->c[2], a->n[19],
              valid[0], a->s[12]);
  conv(&u->layer[13], a->a[19], a->a[20], a->n[20], valid[0], a->s[13]);
  conv(&u->layer[14], a->a[20], a->a[21], a->n[21], valid[0], NULL);
  float *z = a->a[21];
  int64_t usable = len;
#pragma omp parallel for
  for (int64_t p = 0; p < usable; ++p) {
    float top = -INFINITY, sum = 0;
    for (int j = 0; j < ANNO_CLASSES; ++j) {
      if (z[(size_t)j * ANNO_UNET_WINDOW + p] > top)
        top = z[(size_t)j * ANNO_UNET_WINDOW + p];
    }
    for (int j = 0; j < ANNO_CLASSES; ++j) {
      prob[(size_t)p * ANNO_CLASSES + j] =
          expf(z[(size_t)j * ANNO_UNET_WINDOW + p] - top);
      sum += prob[(size_t)p * ANNO_CLASSES + j];
    }
    for (int j = 0; j < ANNO_CLASSES; ++j)
      prob[(size_t)p * ANNO_CLASSES + j] /= sum;
  }
}

/**
 * @brief Evaluate channel-major input features.
 * @param u Model to evaluate.
 * @param x Channel-major input [channel][position].
 * @param n Number of input positions.
 * @param p Position-major output probabilities.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_predict_channels(const Unet *u, const float *x, int64_t n, float *p) {
  if (!u || !x || !p || n < 1 || n > ANNO_UNET_WINDOW)
    return EINVAL;
  Workspace a = {0};
  if (!workspace_init(u, &a))
    return ENOMEM;
  run_forward(u, x, n, &a, p);
  workspace_free(&a);
  return 0;
}

/**
 * @brief Evaluate interleaved input features.
 *
 * @param u Model to evaluate.
 * @param in Interleaved input [position][channel].
 * @param n Number of input positions.
 * @param p Output probabilities [position][class].
 * @return 0 on success, or an errno value otherwise.
 */
int unet_predict(const Unet *u, const float *in, int64_t n, float *p) {
  if (!u || !in || !p || n < 1 || n > ANNO_UNET_WINDOW)
    return EINVAL;
  float *x = (float *)calloc((size_t)ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW,
                             sizeof(float));
  if (!x)
    return ENOMEM;
  for (int64_t i = 0; i < n; ++i)
    for (int c = 0; c < ANNO_INPUT_CHANNELS; ++c)
      x[(size_t)c * ANNO_UNET_WINDOW + i] =
          in[(size_t)i * ANNO_INPUT_CHANNELS + c];
  int rc = unet_predict_channels(u, x, n, p);
  free(x);
  return rc;
}
/**
 * @brief Backpropagate through convolve().
 * @param l Layer parameters.
 * @param a First forward input (split channels).
 * @param b Second forward input, or NULL.
 * @param split Number of channels read from a.
 * @param out Forward output (normalized when stat is set).
 * @param dy Output gradient; overwritten with the pre-ReLU gradient.
 * @param da Receives the gradient of a (overwritten), or NULL to skip.
 * @param db Receives the gradient of b (overwritten), or NULL.
 * @param n Activation length.
 * @param stat Forward (mean, sd) per channel, or NULL for a linear layer.
 * @param gw Weight gradient, accumulated.
 * @param gb Bias gradient, accumulated.
 * @return None.
 */
static void convolve_backward(const Layer *l, const float *a, const float *b,
                              int split, const float *out, float *dy,
                              float *da, float *db, int n, int valid,
                              const float *stat, float *gw, float *gb) {
  int r = l->k / 2;
  for (int o = 0; o < l->out; ++o) {
    float *g = dy + (size_t)o * n;
    for (int p = valid; p < n; ++p)
      g[p] = 0.0f;
    if (!stat)
      continue;
    /* Normalization statistics are over real positions, not zero padding. */
    const float *x = out + (size_t)o * n;
    float mu = stat[2 * o], sd = stat[2 * o + 1], floor = (0.0f - mu) / sd;
    double m1 = 0, m2 = 0;
    for (int p = 0; p < valid; ++p) {
      m1 += g[p];
      m2 += (double)g[p] * x[p];
    }
    m1 /= valid;
    m2 /= valid;
    for (int p = 0; p < valid; ++p)
      g[p] = x[p] > floor ? (float)((g[p] - m1 - x[p] * m2) / sd) : 0.0f;
  }
#pragma omp parallel for schedule(static)
  for (int o = 0; o < l->out; ++o) {
    const float *g = dy + (size_t)o * n;
    float sum = 0;
    for (int p = 0; p < valid; ++p)
      sum += g[p];
    gb[o] += sum;
    for (int i = 0; i < l->in; ++i) {
      const float *xi =
          i < split ? a + (size_t)i * n : b + (size_t)(i - split) * n;
      float *w = gw + ((size_t)o * l->in + i) * l->k;
      for (int q = 0; q < l->k; ++q) {
        int s = q - r, lo = s < 0 ? -s : 0, hi = s > 0 ? valid - s : valid;
        float dw = 0;
        for (int p = lo; p < hi; ++p)
          dw += g[p] * xi[p + s];
        w[q] += dw;
      }
    }
  }
  if (!da)
    return;
#pragma omp parallel for schedule(static)
  for (int i = 0; i < l->in; ++i) {
    float *dx = i < split ? da + (size_t)i * n : db + (size_t)(i - split) * n;
    memset(dx, 0, (size_t)n * sizeof(float));
    for (int o = 0; o < l->out; ++o) {
      const float *g = dy + (size_t)o * n;
      const float *w = l->w + ((size_t)o * l->in + i) * l->k;
      for (int q = 0; q < l->k; ++q) {
        int s = q - r, lo = s < 0 ? -s : 0, hi = s > 0 ? valid - s : valid;
        float wq = w[q];
        for (int p = lo; p < hi; ++p)
          dx[p + s] += wq * g[p];
      }
    }
  }
}

/** @brief Backpropagate through conv(); dx may be NULL. */
static void backward_layer(const Layer *l, const float *x, const float *out,
                           float *y, float *dx, int n, int valid,
                           const float *stat, float *gw, float *gb) {
  convolve_backward(l, x, NULL, l->in, out, y, dx, NULL, n, valid, stat, gw,
                    gb);
}

/** @brief Backpropagate through conv_concat(). */
static void backward_concat_layer(const Layer *l, const float *up,
                                  const float *skip, const float *out,
                                  float *y, float *dup, float *dskip,
                                  int skip_channels, int n, int valid,
                                  const float *stat, float *gw, float *gb) {
  convolve_backward(l, up, skip, l->in - skip_channels, out, y, dup, dskip, n,
                    valid, stat, gw, gb);
}

/* Accumulates into dx: the pooled tensor is also a skip input, whose
 * gradient was written there by backward_concat_layer(). */
static void backward_pool(const float *dy, float *dx, int c, int n, int f,
                          int valid_in, int valid_out) {
  for (int ch = 0; ch < c; ++ch)
    for (int p = 0; p < valid_out; ++p)
      for (int q = 0; q < f && p * f + q < valid_in; ++q)
        dx[(size_t)ch * n + p * f + q] +=
            dy[(size_t)ch * (n / f) + p] / f;
}

static void backward_upsample(const float *dy, float *dx, int c, int n, int f,
                              int valid_in, int valid_out) {
  memset(dx, 0, (size_t)c * n * sizeof(float));
  for (int ch = 0; ch < c; ++ch)
    for (int p = 0; p < valid_in; ++p)
      for (int q = 0; q < f && p * f + q < valid_out; ++q)
        dx[(size_t)ch * n + p] += dy[(size_t)ch * n * f + p * f + q];
}

/* Backpropagate a->d[21] through the whole network of one window. */
static void window_backward(const Unet *u, Workspace *a, float **gw,
                            float **gb, int valid_input) {
  int valid[ANNO_MAX_DEPTH];
  valid_levels(valid_input, u->config.pool, valid);
  backward_layer(&u->layer[14], a->a[20], a->a[21], a->d[21], a->d[20],
                 a->n[21], valid[0], NULL, gw[14], gb[14]);
  backward_layer(&u->layer[13], a->a[19], a->a[20], a->d[20], a->d[19],
                 a->n[20], valid[0], a->s[13], gw[13], gb[13]);
  backward_concat_layer(&u->layer[12], a->a[18], a->a[2], a->a[19], a->d[19],
                        a->d[18], a->d[2], a->c[2], a->n[19], valid[0],
                        a->s[12], gw[12], gb[12]);
  backward_upsample(a->d[18], a->d[17], a->c[17], a->n[17], u->config.pool,
                    valid[1], valid[0]);
  backward_layer(&u->layer[11], a->a[16], a->a[17], a->d[17], a->d[16],
                 a->n[17], valid[1], a->s[11], gw[11], gb[11]);
  backward_concat_layer(&u->layer[10], a->a[15], a->a[5], a->a[16], a->d[16],
                        a->d[15], a->d[5], a->c[5], a->n[16], valid[1],
                        a->s[10], gw[10], gb[10]);
  backward_upsample(a->d[15], a->d[14], a->c[14], a->n[14], u->config.pool,
                    valid[2], valid[1]);
  backward_layer(&u->layer[9], a->a[13], a->a[14], a->d[14], a->d[13],
                 a->n[14], valid[2], a->s[9], gw[9], gb[9]);
  backward_concat_layer(&u->layer[8], a->a[12], a->a[8], a->a[13], a->d[13],
                        a->d[12], a->d[8], a->c[8], a->n[13], valid[2],
                        a->s[8], gw[8], gb[8]);
  backward_upsample(a->d[12], a->d[11], a->c[11], a->n[11], u->config.pool,
                    valid[3], valid[2]);
  backward_layer(&u->layer[7], a->a[10], a->a[11], a->d[11], a->d[10],
                 a->n[11], valid[3], a->s[7], gw[7], gb[7]);
  backward_layer(&u->layer[6], a->a[9], a->a[10], a->d[10], a->d[9], a->n[10],
                 valid[3], a->s[6], gw[6], gb[6]);
  backward_pool(a->d[9], a->d[8], a->c[8], a->n[8], u->config.pool, valid[2],
                valid[3]);
  backward_layer(&u->layer[5], a->a[7], a->a[8], a->d[8], a->d[7], a->n[8],
                 valid[2], a->s[5], gw[5], gb[5]);
  backward_layer(&u->layer[4], a->a[6], a->a[7], a->d[7], a->d[6], a->n[7],
                 valid[2], a->s[4], gw[4], gb[4]);
  backward_pool(a->d[6], a->d[5], a->c[5], a->n[5], u->config.pool, valid[1],
                valid[2]);
  backward_layer(&u->layer[3], a->a[4], a->a[5], a->d[5], a->d[4], a->n[5],
                 valid[1], a->s[3], gw[3], gb[3]);
  backward_layer(&u->layer[2], a->a[3], a->a[4], a->d[4], a->d[3], a->n[4],
                 valid[1], a->s[2], gw[2], gb[2]);
  backward_pool(a->d[3], a->d[2], a->c[2], a->n[2], u->config.pool, valid[0],
                valid[1]);
  backward_layer(&u->layer[1], a->a[1], a->a[2], a->d[2], a->d[1], a->n[2],
                 valid[0], a->s[1], gw[1], gb[1]);
  backward_layer(&u->layer[0], a->a[0], a->a[1], a->d[1], NULL, a->n[1],
                 valid[0], a->s[0], gw[0], gb[0]);
}

/* Independent per-thread scratch: activations, probabilities, gradients
 * and Soft F1 statistics. */
typedef struct {
  Workspace a;
  float *p, *gw[UNET_LAYER_COUNT], *gb[UNET_LAYER_COUNT];
  double sum_p[ANNO_CLASSES], sum_y[ANNO_CLASSES], inter[ANNO_CLASSES], count;
} Worker;

static void worker_free(Worker *w) {
  workspace_free(&w->a);
  free(w->p);
  for (int l = 0; l < UNET_LAYER_COUNT; ++l) {
    free(w->gw[l]);
    free(w->gb[l]);
  }
}

static int worker_init(const Unet *u, Worker *w) {
  memset(w, 0, sizeof(*w));
  int ok = workspace_init(u, &w->a);
  w->p = (float *)malloc((size_t)ANNO_UNET_WINDOW * ANNO_CLASSES *
                         sizeof(float));
  ok = ok && w->p;
  for (int l = 0; ok && l < u->nlayer; ++l) {
    const Layer *y = &u->layer[l];
    w->gw[l] = (float *)calloc((size_t)y->in * y->out * y->k, sizeof(float));
    w->gb[l] = (float *)calloc((size_t)y->out, sizeof(float));
    ok = w->gw[l] && w->gb[l];
  }
  if (!ok)
    worker_free(w);
  return ok;
}

/**
 * @brief Exact gradient of the batch Soft F1 loss.
 *
 * Soft F1 is a ratio of sums, so its statistics are pooled over every valid
 * position of the whole batch before any gradient is formed.  Pass 1 pools
 * the statistics; pass 2 recomputes the forward activations and
 * backpropagates.  Windows are distributed over threads, each with its own
 * workspace and gradient buffers, reduced at the end in thread order.
 * @param u Model to differentiate.
 * @param x Channel-major input batch.
 * @param y Channel-major one-hot target batch.
 * @param m Position mask batch.
 * @param batch Number of batch elements.
 * @param gw Zeroed per-layer weight gradients, accumulated in place.
 * @param gb Zeroed per-layer bias gradients, accumulated in place.
 * @param loss Receives the batch loss.
 * @return 0 on success, EINVAL without valid positions, ENOMEM on failure.
 */
static int window_valid_length(const float *mask) {
  int valid = ANNO_UNET_WINDOW;
  while (valid && mask[valid - 1] <= 0.0f)
    --valid;
  return valid;
}

static int batch_gradient(const Unet *u, const float *x, const float *y,
                          const float *m, size_t batch, float **gw, float **gb,
                          double *loss) {
  int threads = 1;
#ifdef _OPENMP
  threads = omp_get_max_threads();
#endif
  if ((size_t)threads > batch)
    threads = (int)batch;
  Worker *work = (Worker *)calloc((size_t)threads, sizeof(*work));
  int ready = 0;
  while (work && ready < threads && worker_init(u, &work[ready]))
    ++ready;
  if (ready < threads) {
    for (int t = 0; t < ready; ++t)
      worker_free(&work[t]);
    free(work);
    return ENOMEM;
  }
  const size_t xs = (size_t)ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW,
               ys = (size_t)ANNO_CLASSES * ANNO_UNET_WINDOW;
#pragma omp parallel for num_threads(threads) schedule(static)
  for (size_t b = 0; b < batch; ++b) {
    int t = 0;
#ifdef _OPENMP
    t = omp_get_thread_num();
#endif
    Worker *w = &work[t];
    const float *mask = m + b * ANNO_UNET_WINDOW;
    int valid = window_valid_length(mask);
    if (!valid)
      continue;
    run_forward(u, x + b * xs, valid, &w->a, w->p);
    for (int i = 0; i < valid; ++i) {
      if (mask[i] <= 0.0f)
        continue;
      w->count += 1.0;
      for (int j = 0; j < ANNO_CLASSES; ++j) {
        double target = y[b * ys + (size_t)j * ANNO_UNET_WINDOW + i];
        double probability = w->p[(size_t)i * ANNO_CLASSES + j];
        w->sum_p[j] += probability;
        w->sum_y[j] += target;
        w->inter[j] += probability * target;
      }
    }
  }
  double sum_p[ANNO_CLASSES] = {0}, sum_y[ANNO_CLASSES] = {0},
         inter[ANNO_CLASSES] = {0}, count = 0;
  for (int t = 0; t < threads; ++t) {
    count += work[t].count;
    for (int j = 0; j < ANNO_CLASSES; ++j) {
      sum_p[j] += work[t].sum_p[j];
      sum_y[j] += work[t].sum_y[j];
      inter[j] += work[t].inter[j];
    }
  }
  int rc = count ? 0 : EINVAL;
  double numerator[ANNO_CLASSES], denominator[ANNO_CLASSES], total = 0;
  for (int j = 0; j < ANNO_CLASSES; ++j) {
    numerator[j] = 2.0 * inter[j] + ANNO_EPSILON;
    denominator[j] = sum_p[j] + sum_y[j] + ANNO_EPSILON;
    total += 1.0 - numerator[j] / denominator[j];
  }
  *loss = total / ANNO_CLASSES;
  if (!rc) {
#pragma omp parallel for num_threads(threads) schedule(static)
    for (size_t b = 0; b < batch; ++b) {
      int t = 0;
#ifdef _OPENMP
      t = omp_get_thread_num();
#endif
      Worker *w = &work[t];
      const float *mask = m + b * ANNO_UNET_WINDOW;
      int valid = window_valid_length(mask);
      float *d = w->a.d[21];
      memset(d, 0, (size_t)ANNO_CLASSES * ANNO_UNET_WINDOW * sizeof(float));
      if (!valid)
        continue;
      run_forward(u, x + b * xs, valid, &w->a, w->p);
      for (int i = 0; i < valid; ++i) {
        if (mask[i] <= 0.0f)
          continue;
        double dp[ANNO_CLASSES], softmax_dot = 0;
        for (int j = 0; j < ANNO_CLASSES; ++j) {
          double target = y[b * ys + (size_t)j * ANNO_UNET_WINDOW + i];
          /* d/dp_j of [1 - (2I_j + e)/(P_j + Y_j + e)] / C */
          dp[j] = -(2.0 * target * denominator[j] - numerator[j]) /
                  (denominator[j] * denominator[j] * ANNO_CLASSES);
          softmax_dot += dp[j] * w->p[(size_t)i * ANNO_CLASSES + j];
        }
        for (int j = 0; j < ANNO_CLASSES; ++j)
          d[(size_t)j * ANNO_UNET_WINDOW + i] =
              (float)(w->p[(size_t)i * ANNO_CLASSES + j] *
                      (dp[j] - softmax_dot));
      }
      window_backward(u, &w->a, w->gw, w->gb, valid);
    }
    for (int l = 0; l < u->nlayer; ++l) {
      const Layer *q = &u->layer[l];
      size_t weights = (size_t)q->in * q->out * q->k;
      for (int t = 0; t < threads; ++t) {
        for (size_t z = 0; z < weights; ++z)
          gw[l][z] += work[t].gw[l][z];
        for (int z = 0; z < q->out; ++z)
          gb[l][z] += work[t].gb[l][z];
      }
    }
  }
  for (int t = 0; t < threads; ++t)
    worker_free(&work[t]);
  free(work);
  return rc;
}

/**
 * @brief Train one batch using the masked batch Soft F1 objective.
 * @param u Model to update.
 * @param x Channel-major input batch.
 * @param y Channel-major one-hot target batch.
 * @param m Position mask batch.
 * @param batch Number of batch elements.
 * @param lr Adam learning rate.
 * @param loss Receives the batch loss; may be NULL.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_train_batch(Unet *u, const float *x, const float *y, const float *m,
                     size_t batch, float lr, float *loss) {
  if (!u || !x || !y || !m || !batch || !isfinite(lr) || lr <= 0)
    return EINVAL;
  float **gw = (float **)calloc((size_t)u->nlayer, sizeof(*gw));
  float **gb = (float **)calloc((size_t)u->nlayer, sizeof(*gb));
  int rc = gw && gb ? 0 : ENOMEM;
  for (int l = 0; !rc && l < u->nlayer; ++l) {
    size_t weights = (size_t)u->layer[l].in * u->layer[l].out * u->layer[l].k;
    gw[l] = (float *)calloc(weights, sizeof(float));
    gb[l] = (float *)calloc((size_t)u->layer[l].out, sizeof(float));
    if (!gw[l] || !gb[l])
      rc = ENOMEM;
  }
  double total = 0;
  if (!rc)
    rc = batch_gradient(u, x, y, m, batch, gw, gb, &total);
  if (!rc) {
    ++u->optimizer_step;
    double step = (double)u->optimizer_step;
    double bias_correction1 = 1.0 - pow((double)ANNO_ADAM_BETA1, step);
    double bias_correction2 = 1.0 - pow((double)ANNO_ADAM_BETA2, step);
    double step_size = (double)lr * sqrt(bias_correction2) / bias_correction1;
    for (int l = 0; l < u->nlayer; ++l) {
      Layer *layer = &u->layer[l];
      size_t weights = (size_t)layer->in * layer->out * layer->k;
      adam_update(layer->w, gw[l], layer->mw, layer->vw, 1, weights,
                  step_size);
      adam_update(layer->b, gb[l], layer->mb, layer->vb, 1,
                  (size_t)layer->out, step_size);
    }
    if (loss)
      *loss = (float)total;
  }
  if (gw && gb)
    free_gradients(gw, gb, u->nlayer);
  else {
    free(gw);
    free(gb);
  }
  return rc;
}
/**
 * @brief Write a UNet model to a binary file.
 * @param path Model file path.
 * @param u Model to write.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_save(const char *path, const Unet *u) {
  if (!path || !u)
    return EINVAL;
  FILE *f = fopen(path, "wb");
  if (!f)
    return errno;
  int ok = fwrite(&u->config, sizeof(u->config), 1, f) == 1;
  for (int i = 0; ok && i < u->nlayer; ++i) {
    size_t weights = (size_t)u->layer[i].in * u->layer[i].out * u->layer[i].k;
    ok = fwrite(u->layer[i].w, weights * sizeof(float), 1, f) == 1 &&
         fwrite(u->layer[i].b, (size_t)u->layer[i].out * sizeof(float), 1, f) ==
             1;
  }
  int closed = fclose(f) == 0;
  return ok && closed ? 0 : EIO;
}

/**
 * @brief Read a UNet model from a binary file.
 * @param path Model file path.
 * @param out Receives the loaded model.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_load(const char *path, Unet **out) {
  if (!path || !out)
    return EINVAL;
  *out = NULL;
  FILE *f = fopen(path, "rb");
  if (!f)
    return errno;
  UnetConfig config;
  int ok = fread(&config, sizeof(config), 1, f) == 1;
  Unet *u = ok ? unet_create(&config) : NULL;
  for (int i = 0; ok && u && i < u->nlayer; ++i) {
    size_t weights = (size_t)u->layer[i].in * u->layer[i].out * u->layer[i].k;
    ok = fread(u->layer[i].w, weights * sizeof(float), 1, f) == 1 &&
         fread(u->layer[i].b, (size_t)u->layer[i].out * sizeof(float), 1, f) ==
             1;
  }
  fclose(f);
  if (!ok) {
    release(u);
    return EINVAL;
  }
  *out = u;
  return 0;
}
