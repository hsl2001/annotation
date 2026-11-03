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

enum {
  ENCODER_OPERATIONS = 2,
  DECODER_OPERATIONS = 2,
  OUTPUT_OPERATIONS = 1,
  /* Two convolutions, one pool, one skip, one upsample. */
  WORKSPACES_PER_LEVEL = 5,
  /* Input, bottleneck, output, and alignment storage. */
  EXTRA_WORKSPACES = 4,
  CONV_OPERATION = 1,
  POOL_OPERATION = 2,
  CONCAT_OPERATION = 3,
  UPSAMPLE_OPERATION = 4,
  MIN_DIMENSION = 1,
  LINEAR_ACTIVATION = 0,
  RELU_ACTIVATION = 1,
  OUTPUT_KERNEL = 1,
  NO_LAYER = -1,
  NO_SKIP = -1,
  INITIAL_DEPTH = 4,
  INITIAL_KERNEL = 9,
  INITIAL_POOL = 2,
  INITIAL_CHANNEL_0 = 16,
  INITIAL_CHANNEL_1 = 32,
  INITIAL_CHANNEL_2 = 64,
  INITIAL_CHANNEL_3 = 128,
};


#define MAX_LAYERS                                                             \
  (ANNO_MAX_DEPTH * (ENCODER_OPERATIONS + DECODER_OPERATIONS) +                \
   OUTPUT_OPERATIONS)
#define MAX_WORKSPACES                                                         \
  (ANNO_MAX_DEPTH * WORKSPACES_PER_LEVEL + EXTRA_WORKSPACES)

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
  int depth, nlayer;
  uint64_t optimizer_step;
  Layer layer[MAX_LAYERS];
};

/**
 * @brief Validate a UNet configuration.
 * @param c Configuration to validate.
 * @return Nonzero when the configuration is valid, otherwise zero.
 */
static int good_config(const UnetConfig *c) {
  if (!c || c->depth < MIN_DIMENSION || c->depth > ANNO_MAX_DEPTH ||
      c->kernel < MIN_DIMENSION || !(c->kernel & MIN_DIMENSION) ||
      c->pool < MIN_DIMENSION)
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
  for (size_t i = 0; i < weights; ++i) {
    l->w[i] = 0.01 * anno_random_symmetric(seed);
  }
  return 1;
}

/**
 * @brief Apply one Adam update to a parameter vector.
 * @param parameter Parameters to update.
 * @param gradient Accumulated gradients.
 * @param first_moment Adam first-moment state.
 * @param second_moment Adam second-moment state.
 * @param count Number of valid samples contributing to gradient.
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
  u->depth = c->depth;
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
 * @brief Apply one convolution layer.
 * @param l Layer parameters.
 * @param x Channel-major input activation.
 * @param y Channel-major output activation.
 * @param n Number of positions.
 * @param relu Whether to apply ReLU to the output.
 * @return None.
 */
static void conv(const Layer *l, const float *x, float *y, int n, int relu) {
  int r = l->k / 2;
#pragma omp parallel for collapse(2)
  for (int o = 0; o < l->out; ++o)
    for (int p = 0; p < n; ++p) {
      float v = l->b[o];
      for (int i = 0; i < l->in; ++i)
        for (int q = 0; q < l->k; ++q) {
          int z = p + q - r;
          if (z >= 0 && z < n)
            v +=
                l->w[((size_t)o * l->in + i) * l->k + q] * x[(size_t)i * n + z];
        }
      y[(size_t)o * n + p] = relu && v < 0 ? 0 : v;
    }
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

typedef struct {
  float *a[MAX_WORKSPACES], *d[MAX_WORKSPACES];
  int n[MAX_WORKSPACES], c[MAX_WORKSPACES], count;
  float *concat[MAX_WORKSPACES];
  struct {
    int type, layer, in, out, aux, skip;
  } op[MAX_WORKSPACES];
  int n_op;
} Workspace;

/**
 * @brief Free all tensors and operation metadata in a workspace.
 * @param a Workspace to clear.
 * @return None.
 */
static void workspace_free(Workspace *a) {
  for (int i = 0; i < a->count; ++i) {
    free(a->a[i]);
    free(a->d[i]);
  }
  for (int i = 0; i < MAX_WORKSPACES; ++i)
    free(a->concat[i]);
  memset(a, 0, sizeof(*a));
}

/**
 * @brief Allocate one activation and gradient tensor pair.
 * @param a Workspace receiving the tensors.
 * @param c Number of channels.
 * @param n Number of positions.
 * @return Nonzero on success, otherwise zero on allocation failure.
 */
static int workspace_add(Workspace *a, int c, int n) {
  int i = a->count++;
  size_t z = (size_t)c * n;
  a->c[i] = c;
  a->n[i] = n;
  a->a[i] = (float *)calloc(z, sizeof(float));
  a->d[i] = (float *)calloc(z, sizeof(float));
  return a->a[i] && a->d[i];
}

/**
 * @brief Record one forward operation for backpropagation.
 * @param a Workspace receiving the operation.
 * @param type Operation kind.
 * @param layer Associated convolution layer, or NO_LAYER.
 * @param in Input workspace index.
 * @param out Output workspace index.
 * @param aux Operation-specific parameter.
 * @param skip Skip-connection workspace index, or NO_SKIP.
 * @return None.
 */
static void record(Workspace *a, int type, int layer, int in, int out, int aux,
                   int skip) {
  a->op[a->n_op].type = type;
  a->op[a->n_op].layer = layer;
  a->op[a->n_op].in = in;
  a->op[a->n_op].out = out;
  a->op[a->n_op].aux = aux;
  a->op[a->n_op].skip = skip;
  ++a->n_op;
}

/**
 * @brief Run the UNet forward pass and softmax output.
 * @param u Model to evaluate.
 * @param x Channel-major input tensor.
 * @param len Number of valid input positions.
 * @param a Workspace used by the pass.
 * @param prob Position-major class probabilities.
 * @return 0 on success, or ENOMEM on allocation failure.
 */
static int run_forward(const Unet *u, const float *x, int64_t len, Workspace *a,
                       float *prob) {
  int n = ANNO_UNET_WINDOW, ai = 0, li = 0;
  if (!workspace_add(a, ANNO_INPUT_CHANNELS, n))
    return ENOMEM;
  memcpy(a->a[0], x, (size_t)ANNO_INPUT_CHANNELS * n * sizeof(float));
  int skip[MAX_LAYERS], sn = 0;
  for (int lev = 0; lev < u->depth; ++lev) {
    if (!workspace_add(a, u->config.channels[lev], n))
      return ENOMEM;
    conv(&u->layer[li++], a->a[ai], a->a[a->count - 1], n, RELU_ACTIVATION);
    ai = a->count - 1;
    record(a, CONV_OPERATION, li - 1, ai - 1, ai, RELU_ACTIVATION, NO_SKIP);
    if (!workspace_add(a, u->config.channels[lev], n))
      return ENOMEM;
    conv(&u->layer[li++], a->a[ai], a->a[a->count - 1], n, RELU_ACTIVATION);
    ai = a->count - 1;
    record(a, CONV_OPERATION, li - 1, ai - 1, ai, RELU_ACTIVATION, NO_SKIP);
    skip[sn++] = ai;
    if (lev + 1 < u->depth) {
      int pn = n / u->config.pool;
      if (!workspace_add(a, u->config.channels[lev], pn))
        return ENOMEM;
      pool(a->a[ai], a->a[a->count - 1], u->config.channels[lev], n,
           u->config.pool);
      ai = a->count - 1;
      record(a, POOL_OPERATION, NO_LAYER, ai - 1, ai, u->config.pool, NO_SKIP);
      n = pn;
    }
  }
  for (int lev = u->depth - 2; lev >= 0; --lev) {
    int oldn = n, f = u->config.pool;
    n *= f;
    if (!workspace_add(a, u->config.channels[lev + 1], n))
      return ENOMEM;
    float *up = a->a[a->count - 1];
    upsample(a->a[ai], up, u->config.channels[lev + 1], oldn, f);
    record(a, UPSAMPLE_OPERATION, NO_LAYER, ai, a->count - 1, f, NO_SKIP);
    int cc = u->config.channels[lev + 1] + u->config.channels[lev];
    float *cat = (float *)malloc((size_t)cc * n * sizeof(float));
    if (!cat)
      return ENOMEM;
    memcpy(cat, up, (size_t)u->config.channels[lev + 1] * n * sizeof(float));
    memcpy(cat + (size_t)u->config.channels[lev + 1] * n, a->a[skip[lev]],
           (size_t)u->config.channels[lev] * n * sizeof(float));
    if (!workspace_add(a, cc, n)) {
      free(cat);
      return ENOMEM;
    }
    memcpy(a->a[a->count - 1], cat, (size_t)cc * n * sizeof(float));
    free(cat);
    ai = a->count - 1;
    record(a, CONCAT_OPERATION, NO_LAYER, ai - 1, ai,
           u->config.channels[lev + 1], skip[lev]);
    if (!workspace_add(a, u->config.channels[lev], n))
      return ENOMEM;
    conv(&u->layer[li++], a->a[ai], a->a[a->count - 1], n, RELU_ACTIVATION);
    ai = a->count - 1;
    record(a, CONV_OPERATION, li - 1, ai - 1, ai, RELU_ACTIVATION, NO_SKIP);
    if (!workspace_add(a, u->config.channels[lev], n))
      return ENOMEM;
    conv(&u->layer[li++], a->a[ai], a->a[a->count - 1], n, RELU_ACTIVATION);
    ai = a->count - 1;
    record(a, CONV_OPERATION, li - 1, ai - 1, ai, RELU_ACTIVATION, NO_SKIP);
  }
  if (!workspace_add(a, ANNO_CLASSES, n))
    return ENOMEM;
  conv(&u->layer[li], a->a[ai], a->a[a->count - 1], n, 0);
  record(a, CONV_OPERATION, li, ai, a->count - 1, LINEAR_ACTIVATION, NO_SKIP);
  float *z = a->a[a->count - 1];
  int64_t usable = len;
#pragma omp parallel for
  for (int64_t p = 0; p < usable; ++p) {
    float top = -INFINITY, sum = 0;
    for (int j = 0; j < ANNO_CLASSES; ++j) {
      if (z[(size_t)j * n + p] > top)
        top = z[(size_t)j * n + p];
    }
    for (int j = 0; j < ANNO_CLASSES; ++j) {
      prob[(size_t)p * ANNO_CLASSES + j] = expf(z[(size_t)j * n + p] - top);
      sum += prob[(size_t)p * ANNO_CLASSES + j];
    }
    for (int j = 0; j < ANNO_CLASSES; ++j)
      prob[(size_t)p * ANNO_CLASSES + j] /= sum;
  }
  return 0;
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
  int rc = run_forward(u, x, n, &a, p);
  workspace_free(&a);
  return rc;
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
 * @brief Backpropagate through one convolution layer.
 * @param l Layer parameters.
 * @param x Layer input activation.
 * @param y Layer output gradient.
 * @param dx Destination input gradient.
 * @param n Activation length.
 * @param relu Whether the output used ReLU.
 * @param gw Destination weight gradient.
 * @param gb Destination bias gradient.
 * @return None.
 */
static void backward_layer(const Layer *l, const float *x, const float *y,
                           float *dx, int n, int relu, float *gw, float *gb) {
  int r = l->k / 2;
  memset(dx, 0, (size_t)l->in * n * sizeof(float));
#pragma omp parallel for
  for (int o = 0; o < l->out; ++o) {
    float db = 0;
    for (int p = 0; p < n; ++p) {
      float g = y[(size_t)o * n + p];
      if (relu && g < 0)
        g = 0;
      db += g;
    }
    gb[o] += db;
    for (int i = 0; i < l->in; ++i)
      for (int q = 0; q < l->k; ++q) {
        float dw = 0;
        for (int p = 0; p < n; ++p) {
          float g = y[(size_t)o * n + p];
          if (relu && g < 0)
            g = 0;
          int z = p + q - r;
          if (z >= 0 && z < n) {
            dw += g * x[(size_t)i * n + z];
#pragma omp atomic
            dx[(size_t)i * n + z] +=
                g * l->w[((size_t)o * l->in + i) * l->k + q];
          }
        }
        gw[((size_t)o * l->in + i) * l->k + q] += dw;
      }
  }
}
/**
 * @brief Train one batch using the masked Soft F1 objective.
 * @param u Model to update.
 * @param x Channel-major input batch.
 * @param y Channel-major one-hot target batch.
 * @param m Position mask batch.
 * @param batch Number of batch elements.
 * @param lr Adam learning rate.
 * @param loss Receives the mean loss; may be NULL.
 * @return 0 on success, or an errno value otherwise.
 */
int unet_train_batch(Unet *u, const float *x, const float *y, const float *m,
                     size_t batch, float lr, float *loss) {
  if (!u || !x || !y || !m || !batch || !isfinite(lr) || lr <= 0)
    return EINVAL;
  double total = 0;
  double count = 0;
  float **gw = (float **)calloc((size_t)u->nlayer, sizeof(*gw));
  float **gb = (float **)calloc((size_t)u->nlayer, sizeof(*gb));
  if (!gw || !gb) {
    free(gw);
    free(gb);
    return ENOMEM;
  }
  for (int l = 0; l < u->nlayer; ++l) {
    size_t weights = (size_t)u->layer[l].in * u->layer[l].out * u->layer[l].k;
    gw[l] = (float *)calloc(weights, sizeof(float));
    gb[l] = (float *)calloc((size_t)u->layer[l].out, sizeof(float));
    if (!gw[l] || !gb[l]) {
      free_gradients(gw, gb, u->nlayer);
      return ENOMEM;
    }
  }
#pragma omp parallel for schedule(static)
  for (size_t b = 0; b < batch; ++b) {
    float **sgw = (float **)calloc((size_t)u->nlayer, sizeof(*sgw));
    float **sgb = (float **)calloc((size_t)u->nlayer, sizeof(*sgb));
    for (int l = 0; l < u->nlayer; ++l) {
      size_t weights = (size_t)u->layer[l].in * u->layer[l].out * u->layer[l].k;
      sgw[l] = (float *)calloc(weights, sizeof(float));
      sgb[l] = (float *)calloc((size_t)u->layer[l].out, sizeof(float));
    }
    int gradients_ready = 1;
    for (int l = 0; l < u->nlayer; ++l)
      if (!sgw[l] || !sgb[l]) {
        gradients_ready = 0;
        break;
      }
    if (!gradients_ready) {
      free_gradients(sgw, sgb, u->nlayer);
      continue;
    }
    Workspace a = {0};
    float *p = (float *)malloc((size_t)ANNO_UNET_WINDOW * ANNO_CLASSES *
                               sizeof(float));
    if (!p) {
      free_gradients(sgw, sgb, u->nlayer);
      continue;
    }
    int rc =
        run_forward(u, x + b * (size_t)ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW,
                    ANNO_UNET_WINDOW, &a, p);
    if (rc) {
      free(p);
      workspace_free(&a);
      free_gradients(sgw, sgb, u->nlayer);
      continue;
    }
    int out = a.count - 1;
    float sum_p[ANNO_CLASSES] = {0};
    float sum_y[ANNO_CLASSES] = {0};
    float intersection[ANNO_CLASSES] = {0};
    size_t valid = 0;
    for (int64_t i = 0; i < ANNO_UNET_WINDOW; ++i) {
      if (m[b * ANNO_UNET_WINDOW + i] <= 0)
        continue;
      ++valid;
      for (int j = 0; j < ANNO_CLASSES; ++j) {
        float target = y[b * (size_t)ANNO_CLASSES * ANNO_UNET_WINDOW +
                         (size_t)j * ANNO_UNET_WINDOW + i];
        float probability = p[(size_t)i * ANNO_CLASSES + j];
        sum_p[j] += probability;
        sum_y[j] += target;
        intersection[j] += probability * target;
      }
    }
    if (!valid) {
      free(p);
      workspace_free(&a);
      free_gradients(sgw, sgb, u->nlayer);
      continue;
    }
    double sample_loss = 0;
    for (int j = 0; j < ANNO_CLASSES; ++j) {
      float numerator = 2.0f * intersection[j] + ANNO_EPSILON;
      float denominator = sum_p[j] + sum_y[j] + ANNO_EPSILON;
      sample_loss += 1.0 - numerator / denominator;
    }
    sample_loss /= ANNO_CLASSES;
    float **ag = (float **)calloc((size_t)a.count, sizeof(*ag));
    if (!ag) {
      free(p);
      workspace_free(&a);
      free_gradients(sgw, sgb, u->nlayer);
      continue;
    }
    int gradients_allocated = 1;
    for (int z = 0; z < a.count; ++z) {
      ag[z] = (float *)calloc((size_t)a.c[z] * a.n[z], sizeof(float));
      if (!ag[z]) {
        gradients_allocated = 0;
        break;
      }
    }
    if (!gradients_allocated) {
      for (int z = 0; z < a.count; ++z)
        free(ag[z]);
      free(ag);
      free(p);
      workspace_free(&a);
      free_gradients(sgw, sgb, u->nlayer);
      continue;
    }
    for (int64_t i = 0; i < ANNO_UNET_WINDOW; ++i) {
      if (m[b * ANNO_UNET_WINDOW + i] <= 0)
        continue;
      float dp[ANNO_CLASSES];
      float softmax_dot = 0;
      for (int j = 0; j < ANNO_CLASSES; ++j) {
        float target = y[b * (size_t)ANNO_CLASSES * ANNO_UNET_WINDOW +
                         (size_t)j * ANNO_UNET_WINDOW + i];
        float probability = p[(size_t)i * ANNO_CLASSES + j];
        float denominator = sum_p[j] + sum_y[j] + ANNO_EPSILON;
        float numerator = 2.0f * intersection[j] + ANNO_EPSILON;
        dp[j] = -(2.0f * target * denominator - numerator) /
                (denominator * denominator * ANNO_CLASSES);
        softmax_dot += dp[j] * probability;
      }
      for (int j = 0; j < ANNO_CLASSES; ++j) {
        float probability = p[(size_t)i * ANNO_CLASSES + j];
        ag[out][(size_t)j * ANNO_UNET_WINDOW + i] =
            probability * (dp[j] - softmax_dot);
      }
    }
    for (int oi = a.n_op - 1; oi >= 0; --oi) {
      int t = a.op[oi].type, in = a.op[oi].in, outi = a.op[oi].out;
      if (t == 1)
        backward_layer(&u->layer[a.op[oi].layer], a.a[in], ag[outi], ag[in],
                       a.n[outi], a.op[oi].aux, sgw[a.op[oi].layer],
                       sgb[a.op[oi].layer]);
      else if (t == 2) {
        int f = a.op[oi].aux;
        for (int c = 0; c < a.c[in]; ++c)
          for (int q = 0; q < a.n[outi]; ++q)
            for (int z = 0; z < f; ++z)
              ag[in][(size_t)c * a.n[in] + q * f + z] +=
                  ag[outi][(size_t)c * a.n[outi] + q] / f;
      } else if (t == 4) {
        int f = a.op[oi].aux;
        for (int c = 0; c < a.c[in]; ++c)
          for (int q = 0; q < a.n[in]; ++q)
            for (int z = 0; z < f; ++z)
              ag[in][(size_t)c * a.n[in] + q] +=
                  ag[outi][(size_t)c * a.n[outi] + q * f + z];
      } else if (t == 3) {
        int split = a.op[oi].aux, sc = a.c[a.op[oi].skip];
        for (int c = 0; c < split; ++c)
          for (int q = 0; q < a.n[outi]; ++q)
            ag[in][(size_t)c * a.n[in] + q] +=
                ag[outi][(size_t)c * a.n[outi] + q];
        for (int c = 0; c < sc; ++c)
          for (int q = 0; q < a.n[outi]; ++q)
            ag[a.op[oi].skip][(size_t)c * a.n[a.op[oi].skip] + q] +=
                ag[outi][(size_t)(split + c) * a.n[outi] + q];
      }
    }
    for (int z = 0; z < a.count; ++z)
      free(ag[z]);
    free(ag);
    free(p);
    workspace_free(&a);
#pragma omp critical(ann_gradient_reduce)
    {
      total += sample_loss;
      count += 1.0;
      for (int l = 0; l < u->nlayer; ++l) {
        for (size_t z = 0;
             z < (size_t)u->layer[l].in * u->layer[l].out * u->layer[l].k; ++z)
          gw[l][z] += sgw[l][z];
        for (int z = 0; z < u->layer[l].out; ++z)
          gb[l][z] += sgb[l][z];
        free(sgw[l]);
        free(sgb[l]);
      }
    }
    free(sgw);
    free(sgb);
  }
  if (!count) {
    free_gradients(gw, gb, u->nlayer);
    return EINVAL;
  }
  ++u->optimizer_step;
  double step = (double)u->optimizer_step;
  double bias_correction1 = 1.0 - pow((double)ANNO_ADAM_BETA1, step);
  double bias_correction2 = 1.0 - pow((double)ANNO_ADAM_BETA2, step);
  double step_size = (double)lr * sqrt(bias_correction2) / bias_correction1;
  for (int l = 0; l < u->nlayer; ++l) {
    Layer *layer = &u->layer[l];
    size_t weights = (size_t)layer->in * layer->out * layer->k;
    adam_update(layer->w, gw[l], layer->mw, layer->vw, (size_t)count, weights,
                step_size);
    adam_update(layer->b, gb[l], layer->mb, layer->vb, (size_t)count,
                (size_t)layer->out, step_size);
    free(gw[l]);
    free(gb[l]);
  }
  free(gw);
  free(gb);
  if (loss)
    *loss = (float)(total / count);
  return 0;
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
