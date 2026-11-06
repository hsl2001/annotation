/**
 * @file utils.c
 * @brief utility functions
 */

#ifndef UTILS_C
#define UTILS_C

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/**
 * @brief Report an error and terminate.
 * @param[in] format: printf-style diagnostic format.
 * @param[in] ... Values referenced by @p format.
 * @return !!This function does not return!!
 */
_Noreturn void bye(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("surf: error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

/**
 * @brief Write an informational log message to standard error.
 * @param[in] format: printf-style message format.
 * @param[in] ... Values referenced by @p format.
 * @return None.
 */
void cry(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("surf: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
}

/**
 * @brief Multiply two sizes, failing on overflow.
 * @param[in] a First factor.
 * @param[in] b Second factor.
 * @return Product of @p a and @p b.
 */
size_t smul(size_t a, size_t b) {
  if (a && b > SIZE_MAX / a)
    bye("size overflow");
  return a * b;
}

/**
 * @brief Safely allocate zero-initialized memory or terminate on failure.
 * @param[in] count: Number of elements.
 * @param[in] size: Size of each element.
 * @return Allocated memory.
 */
void *scalloc(size_t count, size_t size) {
  void *memory = calloc(count ? count : 1, size ? size : 1);
  if (!memory || (size && count > SIZE_MAX / size))
    bye("out of memory");
  return memory;
}

/**
 * @brief Safely resize an allocation or terminate on failure.
 * @param[in] old: Existing allocation, or NULL.
 * @param[in] count: Number of elements.
 * @param[in] size: Size of each element.
 * @return Resized allocation.
 */
void *srealloc(void *old, size_t count, size_t size) {
  void *memory = realloc(old, smul(count ? count : 1, size));
  if (!memory)
    bye("out of memory");
  return memory;
}

/**
 * @brief Safely duplicate a string or terminate on allocation failure.
 * @param[in] text: String to duplicate.
 * @return Caller-owned duplicate.
 */
char *sstrdup(const char *text) {
  char *copy = strdup(text);
  if (!copy)
    bye("out of memory");
  return copy;
}

/**
 * @brief Read the monotonic clock.
 * @return Clock value in seconds.
 */
double whattime(void) {
  struct timespec time;
  clock_gettime(CLOCK_MONOTONIC, &time);
  return (double)time.tv_sec + 1e-9 * (double)time.tv_nsec;
}

/**
 * @brief Convert A/C/G/T to a base-channel index.
 * @param[in] base: Nucleotide character.
 * @return Index 0..3, or -1 for another character.
 */
int base_index(char base) {
  switch (base) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}

/**
 * @brief Return the DNA complement.
 * @param[in] base: Nucleotide character.
 * @return Complementary nucleotide, or N for an unknown base.
 */
char complement(char base) {
  int index = base_index(base);
  return index < 0 ? 'N' : "TGCA"[index];
}

#endif /* UTILS_C */
