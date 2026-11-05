/**
 * @file util.c
 * @brief Shared diagnostics, allocation, timing, and nucleotide helpers.
 */
#define _POSIX_C_SOURCE 200809L

#include "include/anno_internal.h"

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/**
 * @brief Print an error to stderr and terminate the process.
 * @param[in] format printf-style error format string.
 * @param[in] ... Values referenced by @p format.
 * @return This function does not return.
 */
_Noreturn void anno_fail(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("anno: error: ", stderr);
  vfprintf(stderr, format, arguments);
  fputc('\n', stderr);
  va_end(arguments);
  exit(EXIT_FAILURE);
}

/**
 * @brief Print a log line to stderr.
 * @param[in] format printf-style log format string.
 * @param[in] ... Values referenced by @p format.
 * @return None.
 */
void anno_log(const char *format, ...) {
  va_list arguments;
  va_start(arguments, format);
  fputs("anno: ", stderr);
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
size_t anno_mul_size(size_t a, size_t b) {
  if (a && b > SIZE_MAX / a)
    anno_fail("size overflow");
  return a * b;
}

/**
 * @brief Allocate zero-initialized memory or terminate on failure.
 * @param[in] count Number of elements.
 * @param[in] size Size of each element.
 * @return Allocated zero-initialized memory.
 */
void *anno_xcalloc(size_t count, size_t size) {
  void *memory = calloc(count ? count : 1, size ? size : 1);
  if (!memory || (size && count > SIZE_MAX / size))
    anno_fail("out of memory");
  return memory;
}

/**
 * @brief Resize an allocation or terminate on failure.
 * @param[in] old Existing allocation, or NULL.
 * @param[in] count Number of elements.
 * @param[in] size Size of each element.
 * @return Resized allocation.
 */
void *anno_xrealloc(void *old, size_t count, size_t size) {
  void *memory = realloc(old, anno_mul_size(count ? count : 1, size));
  if (!memory)
    anno_fail("out of memory");
  return memory;
}

/**
 * @brief Duplicate a string or terminate on allocation failure.
 * @param[in] text String to duplicate.
 * @return Newly allocated string.
 */
char *anno_xstrdup(const char *text) {
  char *copy = strdup(text);
  if (!copy)
    anno_fail("out of memory");
  return copy;
}

/**
 * @brief Return monotonic elapsed time in seconds.
 * @return Monotonic clock value in seconds.
 */
double anno_now_seconds(void) {
  struct timespec time;
  clock_gettime(CLOCK_MONOTONIC, &time);
  return (double)time.tv_sec + 1e-9 * (double)time.tv_nsec;
}

/**
 * @brief Convert a nucleotide to its one-hot channel index.
 * @param[in] base Nucleotide character.
 * @return A/C/G/T index in the range 0..3, or -1 for another character.
 */
int anno_base_index(char base) {
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
 * @brief Return the complement of a nucleotide.
 * @param[in] base Nucleotide character.
 * @return Complementary nucleotide, or N for an unknown base.
 */
char anno_complement(char base) {
  int index = anno_base_index(base);
  return index < 0 ? 'N' : "TGCA"[index];
}
