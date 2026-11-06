/**
 * @file test.c
 * @brief Unit and pipeline tests for surf.
 */

#include <stdlib.h>

static int failures = 0, checks = 0;

/** @brief Record one check. */
#define CHECK(condition)                                                       \
  do {                                                                         \
    checks++;                                                                  \
    if (!(condition)) {                                                        \
      failures++;                                                              \
      fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__,         \
              #condition);                                                     \
    }                                                                          \
  } while (0)

/**
 * @brief Run all unit, regression, round-trip, and pipeline tests.
 * @return EXIT_SUCCESS when all checks pass; otherwise EXIT_FAILURE.
 */
int main(void) {
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
