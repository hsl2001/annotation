/**
 * @file test.c
 * @brief Unit and integration tests
 */

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/wait.h>
#include <unistd.h>

static int failures = 0;

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
// TODO!!!
return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
