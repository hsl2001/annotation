#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>

_Noreturn void bye(const char *format, ...);
void cry(const char *format, ...);
size_t smul(size_t a, size_t b);
void *scalloc(size_t count, size_t size);
void *srealloc(void *old, size_t count, size_t size);
char *sstrdup(const char *text);
double whattime(void);
int base_index(char base);
char complement(char base);

#endif
