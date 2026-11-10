/* Minimal gb_flip.c implementation for ssxcc */
#include <stdlib.h>
#include "gb_flip.h"

static unsigned long long seed_value = 0;

void gb_init_rand(int seed) {
    seed_value = (unsigned long long)seed;
    if (seed_value == 0) seed_value = 1;
}

unsigned int gb_unif_rand(unsigned int n) {
    if (n == 0) return 0;
    seed_value = seed_value * 6364136223846793005ULL + 1;
    return (unsigned int)((seed_value >> 32) % n);
}
