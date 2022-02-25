#ifndef Rand_h
#define Rand_h


#include <cstdint>


struct Rand {
    uint64_t u, v, w;
    Rand(uint64_t j) : v(4101842887655102017), w(1) {
        u = j ^ v;
        uint64();
        v = u;
        uint64();
        w = v;
        uint64();
    }
    inline uint64_t uint64() {
        u = u * 2862933555777941757 + 7046029254386353087;
        v ^= v >> 17;
        v ^= v << 31;
        v ^= v >> 8;
        w = 4294957665 * (w & 0xffffffff) + (w >> 32);
        uint64_t x = u ^ (u << 21);
        x ^= x >> 35;
        x ^= x << 4;
        return (x + v) ^ w;
    }
    inline double float64() { return 5.42101086242752217E-20 * uint64(); }
    inline uint32_t uint32() { return (uint32_t)uint64(); }
    inline uint64_t roll100() { return ((uint64() >> 7) * 100) >> (64-7); }
    inline uint64_t rollN(unsigned nLessThan65536) {
        return ((uint64() >> 16) * nLessThan65536) >> (64-16);
    }
};


#endif
