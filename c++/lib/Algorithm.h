#include "Bint.h"

namespace Algorithm {

Bint fact(uint32_t n);

// 2-adic valuation of n.
// Returns a pair containing \\nu_2(n) and n / 2 ** \\nu_2(n)
std::pair<uint32_t, Bint> two_adic_val(Bint n);

// If n = a ** b where a is minimal, returns {a, b}.
// If n is not a perfect power, then b = 1, so this returns {n, 1}.
std::pair<Bint, uint32_t> factor_perfect_power(Bint n);

// Miller Rabin primality test
bool is_probable_prime(Bint n, uint64_t rounds = 1000);

// Pollard's rho prime factorization algorithm.
// If n is composite, returns a proper factor > 1 of n.
// Otherwise, returns 1.
Bint pollard_rho(Bint n);

}