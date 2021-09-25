#include "Bint.h"
#include <functional>
#include <map>

namespace Algorithm {

Bint factorial(uint32_t n);

// 2-adic valuation of n.
// Returns a pair containing \\nu_2(n) and n / 2 ** \\nu_2(n)
std::pair<uint32_t, Bint> two_adic_val(Bint n);

// If n = a ** b where a is minimal, returns {a, b}.
// If n is not a perfect power, then b = 1, so this returns {n, 1}.
std::pair<Bint, uint32_t> factor_perfect_power(Bint n);

// Miller Rabin primality test
bool is_probable_prime(Bint n, uint64_t rounds = 40);

// Given an n by m matrix, find a subset of rows that sum to 0 in GF(2)
// outputs list of linearly independent subsets (hopefully the maximum amount)
// Algorithm from:
// https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
std::vector<std::vector<int>>
gf2_gaussian_elimination(const std::vector<std::vector<uint32_t>> &matrix);

// Wheel augmented sieve of Eratosthenes with prime basis 2, 3, 5, 7.
//
// Returns a vector of all primes less than n.
std::vector<uint32_t> primes_less_than(const uint32_t &n);

// Implementation of Euler's Criterion
int legendre_symbol(const Bint &a, const Bint &p);

// Power of two order, i.e. the least n such that t ** (2 ** n) is 1 mod p
uint32_t pow2_ord(const Bint &t, const Bint &p);

// Tonelli Shanks Algorithm:
// Finds a square root R of n mod p; the other root will be -R.
//
// Assumes that p is prime; throws std::domain_error if n is not a
// quadratic residue modulo p.
Bint square_root_modulo_prime(const Bint &n, const Bint &p);

// Calculates a square root R of n mod p ** exp; the other root will be -R.
// Optionally pass a precomputed value of pre_comp that is a square root
// of n more p ** (exp - 1) when exp > 1.
//
// Assumes that p is prime; throws std::domain_error if n is not a
// quadratic residue modulo p.
Bint square_root_modulo_prime_power(const Bint &n, const Bint &p,
                                    const uint32_t &exp,
                                    const Bint &pre_comp = 0);

// Returns a map (p -> \\nu_p(n)) of the prime factorization of n
// given a function to find nontrivial factors of n. See pollard_rho for a
// working example.
std::map<Bint, uint32_t>
prime_factors(const Bint &n,
              std::function<Bint(const Bint &)> &&find_nontrivial_factor);

// Pollard's rho prime factorization algorithm.
// If n is composite, returns a proper factor > 1 of n.
// Otherwise, returns 1.
Bint pollard_rho(const Bint &n);

// Quadratic sieve factorization algorithm.
Bint quadratic_sieve(const Bint &n);

} // namespace Algorithm