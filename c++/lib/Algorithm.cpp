#include "Bint.h"
#include <functional>
#include <map>
#include <queue>

namespace Algorithm {

Bint fact(uint32_t n) {
  Bint res = 1;
  for (uint32_t i = 2; i <= n; i++) {
    res *= i;
  }
  return res;
}

std::pair<uint32_t, Bint> two_adic_val(Bint n) {
  for (uint32_t r = 0;; r++) {
    if ((n & 1).to_bool())
      return {r, n};
    n >>= 1;
  }
}

std::pair<Bint, uint32_t> factor_perfect_power(Bint n) {
  uint32_t rmax = n.bit_length() - 1;
  for (int r = rmax; r >= 2; r--) {
    Bint root = Bint::integral_rth_root(n, r);

    if (Bint::pow(root, r) == n) {
      return {root, r};
    }
  }
  return {n, 1};
}

bool is_probable_prime(Bint n, uint64_t rounds = 40) {
  if (!(n & 1).to_bool())
    return n == 2;
  if (n % 3 == 0)
    return n == 3;
  if (n % 5 == 0)
    return n == 5;
  if (n % 7 == 0)
    return n == 7;

  Bint n_minus_one = n - 1;

  std::pair<uint32_t, Bint> rd = two_adic_val(n_minus_one);

  for (int round = 0; round < rounds; round++) {
    Bint a = Bint::rand(n - 3) + 2;

    Bint x = Bint::pow(a, rd.second, n);
    if (x != 1 && x != n_minus_one) {
      bool works = false;

      for (int j = 0; j < rd.first - 1; j++) {
        x = x * x % n;
        if (x == n_minus_one) {
          works = true;
          break;
        }
      }
      if (!works) {
        return false;
      }
    }
  }
  return true;
}

Bint pollard_rho(const Bint &n) {
  if (is_probable_prime(n)) {
    return 1;
  }

  auto poly = [&n](Bint x) { return (x * x + 1) % n; };

  std::pair<Bint, uint32_t> p = factor_perfect_power(n);
  if (p.second != 1) {
    return Bint::pow(p.first, p.second / 2);
  }

  Bint g = n;
  Bint n_minus_two = n - 2;
  while (g == n) {
    Bint turtle = Bint::rand(n_minus_two) + 1;
    Bint hare = turtle;
    g = 1;

    while (g == 1) {
      turtle = poly(turtle);
      hare = poly(poly(hare));
      g = Bint::gcd(turtle - hare, n);
    }
  }
  return Bint::gcd(g, n);
}

std::map<Bint, uint32_t>
prime_factors(const Bint &n,
              std::function<Bint(const Bint &)> &&find_nontrivial_factor) {
  std::priority_queue<Bint> pq;
  pq.push(n);

  Bint one = 1;

  while (true) {
    Bint t = pq.top();
    Bint f = find_nontrivial_factor(t);

    if (f == one) {
      break;
    }
    pq.pop();

    pq.push(f);
    pq.push(t / f);
  }

  std::map<Bint, uint32_t> prime_fact;

  while (!pq.empty()) {
    Bint p = pq.top();
    pq.pop();

    prime_fact[p]++;
  }
  return prime_fact;
}

} // namespace Algorithm