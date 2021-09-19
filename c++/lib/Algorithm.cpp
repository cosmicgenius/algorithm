#include "Bint.h"
#include <functional>
#include <map>
#include <math.h>
#include <queue>
#include <string.h>

namespace Algorithm {

Bint factorial(uint32_t n) {
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

std::vector<std::vector<int>>
gf2_gaussian_elimination(std::vector<std::vector<bool>> matrix) {
  std::vector<std::vector<int>> res;

  size_t n, m;
  n = matrix.size();
  if (n == 0) {
    return res;
  }
  m = matrix[0].size();
  if (m == 0) {
    return res;
  }

  for (std::vector<bool> v : matrix) {
    if (v.size() != m) {
      return res;
    }
  }

  bool marked[n];
  memset(marked, false, sizeof(marked));

  int row_marked[m];

  for (int c = 0; c < m; c++) {
    int pivot = -1;
    for (int r = 0; r < n; r++) {
      if (matrix[r][c]) {
        pivot = r;
        break;
      }
    }

    row_marked[c] = pivot;

    if (pivot != -1) {
      marked[pivot] = true;

      for (int c2 = 0; c2 < m; c2++) {
        if (matrix[pivot][c2] && c2 != c) {
          for (int r = 0; r < n; r++) {
            matrix[r][c2] = matrix[r][c2] ^ matrix[r][c];
          }
        }
      }
    }
  }

  for (int r = 0; r < n; r++) {
    if (!marked[r]) {
      std::vector<int> needed_rows(1, r);

      for (int c = 0; c < m; c++) {
        if (matrix[r][c]) {
          needed_rows.push_back(row_marked[c]);
        }
      }
      res.push_back(needed_rows);
    }
  }
  return res;
}

std::vector<uint32_t> primes_less_than(const uint32_t &n) {
  if (n <= 2) {
    return {};
  }
  if (n <= 3) {
    return {2};
  }
  if (n <= 5) {
    return {2, 3};
  }
  if (n <= 7) {
    return {2, 3, 5};
  }

  uint32_t wheel_residues[48];       // totient(210) = 48
  std::map<uint32_t, uint32_t> back; // opposite direction of wheel_residues
  int min_res = -1;

  int residue = 0;

  for (int i = 1; i < 210; i++) {
    if (i % 2 != 0 && i % 3 != 0 && i % 5 != 0 && i % 7 != 0) {
      back[i] = residue;

      if (min_res == -1 && i > n % 210) {
        min_res = i;
      }
      wheel_residues[residue++] = i;
    }
  }

  if (min_res == -1) {
    min_res = 1;
  }

  uint32_t rootn = (uint32_t)sqrt(n);
  uint32_t q = rootn / 210;

  uint32_t P = 48 * (n / 210) + back[min_res];

  bool is_prime[P];
  memset(is_prime, true, sizeof(is_prime));

  // 1 is not prime
  is_prime[0] = false;

  auto get_is_prime = [&is_prime, &back](uint32_t x) {
    return is_prime[48 * (x / 210) + back[x % 210]];
  };

  auto set_is_prime = [&is_prime, &back](uint32_t x, uint32_t v) {
    is_prime[48 * (x / 210) + back[x % 210]] = v;
  };

  uint32_t qn = n / 210;

  for (uint32_t i = 0; i <= q; i++) {
    for (uint32_t j_idx1 = 0; j_idx1 < 48; j_idx1++) {

      uint32_t c = i * 210 + wheel_residues[j_idx1];
      if (c > rootn) {
        break;
      }
      if (c == 1) {
        continue;
      }

      if (get_is_prime(c)) {
        // c2 should be at least c
        bool done = false;
        for (uint32_t j_idx2 = j_idx1; j_idx2 < 48; j_idx2++) {
          uint32_t c2 = i * 210 + wheel_residues[j_idx2];
          if (c2 * c >= n) {
            done = true;
            break;
          }
          set_is_prime(c * c2, false);
        }

        if (done) {
          continue;
        }

        for (uint32_t i2 = i + 1; i2 < qn / c; i2++) {
          for (uint32_t j_idx2 = 0; j_idx2 < 48; j_idx2++) {
            set_is_prime((i2 * 210 + wheel_residues[j_idx2]) * c, false);
          }
        }

        for (uint32_t j_idx2 = 0; j_idx2 < 48; j_idx2++) {
          uint32_t c2 = qn / c * 210 + wheel_residues[j_idx2];
          if (c2 * c >= n) {
            break;
          }
          set_is_prime(c2 * c, false);
        }
      }
    }
  }

  std::vector<uint32_t> primes = {2, 3, 5, 7};
  for (uint32_t x = 0; x < P; x++) {
    if (is_prime[x]) {
      uint32_t p = 210 * (x / 48) + wheel_residues[x % 48];

      if (p >= n) {
        break;
      }

      primes.push_back(p);
    }
  }
  return primes;
}

int legendre_symbol(const Bint &a, const Bint &p) {
  Bint a_mod_p = a % p;

  if (p == 2) {
    if (a_mod_p == 1) {
      return 1;
    }
    return 0;
  }

  if (a_mod_p == 0) {
    return 0;
  }

  Bint power = Bint::pow(a_mod_p, (p - 1) >> 1, p);
  if (power == 1) {
    return 1;
  }
  return -1;
}

uint32_t pow2_ord(const Bint &t, const Bint &p) {
  Bint t_mod_p = t % p;
  if (t_mod_p == 0) {
    return 0;
  }

  Bint v = t_mod_p;
  for (uint32_t o = 0;; o++) {
    if (v == 1) {
      return o;
    }
    v = v * v % p;
  }
}

Bint square_root_modulo_prime(const Bint &n, const Bint &p) {
  Bint n_mod_p = (n % p + p) % p;
  if (n_mod_p == 0) {
    return 0;
  }
  if (p == 2) {
    return n_mod_p;
  }

  if (legendre_symbol(n, p) == -1) {
    throw std::domain_error(
        n.to_string() + " is not a quadratic residue modulo " + p.to_string());
  }

  if (!((p + 1) & 3).to_bool()) {
    return Bint::pow(n, (p + 1) >> 2, p);
  }

  std::pair<uint32_t, Bint> res = two_adic_val(p - 1);
  uint32_t S = res.first;
  Bint Q = res.second;

  // std::cout << S << " " << Q << std::endl;

  Bint z = 2;
  while (legendre_symbol(z, p) == 1) {
    z = Bint::rand(p - 3) + 2;
  }
  // std::cout << z << std::endl;

  Bint L = Bint::pow(n, (Q - 1) >> 1, p);

  uint32_t M = S;
  Bint c = Bint::pow(z, Q, p);
  Bint R = L * n % p;
  Bint t = L * R % p;

  // std::cout << M << " " << c << " " << R << " " << t << std::endl;

  while (t != 1) {
    uint32_t o = pow2_ord(t, p);
    // std::cout << o << std::endl;
    Bint b = Bint::pow(c, Bint::pow(2, M - o - 1, p - 1), p);

    M = o;
    c = b * b % p;
    t = t * c % p;
    R = R * b % p;
  }

  return R;
}

Bint square_root_modulo_prime_power(const Bint &n, const Bint &p,
                                    const uint32_t &exp,
                                    const Bint &pre_comp = 0) {
  if (exp == 1) {
    return square_root_modulo_prime(n, p);
  }

  Bint n_mod_p = (n % p + p) % p;
  if (n_mod_p == 0) {
    return 0;
  }

  Bint R = (pre_comp == 0) ? square_root_modulo_prime_power(n, p, exp - 1)
                           : pre_comp;

  Bint pow_p = Bint::pow(p, exp - 1);
  std::cout << exp - 1 << std::endl;
  if (p == 2) {
    if (((n - R * R) / pow_p & 1).to_bool()) {
      throw std::domain_error(n.to_string() +
                              " is not a quadratic residue modulo " +
                              p.to_string() + " ** " + std::to_string(exp));
    }
    R += pow_p * ((((n - R * R) / pow_p) >> 1) % p);
  } else {
    R += pow_p *
         ((((n - R * R) / pow_p) * Bint::pow(R * 2, p - 2, p) % p + p) % p);
  }
  return R;
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