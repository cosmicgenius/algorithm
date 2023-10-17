#include "Bint.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <map>
#include <stack>
#include <vector>

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
  if (n == 1)
    return false;
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
gf2_gaussian_elimination(const std::vector<std::vector<uint32_t>> &matrix) {
  std::vector<std::vector<int>> res;

  const int n = matrix.size();
  const int N = n / 32 + (n % 32 != 0);
  if (n == 0) {
    throw std::invalid_argument("The given matrix has no rows.");
  }
  const int m = matrix[0].size();
  if (m == 0) {
    throw std::invalid_argument("The given matrix has no column.");
  }

  for (std::vector<uint32_t> v : matrix) {
    if (v.size() != m) {
      throw std::invalid_argument(
          "The given matrix has inconsistent row lengths.");
    }
  }

  std::vector<std::vector<uint32_t>> gf2_matrix_transpose;
  gf2_matrix_transpose.reserve(m);

  for (int c = 0; c < m; c++) {
    std::vector<uint32_t> gf2_vector;
    gf2_vector.reserve(N);

    for (int rq = 0; rq <= n - 32; rq += 32) {
      uint32_t t = 0;
      for (int rr = 0; rr < 32; rr++) {
        if (matrix[rq + rr][c] % 2 == 1) {
          t |= (1 << rr);
        }
      }
      gf2_vector.push_back(t);
    }
    if (n % 32) {
      uint32_t t = 0;
      for (size_t rr = 0; rr < n % 32; rr++) {
        if (matrix[(n / 32) * 32 + rr][c] % 2 == 1) {
          t |= (1 << rr);
        }
      }
      gf2_vector.push_back(t);
    }
    gf2_matrix_transpose.push_back(gf2_vector);
  }

  bool marked[n];
  memset(marked, false, sizeof(marked));

  int row_marked[m];

  for (int c = 0; c < m; c++) {
    int pivot = -1;
    for (int r = 0; r < n; r++) {
      if (gf2_matrix_transpose[c][r / 32] & (1 << (r % 32))) {
        pivot = r;
        break;
      }
    }

    row_marked[c] = pivot;

    if (pivot != -1) {
      marked[pivot] = true;

      int pivot_q = pivot / 32;
      int pow2_pivot_r = (1 << (pivot % 32));

      // for (int c2 = 0; c2 < m; c2++) {
      //   if (gf2_matrix[pivot][c2] && c2 != c) {
      //     for (int r = 0; r < n; r++) {
      //       gf2_matrix[r][c2] = gf2_matrix[r][c2] ^ gf2_matrix[r][c];
      //     }
      //   }
      // }
      for (int c2 = 0; c2 < c; c2++) {
        if (gf2_matrix_transpose[c2][pivot_q] & pow2_pivot_r) {
          for (int R = 0; R < N; R++) {
            gf2_matrix_transpose[c2][R] ^= gf2_matrix_transpose[c][R];
          }
        }
      }
      for (int c2 = c + 1; c2 < m; c2++) {
        if (gf2_matrix_transpose[c2][pivot_q] & pow2_pivot_r) {
          for (int R = 0; R < N; R++) {
            gf2_matrix_transpose[c2][R] ^= gf2_matrix_transpose[c][R];
          }
        }
      }
    }
  }

  for (int r = 0; r < n; r++) {
    if (!marked[r]) {
      std::vector<int> needed_rows(1, r);

      int rq = r / 32, pow2_rr = (1 << (r % 32));

      for (int c = 0; c < m; c++) {
        if (gf2_matrix_transpose[c][rq] & pow2_rr) {
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

  // Euler Criterion is not good

  // Bint power = Bint::pow(a_mod_p, (p - 1) >> 1, p);
  // if (power == 1) {
  //   return 1;
  // }
  // return -1;

  Bint x = a_mod_p, y = p;
  int answer = 1;
  // We recursively calculate the Jacobi Symbol (x/y)
  while (x != 0) {
    // std::cout << x << " " << y << " " << answer << std::endl;
    if ((y & 7) == 1 || (y & 7) == 7) { // i.e. if (2/y) = 1
      while (!(x & 1).to_bool()) {
        x >>= 1;
      }
    } else {
      while (!(x & 1).to_bool()) {
        x >>= 1;
        answer = -answer;
      }
    }

    Bint temp = y;
    y = x;
    x = temp;

    if ((x & 3) == 3 && (y & 3) == 3) { // Quadratic Reciprocity
      answer = -answer;
    }
    x %= y;
  }
  return answer;
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

std::map<Bint, uint32_t>
prime_factors(const Bint &n,
              std::function<Bint(const Bint &)> &&find_nontrivial_factor) {

  std::stack<Bint> st;
  std::map<Bint, uint32_t> prime_fact;

  if (!(n & 1).to_bool()) {
    std::pair<uint32_t, Bint> p = two_adic_val(n);
    prime_fact[2] = p.first;
    st.push(p.second);
  } else {
    st.push(n);
  }

  Bint one = 1;

  while (!st.empty()) {
    Bint t = st.top();
    st.pop();

    std::pair<Bint, uint32_t> p = factor_perfect_power(t);
    if (p.second == 1) {
      Bint f = find_nontrivial_factor(t);

      if (f == one) {
        prime_fact[t]++;
      } else {
        st.push(f);
        st.push(t / f);
      }
    } else {
      prime_fact[p.first] += p.second;
    }
  }

  return prime_fact;
}

Bint pollard_rho(const Bint &n) {
  // std::cout << "[ " << n << std::flush;
  if (is_probable_prime(n)) {
    // std::cout << " -> 1 ]" << std::flush;
    return 1;
  }

  auto poly = [&n](Bint x) { return (x * x + 1) % n; };

  std::pair<Bint, uint32_t> p = factor_perfect_power(n);
  if (p.second != 1) {
    // std::cout << " -> " << Bint::pow(p.first, p.second / 2) << "]"
    //           << std::flush;
    return Bint::pow(p.first, p.second / 2);
  }

  Bint g = n;
  Bint n_minus_two = n - 2;
  while (g == n) {
    Bint turtle = Bint::rand(n_minus_two) + 1;
    Bint hare = turtle;
    // std::cout << n_minus_two << " " << turtle << " " << hare << std::endl;
    g = 1;

    while (g == 1) {
      turtle = poly(turtle);
      hare = poly(poly(hare));
      g = Bint::gcd(turtle - hare, n);
    }
  }

  // std::cout << " -> " << Bint::gcd(g, n) << "]" << std::flush;

  return Bint::gcd(g, n);
}

std::tuple<Bint, Bint, Bint> extended_euclidean(const Bint &a, const Bint &b) {
  Bint u1 = 1, u2 = 0, v1 = 0, v2 = 1;
  Bint u3 = a * u1 + b * u2, v3 = a * v1 + b * v2;

  while (v3 != 0) {
    Bint q = u3 / v3;
    Bint t1 = u1 - q * v1;
    Bint t2 = u2 - q * v2;

    u1 = v1;
    u2 = v2;

    v1 = t1;
    v2 = t2;

    u3 = a * u1 + b * u2;
    v3 = a * v1 + b * v2;
  }
  return std::make_tuple(u1, u2, a * u1 + b * u2);
}

Bint modular_inv(const Bint &a, const Bint &m) {
  std::tuple<Bint, Bint, Bint> t = extended_euclidean(a % m, m);

  if (std::get<2>(t) != 1) {
    throw std::domain_error("a = " + a.to_string() +
                            " and m = " + m.to_string() + " are not coprime.");
  }

  return (std::get<0>(t) % m + m) % m;
}

Bint modular_inv_mod_prime(const Bint &a, const Bint &p) {
  if (p.bit_length() <= 31) {
    int p_int = p.to_uint32_t();

    int u = Bint::div_m(a, p).second.to_uint32_t(), v = p_int;
    int x1 = 1, x2 = 0;

    while (u != 1) {
      int q = v / u, r = v % u;
      int x = x2 - q * x1;

      v = u;
      u = r;
      x2 = x1;
      x1 = x;
    }
    return (x1 + p_int) % p_int;
  }

  Bint u = Bint::div_m(a, p).second, v = p;
  Bint x1 = 1, x2 = 0;

  while (u != 1) {
    std::pair<Bint, Bint> d = Bint::div(v, u);
    Bint x = x2 - d.first * x1;

    v = u;
    u = d.second;
    x2 = x1;
    x1 = x;
  }
  return Bint::div_m(x1, p).second;
}

} // namespace Algorithm
