#include "Bint.h"

#include <algorithm>
#include <functional>
#include <map>
#include <math.h>
#include <set>
#include <stack>
#include <string.h>
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

Bint quadratic_sieve(const Bint &n) {
  if (is_probable_prime(n)) {
    return 1;
  }

  std::cout << "factoring " << n << std::endl;

  double approx_B;
  uint32_t B, M, V;
  std::set<uint32_t> prime_base;
  // std::map<uint32_t, uint32_t> prime_base_back;

  Bint sqrt_n = Bint::integral_rth_root(n, 2);
  int sqrt_n_bit_length = sqrt_n.bit_length();

  double ln_n = Bint::log_2(n) / 1.4426950408889634;
  approx_B = exp(sqrt(ln_n * log(ln_n) / 8));

  V = (uint32_t)(approx_B * (log(approx_B) + log(log(approx_B))) * 2);

  std::vector<uint32_t> possible_primes = primes_less_than(V);

  // size_t sz = 0;
  for (uint32_t p : possible_primes) {
    if (legendre_symbol(n, p) == 1) {
      prime_base.insert(p);
      // prime_base_back[p] = sz++;
    }
  }

  B = prime_base.size();
  M = std::min(std::min((uint32_t)100000000, B * B * B * 2),
               sqrt_n.to_uint32_t()); // conserve memory by sieving max
                                      // 100 million at a time

  std::cout << "approx_B = " << approx_B << std::endl;
  std::cout << "B = " << B << std::endl;
  std::cout << "M = " << M << std::endl;
  std::cout << "V = " << V << std::endl;
  std::cout << "prime_base = [ ";
  for (uint32_t p : prime_base) {
    std::cout << p << " ";
  }
  std::cout << "]" << std::endl;

  uint32_t sieving_idx = 1;

  std::vector<Bint> sieve_results, sieve_result_root;
  uint32_t verified_sieve_results = 0;
  int *vals = new int[M];

  std::map<uint32_t, std::map<uint32_t, Bint>> nsqrt;
  std::map<uint32_t, uint32_t> nlog;

  std::cout << "initialization step" << std::endl;
  for (uint32_t p : prime_base) {
    uint32_t q = p;
    int logp = Bint(p).bit_length();
    nlog[p] = logp;

    Bint rt = 0;
    for (uint32_t e = 1; e <= Bint::log_2(M) / logp; e++) {
      try {
        rt = square_root_modulo_prime_power(n, p, e, rt);
        nsqrt[p][e] = rt;
      } catch (const std::domain_error &err) {
        break;
      }
    }
  }

  std::vector<std::vector<uint32_t>> gf2_smooth_matrix;
  std::vector<uint32_t> gf2_smooth_vector;
  gf2_smooth_vector.resize(B + 1);

  do {
    uint32_t offset = sieving_idx / 2;

    bool done = false;
    Bint a = sieving_idx % 2 == 0 ? (sqrt_n + (Bint)offset * M)
                                  : (sqrt_n - (Bint)offset * M);
    Bint b = n;
    Bint max_val = (sqrt_n * (offset + 1) * M) << 1;
    if (a < 0) {
      continue;
    }
    std::cout << "sieving with new a: " << a << std::endl;

    for (uint32_t i = 0; i < M; i++) {
      vals[i] = -sqrt_n_bit_length;
    }

    for (uint32_t p : prime_base) {
      // std::cout << p << " started" << std::endl;

      uint32_t q = p;
      int logp = Bint(p).bit_length();

      for (const std::pair<uint32_t, Bint> &p_nsqrt : nsqrt[p]) {
        const uint32_t &e = p_nsqrt.first;
        const Bint &rt = p_nsqrt.second;

        uint32_t st1 = 0, st2 = 0;

        if (p == 2 && e == 1) {
          st2 = ((a & 1).to_bool() + (b & 1).to_bool()) % 2;
        } else {
          st1 = (((-a - rt) % q + q) % q).to_uint32_t();
          st2 = (((-a + rt) % q + q) % q).to_uint32_t();

          // std::cout << st1 << " " << st2 << std::endl;

          for (uint32_t x = st1; x < M; x += q) {
            vals[x] += logp;
          }
        }
        for (uint32_t x = st2; x < M; x += q) {
          vals[x] += logp;
        }
        q *= p;
      }
    }

    // std::cout << "checking sieved results" << std::endl;
    // std::cout << "vals array = [ ";
    // for (int i = 0; i < M; i++) {
    //   std::cout << vals[i] << " ";
    // }
    // std::cout << "]" << std::endl;

    Bint t = a * a - b;
    Bint a2 = a << 1;

    size_t prev_results = sieve_results.size();

    int lazy_bit_length =
        ((Bint)offset * M + ((sieving_idx % 2 == 0) ? 0 : M)).bit_length();

    // std::cout << "lazy_bit_length = " << lazy_bit_length << std::endl;
    for (uint32_t x = 0; x < M; x++) {
      // if (vals[x] >= ((a2 + x) * x + t).bit_length()) {
      // Bint sqrt_n_offset = sieving_idx % 2 == 0 ? ((Bint)offset * M + x)
      //                                           : (-(Bint)offset * M + x);
      // if ((vals[x] >= (int)sqrt_n_offset.bit_length())) {
      if ((vals[x] >= lazy_bit_length)) {
        Bint v = (a + x) * (a + x) - b;

        // we want to multiply things to make a square,
        // if it's already a square, technically we could just immediately try
        // to factor it, but let's just ignore this for now
        //
        // edit: this makes no sense, removed for now
        //
        // Bint isqrt = Bint::integral_rth_root(v, 2);
        // if (isqrt * isqrt == v) {
        //   continue;
        // }

        // std::cout << "found " << x << " " << v << std::endl;

        sieve_results.push_back(v);
        sieve_result_root.push_back(a + x);
      }
    }

    uint32_t S = sieve_results.size();

    std::cout << S - prev_results << " results found, " << S << " / "
              << (uint32_t)(B + ln_n) << " results in total" << std::endl;

    while (S > B + ln_n) {
      std::cout << "generating matrix from sieved results" << std::endl;

      std::vector<Bint>::iterator res_it = sieve_results.begin() +
                                           verified_sieve_results,
                                  res_rt_it = sieve_result_root.begin() +
                                              verified_sieve_results;
      std::vector<std::vector<uint32_t>>::iterator matrix_it;

      uint32_t new_res = 0, tot = 0;

      while (res_it != sieve_results.end()) {
        tot++;
        // std::cout << "sieve results: [ ";
        // for (Bint b : sieve_results) {
        //   std::cout << b << " ";
        // }
        // std::cout << "]" << std::endl;

        // we will use the first value as the "exponent" of -1
        gf2_smooth_vector[0] = (*res_it < 0) ? 1 : 0;

        // std::cout << ">" << std::flush;
        // std::map<Bint, uint32_t> pf =
        //     prime_factors((*res_it).abs(), [](const Bint &b) {
        //       return Algorithm::pollard_rho(b);
        //     });
        // bool actually_smooth = true;
        // for (const std::pair<Bint, uint32_t> &prime_power : pf) {
        //   if (prime_base.find(prime_power.first.to_uint32_t()) ==
        //       prime_base.end()) {

        //     // std::cout << *res_it << " fails; divisible by "
        //     //           << prime_power.first.to_uint32_t() << std::endl;

        //     // std::cout << ".";
        //     actually_smooth = false;
        //     break;
        //   }
        // }

        // Trial division is actually faster
        Bint v = (*res_it).abs();
        int idx = 0;
        for (uint32_t p : prime_base) {
          gf2_smooth_vector[++idx] = 0;
          while (v % p == 0) {
            gf2_smooth_vector[idx]++;
            v /= p;
          }
        }

        if (v != 1) {
          res_it = sieve_results.erase(res_it);
          res_rt_it = sieve_result_root.erase(res_rt_it);

          // std::cout << "deleted one" << std::endl;

          continue;
        }

        // std::cout << *res_it << " is verified" << std::endl;
        std::cout << "verified smooth numbers: " << new_res + 1 << " / " << tot
                  << " (" << 1000.0 * (new_res + 1) / tot / 10.0 << "%)" << '\r'
                  << std::flush;

        // int idx = 0;
        // for (uint32_t p : prime_base) {
        //   gf2_smooth_vector[++idx] = pf[p];
        // }

        gf2_smooth_matrix.push_back(gf2_smooth_vector);

        res_it++;
        res_rt_it++;
        // std::cout << "./";
        new_res++;
      }
      S = sieve_results.size();
      verified_sieve_results = S;

      if (S <= B + ln_n) {
        break;
      }

      // std::vector<std::vector<bool>> gf2_smooth_matrix_bool =
      // gf2_smooth_matrix;

      std::cout << std::endl;
      std::cout << "starting gaussian elimination after finding " << S
                << " > B + ln(n) = " << B + (int)ln_n << " smooth numbers"
                << std::endl;

      std::vector<std::vector<int>> all_subsets =
          gf2_gaussian_elimination(gf2_smooth_matrix);

      std::cout << "found " << all_subsets.size()
                << " linear dependencies, extracting g" << std::endl;

      bool marked_for_deletion[S];

      memset(marked_for_deletion, 0, sizeof(marked_for_deletion));

      for (const std::vector<int> &subset : all_subsets) {
        Bint square1_rt = 1;
        Bint square2_rt = 1;

        uint32_t prime_powers[B];
        memset(prime_powers, 0, sizeof(prime_powers));

        // std::cout << "subset = [ ";
        for (int s : subset) {
          // std::cout << s << " ";
          square1_rt *= sieve_result_root[s];
          for (uint32_t i = 0; i < B; i++) {
            prime_powers[i] += gf2_smooth_matrix[s][i + 1];
          }
        }
        // std::cout << "]" << std::endl;
        std::set<uint32_t>::iterator it = prime_base.begin();
        for (uint32_t i = 0; i < B; i++) {
          // std::cout << prime_powers[i] << std::endl;
          square2_rt *= Bint::pow(*it++, prime_powers[i] / 2);
        }
        // std::cout << square1_rt << " " << square2_rt << std::endl;
        Bint g = Bint::gcd((square1_rt - square2_rt).abs(), n);

        std::cout << "g = " << g << std::endl;

        if (g > 1 && g < n) {
          delete[] vals;
          return g;
        }

        marked_for_deletion[subset[0]] = true;
      }

      // It failed
      res_it = sieve_results.begin();
      res_rt_it = sieve_result_root.begin();
      matrix_it = gf2_smooth_matrix.begin();
      for (uint32_t idx = 0; idx < S; idx++) {
        if (marked_for_deletion[idx]) {
          res_it = sieve_results.erase(res_it);
          res_rt_it = sieve_result_root.erase(res_rt_it);
          matrix_it = gf2_smooth_matrix.erase(matrix_it);
        } else {
          res_it++;
          res_rt_it++;
          matrix_it++;
        }
      }

      S = sieve_results.size();
      verified_sieve_results = S;
    }
    // delete[] vals;
    // return 1;
  } while (++sieving_idx);

  return 1;
}

} // namespace Algorithm