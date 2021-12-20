#include "lib/Algorithm.h"
#include "lib/Bint.h"
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <set>

Bint quadratic_sieve(const Bint &n) {
  if (Algorithm::is_probable_prime(n)) {
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

  V = (uint32_t)(2 * approx_B *
                 (log(2 * approx_B) + log(log(2 * approx_B)) - 1));

  std::vector<uint32_t> possible_primes = Algorithm::primes_less_than(V);

  // size_t sz = 0;
  for (uint32_t p : possible_primes) {
    if (Algorithm::legendre_symbol(n, p) == 1) {
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
        rt = Algorithm::square_root_modulo_prime_power(n, p, e, rt);
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
          Algorithm::gf2_gaussian_elimination(gf2_smooth_matrix);

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

void print_prime_fact(Bint n) {
  if (n == 1) {
    std::cout << "n = 1" << std::endl;
    return;
  }

  std::map<Bint, uint32_t> pf = Algorithm::prime_factors(
      n, [](const Bint &b) { return quadratic_sieve(b); });

  auto power_text = [](std::pair<Bint, uint32_t> p) {
    if (p.second == 1) {
      return p.first.to_string();
    }
    return p.first.to_string() + " ^ " + std::to_string(p.second);
  };

  std::cout << "the prime factorization of n is n = ";

  for (auto it = pf.begin(); it != pf.end(); it++) {
    if (it == pf.begin()) {
      std::cout << power_text(*it);
    } else {
      std::cout << " * " << power_text(*it);
    }
  }
  std::cout << std::endl;
}

int main() {
  while (true) {
    std::cout << "n: ";
    Bint n;
    std::cin >> n;

    clock_t tStart = clock();
    print_prime_fact(n);

    std::cout << "Time taken: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;
  }
}
