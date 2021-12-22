#include "lib/Algorithm.h"
#include "lib/Bint.h"
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <set>

Bint mpqs(const Bint &n) {
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
  M = floor(pow(B, 1.2)) * 60;

  std::cout << "approx_B = " << approx_B << std::endl;
  std::cout << "B = " << B << std::endl;
  std::cout << "M = " << M << std::endl;
  std::cout << "V = " << V << std::endl;
  // std::cout << "prime_base = [ ";
  // for (uint32_t p : prime_base) {
  //   std::cout << p << " ";
  // }
  // std::cout << "]" << std::endl;

  uint32_t poly = 1;

  std::vector<Bint> sieve_results, sieve_result_root;
  uint32_t verified_sieve_results = 0;
  int *vals = new int[2 * M + 1];

  std::map<uint32_t, Bint> nsqrt;
  std::map<uint32_t, uint32_t> nlog;

  std::cout << "initialization step" << std::endl;
  for (uint32_t p : prime_base) {
    uint32_t q = p;
    int logp = Bint(p).bit_length();
    nlog[p] = logp;

    Bint rt = 0;
    try {
      rt = Algorithm::square_root_modulo_prime(n, p);
    } catch (const std::domain_error &err) {
      break;
    }
    nsqrt[p] = rt;
  }

  std::vector<std::vector<uint32_t>> gf2_smooth_matrix;
  std::vector<uint32_t> gf2_smooth_vector;
  gf2_smooth_vector.resize(B + 1);

  Bint q_min =
      Bint::integral_rth_root(Bint::integral_rth_root(n << 1, 2) / M, 2);
  if (q_min < 3) {
    q_min = 3;
  }

  do {
    bool done = false;
    // std::cout << q_min << std::endl;

    Bint q = q_min;

    while (!Algorithm::is_probable_prime(q) ||
           Algorithm::legendre_symbol(n, q) != 1) {
      // std::cout << q << std::endl;
      q++;
    }

    Bint q_inv = Algorithm::modular_inv(q, n);

    q_min = q + 1;
    Bint a = q * q;
    Bint b = Algorithm::square_root_modulo_prime_power(n, q, 2);

    auto find_q =
        std::find(prime_base.begin(), prime_base.end(), q.to_uint32_t());
    if (find_q != prime_base.end()) {
      prime_base.erase(find_q);
    }

    // std::cout << "sieving with new a, b: (" << a << " " << b << ")"
    // << std::endl;

    for (uint32_t i = 0; i <= 2 * M; i++) {
      vals[i] = -sqrt_n_bit_length;
    }

    for (uint32_t p : prime_base) {
      // std::cout << "p = " << p << " started" << std::endl;

      uint32_t a_inv = Bint::pow(a, p - 2, p).to_uint32_t();
      // std::cout << "p, a_inv: " << p << ", " << a_inv << std::endl;
      int logp = Bint(p).bit_length();

      const Bint &rt = nsqrt[p];

      uint32_t st1 = 0, st2 = 0;

      if (p == 2) {
        st2 = ((a & 1).to_bool() + (b & 1).to_bool()) % 2;
      } else {
        st1 = ((((-rt - b) * a_inv + M) % p + p) % p).to_uint32_t();
        st2 = ((((rt - b) * a_inv + M) % p + p) % p).to_uint32_t();

        // std::cout << "rt: " << rt << std::endl;
        // std::cout << "(st1, st2): (" << st1 << ", " << st2 << ")" <<
        // std::endl;

        for (uint32_t x = st1; x <= 2 * M; x += p) {
          vals[x] += logp;
        }
      }
      for (uint32_t x = st2; x <= 2 * M; x += p) {
        vals[x] += logp;
      }
    }

    // std::cout << "checking sieved results" << std::endl;
    // std::cout << "vals array = [ ";
    // for (int i = 0; i <= 2 * M; i++) {
    //   std::cout << vals[i] << " ";
    // }
    // std::cout << "]" << std::endl;

    size_t prev_results = sieve_results.size();
    int lazy_bit_length = Bint(M).bit_length();

    // std::cout << "lazy_bit_length = " << lazy_bit_length << std::endl;
    for (uint32_t x = 0; x < 2 * M + 1; x++) {
      // if (vals[x] >= ((a2 + x) * x + t).bit_length()) {
      // Bint sqrt_n_offset = sieving_idx % 2 == 0 ? ((Bint)offset * M + x)
      //                                           : (-(Bint)offset * M + x);
      // if ((vals[x] >= (int)sqrt_n_offset.bit_length())) {
      if (vals[x] >= lazy_bit_length - 1) {
        Bint u = a * (Bint(x) - M) + b;
        Bint v = (u * u - n) / a;
        // std::cout << v << " " << x << std::endl;

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
        sieve_result_root.push_back(u * q_inv % n);
      }
    }

    uint32_t S = sieve_results.size();

    std::cout << "Sieved " << poly << " polynomials. " << S - prev_results
              << " results found in the last polynomial. " << S << " / "
              << (uint32_t)(B + ln_n) << " results in total"
              << " (" << 1000.0 * S / ((uint32_t)(B + ln_n)) / 10.0 << "%). "
              << "Around " << round(1000.0 * S / poly) / 1000.0
              << " per polynomial    " << '\r' << std::flush;

    while (S > B + ln_n) {
      std::cout << std::endl;
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
    // return 71;
  } while (++poly);

  return 1;
}

void print_prime_fact(Bint n) {
  if (n == 1) {
    std::cout << "n = 1" << std::endl;
    return;
  }

  std::map<Bint, uint32_t> pf =
      Algorithm::prime_factors(n, [](const Bint &b) { return mpqs(b); });

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
