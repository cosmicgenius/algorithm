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
  uint32_t B, M, V, largest_partial_prime;
  std::set<uint32_t> prime_base;

  // std::map<uint32_t, uint32_t> prime_base_back;

  Bint sqrt_n = Bint::integral_rth_root(n, 2);
  int sqrt_n_bit_length = sqrt_n.bit_length();

  double ln_n = Bint::log_2(n) / 1.4426950408889634;
  approx_B = exp(sqrt(ln_n * log(ln_n) / 9) + 0.1);

  V = (uint32_t)(2 * approx_B *
                 (log(2 * approx_B) + log(log(2 * approx_B)) - 1));

  std::vector<uint32_t> possible_primes = Algorithm::primes_less_than(V);
  largest_partial_prime = V * 100;

  // size_t sz = 0;
  for (const uint32_t &p : possible_primes) {
    if (Algorithm::legendre_symbol(n, p) == 1) {
      prime_base.insert(p);
      // prime_base_back[p] = sz++;
    }
  }

  B = prime_base.size();
  M = floor(pow(B, 1.5) * 0.8);

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

  std::vector<Bint> /*sieve_results,*/ sieve_result_root;
  int *vals = new int[2 * M + 1];

  std::map<uint32_t, Bint> nsqrt;
  std::map<uint32_t, uint32_t> nlog;

  std::cout << "initialization step" << std::endl;
  for (const uint32_t &p : prime_base) {
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

  std::map<uint32_t, std::vector<std::pair<Bint, Bint>>> partial_cong;

  std::vector<std::vector<uint32_t>> gf2_smooth_matrix;
  std::vector<uint32_t> gf2_smooth_vector;
  gf2_smooth_vector.resize(B + 1);

  int lazy_bit_length = sqrt_n_bit_length + Bint(M).bit_length() -
                        Bint(largest_partial_prime).bit_length();

  Bint q_min =
      Bint::integral_rth_root(Bint::integral_rth_root(n << 1, 2) / M, 2);
  if (q_min < 3) {
    q_min = 3;
  }

  // clock_t q_find_time = 0, init_time = 0, calc_time = 0, assign_time = 0;
  // clock_t last = clock();

  size_t prev_results = 0;
  int c = 0;

  do {
    // last = clock();
    bool done = false;
    // std::cout << q_min << std::endl;

    Bint q = q_min;

    while (!Algorithm::is_probable_prime(q) ||
           Algorithm::legendre_symbol(n, q) != 1) {
      // std::cout << q << std::endl;
      q++;
    }
    // q_find_time += clock() - last;
    // last = clock();

    Bint q_inv = Algorithm::modular_inv(q, n);

    q_min = q + 1;
    Bint a = q * q;
    Bint b = Algorithm::square_root_modulo_prime_power(n, q, 2);

    // Usually q >> any prime in the prime base, but just to be safe, remove it
    prime_base.erase(q.to_uint32_t());

    // init_time += clock() - last;
    // last = clock();

    // std::cout << "sieving with new a, b: (" << a << " " << b << ")"
    //           << std::endl;

    for (uint32_t i = 0; i <= 2 * M; i++) {
      vals[i] = -lazy_bit_length;
    }

    // assign_time += clock() - last;
    // last = clock();

    for (uint32_t p : prime_base) {
      // std::cout << "p = " << p << " started" << std::endl;

      // When p is small (and also prime), this is faster than the extended
      // euclidean algorithm
      // uint32_t a_inv = Bint::pow(a, p - 2, p).to_uint32_t();
      uint32_t a_inv = Algorithm::modular_inv_mod_prime(a, p).to_uint32_t();
      // std::cout << "p, a_inv: " << p << ", " << a_inv << std::endl;
      int logp = Bint(p).bit_length();

      const Bint &rt = nsqrt[p];

      uint32_t st1 = 0, st2 = 0;

      // init_time += clock() - last;
      // last = clock();

      if (p == 2) {
        st2 = ((a & 1).to_bool() + (b & 1).to_bool()) % 2;
      } else {
        st1 = Bint::div_m((-rt - b) * a_inv + M, p).second;
        st2 = Bint::div_m((rt - b) * a_inv + M, p).second;
        // std::cout << "rt: " << rt << std::endl;
        // std::cout << "(st1, st2): (" << st1 << ", " << st2 << ")" <<
        // std::endl;
        // calc_time += clock() - last;
        // last = clock();

        for (uint32_t x = st1; x <= 2 * M; x += p) {
          vals[x] += logp;
        }
      }
      for (uint32_t x = st2; x <= 2 * M; x += p) {
        vals[x] += logp;
      }

      // assign_time += clock() - last;
      // last = clock();
    }

    // std::cout << "checking sieved results" << std::endl;
    // std::cout << "vals array = [ ";
    // for (int i = 0; i <= 2 * M; i++) {
    //   std::cout << vals[i] << " ";
    // }
    // std::cout << "]" << std::endl;

    // std::cout << "lazy_bit_length = " << lazy_bit_length << std::endl;
    for (uint32_t x = 0; x < 2 * M + 1; x++) {
      // if (vals[x] >= ((a2 + x) * x + t).bit_length()) {
      // Bint sqrt_n_offset = sieving_idx % 2 == 0 ? ((Bint)offset * M + x)
      //                                           : (-(Bint)offset * M + x);
      // if ((vals[x] >= (int)sqrt_n_offset.bit_length())) {
      if (vals[x] >= 0) {
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

        // std::vector<Bint>::iterator res_it = sieve_results.begin() +
        //                                      verified_sieve_results,
        //                             res_rt_it = sieve_result_root.begin() +
        //                                         verified_sieve_results;
        // std::vector<std::vector<uint32_t>>::iterator matrix_it;
        // we will use the first value as the "exponent" of -1

        Bint w = v.abs();
        int idx = 0;
        // for (uint32_t p : prime_base) {
        //   if (w % p == 0) {
        //     c++;
        //   }
        // }

        for (const uint32_t &p : prime_base) {
          std::pair<Bint, uint32_t> res = Bint::div_m(w, p);
          while (res.second == 0) {
            w = res.first;
            if (w == 1) {
              goto done;
            }
            res = Bint::div_m(w, p);
          }
        }
      done:

        if (w == 1) {
          // sieve_results.push_back(v);
          sieve_result_root.push_back(u * q_inv % n);

          gf2_smooth_vector[0] = (v < 0) ? 1 : 0;
          w = v.abs();
          for (const uint32_t &p : prime_base) {
            gf2_smooth_vector[++idx] = 0;
            while (w % p == 0) {
              gf2_smooth_vector[idx]++;
              w /= p;
            }
          }
          gf2_smooth_matrix.push_back(gf2_smooth_vector);
        } else if (w.bit_length() <= 32 &&
                   w.to_uint32_t() < largest_partial_prime) {
          // might not have to be a prime, but let's keep this for now
          partial_cong[w.to_uint32_t()].emplace_back(u * q_inv % n, v);
        }
      }
    }

    uint32_t partials = 0;
    for (const auto &partial_list : partial_cong) {
      partials += (uint32_t)(partial_list.second.size()) / 2;
    }
    uint32_t S = sieve_result_root.size();

    std::cout << "Polynomial " << poly << ": " << S + partials - prev_results
              << " congruences found. " << S + partials << " / "
              << (uint32_t)(B + ln_n) << " congruences in total"
              << " (" << 1000.0 * (S + partials) / ((uint32_t)(B + ln_n)) / 10.0
              << "%). "
              << "Around " << round(1000.0 * (S + partials) / poly) / 1000.0
              << " per polynomial    " << '\r' << std::flush;
    prev_results = S + partials;

    while (S + partials > B + ln_n) {
      std::cout << std::endl;
      std::cout << "Generating matrix from sieved results" << std::endl;
      // std::cout << "Sieve time estimates: q_find=" << q_find_time
      //           << " init=" << init_time << " calc=" << calc_time
      //           << " assign=" << assign_time << ", efficiency="
      //           << (double)(q_find_time + init_time + calc_time) /
      //           assign_time
      //           << std::endl;

      std::cout << "Checking partial congruences..." << std::endl;
      // Turn partials into real congruences
      for (const auto &partial_list : partial_cong) {
        if (partial_list.second.size() > 1) {
          for (size_t i = 0; i < partial_list.second.size() - 1; i += 2) {
            Bint u1 = partial_list.second[i].first,
                 u2 = partial_list.second[i + 1].first;
            Bint v1 = partial_list.second[i].second,
                 v2 = partial_list.second[i + 1].second;

            gf2_smooth_vector[0] = (v1 < 0) ? 1 : 0;

            int idx = 0;
            Bint w1 = v1, w2 = v2;
            for (uint32_t p : prime_base) {
              gf2_smooth_vector[++idx] = 0;
              while (w1 % p == 0) {
                gf2_smooth_vector[idx]++;
                w1 /= p;
              }
              while (w2 % p == 0) {
                gf2_smooth_vector[idx]++;
                w2 /= p;
              }
            }

            const Bint &k = partial_list.first;
            Bint k_inv = Algorithm::modular_inv(k, n);
            Bint u = (u1 * u2 % n) * k_inv % n;

            sieve_result_root.push_back(u);
            gf2_smooth_matrix.push_back(gf2_smooth_vector);
          }
        }
      }

      partials = 0;
      S = sieve_result_root.size();

      if (S <= B + ln_n) {
        break;
      }

      std::cout << "Starting gaussian elimination after finding " << S
                << " > B + ln(n) = " << B + (int)ln_n << " smooth numbers"
                << std::endl;

      std::vector<std::vector<int>> all_subsets =
          Algorithm::gf2_gaussian_elimination(gf2_smooth_matrix);

      std::cout << "Found " << all_subsets.size()
                << " linear dependencies, extracting g" << std::endl;

      bool marked_for_deletion[S];

      memset(marked_for_deletion, 0, sizeof(marked_for_deletion));

      for (const std::vector<int> &subset : all_subsets) {
        Bint square1_rt = 1;
        Bint square2_rt = 1;

        uint32_t prime_powers[B];
        memset(prime_powers, 0, sizeof(prime_powers));

        // std::cout << "subset = [ ";
        for (const int &s : subset) {
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
      std::vector<Bint>::iterator /*res_it = sieve_results.begin(),*/
          res_rt_it = sieve_result_root.begin();
      std::vector<std::vector<uint32_t>>::iterator matrix_it =
          gf2_smooth_matrix.begin();
      for (uint32_t idx = 0; idx < S; idx++) {
        if (marked_for_deletion[idx]) {
          // res_it = sieve_results.erase(res_it);
          res_rt_it = sieve_result_root.erase(res_rt_it);
          matrix_it = gf2_smooth_matrix.erase(matrix_it);
        } else {
          // res_it++;
          res_rt_it++;
          matrix_it++;
        }
      }

      S = sieve_result_root.size();
    }
    // delete[] vals;
    // return 71;
  } while (++poly);
  std::cout << c;
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

  std::cout << "The prime factorization of n is n = ";

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
