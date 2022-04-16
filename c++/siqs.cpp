#include "lib/Algorithm.h"
#include "lib/Bint.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <set>

Bint siqs(const Bint &n) {
  if (Algorithm::is_probable_prime(n)) {
    return 1;
  }

  std::cout << "factoring " << n << std::endl;

  double approx_B;
  uint32_t B, M, V, largest_partial_prime;
  std::vector<uint32_t> prime_base;

  // std::map<uint32_t, uint32_t> prime_base_back;

  Bint sqrt_n = Bint::integral_rth_root(n, 2);
  int sqrt_n_bit_length = sqrt_n.bit_length();

  double ln_n = Bint::log_2(n) / 1.4426950408889634;
  approx_B = exp(sqrt(ln_n * log(ln_n) / 9.6));

  V = (uint32_t)(2 * approx_B *
                 (log(2 * approx_B) + log(log(2 * approx_B)) - 1));

  std::vector<uint32_t> possible_primes = Algorithm::primes_less_than(V);
  largest_partial_prime = V * 100;

  // size_t sz = 0;
  for (const uint32_t &p : possible_primes) {
    if (Algorithm::legendre_symbol(n, p) == 1) {
      prime_base.push_back(p);
      // prime_base_back[p] = sz++;
    }
  }

  B = prime_base.size();
  M = floor(pow(B, 1.5) * 0.7);

  std::cout << "approx_B = " << approx_B << std::endl;
  std::cout << "B = " << B << std::endl;
  std::cout << "M = " << M << std::endl;
  std::cout << "V = " << V << std::endl;
  // std::cout << "prime_base = [ ";
  // for (uint32_t p : prime_base) {
  //   std::cout << p << " ";
  // }
  // std::cout << "]" << std::endl;

  uint32_t poly_group = 1;
  uint32_t total_poly = 1;

  std::vector<Bint> /*sieve_results,*/ sieve_result_root;
  int *vals = new int[2 * M + 1];

  std::map<uint32_t, Bint> nsqrt;
  std::map<uint32_t, int16_t> nlog;

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

  std::map<uint32_t, std::vector<std::pair<Bint, Bint>>> partial_cong;

  std::vector<std::vector<uint32_t>> gf2_smooth_matrix;
  std::vector<uint32_t> gf2_smooth_vector;
  gf2_smooth_vector.resize(B + 1);

  int lazy_bit_length = sqrt_n_bit_length + Bint(M).bit_length() -
                        Bint(largest_partial_prime).bit_length();

  // uint32_t shuffled_prime_base[B];
  // int prime_base_idx = 0;
  // for (const uint32_t &p : prime_base) {
  //   shuffled_prime_base[prime_base_idx++] = p;
  // }
  // std::shuffle(shuffled_prime_base, shuffled_prime_base + B,
  //              std::default_random_engine(4));
  // int critical_min_idx = 0;

  size_t target_a_approx_len = (n.bit_length() + 1) / 2 - Bint(M).bit_length();
  Bint target_a = Bint::integral_rth_root(n << 1, 2) / M;

  // clock_t q_find_time = 0, init_time = 0, calc_time = 0, assign_time = 0;
  // clock_t last = clock();

  uint32_t S = round(
      (double)target_a_approx_len /
      11); // we pick critical primes of size ~2000, and 2^11 ~ 2000, hence 11

  uint32_t T = B - S; // number of non-critical primes

  // critical primes should be around target_a^(1/S)
  uint32_t critical_target = Bint::integral_rth_root(target_a, S).to_uint32_t();
  int critical_lower = 0, critical_upper = 0;
  for (int i = 0; i < B; i++) {
    if (prime_base[i] > critical_target / 2) {
      critical_lower = i;
      break;
    }
  }

  for (int i = B - 1; i >= 0; i--) {
    if (prime_base[i] < critical_target * 2) {
      critical_upper = i + 1;
      break;
    }
  }

  size_t prev_results = 0;

  do {
    // last = clock();
    bool done = false;

    Bint a = 1;
    size_t a_len = a.bit_length();
    Bint b = 0;

    std::set<uint32_t> critical_primes;
    std::vector<uint32_t> non_critical_primes;

    // for (int i = critical_min_idx; i < B; i++) {
    //   if (shuffled_prime_base[i] > 2000) {
    //     critical_primes.insert(shuffled_prime_base[i]);
    //     a *= shuffled_prime_base[i];
    //     if (abs(a.bit_length() + Bint(shuffled_prime_base[i +
    //     1]).bit_length() -
    //             target_a_approx_len) >
    //         abs(a.bit_length() - target_a_approx_len)) {
    //       // heuristic, if adding the next prime gets us further away, stop
    //       critical_min_idx += 2;
    //       break;
    //     }
    //   }
    // }
    // for (int i = critical_min; i < B; i++) {
    //   if (prime_base[i] > 1000) {
    //     critical_primes.insert(prime_base[i]);
    //     a *= prime_base[i];
    //     if (a.bit_length() >
    //         target_a_approx_len - (Bint(prime_base[i]).bit_length()) /
    //         3) {
    //       critical_min += 2;
    //       break;
    //     }
    //   }
    // }

    // We find critical primes through a method described here:
    // https://www.mersenneforum.org/showthread.php?p=535652
    //
    // Essentially we randomly pick S - 1 primes in the prime base from 1000 to
    // 4000, then hand pick the last prime to be as close as possible
    //
    // Technically, this method ever so slightly skews the last prime to be
    // larger due to 1. the distribution of primes, and 2. the distribution of
    // rand, but it probably doesn't matter.

    for (int i = 0; i < S - 1; i++) {
      int picked_idx = 0;
      do {
        picked_idx =
            rand() % (critical_upper - critical_lower) + critical_lower;
      } while (critical_primes.find(prime_base[picked_idx]) !=
               critical_primes.end());

      critical_primes.insert(prime_base[picked_idx]);
      a *= prime_base[picked_idx];
    }

    uint32_t last_prime_predict = (target_a / a).to_uint32_t();
    std::vector<uint32_t>::iterator last_prime_pick = std::lower_bound(
        prime_base.begin(), prime_base.end(), last_prime_predict);

    for (auto it = last_prime_pick; it != prime_base.end(); it++) {
      if (critical_primes.find(*it) == critical_primes.end()) {
        critical_primes.insert(*it);
        a *= *it;
        break;
      }
    }

    for (const uint32_t &p : prime_base) {
      if (critical_primes.find(p) == critical_primes.end()) {
        non_critical_primes.push_back(p);
      }
    }

    Bint crt_partial[S];
    uint32_t a_inv[T];
    std::vector<uint32_t> B_a_inv2[T];
    std::pair<uint32_t, uint32_t> soln[T];

    uint32_t prime_idx = 0;
    for (const uint32_t &p : critical_primes) {
      Bint product_others = a / p;

      uint32_t gamma =
          Bint::div_m(
              nsqrt[p] * Algorithm::modular_inv_mod_prime(product_others, p), p)
              .second;
      if (gamma > p / 2) {
        gamma = p - gamma;
      }

      Bint B = product_others * gamma;
      crt_partial[prime_idx++] = B;
      b += B;
    }
    for (int i = 0; i < T; i++) {
      const uint32_t &p = non_critical_primes[i];
      a_inv[i] = Algorithm::modular_inv_mod_prime(a, p).to_uint32_t();

      B_a_inv2[i].resize(S);
      for (int j = 0; j < S; j++) {
        B_a_inv2[i][j] =
            Bint::div_m(Bint(a_inv[i]) * crt_partial[j] << 1, p).second;
      }
    }
    uint32_t num_polys = (1 << (S - 1));

    for (uint32_t poly = 0; poly < num_polys; poly++, total_poly++) {
      for (uint32_t i = 0; i <= 2 * M; i++) {
        vals[i] = -lazy_bit_length;
      }

      if (poly == 0) {
        for (int i = 0; i < T; i++) {
          const uint32_t &p = non_critical_primes[i];
          const Bint &rt = nsqrt[p];

          soln[i] = {Bint::div_m((-rt - b) * a_inv[i] + M, p).second,
                     Bint::div_m((rt - b) * a_inv[i] + M, p).second};
        }
      } else {
        int nu = 0;
#ifdef __GNUC__
        nu = __builtin_ctz(poly);
#else
        uint32_t poly2 = poly;
        poly2 = poly2 | (poly2 << 1);
        poly2 = poly2 | (poly2 << 2);
        poly2 = poly2 | (poly2 << 4);
        poly2 = poly2 | (poly2 << 8);
        poly2 = poly2 | (poly2 << 16);
        nu = std::bitset<32>(~poly2).count();
#endif
        int gray_code_indicator = (poly >> (nu + 1));
        if (gray_code_indicator % 2 == 0) {
          b -= (crt_partial[nu] << 1);
        } else {
          b += (crt_partial[nu] << 1);
        }
        // std::cout << "nu: " << nu << std::endl;
        // std::cout << "b: " << b << std::endl;
        // std::cout << "B_" << nu << "=" << crt_partial[nu] << std::endl;

        for (int i = 0; i < T; i++) {
          const uint32_t &p = non_critical_primes[i];
          if (gray_code_indicator % 2 == 0) {
            soln[i].first = (soln[i].first + B_a_inv2[i][nu]) % p;
            soln[i].second = (soln[i].second + B_a_inv2[i][nu]) % p;
          } else {
            soln[i].first = (soln[i].first + p - B_a_inv2[i][nu]) % p;
            soln[i].second = (soln[i].second + p - B_a_inv2[i][nu]) % p;
          }
        }
      }

      for (int i = 0; i < T; i++) {
        const uint32_t &p = non_critical_primes[i];
        const int &logp = nlog[p];

        if (p == 2) {
          uint32_t soln = ((a & 1).to_bool() + (b & 1).to_bool()) % 2;

          for (uint32_t x = soln; x <= 2 * M; x += p) {
            vals[x] += logp;
          }
        } else {
          const std::pair<uint32_t, uint32_t> &solns = soln[i];

          for (uint32_t x = solns.first; x <= 2 * M; x += p) {
            vals[x] += logp;
          }
          for (uint32_t x = solns.second; x <= 2 * M; x += p) {
            vals[x] += logp;
          }
        }
      }

      for (uint32_t x = 0; x < 2 * M + 1; x++) {
        if (vals[x] >= 0) {
          Bint u = a * (Bint(x) - M) + b;
          Bint v = u * u - n;

          Bint w = v.abs();
          int idx = 0;

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
            sieve_result_root.push_back(u);

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
            // std::cout << u << " " << v << std::endl;
          } else if (w.bit_length() <= 32 &&
                     w.to_uint32_t() < largest_partial_prime) {
            // might not have to be a prime, but let's keep this for now
            partial_cong[w.to_uint32_t()].emplace_back(u, v);
          }
        }
      }
    }

    uint32_t real_cong = sieve_result_root.size();
    uint32_t partials = 0;
    for (const auto &partial_list : partial_cong) {
      partials += (uint32_t)(partial_list.second.size()) / 2;
    }

    std::cout << "Polynomial group " << poly_group << " with " << S
              << " primes: " << real_cong + partials - prev_results
              << " congruences found. " << real_cong + partials << " / "
              << (uint32_t)(B + ln_n) << " congruences in total"
              << " ("
              << 1000.0 * (real_cong + partials) / ((uint32_t)(B + ln_n)) / 10.0
              << "%). Around "
              << round(1000.0 * (real_cong + partials) / poly_group) / 1000.0
              << " per polynomial group.    " << '\r' << std::flush;
    prev_results = real_cong + partials;

    while (real_cong + partials > B + ln_n) {
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

            int idx = 0;
            for (const uint32_t &p : prime_base) {
              gf2_smooth_vector[++idx] = 0;
              while (v1 % p == 0) {
                gf2_smooth_vector[idx]++;
                v1 /= p;
              }
              while (v2 % p == 0) {
                gf2_smooth_vector[idx]++;
                v2 /= p;
              }
            }

            const Bint &k = partial_list.first;

            gf2_smooth_vector[0] = (v1 * v2 < 0) ? 1 : 0;

            try {
              Bint k_inv = Algorithm::modular_inv(k, n);
              Bint u = (u1 * u2 % n) * k_inv % n;

              sieve_result_root.push_back(u);
              gf2_smooth_matrix.push_back(gf2_smooth_vector);
            } catch (const std::domain_error &err) {
              delete[] vals;
              return Bint::gcd(k, n);
            }
          }
        }
      }

      partials = 0;
      real_cong = sieve_result_root.size();

      std::cout << "Starting gaussian elimination after finding " << real_cong
                << " > B + ln(n) = " << B + (int)ln_n << " smooth numbers"
                << std::endl;

      std::vector<std::vector<int>> all_subsets =
          Algorithm::gf2_gaussian_elimination(gf2_smooth_matrix);

      std::cout << "Found " << all_subsets.size()
                << " linear dependencies, extracting g" << std::endl;

      bool marked_for_deletion[real_cong];

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
        std::vector<uint32_t>::iterator it = prime_base.begin();
        for (uint32_t i = 0; i < B; i++) {
          // std::cout << prime_base[i] << " " << prime_powers[i] << std::endl;
          square2_rt *= Bint::pow(prime_base[i], prime_powers[i] / 2);
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
      for (uint32_t idx = 0; idx < real_cong; idx++) {
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

      real_cong = sieve_result_root.size();
    }
  } while (++poly_group);
  return 1;
}

void print_prime_fact(Bint n) {
  if (n == 1) {
    std::cout << "n = 1" << std::endl;
    return;
  }

  std::map<Bint, uint32_t> pf =
      Algorithm::prime_factors(n, [](const Bint &b) { return siqs(b); });

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
