#include "lib/Algorithm.h"
#include "lib/Bint.h"
#include <iostream>
#include <map>
#include <time.h>

void print_prime_fact(Bint n) {
  if (n == 1) {
    std::cout << "n = 1" << std::endl;
    return;
  }

  std::map<Bint, uint32_t> pf = Algorithm::prime_factors(
      n, [](const Bint &b) { return Algorithm::pollard_rho(b); });

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
