#include "lib/Algorithm.h"
#include "lib/Bint.h"
#include <iostream>
#include <map>
#include <time.h>

using namespace std;

void print_prime_fact(Bint n) {
  if (n == 1) {
    cout << "n = 1" << endl;
    return;
  }

  map<Bint, uint32_t> pf = Algorithm::prime_factors(
      n, [](const Bint &b) { return Algorithm::pollard_rho(b); });

  auto power_text = [](pair<Bint, uint32_t> p) {
    if (p.second == 1) {
      return p.first.to_string();
    }
    return p.first.to_string() + " ^ " + to_string(p.second);
  };

  cout << "the prime factorization of n is n = ";

  for (auto it = pf.begin(); it != pf.end(); it++) {
    if (it == pf.begin()) {
      cout << power_text(*it);
    } else {
      cout << " * " << power_text(*it);
    }
  }
  cout << endl;
}

int main() {
  while (true) {
    cout << "n: ";
    Bint n;
    cin >> n;

    clock_t tStart = clock();
    print_prime_fact(n);

    cout << "Time taken: " << fixed << setprecision(3)
         << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s" << endl;
  }
}
