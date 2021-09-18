#include "lib/Algorithm.h"
#include "lib/Bint.h"
#include <iostream>
#include <map>
#include <time.h>

using namespace std;

vector<Bint> prime_fact;

void populate_prime_fact(Bint n) {
  // cout << "Factoring: " << n << endl;
  Bint f1 = Algorithm::pollard_rho(n);

  if (f1 == 1) {
    prime_fact.push_back(n);
    return;
  }

  Bint f2 = n / f1;

  populate_prime_fact(f1);
  populate_prime_fact(f2);
}

void print_prime_fact(Bint n) {
  if (n == 1) {
    cout << "n = 1" << endl;
    return;
  }

  populate_prime_fact(n);
  sort(prime_fact.begin(), prime_fact.end());

  vector<pair<Bint, uint32_t>> prime_power(1, {prime_fact[0], 1});

  for (int i = 1; i < prime_fact.size(); i++) {
    if (prime_fact[i] == prime_power.back().first) {
      prime_power.back().second++;
    } else {
      prime_power.push_back({prime_fact[i], 1});
    }
  }
  prime_fact.clear();

  cout << "the prime factorization of n is n = ";

  auto power_text = [](pair<Bint, uint32_t> p) {
    if (p.second == 1) {
      return p.first.to_string();
    }
    return p.first.to_string() + " ^ " + to_string(p.second);
  };

  cout << power_text(prime_power[0]);

  for (int i = 1; i < prime_power.size(); i++) {
    cout << " * " << power_text(prime_power[i]);
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
