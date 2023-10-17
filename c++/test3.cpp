#include <bits/stdc++.h>
using namespace std;

int primes_less_than(const uint32_t &n) {
  if (n <= 2) {
      return 0;
  }
  if (n <= 3) {
    return 1;
  }
  if (n <= 5) {
    return 2;
  }
  if (n <= 7) {
    return 3;
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

  bool* is_prime = new bool[P];
  memset(is_prime, true, P * sizeof(bool));

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

  int nums = 4;
  for (uint32_t x = 0; x < P; x++) {
    if (is_prime[x]) {
      uint32_t p = 210 * (x / 48) + wheel_residues[x % 48];

      if (p >= n) {
        break;
      }

      nums++;
    }
  }

  delete[] is_prime;

  return nums;
}

int main() {
    int p = primes_less_than(1000000000);
    std::cout << p;
}
