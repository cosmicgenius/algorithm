#include <bits/stdc++.h>
using namespace std;

vector<int> v() {
  vector<int> vec;
  vec.reserve(10000);
  for (int i = 0; i < 10000; i++) {
    vec[i] = ((i >> 3) ^ i) * 129837192;
  }
  return vec;
}

int main() {
  clock_t tStart = clock();
  for (int i = 0; i < 10000; i++) {
    vector<int> vec = v();
  }
  cout << "Time taken: " << fixed << setprecision(3)
       << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s" << endl;
}