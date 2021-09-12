#include "Bint.h"

#include <bitset>
#include <chrono>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

// Blocks when parsing decimal
const size_t DBLOCK_LEN = 9;
const BLOCK DBLOCK_SIZE = 1e9;

const size_t KARATSUBA_CUTOFF = 32;

void Bint::assign_decimal(const std::string &str) {
  int n = str.length();
  data.clear();

  if (n == 0) {
    sgn = 0;
    return;
  }

  sgn = 1;
  int8_t new_sgn = 1;
  int st_idx = 0;
  if (str[0] == '-') {
    new_sgn = -1;
    st_idx = 1;
    n--;
  }

  for (int i = st_idx; i < n; i++) {
    if (str[i] < '0' || str[i] > '9') {
      throw std::invalid_argument("Invalid character in decimal '" +
                                  std::string(1, str[i]) + "'.");
    }
  }

  BLOCK block = 0;
  size_t fblock_len = (n - 1) % DBLOCK_LEN + 1;
  for (int i = 0; i < fblock_len; i++) {
    block *= 10;
    block += str[st_idx + i] - '0';
  }
  data.push_back(block);

  for (int b = 0; b * DBLOCK_LEN + fblock_len < n; b++) {
    block = 0;
    for (int i = 0; i < DBLOCK_LEN; i++) {
      block *= 10;
      block += str[st_idx + b * DBLOCK_LEN + fblock_len + i] - '0';
    }
    *this *= DBLOCK_SIZE;
    *this += block;
  }

  sgn = new_sgn; // set this at the end, otherwise we will have problems when
                 // trying to create the number by adding
}

void Bint::resize() {
  size_t most_sigfig = data.size();
  while (most_sigfig > 0 && data[most_sigfig - 1] == 0) {
    most_sigfig--;
  }

  if (most_sigfig == 0) {
    data.clear();
    sgn = 0;
  } else {
    data.resize(most_sigfig);
  }
}

uint8_t Bint::clz(const BLOCK &n) {
  uint8_t leading_zeroes = 0;
#ifdef __GNUC__
  leading_zeroes = __builtin_clz(n);
#else
  BLOCK leading_digit = n;
  leading_digit = leading_digit | (leading_digit >> 1);
  leading_digit = leading_digit | (leading_digit >> 2);
  leading_digit = leading_digit | (leading_digit >> 4);
  leading_digit = leading_digit | (leading_digit >> 8);
  leading_digit = leading_digit | (leading_digit >> 16);
  leading_zeroes = std::bitset<32>(~leading_digit).count();
#endif
  return leading_zeroes;
}

// Constructors

Bint::Bint() : sgn(0) {}

Bint::Bint(int32_t n) : sgn((n > 0) - (n < 0)) {
  if (n != 0) {
    data.push_back(std::abs(n));
  }
}
Bint::Bint(int64_t n) : sgn((n > 0) - (n < 0)) {
  data.push_back(llabs(n));
  data.push_back(llabs(n) >> 32);
  resize();
}
Bint::Bint(uint32_t n) : sgn((n > 0) - (n < 0)) {
  if (n != 0) {
    data.push_back(n);
  }
}
Bint::Bint(uint64_t n) : sgn((n > 0) - (n < 0)) {
  data.push_back(n);
  data.push_back(n >> 32);
  resize();
}

Bint::Bint(const char *str) { assign_decimal(str); }
Bint::Bint(const std::string &str) { assign_decimal(str); }

Bint::Bint(const Bint &b) : sgn(b.sgn) { data = b.data; }

// Arithmetic

Bint Bint::operator+(const Bint &rhs) const {
  Bint res = *this;
  res += rhs;
  return res;
}

Bint Bint::operator-(const Bint &rhs) const {
  Bint res = *this;
  res -= rhs;
  return res;
}

Bint Bint::operator*(const Bint &rhs) const { return mul(*this, rhs); }
Bint Bint::operator/(const Bint &rhs) const { return div(*this, rhs).first; }
Bint Bint::operator%(const Bint &rhs) const { return div(*this, rhs).second; }

Bint Bint::operator*(const BLOCK &rhs) const {
  Bint res = *this;
  res *= rhs;
  return res;
}

Bint Bint::operator/(const BLOCK &rhs) const { return div(*this, rhs).first; }
Bint Bint::operator%(const BLOCK &rhs) const { return div(*this, rhs).second; }

// Assignment with Arithmetic

const Bint &Bint::operator+=(const Bint &rhs) {
  // std::cout << "Adding " << *this << " and " << rhs << std::endl;

  if (rhs.sgn == 0) {
    return *this;
  }
  if (this->sgn == 0) {
    this->sgn = rhs.sgn;
    this->data = rhs.data;
    return *this;
  }
  size_t ls = this->data.size(), rs = rhs.data.size();
  size_t ms = std::max(ls, rs);

  this->data.resize(ms);

  // Normal addition
  if (this->sgn == rhs.sgn) {
    LBLOCK carry = 0;
    for (int i = 0; i < rs; i++) {
      carry += data[i];
      carry += rhs.data[i];
      data[i] = carry; // this will keep only the last 32 bits
      carry >>= 32;
    }
    for (int i = rs; i < ms; i++) {
      if (carry == 0)
        break;
      carry += data[i];
      data[i] = carry; // this will keep only the last 32 bits
      carry >>= 32;
    }
    if (carry > 0) {
      data.push_back(carry);
    }
  }
  // 2^32's complement subtraction
  else {
    if (this->data == rhs.data) {
      data.clear();
      sgn = 0;
    } else {
      // what will the sign of the result be?
      int8_t new_sgn = (-*this > rhs) ? rhs.sgn : this->sgn;
      Bint rhs2 = rhs;
      rhs2.data.resize(ms);

      if (sgn < 0) {
        *this = ~*this;
        LBLOCK carry = 1;
        for (int i = 0; i < data.size(); i++) {
          if (carry == 0)
            break;
          carry += data[i];
          data[i] = carry;
          carry >>= 32;
        }
      } else if (rhs2.sgn < 0) {
        rhs2 = ~rhs2;
        LBLOCK carry = 1;
        for (int i = 0; i < rhs2.data.size(); i++) {
          if (carry == 0)
            break;
          carry += rhs2.data[i];
          rhs2.data[i] = carry;
          carry >>= 32;
        }
      }

      LBLOCK carry = 0;
      for (int i = 0; i < ms; i++) {
        carry += data[i];
        carry += rhs2.data[i];
        data[i] = carry; // this will keep only the last 32 bits
        carry >>= 32;
      }
      // don't add the carry at the very end

      // If it's negative, convert back from 2^32's complement
      if (new_sgn < 0) {
        *this = ~*this;
        LBLOCK carry = 1;
        for (int i = 0; i < data.size(); i++) {
          if (carry == 0)
            break;
          carry += data[i];
          data[i] = carry;
          carry >>= 32;
        }
      }
      sgn = new_sgn;
      resize();
    }
  }
  return *this;
}

const Bint &Bint::operator-=(const Bint &rhs) {
  *this += -rhs;
  return *this;
}

const Bint &Bint::operator*=(const BLOCK &rhs) {
  LBLOCK carry = 0;
  LBLOCK rhs2 = rhs;

  for (int i = 0; i < data.size(); i++) {
    carry += data[i] * rhs2;
    data[i] = carry; // this will keep only the last 32 bits
    carry >>= 32;
  }
  if (carry > 0) {
    data.push_back(carry);
  }
  return *this;
}

const Bint &Bint::operator/=(const Bint &rhs) {
  *this = *this / rhs;
  return *this;
}

const Bint &Bint::operator%=(const Bint &rhs) {
  *this = *this % rhs;
  return *this;
}

const Bint &Bint::operator*=(const Bint &rhs) {
  *this = *this * rhs;
  return *this;
}

const Bint &Bint::operator/=(const BLOCK &rhs) {
  *this = *this / rhs;
  return *this;
}

const Bint &Bint::operator%=(const BLOCK &rhs) {
  *this = *this % rhs;
  return *this;
}

// const Bint &Bint::operator/=(const BLOCK &rhs) {
//   LBLOCK carry = 0;
//   LBLOCK rhs2 = rhs;

//   for (int i = 0; i < data.size(); i++) {
//     carry += data[i] * rhs2;
//     data[i] = carry; // this will keep only the last 32 bits
//     carry >>= 32;
//   }
//   if (carry > 0) {
//     data.push_back(carry);
//   }
//   return *this;
// }

// Bitwise operators
Bint Bint::operator&(const Bint &rhs) const {
  Bint res = *this;
  res &= rhs;
  return res;
}
Bint Bint::operator|(const Bint &rhs) const {
  Bint res = *this;
  res |= rhs;
  return res;
}
Bint Bint::operator^(const Bint &rhs) const {
  Bint res = *this;
  res ^= rhs;
  return res;
}

// Assignment with Bitwise operators
Bint &Bint::operator&=(const Bint &rhs) {
  size_t ms = std::min(this->data.size(), rhs.data.size());
  data.resize(ms);

  // only up to the point where both numbers have digits
  // since everything else will always be identically zero
  for (int i = 0; i < ms; i++) {
    data[i] &= rhs.data[i];
  }
  resize();
  return *this;
}
Bint &Bint::operator|=(const Bint &rhs) {
  size_t ms = std::max(this->data.size(), rhs.data.size());
  data.resize(ms);

  // only up to the point where the rhs has digits
  // since everything else is just copied from the lhs.
  // we could do this for the lhs as well, but
  // that would require copying over the rhs which is
  // basically the same as | as well
  // similar thing for ^=
  for (int i = 0; i < rhs.data.size(); i++) {
    data[i] |= rhs.data[i];
  }

  return *this;
}

Bint &Bint::operator^=(const Bint &rhs) {
  size_t ms = std::max(this->data.size(), rhs.data.size());
  data.resize(ms);

  for (int i = 0; i < rhs.data.size(); i++) {
    data[i] ^= rhs.data[i];
  }
  resize();
  return *this;
}

// Bit shifts

Bint &Bint::operator<<=(const size_t shift) {
  size_t n = data.size();

  size_t block_shifts = shift / 32;
  size_t inside_shift = shift % 32;
  if (inside_shift > 0) {
    LBLOCK carry = 0;
    size_t anti_shift = 32 - inside_shift;
    for (int i = 0; i < n; i++) {
      carry += ((uint64_t)data[i] << 32);
      data[i] = (carry >> anti_shift);
      carry >>= 32;
    }
    if ((carry >> anti_shift) > 0) {
      data.push_back(carry >> anti_shift);
    }
  }
  data.insert(data.begin(), block_shifts, 0);
  return *this;
}

Bint &Bint::operator>>=(const size_t shift) {
  size_t n = data.size();

  size_t block_shifts = shift / 32;
  size_t inside_shift = shift % 32;
  // if the leading digit became 0 after the shift inside the blocks
  bool leading_zero = false;

  if (inside_shift > 0) {
    LBLOCK carry = data[n - 1];
    data[n - 1] = (carry >> inside_shift);
    leading_zero = (data[n - 1] == 0);
    carry <<= 32;

    for (int i = data.size() - 2; i >= 0; i--) {
      carry += data[i];
      data[i] = (carry >> inside_shift);
      carry <<= 32;
    }
  }

  // literally became zero
  if (block_shifts >= n || (leading_zero && block_shifts == n - 1)) {
    sgn = 0;
    data.clear();
    return *this;
  }

  if (block_shifts > 0 || leading_zero) {
    data = std::vector<BLOCK>(data.begin() + block_shifts,
                              data.end() - (leading_zero ? 1 : 0));
  }

  return *this;
}

Bint Bint::operator<<(const size_t shift) const {
  Bint res = *this;
  res <<= shift;
  return res;
}
Bint Bint::operator>>(const size_t shift) const {
  Bint res = *this;
  res >>= shift;
  return res;
}

// Increment and Decrement

Bint &Bint::operator++() {
  *this += 1;
  return *this;
}
Bint Bint::operator++(int) {
  Bint temp = *this;
  ++*this;
  return temp;
}
Bint &Bint::operator--() {
  *this += -1;
  return *this;
}
Bint Bint::operator--(int) {
  Bint temp = *this;
  --*this;
  return temp;
}

// Comparison

bool Bint::operator==(const Bint &rhs) const {
  if (this->sgn != rhs.sgn) {
    return false;
  }
  return this->data == rhs.data;
}

bool Bint::operator<(const Bint &rhs) const {
  if (this->sgn > rhs.sgn) {
    return false;
  }
  if (this->sgn < rhs.sgn) {
    return true;
  }
  if (this->sgn == rhs.sgn && this->sgn == 0) {
    return false;
  }

  bool reversed = sgn == -1;

  if (this->data.size() != rhs.data.size()) {
    return (this->data.size() < rhs.data.size()) ^ reversed;
  }

  for (int i = data.size() - 1; i >= 0; i--) {
    if (this->data[i] != rhs.data[i]) {
      return (this->data[i] < rhs.data[i]) ^ reversed;
    }
  }
  return false;
}

bool Bint::operator!=(const Bint &rhs) const { return !(*this == rhs); }
bool Bint::operator<=(const Bint &rhs) const {
  return (*this < rhs || *this == rhs);
} // a bit slow, but I don't care
bool Bint::operator>(const Bint &rhs) const { return !(*this <= rhs); }
bool Bint::operator>=(const Bint &rhs) const { return !(*this < rhs); }

// Printing
std::istream &operator>>(std::istream &strm, Bint &n) {
  std::string str;
  strm >> str;
  n.assign_decimal(str);
  return strm;
}

std::ostream &operator<<(std::ostream &strm, const Bint &n) {
  strm << n.to_string();
  return strm;
}

// Explicit conversion

bool Bint::to_bool() const { return sgn != 0; }

std::string Bint::to_string() const {
  if (this->sgn == 0) {
    return "0";
  }

  std::string res = "";

  Bint r = *this;
  int8_t sgn = r.sgn;

  while (r.sgn != 0) {
    std::pair<Bint, BLOCK> p = div_m(r, DBLOCK_SIZE);
    r = p.first;

    std::string dblock = std::to_string(p.second);
    dblock.insert(0, 9 - dblock.length(), '0');
    res = dblock + res;
  }

  res.erase(0, res.find_first_not_of('0'));

  // res += "{ ";
  // res += (sgn < 0 ? "- " : "");
  // for (BLOCK b : data) {
  //   res += std::to_string(b);
  //   res += " ";
  // }
  // res += "}";

  return ((sgn == -1) ? "-" : "") + res;
}

// Other

Bint Bint::operator-() const {
  Bint res = *this;
  res.sgn = -res.sgn;
  return res;
}
Bint Bint::operator~() const {
  Bint res = *this;
  for (int i = 0; i < res.data.size(); i++) {
    res.data[i] = ~res.data[i];
  }
  return res;
}

Bint Bint::abs() const {
  Bint res = *this;
  res.sgn = 1;
  return res;
}

size_t Bint::bit_length() const {
  if (sgn == 0) {
    return 32;
  }

  return 32 * data.size() - (int)clz(data.back());
}

std::pair<Bint, BLOCK> Bint::div_m(const Bint &dividend, const BLOCK &divisor) {
  if (divisor == 0) {
    throw std::invalid_argument("Cannot divide by 0");
  }
  if (dividend.sgn == 0) {
    return {0, 0};
  }
  if (dividend.abs() < divisor) {
    return {0, dividend.data[0]};
  }

  LBLOCK d = divisor, r = 0;

  Bint q;
  q.data.reserve(dividend.data.size());

  for (int i = dividend.data.size() - 1; i >= 0; i--) {
    r <<= 32;
    r += dividend.data[i];

    q.data.push_back(r / d);

    r %= d;
  }

  std::reverse(q.data.begin(), q.data.end());
  q.resize();
  q.sgn = 1;

  if (dividend.sgn == -1) {
    if (r == 0) {
      return {-q, 0};
    }
    return {-(q + 1), d - r};
  }
  return {q, r};
}

std::pair<Bint, Bint> Bint::div(const Bint &dividend, const Bint &divisor) {
  if (divisor.sgn == 0) {
    throw std::invalid_argument("Cannot divide by 0");
  }

  if (dividend.abs() < divisor.abs()) {
    return {0, dividend};
  }
  std::pair<Bint, Bint> res;

  if (divisor.data.size() == 1) {
    res = div_m(dividend.abs(), divisor.data[0]);
  } else {
    // Knuth's algorithm D

    Bint u = dividend, v = divisor, q;
    u.sgn = 1;
    v.sgn = 1;
    q.sgn = 1;

    size_t leading_zeroes = 32 * v.data.size() - (int)v.bit_length();

    u <<= leading_zeroes;
    v <<= leading_zeroes;

    size_t n = v.data.size();
    size_t m = u.data.size() - n;

    // std::cout << u << " " << v << " " << n << " " << m << std::endl;

    q.data.resize(m + 1);
    u.data.push_back(0);

    for (int j = m; j >= 0; j--) {
      // std::cout << u << std::endl;

      LBLOCK b = (LBLOCK)1 << 32;
      LBLOCK rhat = (((LBLOCK)u.data[j + n] << 32) + u.data[j + n - 1]);
      LBLOCK qhat = rhat / v.data[n - 1];
      rhat %= v.data[n - 1];

      // std::cout << "qhat before: " << qhat << std::endl;

      if (qhat == 0) {
        continue;
      }

      if (qhat == b) {
        qhat--;
        rhat += v.data[n - 1];
      }

      while (rhat < b &&
             qhat * v.data[n - 2] > (rhat << 32) + u.data[j + n - 2]) {
        qhat--;
        rhat += v.data[n - 1];
      }

      // std::cout << "qhat after: " << qhat << std::endl;

      Bint rm = v * (BLOCK)qhat;
      rm.data.resize(n + 1);
      rm = ~rm + 1;
      rm.data.resize(n + 1);

      // std::cout << j << " " << rm << std::endl;

      LBLOCK carry = 0;
      for (int i = 0; i <= n; i++) {
        carry += u.data[j + i];
        carry += rm.data[i];
        u.data[j + i] = carry; // this will keep only the last 32 bits
        carry >>= 32;
      }

      // if 2 ** 32's complement subtraction didn't carry
      // then the result was negative (i.e. qhat is too large)
      if (carry == 0) {
        q.data[j] = qhat - 1;
        carry = 0;

        for (int i = 0; i < n; i++) {
          carry += u.data[j + i];
          carry += v.data[i];
          u.data[i] = carry; // this will keep only the last 32 bits
          carry >>= 32;
        }
        carry += u.data[j + n];
        u.data[j + n] = carry;
      } else {
        q.data[j] = qhat;
      }
    }
    q.resize();
    u.resize();

    u >>= leading_zeroes;

    res = {q, u};
  }

  res.first.sgn = (dividend.sgn == divisor.sgn) ? 1 : -1;
  res.second.sgn = res.second.sgn == 0 ? 0 : dividend.sgn;
  return res;
}

Bint Bint::mul(const Bint &a, const Bint &b) {
  if (a.sgn == 0 || b.sgn == 0) {
    return 0;
  }

  size_t n = a.data.size(), m = b.data.size();
  // std::cout << n << " " << m << std::endl;
  if (m == 1) {
    return a * b.data[0];
  }
  if (n == 1) {
    return b * a.data[0];
  }

  size_t O = std::max(n, m) / 2;

  if (O < KARATSUBA_CUTOFF) {
    Bint res;
    res.sgn = 1;
    res.data.resize(n + m + 1);

    for (int i = 0; i < m; i++) {
      LBLOCK carry = 0;

      for (int j = 0; j < n; j++) {
        carry += res.data[i + j] + (LBLOCK)a.data[j] * b.data[i];
        res.data[i + j] = carry; // this will keep only the last 32 bits
        carry >>= 32;
      }
      res.data[i + n] = carry;
    }
    res.resize();

    // std::cout << "Found: " << a << " * " << b << " = " << res << std::endl;

    return res;
  }

  // Karatsuba, TODO: optimize further

  Bint x1, x2;
  Bint y1, y2;

  x1.data = std::vector<BLOCK>(a.data.begin() + O, a.data.end());
  y1.data = std::vector<BLOCK>(b.data.begin() + O, b.data.end());
  x2.data = std::vector<BLOCK>(a.data.begin(), a.data.begin() + O);
  y2.data = std::vector<BLOCK>(b.data.begin(), b.data.begin() + O);

  x1.sgn = 1;
  x2.sgn = 1;
  y1.sgn = 1;
  y2.sgn = 1;

  // for (int i = 0; i < O; i++) {
  //   x2.data.push_back(a.data[i]);
  //   y2.data.push_back(b.data[i]);
  // }
  // for (int i = 0; i < n - O; i++) {
  //   x1.data.push_back(a.data[O + i]);
  // }
  // for (int i = 0; i < m - O; i++) {
  //   y1.data.push_back(b.data[O + i]);
  // }

  // std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;

  // Bint x1 = a, x2 = a;
  // Bint y1 = b, y2 = b;

  // x2.data.resize(O);
  // x1 >>= O * 32;
  // y2.data.resize(O);
  // y1 >>= O * 32;

  Bint z1 = mul(x1 + x2, y1 + y2);
  z1 -= x1 + x2;

  Bint x1y1 = mul(x1, y1);
  Bint res = mul(x2, y2);

  res.data.resize(n + m + 1);
  std::copy(x1.data.begin(), x1.data.end(), res.data.begin() + O * 2);

  // size_t offset_2 = O * 2;
  // for (int i = 0; i < x1.data.size(); i++) {
  //   res.data[i + offset_2] = x1.data[i];
  // }

  LBLOCK carry = 0;
  for (int i = O; i < z1.data.size() + O; i++) {
    carry += res.data[i];
    carry += z1.data[i - O];
    res.data[i] = carry; // this will keep only the last 32 bits
    carry >>= 32;
  }
  for (int i = z1.data.size() + O; i < n + m + 1; i++) {
    if (carry == 0)
      break;

    carry += res.data[i];
    res.data[i] = carry; // this will keep only the last 32 bits
    carry >>= 32;
  }
  if (carry > 0) {
    res.data.push_back(carry);
  }

  res.resize();
  // std::cout << x1 + x2 << " " << z1 << std::endl;
  // std::cout << "Found: " << a << " * " << b << " = " << res << std::endl;

  res.sgn = (a.sgn == b.sgn) ? 1 : -1;

  return res;
}

Bint Bint::pow(const Bint &a, const BLOCK &b) {
  Bint res = 1;
  Bint a2 = a;

  for (int j = 0; j < 32 - clz(b); j++) {
    if (b & (1 << j)) {
      res *= a2;
    }

    a2 *= a2;
  }
  return res;
}

Bint Bint::pow(const Bint &a, const Bint &b, const Bint &M) {
  Bint res = 1;
  Bint a2 = a;

  for (int i = 0; i < b.data.size() - 1; i++) {
    for (int j = 0; j < 32; j++) {
      if (b.data[i] & (1 << j)) {
        res *= a2;
        res %= M;
      }
      a2 *= a2;
      a2 %= M;
    }
  }

  for (int j = 0; j < 32 - clz(b.data.back()); j++) {
    if (b.data.back() & (1 << j)) {
      res *= a2;
      res %= M;
    }
    a2 *= a2;
    a2 %= M;
  }
  return res;
}

Bint Bint::integral_rth_root(const Bint &n, const BLOCK &r) {
  Bint lo = (1 << ((n.bit_length() - 1) / r));
  Bint hi = lo << 1;

  while (hi > lo) {
    Bint mid = (hi + lo) >> 1;

    if (pow(mid, r) < n) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

Bint Bint::rand(const Bint &n) {
  if (n.sgn == 0) {
    return 0;
  }
  Bint r = n;

  uint32_t seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 twister(seed);

  std::uniform_int_distribution<BLOCK> first_digit(
      0, (1 << (32 - clz(n.data.back()))));

  do {
    r.data.back() = first_digit(twister);
    for (int i = 0; i < r.data.size() - 1; i++) {
      r.data[i] = twister();
    }
  } while (n.sgn == -1 ? r <= n : r >= n);
  r.resize();
  return r;
}

Bint Bint::gcd(const Bint &a, const Bint &b) {
  if (b == 0) {
    return a.abs();
  }
  return gcd(b, a % b);
}
