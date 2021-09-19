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

void Bint::set_0() {
  sgn = 0;
  data.clear();
}

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
  data.emplace_back(block);

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

bool Bint::magnitude_comp(const Bint &a, const Bint &b) {
  if (b.sgn == 0) {
    return false;
  }
  if (a.sgn == 0) {
    return true;
  }

  int adsz = a.data.size(), bdsz = b.data.size();

  if (adsz != bdsz) {
    return adsz < bdsz;
  }

  for (int i = adsz - 1; i >= 0; i--) {
    if (a.data[i] != b.data[i]) {
      return a.data[i] < b.data[i];
    }
  }
  return false;
}

void Bint::assign_mul(const Bint &a, const Bint &b, Bint &res) {
  if (a.sgn == 0 || b.sgn == 0) {
    return;
  }

  size_t n = a.data.size(), m = b.data.size();
  // std::cout << n << " " << m << std::endl;
  if (m == 1) {
    res = a * b.data[0];
    return;
  }
  if (n == 1) {
    res = b * a.data[0];
    return;
  }

  size_t O = std::max(n, m) / 2;

  if (O < KARATSUBA_CUTOFF) {
    if (a == b) {
      res.sgn = 1;
      res.data.resize(2 * n + 1);

      for (int i = 0; i < n; i++) {
        LBLOCK carry = 0;

        for (int j = 0; j < i; j++) {
          carry += res.data[i + j] + (LBLOCK)a.data[i] * a.data[j];
          res.data[i + j] = carry; // this will keep only the last 32 bits
          carry >>= 32;
        }
        res.data[i * 2] = carry;
      }
      res <<= 1;

      LBLOCK carry = 0;
      for (int i = 0; i < n; i++) {
        carry += res.data[i * 2] + (LBLOCK)a.data[i] * a.data[i];
        res.data[i * 2] = carry; // this will keep only the last 32 bits
        carry >>= 32;

        carry += res.data[i * 2 + 1];
        res.data[i * 2 + 1] = carry;
        carry >>= 32;
      }

      res.resize();

      // std::cout << "Found: " << a << " * " << b << " = " << res << std::endl;

      return;
    } else {
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

      return;
    }
  }

  // Karatsuba, TODO: optimize further

  if (a == b) {
    Bint x1, x2;

    x1.data = std::vector<BLOCK>(a.data.begin() + O, a.data.end());
    x2.data = std::vector<BLOCK>(a.data.begin(), a.data.begin() + O);

    x1.sgn = 1;
    x2.sgn = 1;

    Bint z1;
    assign_mul(x1 + x2, x1 + x2, z1);

    Bint x1x1;
    assign_mul(x1, x1, x1x1);
    assign_mul(x2, x2, res);
    z1 -= x1x1 + res;

    res.data.resize(n * 2 + 1);
    std::copy(x1x1.data.begin(), x1x1.data.end(), res.data.begin() + O * 2);

    LBLOCK carry = 0;
    int sz1 = z1.data.size();
    for (int i = O; i < sz1 + O; i++) {
      carry += res.data[i];
      carry += z1.data[i - O];
      res.data[i] = carry; // this will keep only the last 32 bits
      carry >>= 32;
    }
    for (int i = sz1 + O; i < n + m + 1; i++) {
      if (carry == 0)
        break;

      carry += res.data[i];
      res.data[i] = carry; // this will keep only the last 32 bits
      carry >>= 32;
    }
    if (carry > 0) {
      res.data.emplace_back(carry);
    }

    res.resize();

    return;
  } else {
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
    //   x2.data.emplace_back(a.data[i]);
    //   y2.data.emplace_back(b.data[i]);
    // }
    // for (int i = 0; i < n - O; i++) {
    //   x1.data.emplace_back(a.data[O + i]);
    // }
    // for (int i = 0; i < m - O; i++) {
    //   y1.data.emplace_back(b.data[O + i]);
    // }

    // std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;
    // std::cout << O << std::endl;

    // Bint x1 = a, x2 = a;
    // Bint y1 = b, y2 = b;

    // x2.data.resize(O);
    // x1 >>= O * 32;
    // y2.data.resize(O);
    // y1 >>= O * 32;

    Bint z1;
    assign_mul(x1 + x2, y1 + y2, z1);

    Bint x1y1;
    assign_mul(x1, y1, x1y1);
    assign_mul(x2, y2, res);
    z1 -= x1y1 + res;

    res.data.resize(n + m + 1);
    std::copy(x1y1.data.begin(), x1y1.data.end(), res.data.begin() + O * 2);

    // std::cout << res << std::endl;

    // size_t offset_2 = O * 2;
    // for (int i = 0; i < x1.data.size(); i++) {
    //   res.data[i + offset_2] = x1.data[i];
    // }

    LBLOCK carry = 0;
    int sz1 = z1.data.size();
    for (int i = O; i < sz1 + O; i++) {
      carry += res.data[i];
      carry += z1.data[i - O];
      res.data[i] = carry; // this will keep only the last 32 bits
      carry >>= 32;
    }
    for (int i = sz1 + O; i < n + m + 1; i++) {
      if (carry == 0)
        break;

      carry += res.data[i];
      res.data[i] = carry; // this will keep only the last 32 bits
      carry >>= 32;
    }
    if (carry > 0) {
      res.data.emplace_back(carry);
    }

    res.resize();
    // std::cout << x1 + x2 << " " << z1 << std::endl;
    // std::cout << "Found: " << a << " * " << b << " = " << res << std::endl;

    res.sgn = (a.sgn == b.sgn) ? 1 : -1;

    return;
  }
}

void Bint::assign_div_m(const Bint &dividend, const BLOCK &divisor,
                        Bint &quotient, BLOCK &remainder) {
  if (divisor == 0) {
    throw std::domain_error("Cannot divide by 0");
  }
  if (dividend.sgn == 0) {
    return;
  }
  if (dividend.data.size() == 1 && dividend.data[0] < divisor) {
    remainder = dividend.data[0];
    return;
  }

  LBLOCK d = divisor, r = 0;

  quotient.data.resize(dividend.data.size());

  for (int i = dividend.data.size() - 1; i >= 0; i--) {
    r <<= 32;
    r += dividend.data[i];

    quotient.data[i] = r / d;

    r %= d;
  }

  quotient.resize();
  quotient.sgn = 1;

  if (dividend.sgn == -1) {
    quotient.sgn *= -1;
    if (r == 0) {
      return;
    }

    quotient += -1;
    remainder = d - r;
    return;
  }

  remainder = r;
  return;
}

void Bint::assign_div(const Bint &dividend, const Bint &divisor, Bint &quotient,
                      Bint &remainder) {
  if (divisor.sgn == 0) {
    throw std::domain_error("Cannot divide by 0");
  }

  if (magnitude_comp(dividend, divisor)) {
    remainder = dividend;
    return;
  }

  if (divisor.data.size() == 1) {
    BLOCK r = 0;
    assign_div_m(dividend.abs(), divisor.data[0], quotient, r);
    remainder = r;
  } else {
    // Knuth's algorithm D

    Bint &u = remainder, v = divisor,
         &q = quotient; // imagine actually converting code
    u = dividend;

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
    u.data.emplace_back(0);

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

      Bint rm = v;
      rm *= (BLOCK)qhat;
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
  }

  quotient.sgn = (dividend.sgn == divisor.sgn) ? 1 : -1;
  remainder.sgn = remainder.sgn == 0 ? 0 : dividend.sgn;
}

// Constructors

Bint::Bint() : sgn(0) {}

Bint::Bint(int32_t n) : sgn((n > 0) - (n < 0)) {
  if (n != 0) {
    data.emplace_back(std::abs(n));
  }
}
Bint::Bint(int64_t n) : sgn((n > 0) - (n < 0)) {
  data.emplace_back(llabs(n));
  data.emplace_back(llabs(n) >> 32);
  resize();
}
Bint::Bint(uint32_t n) : sgn((n > 0) - (n < 0)) {
  if (n != 0) {
    data.emplace_back(n);
  }
}
Bint::Bint(uint64_t n) : sgn((n > 0) - (n < 0)) {
  data.emplace_back(n);
  data.emplace_back(n >> 32);
  resize();
}

Bint::Bint(const char *str) { assign_decimal(str); }
Bint::Bint(const std::string &str) { assign_decimal(str); }

Bint::Bint(const Bint &b) : sgn(b.sgn) { data = b.data; }
Bint Bint::operator=(Bint const &b) {
  sgn = b.sgn;
  data = b.data;
  return *this;
}

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

Bint Bint::operator*(const Bint &rhs) const {
  Bint res;
  assign_mul(*this, rhs, res);
  return res;
}
Bint Bint::operator/(const Bint &rhs) const {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  return p.first;
}
Bint Bint::operator%(const Bint &rhs) const {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  return p.second;
}

Bint Bint::operator*(const BLOCK &rhs) const {
  Bint res = *this;
  res *= rhs;
  return res;
}

Bint Bint::operator/(const BLOCK &rhs) const {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  return p.first;
}
Bint Bint::operator%(const BLOCK &rhs) const {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  return p.second;
}

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
      data.emplace_back(carry);
    }
  }
  // 2^32's complement subtraction
  else {
    if (this->data == rhs.data) {
      data.clear();
      sgn = 0;
    } else {
      // what will the sign of the result be?
      // int8_t new_sgn = (-*this > rhs) ? rhs.sgn : this->sgn;
      // std::cout << (int)new_sgn << " " << -*this << " " << rhs << " "
      // << (-*this > rhs) << std::endl;
      Bint rhs2 = rhs;
      rhs2.data.resize(ms);

      LBLOCK b = ((LBLOCK)1 << 32);
      LBLOCK borrow = 0;

      if (magnitude_comp(rhs, *this)) {
        // this dominates
        for (int i = 0; i < ms; i++) {
          borrow += data[i];
          borrow -= rhs2.data[i];
          data[i] = borrow; // this will keep only the last 32 bits
          if (borrow <= b) {
            borrow = 0;
          } else {
            borrow = -1;
          }
        }
        // sgn stays the same
      } else {
        // rhs dominates
        for (int i = 0; i < ms; i++) {
          borrow += rhs2.data[i];
          borrow -= data[i];
          data[i] = borrow; // this will keep only the last 32 bits
          if (borrow <= b) {
            borrow = 0;
          } else {
            borrow = -1;
          }
        }
        sgn = -sgn;
      }

      // if (sgn < 0) {
      //   *this = ~*this;
      //   LBLOCK carry = 1;
      //   for (int i = 0; i < data.size(); i++) {
      //     if (carry == 0)
      //       break;
      //     carry += data[i];
      //     data[i] = carry;
      //     carry >>= 32;
      //   }
      // } else if (rhs2.sgn < 0) {
      //   rhs2 = ~rhs2;
      //   LBLOCK carry = 1;
      //   for (int i = 0; i < rhs2.data.size(); i++) {
      //     if (carry == 0)
      //       break;
      //     carry += rhs2.data[i];
      //     rhs2.data[i] = carry;
      //     carry >>= 32;
      //   }
      // }

      // LBLOCK carry = 0;
      // for (int i = 0; i < ms; i++) {
      //   carry += data[i];
      //   carry += rhs2.data[i];
      //   data[i] = carry; // this will keep only the last 32 bits
      //   carry >>= 32;
      // }
      // // don't add the carry at the very end

      // // If it's negative, convert back from 2^32's complement
      // if (new_sgn < 0) {
      //   *this = ~*this;
      //   LBLOCK carry = 1;
      //   for (int i = 0; i < data.size(); i++) {
      //     if (carry == 0)
      //       break;
      //     carry += data[i];
      //     data[i] = carry;
      //     carry >>= 32;
      //   }
      // }
      // sgn = new_sgn;
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

  int s = data.size();
  for (int i = 0; i < s; i++) {
    carry += data[i] * rhs2;
    data[i] = carry; // this will keep only the last 32 bits
    carry >>= 32;
  }
  if (carry > 0) {
    data.emplace_back(carry);
  }
  return *this;
}

const Bint &Bint::operator/=(const Bint &rhs) {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  *this = p.first;
  return *this;
}

const Bint &Bint::operator%=(const Bint &rhs) {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  *this = p.second;
  return *this;
}

const Bint &Bint::operator*=(const Bint &rhs) {
  *this = *this * rhs;
  return *this;
}

const Bint &Bint::operator/=(const BLOCK &rhs) {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  *this = p.first;
  return *this;
}

const Bint &Bint::operator%=(const BLOCK &rhs) {
  std::pair<Bint, Bint> p;
  assign_div(*this, rhs, p.first, p.second);
  *this = p.second;
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
//     data.emplace_back(carry);
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
  int s = rhs.data.size();
  for (int i = 0; i < s; i++) {
    data[i] |= rhs.data[i];
  }

  return *this;
}

Bint &Bint::operator^=(const Bint &rhs) {
  int s = rhs.data.size();
  data.resize(std::max(this->data.size(), (size_t)s));

  for (int i = 0; i < s; i++) {
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
      data.emplace_back(carry >> anti_shift);
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

  int ldsz = this->data.size(), rdsz = rhs.data.size();

  bool reversed = sgn == -1;

  if (ldsz != rdsz) {
    return (ldsz < rdsz) ^ reversed;
  }

  for (int i = ldsz - 1; i >= 0; i--) {
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
    std::pair<Bint, BLOCK> p;
    assign_div_m(r, DBLOCK_SIZE, p.first, p.second);
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
  int s = res.data.size();
  for (int i = 0; i < s; i++) {
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

std::pair<Bint, Bint> Bint::div(const Bint &dividend, const Bint &divisor) {
  std::pair<Bint, Bint> p;
  assign_div(dividend, divisor, p.first, p.second);
  return p;
}

Bint Bint::pow(const Bint &a, const BLOCK &b) {
  Bint res = 1;
  Bint a2 = a;

  int clzb = clz(b);
  for (int j = 0; j < 32 - clzb; j++) {
    if (b & (1 << j)) {
      res *= a2;
    }

    a2 *= a2;
  }
  return res;
}

Bint Bint::pow(const Bint &a, const Bint &b, const Bint &M) {
  if (b == 0) {
    return 1;
  }
  if (b.sgn == -1) {
    return 0;
  }

  Bint res = 1;
  Bint a2 = a;

  int bsz = b.data.size();
  for (int i = 0; i < bsz - 1; i++) {
    for (int j = 0; j < 32; j++) {
      if (b.data[i] & (1 << j)) {
        res *= a2;
        res %= M;
      }
      a2 *= a2;
      a2 %= M;
    }
  }

  int clzbb = clz(b.data.back());
  for (int j = 0; j < 32 - clzbb; j++) {
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

  int rsz_minus_1 = r.data.size() - 1;
  do {
    r.data.back() = first_digit(twister);
    for (int i = 0; i < rsz_minus_1; i++) {
      r.data[i] = twister();
    }
  } while (n.sgn == -1 ? r <= n : r >= n);
  r.resize();
  return r;
}

Bint Bint::gcd(const Bint &a, const Bint &b) {
  Bint A = a, B = b;
  bool parity = false;
  while (A.sgn != 0 && B.sgn != 0) {
    if (parity) {
      A %= B;
    } else {
      B %= A;
    }
    parity ^= 1;
  }
  if (parity) {
    // std::cout << "gcd(" << a << ", " << b << ") = " << A << std::endl;
    return A.abs();
  }
  // std::cout << "gcd(" << a << ", " << b << ") = " << B << std::endl;
  return B.abs();
}
