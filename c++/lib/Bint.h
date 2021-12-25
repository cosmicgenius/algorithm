#ifndef BINT_H_
#define BINT_H_

#include <algorithm>
#include <bitset>
#include <climits>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

typedef uint32_t BLOCK;
typedef uint64_t LBLOCK;

class Bint {
private:
  friend std::ostream &operator<<(std::ostream &strm, const Bint &n);
  friend std::istream &operator>>(std::istream &strm, Bint &n);

  // Sets to 0
  void set_0();

  // Assign number to a string in decimal
  void assign_decimal(const std::string &str);

  // Resize data vector to remove leading zeroes
  void resize();

  // Counts leading zeroes
  static uint8_t clz(const BLOCK &n);

  // Checks if a has a smaller absolute value than b
  // Returns abs(a) < abs(b)
  static bool magnitude_comp(const Bint &a, const Bint &b);

  // Multiplication
  // Assigns answer to res (which is assumed to start at 0)
  static void assign_mul(const Bint &a, const Bint &b, Bint &res);

  /*
    Returns a pair containing floor(|divident| / |divisor|) and |dividend| %
    |divisor|.

    Assigns the answer to {quotient, remainder} (which are assumed to start 0)
  */
  static void assign_div_abs(const Bint &dividend, const BLOCK &divisor,
                             Bint &quotient, BLOCK &remainder);
  static void assign_div_abs(const Bint &dividend, const Bint &divisor,
                             Bint &quotient, Bint &remainder);

  /*
    Returns a pair containing floor(divident / divisor) and
    dividend % divisor, where 0 <= dividend % divisor < |divisor|.
    Hence, -7 / 3 = -3 and -7 % 3 = 2.

    Assigns the answer to {quotient, remainder} (which are assumed to start 0)
  */
  static void assign_div_m(const Bint &dividend, const Bint &divisor,
                           Bint &quotient, Bint &remainder);
  static void assign_div_m(const Bint &dividend, const BLOCK &divisor,
                           Bint &quotient, BLOCK &remainder);
  /*
    Division, see the public div.

    Assigns answer to {quotient, remainder} (which are assumed to start 0)
  */
  static void assign_div(const Bint &dividend, const Bint &divisor,
                         Bint &quotient, Bint &remainder);

  // Actual data for the big int
  // stored in base 2 ** 32 from smallest to largest, i.e.
  // actual value is data[0] + 2 ** 32 * data[1] + ...
  std::vector<BLOCK> data;

  // -1, 0, or 1 depending on the sign of the number
  int8_t sgn;

public:
  // Empty constructor

  Bint();

  // Numerical constructors

  Bint(int32_t n);
  Bint(int64_t n);
  Bint(uint32_t n);
  Bint(uint64_t n);

  Bint(const Bint &b);

  Bint operator=(Bint const &b);

  // String constructors

  Bint(const char *str);
  Bint(const std::string &str);

  // Arithmetic

  Bint operator+(const Bint &rhs) const;
  Bint operator-(const Bint &rhs) const;
  Bint operator*(const Bint &rhs) const;
  Bint operator/(const Bint &rhs) const;
  Bint operator%(const Bint &rhs) const;

  Bint operator*(const BLOCK &rhs) const;
  Bint operator/(const BLOCK &rhs) const;
  Bint operator%(const BLOCK &rhs) const;

  // Assignment with Arithmetic

  const Bint &operator+=(const Bint &rhs);
  const Bint &operator-=(const Bint &rhs);
  const Bint &operator*=(const Bint &rhs);
  const Bint &operator/=(const Bint &rhs);
  const Bint &operator%=(const Bint &rhs);

  const Bint &operator*=(const BLOCK &rhs);
  const Bint &operator/=(const BLOCK &rhs);
  const Bint &operator%=(const BLOCK &rhs);

  // Bitwise operators

  Bint operator&(const Bint &rhs) const;
  Bint operator|(const Bint &rhs) const;
  Bint operator^(const Bint &rhs) const;

  // Assignment with Bitwise operators

  Bint &operator&=(const Bint &rhs);
  Bint &operator|=(const Bint &rhs);
  Bint &operator^=(const Bint &rhs);

  // Bit shifts

  Bint &operator<<=(const size_t shift);
  Bint &operator>>=(const size_t shift);
  Bint operator<<(const size_t shift) const;
  Bint operator>>(const size_t shift) const;

  // Increment and Decrement

  Bint &operator++();
  Bint operator++(int);

  Bint &operator--();
  Bint operator--(int);

  // Comparison

  bool operator==(const Bint &rhs) const;
  bool operator<(const Bint &rhs) const;
  bool operator!=(const Bint &rhs) const;
  bool operator<=(const Bint &rhs) const;
  bool operator>(const Bint &rhs) const;
  bool operator>=(const Bint &rhs) const;

  // Explicit conversion

  bool to_bool() const;
  std::string to_string() const;
  uint32_t to_uint32_t() const;

  // Other

  Bint operator-() const;
  Bint operator~() const;
  Bint abs() const;

  // The number of bits needed to represent the number
  // Specifically, b such that 2 ** (b - 1) <= abs(n) < 2 ** b for
  // all nonzero n, and 0 for n = 0
  size_t bit_length() const;

  /*
    Returns a pair containing divident / divisor and dividend % divisor

    Division is done in the normal c++ way, i.e. -7 / 3 = -2 and -7 % 3 = -1.
  */
  static std::pair<Bint, Bint> div(const Bint &dividend, const Bint &divisor);

  /*
    Returns a pair containing floor(divident / divisor) and
    dividend % divisor, where 0 <= dividend % divisor < |divisor|.
    Hence, -7 / 3 = -3 and -7 % 3 = 2.
  */
  static std::pair<Bint, Bint> div_m(const Bint &dividend, const Bint &divisor);
  static std::pair<Bint, BLOCK> div_m(const Bint &dividend,
                                      const BLOCK &divisor);

  // Finds a ** b.
  //
  // Note: b is a 32 bit integer here since even 2 ** (2 ** 32) > 10 **
  // 1000000000
  static Bint pow(const Bint &a, const BLOCK &b);

  // Finds a ** b % M.
  //
  // Note: b does NOT have to be small here since we are modding.
  static Bint pow(const Bint &a, const Bint &b, const Bint &M);

  // For postive integers n, r, finds floor(n ** (1/r)) through binary search
  // Very slow, but probably fine to use once
  static Bint integral_rth_root(const Bint &n, const BLOCK &r);

  // Approximate binary logarithm using the first two digits
  static double log_2(const Bint &n);

  // (Pseudo)random integer from 0 to n - 1 inclusive if n is positive,
  // and n + 1 to 0 inclusive if n is negative. If n is 0, returns 0.
  //
  // Uses Mersenne Twister.
  static Bint rand(const Bint &n);

  // Self explanatory
  static Bint gcd(const Bint &a, const Bint &b);
};
#endif