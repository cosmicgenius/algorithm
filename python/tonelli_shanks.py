import math
import itertools
import random
from miller_rabin import probably_prime
from perfect_power import perfect_power
from nu_two import nu_two

def legendre_symbol(x, p):
    if p == 2:
        return (x % p + p) % p
    power = pow(x % p, (p - 1) // 2, p)
    if power == p - 1:
        return -1
    return power

def pow2_ord(t, p):
    t %= p
    if t == 0:
        return 0

    v = t
    for o in itertools.count(0):
        if v == 1:
            return o
        v = v * v % p


def square_root_prime(n, p):
    n = ((n % p) + p) % p
    if n == 0:
        return 0
    if p == 2:
        return n

    if legendre_symbol(n, p) == -1:
        return None
    
    if (p + 1) % 4 == 0:
        return pow(n, (p + 1) // 4, p)

    s, q = nu_two(p - 1)
    z = 2

    while legendre_symbol(z, p) != -1:
        z = random.randint(2, p - 2)

    L = pow(n, (q - 1) // 2, p)

    M = s
    c = pow(z, q, p)
    R = L * n % p
    t = L * R % p

    while t != 1:
        o = pow2_ord(t, p)
        b = pow(c, pow(2, M - o - 1, p - 1), p)

        M = o
        c = b * b % p
        t = t * c % p
        R = R * b % p

    return R

# if already computed for a smaller power of p
def square_root(n, m, r_precomp=None):
    if probably_prime(m):
        return square_root_prime(n, m)

    res = perfect_power(m)

    if res == None:
        return None
    else: 
        p, e = res
        if not probably_prime(p):
            return None
        
        if r_precomp == None:
            r = square_root_prime(n, p)
            if r == None:
                return None
            
            if p == 2:
                pow_p = p
                for k in range(1, e):
                    if (n - r * r) // pow_p % 2 == 1:
                        return None
                    r += pow_p * (((n - r * r) // pow_p // 2) % p)
                    pow_p *= p
            else:
                pow_p = p
                for k in range(1, e):
                    r += pow_p * ((((n - r * r) // pow_p) * pow(2 * r, -1, p) % p + p) % p)
                    pow_p *= p
            
            return r
        else:
            pow_p = m // p
            r = r_precomp
            if p == 2:
                if (n - r * r) // pow_p % 2 == 1:
                    return None
                r += pow_p * (((n - r * r) // pow_p // 2) % p)
            else:
                r += pow_p * ((((n - r * r) // pow_p) * pow(2 * r, -1, p) % p + p) % p)
            return r

if __name__ == '__main__':
    while True:
        n = int(input('n: '))
        m = int(input('m: '))
        
        r = square_root(n, m)
        if r == None:
            print('m is not a power of a prime, or n is not a quadratic residue mod m')
        else:
            r1 = min(r, m - r)
            r2 = max(r, m - r)
            if r1 == r2:
                print(f'the root of x^2 = n mod m is {r1}')
            else:
                print(f'the roots of x^2 = n mod m are {r1} and {r2}')
