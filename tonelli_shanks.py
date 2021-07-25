import math
import itertools
import random
from miller_rabin import probably_prime
from nu_two import nu_two

def legendre_symbol(x, p):
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


def square_root(n, p):
    n = ((n % p) + p) % p
    if n == 0:
        return 0

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

if __name__ == '__main__':
    while True:
        n = int(input('n: '))
        p = int(input('p: '))
        
        if p < 2 or not probably_prime(p):
            print('p is not prime')
        else:
            r = square_root(n, p)
            if r == None:
                print('n is not a qr mod p')
            elif r == 0:
                print('the root of x^2 = n mod p is 0')
            else:
                print(f'the roots of x^2 = n mod p are {min(r, p - r)} and {max(r, p - r)}')
