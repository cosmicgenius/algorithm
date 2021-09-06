import random
import math
import time
import functools
from miller_rabin import probably_prime
from perfect_power import perfect_power

B = 10**5

# basically useless
@functools.lru_cache(maxsize=65536)
def invert(x, n):
    return pow(x, -1, n)

def ec_group_add(a, b, x1, y1, x2, y2, n):
    #print(x1, y1, x2, y2)
    if x1 == None and y1 == None:
        return (x2, y2)

    if math.gcd(abs(x1 - x2), n) > 1 and not (x1 == x2 and y1 == y2):
        return (None, None)
    if math.gcd(y1, n) > 1:
        return (None, None)

    m = 1
    yint = 0

    if x1 == x2 and y1 == y2:
        if math.gcd(2 * y1 % n, n) != 1:
            print('-', 2 * y1 % n)
        m = ((3 * x1 * x1 + a) % n) * invert(2 * y1 % n, n) % n
    else:
        if math.gcd((x2 - x1 + n) % n, n) != 1:
            print('+', (x2 - x1 + n) % n)
        m = abs(y2 - y1) * invert((x2 - x1 + n) % n, n) % n

    yint = ((y1 - x1 * m) % n + n) % n

    x3 = ((m * m - x1 - x2) % n + n) % n
    y3 = ((m * x3 + yint) % n + n) % n

    return (x3, y3)

# computes k(x1, y1) under group operation and returns 
# (k, val) where val is the answer if it never hits the identity
# (l, (None, None)) where l is a multiple of the order of (x1, y1) if it does
def fast_power(k, a, b, x1, y1, n):
    power_2 = (x1, y1)
    val = (None, None)

    while k > 0:
        if k % 2 == 1:
            s = ec_group_add(a, b, val[0], val[1], power_2[0], power_2[1], n)
            if s == (None, None):
                return (abs(val[0] - power_2[0]), (None, None))
            val = s

        s = ec_group_add(a, b, power_2[0], power_2[1], power_2[0], power_2[1], n)
        if s == (None, None):
            return (0, power_2)
        power_2 = s

        k //= 2
    return (k, val)

def proper_factor(n):
    res = perfect_power(n)
    if res != None:
        a, b = res
        return a ** (b // 2)

    if n % 2 == 0:
        return 2

    g = n
    while g == n:
        x0 = random.randint(0, n - 1)
        y0 = random.randint(0, n - 1)
        a = random.randint(0, n - 1)
        b = ((y0 * y0 - x0 * x0 * x0 - a * x0) % n + n) % n

        cur = (x0, y0)
        for k in range(1, B+1):
            #if k % 1000 == 0:
            #    print(k)
            res = fast_power(k, a, b, cur[0], cur[1], n)
            if res[1] == (None, None):
                g = math.gcd(n, res[0])
                break
            else:
                cur = res[1]
        #print(invert.cache_info())
    return g

def prime_factors(n):
    if probably_prime(n):
        return [n]

    f1 = proper_factor(n)
    f2 = n // f1

    primes1 = prime_factors(f1)
    primes2 = prime_factors(f2)

    return [*primes1, *primes2]

if __name__ == '__main__':
    while True:
        n = int(input('n: '))

        if n == 1:
            print('n is 1')
            continue

        start = time.time_ns()
    
        prime_power = {}
        prime_fact = prime_factors(n)
        prime_fact.sort()

        for p in prime_fact:
            if p not in prime_power:
                prime_power[p] = 0
            prime_power[p] += 1

        def power_text(p, e):
            if e > 1:
                return f'{p} ^ {e}'
            return str(p)

        print(f'the prime factorization of n is n = {" * ".join([power_text(p, e) for p, e in prime_power.items()])}')
        print(f'time taken: {round((time.time_ns() - start) / 10**9, 3)}s')
