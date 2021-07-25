import math
import random
from miller_rabin import probably_prime
from perfect_power import perfect_power

def poly(x, n):
    return (x * x + 1) % n

def proper_factor(n):
    res = perfect_power(n)
    if res != None:
        a, b = res
        return a ** (b // 2)

    g = n
    while g == n:
        turtle = random.randint(1, n - 1)
        hare = turtle
        g = 1

        while g == 1:
            turtle = poly(turtle, n)
            hare = poly(poly(hare, n), n)
            g = math.gcd(abs(turtle - hare), n)
    return math.gcd(g, n)

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

        
