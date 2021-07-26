import math
import time
import itertools
from miller_rabin import probably_prime
from perfect_power import integral_rth_root, perfect_power
from prime_sieve import primes_less_than
from tonelli_shanks import legendre_symbol, square_root
from nu_two import nu_two, nu_p
from f2_gaussian_elimination import generate_subset
from pollard_rho import prime_factors as prime_factors_small

def proper_factor(n):
    res = perfect_power(n)
    if res != None:
        a, b = res
        return a ** (b // 2)

    if n % 2 == 0:
        nu, _ = nu_two(n)
        return 2 ** nu

    approx_B = math.exp(math.sqrt(math.log(n) * math.log(math.log(n)) / 8))
    V = int(approx_B * (math.log(approx_B) + 1.25) * 2.25)

    prime_base = primes_less_than(V)
    prime_base = [p for p in prime_base if legendre_symbol(n, p) == 1]
    prime_base_back = { p: i for i, p in enumerate(prime_base) }
    B = len(prime_base)
    M = int(B ** 3 * 1.8)

    print(approx_B, B, M, V)
    print(prime_base)

    for i in itertools.count(0):
        # evaluate (x+a)^2 - b
        a = math.ceil(math.sqrt(n)) + i * M * 2
        b = n
        print('new a:', a, 'b:', b)

        #vals = [(x + a - M) ** 2 - b for x in range(2 * M)]
        vals = [0] * (2 * M)
        for x in range(2 * M):
            if x % 1000000 == 0:
                print(x)
            vals[x] = (x + a - M) ** 2 - b
        #val_pf = [{} for x in range(2 * M)]
        #print(vals)
        
        for p in prime_base:
            print(p, 'started')
            st1, st2 = 0, 0
            if p == 2:
                st2 = (a + b + M) % 2
            else:
                rt = square_root(b, p)
                st1 = ((-a - rt + M) % p + p) % p
                st2 = ((-a + rt + M) % p + p) % p
                #print(st1, st2)
                
                for i in range(st1, 2 * M, p):
                    #print('-', i, vals[i], vals[i] % p)
                    r, d = nu_p(vals[i], p)
                    vals[i] = d
                    #val_pf[i][p] = r
                
            for i in range(st2, 2 * M, p):
                #print('+', i, vals[i])
                r, d = nu_p(vals[i], p)
                vals[i] = d
                #val_pf[i][p] = r

        #print(vals)
        print('a')
        smooth_matrix = []
        matrix_indices = []
        for i in range(2 * M):
            if vals[i] == 1:
                v = (i + a - M) ** 2 - b
                print('===', i, v)

                # we want to multiply things to make a square,
                # if it's already a square, it's not useful to us
                isqrt = integral_rth_root(v, 2)
                #print(isqrt)
                if isqrt ** 2 == v:
                    continue

                smooth_vector = [0] * B

                prime_fac = prime_factors_small(v)
                #print(prime_fac)
                for p in prime_fac:
                    smooth_vector[prime_base_back[p]] += 1
                
                #print(smooth_vector)
        
                smooth_matrix.append(smooth_vector)
                matrix_indices.append(i)
                if len(smooth_matrix) == B + 1:
                    break

        # deep copy
        smooth_matrix2 = [vec[:] for vec in smooth_matrix]
        print(B, len(smooth_matrix))
        #print(smooth_matrix)

        subset = generate_subset(smooth_matrix)
        print(subset)
        
        if subset == None:
            continue

        square1_rt = 1
        square2_rt = 1
        prime_powers = [0] * B

        for s in subset:
            square1_rt *= abs(matrix_indices[s] + a - M)
            for i, e in enumerate(smooth_matrix2[s]):
                prime_powers[i] += e
        #print(prime_powers)

        for p, e in zip(prime_base, prime_powers):
            square2_rt *= pow(p, e // 2)

        print(f'{square1_rt=}, {square2_rt=}')
        
        g = math.gcd(square1_rt - square2_rt, n)
        print(f'{g=}')
        if g > 1 and g < n:
            return g

        g = math.gcd(square1_rt + square2_rt, n)
        print(f'{g=}')
        if g > 1 and g < n:
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
