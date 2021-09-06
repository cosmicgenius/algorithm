import math
import time
import itertools
import random
from miller_rabin import probably_prime
from perfect_power import integral_rth_root, perfect_power
from prime_sieve import primes_less_than
from tonelli_shanks import legendre_symbol, square_root
from nu_two import nu_two, nu_p
from f2_gaussian_elimination import generate_subsets
from pollard_rho import prime_factors as prime_factors_small

EPSILON = 10**(-6)

def proper_factor(n):
    res = perfect_power(n)
    print(f'factoring {n}')
    if res != None:
        a, b = res
        return a ** (b // 2)

    if n % 2 == 0:
        nu, _ = nu_two(n)
        return 2 ** nu

    approx_B = math.exp(math.sqrt(math.log(n) * math.log(math.log(n)) / 8))
    V = int(approx_B * (math.log(approx_B) + 1.2) * 2.25)

    prime_base = primes_less_than(V)
    prime_base = [p for p in prime_base if legendre_symbol(n, p) == 1]
    prime_base_back = { p: i for i, p in enumerate(prime_base) }
    B = len(prime_base)
    M = min(10 ** 7 * 2, B ** 3 * 2, int(math.sqrt(n))) # conserve memory by sieving max 20 million at a time

    print(f'{approx_B = }')
    print(f'{B = }')
    print(f'{M = }')
    print(f'{V = }')
    print(f'{prime_base = }')

    sieve_results = []
    sieve_result_root = []

    # sieving for smooth numbers
    for i in itertools.count(0):
        done = False
        # evaluate (x+a)^2 - b
        a = integral_rth_root(n, 2) + i * M
        #print('a', a * a)
        b = n
        #print('b', b)
        print('sieving with new a:', a)

        #logrtn = 1 + math.log(n, 2) / 2
        vals = [0] * M
        #val_pf = [{} for x in range(2 * M)]
        #print(vals)

        max_val = (M + a - 1) ** 2 - b
        #print(max_val)
        #print(M, a, b)
        
        for p in prime_base:
            print(p, 'started')
            q = p
            logp = p.bit_length()
            rt = None
            for e in range(1, math.floor(math.log(max_val, p))):
                st1, st2 = 0, 0
                if p == 2 and e == 1:
                    st2 = (a + b) % 2
                else:
                    #print(b, q, rt)
                    rt = square_root(b, q, rt)
                    #print(rt)

                    # if b stops being a qr, then no multiple will
                    # ever have b as a qr 
                    if rt == None:
                        break

                    st1 = ((-a - rt) % q + q) % q
                    st2 = ((-a + rt) % q + q) % q

                    #print(f'{q}, {st1=}', a, rt)
                    
                    for x in range(st1, M, q):
                        vals[x] += logp
                    
                #print(f'{q}, {st2=}')
                for x in range(st2, M, q):
                    vals[x] += logp

                q *= p

        #print(vals)
        print('checking sieved results')
        #print(vals)
        t = a ** 2 - b
        #print([(x * (x + 2 * a) + t).bit_length() for x in range(M)])

        prev_results = len(sieve_results)
        for x in range(M):
            if vals[x] >= (x * (x + 2 * a) + t).bit_length():
                v = (x + a) ** 2 - b

                # we want to multiply things to make a square,
                # if it's already a square, it's not useful to us
                #print(v)
                isqrt = integral_rth_root(v, 2)
                #print(isqrt)
                if isqrt ** 2 == v:
                    continue
                print('found', x, v)
                
                sieve_results.append(v)
                sieve_result_root.append(x + a)
                #if len(sieve_results) > B:
                    #done = True
                    # don't break here, we need a bit over our prime base size 
                    # because it's possible that some of these are false positives
                    # or maybe we try breaking

        print(f'{len(sieve_results) - prev_results} results found')
        print(f'{len(sieve_results)} results in total')

        while len(sieve_results) > B:
            # evaluate results
            print('generating matrix from sieved results')
            smooth_matrix = []

            # classical trick of iterating through list backwards when 
            # you need to conditionally remove elements
            for idx, v in reversed(list(enumerate(sieve_results))):
                smooth_vector = [0] * B

                prime_fac = prime_factors_small(v)
                not_actually_smooth = False
                #print(prime_fac)
                for p in prime_fac:
                    if p not in prime_base_back:
                        not_actually_smooth = True
                        break

                if not_actually_smooth:
                    del sieve_results[idx]
                    del sieve_result_root[idx]
                    continue

                print(f'{v} is verified')
                
                for p in prime_fac:
                    smooth_vector[prime_base_back[p]] += 1
                
                #print(smooth_vector)

                smooth_matrix.append(smooth_vector)

            if len(sieve_results) <= B:
                break

            # since this list was created backwards, we need to reverse it
            smooth_matrix = list(reversed(smooth_matrix))

            # deep copy
            smooth_matrix2 = [vec[:] for vec in smooth_matrix]
            print(B, len(smooth_matrix), len(sieve_results))
            #print(smooth_matrix)

            all_subsets = generate_subsets(smooth_matrix)
            #print(all_subsets)
                

            
            for subset in all_subsets:
                square1_rt = 1
                square2_rt = 1
                prime_powers = [0] * B

                marked_for_deletion = [False] * len(sieve_results)

                for s in subset:
                    square1_rt *= sieve_result_root[s]
                    for i, e in enumerate(smooth_matrix2[s]):
                        prime_powers[i] += e
                #print(prime_powers)

                for p, e in zip(prime_base, prime_powers):
                    square2_rt *= pow(p, e // 2)

                #print(f'{square1_rt=}, {square2_rt=}')
                
                g = math.gcd(square1_rt - square2_rt, n)
                print(f'{g=}')
                if g > 1 and g < n:
                    return g

                #g = math.gcd(square1_rt + square2_rt, n)
                #print(f'{g=}')
                #if g > 1 and g < n:
                #    return g
                
                # otherwise, results are useless, throw some of them out
                #new_sieve_results = []
                #new_sieve_result_root = []
                #
                #for i in range(len(sieve_results)):
                #    if random.randint(0, 1) == 0:
                #        new_sieve_results.append(sieve_results[i])
                #        new_sieve_result_root.append(sieve_result_root[i])

                #sieve_results = new_sieve_results
                #sieve_result_root = new_sieve_result_root

                marked_for_deletion[subset[0]] = True

            for idx in range(len(sieve_results) - 1, -1, -1):
                if marked_for_deletion[idx]:
                    del sieve_results[idx]
                    del sieve_result_root[idx]

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
