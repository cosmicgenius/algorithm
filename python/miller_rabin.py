import random
from nu_two import nu_two

def probably_prime(n, rounds=1000):
    if n % 2 == 0: return n == 2
    if n % 3 == 0: return n == 3
    if n % 5 == 0: return n == 5
    if n % 7 == 0: return n == 7
    r, d = nu_two(n - 1)
    for i in range(rounds):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        
        if x != 1 and x != n - 1:
            works = False
            for j in range(r - 1):
                x = x * x % n
                if x == n - 1:
                    works = True
                    break
            if not works: 
                return False
    return True

def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    if n % 3 == 0: return n == 3
    if n % 5 == 0: return n == 5
    if n % 7 == 0: return n == 7

    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d >>= 1

    for i in range(2, min(n-1, int(2*(n.bit_length() / 1.44269504)**2))):
        x = pow(i, d, n)
        if x != 1 and x != n - 1:
            works = False
            for j in range(r - 1):
                x = x * x % n
                if x == n - 1:
                    works = True
                    break
            if not works:
                return False
    return True

if __name__ == '__main__':
    for i in range(1000000, 2000000):
        if i % 10000 == 0:
            print(i)
        if probably_prime(i) != is_prime(i):
            print(i, 'BAD')
    #while True:
    #    n = int(input('n: '))
    #    if n == 1:
    #        print('n is 1')
    #    elif n < 0:
    #        print('n is negative')
    #    else:
    #        print(f'n is {"probably prime" if probably_prime(n) else "composite"}')
