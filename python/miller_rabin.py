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

if __name__ == '__main__':
    while True:
        n = int(input('n: '))
        if n == 1:
            print('n is 1')
        elif n < 0:
            print('n is negative')
        else:
            print(f'n is {"probably prime" if probably_prime(n) else "composite"}')
