import itertools

def nu_two(n):
    for r in itertools.count(0):
        if n % 2 == 1:
            return (r, n)
        n //= 2
        
def nu_p(n, p):
    for r in itertools.count(0):
        if n % p != 0:
            return (r, n)
        n //= p
