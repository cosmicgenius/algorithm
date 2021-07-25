import itertools

def nu_two(n):
    for r in itertools.count(0):
        if n % 2 == 1:
            return (r, n)
        n //= 2
