from timeit import timeit
import itertools

def rd(n):
    for r in itertools.count(0):
        if n % 2 == 1:
            return (r, n)
        n //= 2

def fact(n):
    ans = 1
    for i in range(2, n + 1):
        ans *= i
    return ans

def a():
    r, d = (rd(fact(52300)))
    print(r)

print(timeit(lambda: a(), number=1))
