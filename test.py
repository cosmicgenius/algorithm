from timeit import timeit
from random import randrange
from memory_profiler import profile

@profile
def run():
    N = 100000000

    a1 = [False] * N
    a2 = 0

    def flip_random1():
        r = randrange(0, N)
        a1[r] = not a1[r]
        
    def flip_random2():
        nonlocal a2
        a2 ^= 1 << randrange(0, N)

    print(timeit(lambda: flip_random1(), number=100))
    print(timeit(lambda: flip_random2(), number=100))

if __name__ == '__main__':
    run()
