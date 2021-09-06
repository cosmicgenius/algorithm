from timeit import timeit
from memory_profiler import profile
import math
import random

def run():
    def fast(n):
        a = n.bit_length()
        
    def slow(n):
        a = math.floor(math.log(n, 2))

    print(timeit(lambda: fast(random.randint(1, 10**50)), number=1000000))
    print(timeit(lambda: slow(random.randint(1, 10**50)), number=1000000))

if __name__ == '__main__':
    run()
    run()
