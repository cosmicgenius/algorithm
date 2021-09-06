# very slow, but it is alright

import math

# equivalent ceil(n**(1/r)) but doesn't need to convert to float
def integral_rth_root(n, r):
    lo, hi = 2**(int(math.log(n, 2) // r)), 2**(int(math.log(n, 2) // r + 1))

    while hi > lo:
        mid = (hi + lo) // 2

        if mid ** r < n:
            lo = mid + 1
        else:
            hi = mid
    return lo

def perfect_power(n):
    for r in range(int(math.log(n, 2)), 1, -1):
        root = integral_rth_root(n, r)

        if root ** r == n:
            return (root, r)

    return None

if __name__ == '__main__':
    while True:
        n = int(input('n: '))
        if n == 1:
            print('n is 1')
        elif n < 0:
            print('n is negative')
        else:
            res = perfect_power(n)

            if res == None:
                print('n is not a perfect power')
            else:
                print(f'n = {res[0]} ^ {res[1]}')

