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
        
    for i in range(2, min(n-1, int(2*(n.bit_length() / 1.44269504)**2))): # oops import math is allowed
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
print(is_prime(1128736817236879))
