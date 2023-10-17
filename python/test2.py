from miller_rabin import probably_prime

b = [2, 3, 5, 7]
def works(n):
    going = True
    v = [0, 0, 0]
    while going:
        going = False
        if n % 2 == 0:
            going = True
            v[0] += 1
            n //= 2
        if n % 3 == 0:
            going = True
            v[1] += 1
            n //= 3
        if n % 5 == 0:
            going = True
            v[2] += 1
            n //= 5
    if n == 1 or n == -1:
        return v
    return None
        
#for N in range(12983, 1000000, 6):
N = 17111
possible_v = set()
for i in range(1, N):
    m = i * i  - N * 4
    if m == 0:
        continue
    v = works(m)
    if v != None:
        possible_v.add((v[0] % 2, v[1] % 2, v[2] % 2))
        print(i, m, v)

print('P' if probably_prime(N) else 'C', N, list(possible_v))
