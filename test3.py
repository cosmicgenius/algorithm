import functools

primes = [2, 3, 5, 7]
pow2 = [2 ** n for n in range(20)]

for i in range(11, 500000):
    works = True
    for p in primes:
        if i % p == 0:
            works = False
            break
    if works:
        primes.append(i)

#print(primes)

@functools.lru_cache(maxsize=None)
def varphi(n):
    ans = n
    
    for p in primes:
        if p > n:
            break
        if n % p == 0:
            ans *= (p - 1)
            ans //= p
    return ans
    
N = 50000
best = 0
for x in range(1, 2):
    for y in primes:
        a = [y, x, y]
        for i in range(3, 20):
            a2 = varphi(a[i - 1]) + varphi(a[i - 2])
            a.append(a2)
        works = False

        #print(a)

        for i in range(1, 19):
            if a[i + 1] == a[i] and 2 * varphi(a[i]) >= a[i]:
                works = True

                if a[i] in pow2:
                    works = False

                break

        if works:
            pass
            #print(a)
        else:
            continue

        for i in range(19):
            if a[i] == a[i + 1] and a[i] in pow2:
                works = False
                break
        if works:
            if a[19] < 1000:
                print('!!!!!!!!!!!!!!!!!!', a, a[19])
                best = max(best, a[19])

print(best)
