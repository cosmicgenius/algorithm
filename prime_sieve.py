import math

def primes_less_than(n):
    if n <= 2:
        return []
    if n <= 3:
        return [2]
    if n <= 5:
        return [2, 3]
    if n <= 7:
        return [2, 3, 5]

    wheel_residues = []
    back = {}
    min_res = -1

    for i in range(1, 210):
        if i % 2 != 0 and i % 3 != 0 and i % 5 != 0 and i % 7 != 0:
            back[i] = len(wheel_residues)
            if min_res == -1 and i > n % 210:
                min_res = i
            wheel_residues.append(i)

    if min_res == -1:
        min_res = 1
    
    rtn = int(math.sqrt(n))
    q = rtn // 210
    is_prime = [True] * (48 * (n // 210) + back[min_res])

    # 1 is not prime
    is_prime[0] = False

    def get_is_prime(x):
        return is_prime[48 * (x // 210) + back[x % 210]]

    def set_is_prime(x, v):
        is_prime[48 * (x // 210) + back[x % 210]] = v

    qn = n // 210

    for i in range(q + 1):
        for j_idx1 in range(48):
            c = i * 210 + wheel_residues[j_idx1]
            if c >= rtn:
                break
            if c == 1:
                continue
            
            if get_is_prime(c):
                #print(c)
                # c2 should be at least c
                done = False
                for j_idx2 in range(j_idx1, 48):
                    c2 = i * 210 + wheel_residues[j_idx2]
                    if c2 * c >= n:
                        break
                        done = True
                    #print(1, f'set {c * c2}')
                    set_is_prime(c * c2, False)
                
                if done:
                    continue

                for i2 in range(i + 1, qn // c):
                    for j2 in wheel_residues:
                        #print(2, f'set {c * (i2 * 210 + j2)}')
                        set_is_prime((i2 * 210 + j2) * c, False)
                
                for j2 in wheel_residues:
                    c2 = qn // c * 210 + j2
                    if c2 * c >= n:
                        break
                    #print(3, f'set {c * c2}')
                    set_is_prime(c2 * c, False)
    primes = [210 * (x // 48) + wheel_residues[x % 48] for x, v in enumerate(is_prime) if v]
    primes = [x for x in primes if x < n]
    return [2, 3, 5, 7, *primes]


if __name__ == '__main__':
    while True:
        n = int(input('n: '))
        if n == 1:
            print('n is 1')
        elif n < 0:
            print('n is negative')
        else:
            #print(f'the primes less than n are {primes_less_than(n)}')
            print(f'there are {len(primes_less_than(n))} primes less than n')
