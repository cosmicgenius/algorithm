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
    if n == 1:
        return v
    return None

for N in range(192837, 10000000):        
    diff = False
    mod2 = -1
    for i in range(1, N):
        m = i * i % N
        if i * i > N: #int(m ** (1/2)) ** 2 != m
            if m == 0:
                continue
            v = works(m)
            if v != None:
                if mod2 == -1:
                    mod2 = v[2] % 2
                if v[2] % 2 != mod2:
                    print(N, i, m, v, mod2)
                    diff = True
                    break

    if not diff:
        print(N)
