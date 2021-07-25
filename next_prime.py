from miller_rabin import probably_prime

if __name__ == '__main__':
    while True:
        n = int(input('n: '))
        while not probably_prime(n):
            n += 1
        print(f'the next prime after n is {n}')
    
