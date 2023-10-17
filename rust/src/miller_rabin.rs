use num_bigint::BigUint;
use num_traits::{sa, FromPrimitive, Zero};
use rand::Rng;

const ROUNDS: i32 = 40;
const ONE: BigUint = BigUint::from(1 as u32);

fn two_adic_val(n: BigUint) -> (i32, BigUint) {
    let mut d = n;
    let mut r = 1;

    while (d & ONE).is_zero() {
        d >>= 1;
        r += 1;
    }
    return (r, d);
}

fn probable_prime(n: BigUint) -> bool {
    if n == BigUint::from(1 as u32) {
        return false;
    }
    if n % 2 == 0 {
        return n == BigUint::from(2 as u32);
    }
    if n % 3 == 0 {
        return n == BigUint::from(3 as u32);
    }
    if n % 5 == 0 {
        return n == BigUint::from(5 as u32);
    }
    if n % 7 == 0 {
        return n == BigUint::from(7 as u32);
    }

    let (r, d) = two_adic_val(n - 1);

    for _ in 1..ROUNDS {
        let a = rand::thread_rng().gen_range(BigUint::from(2 as u32)..n - 1);

        let mut x = BigUint::modpow(a, d, n);
    }

    return true;
}
