use std::io;

pub mod miller_rabin;

fn main() {
    println!("n:");
    let mut n_str = String::new();

    io::stdin().read_line(&mut n_str).expect("bad input");
    let n: i32 = n_str.trim().parse().expect("bad input");
}
