mod ntru;
mod polynomial;
use ntru::*;
use polynomial::*;

fn main() {
    let n = 7;
    let p = 3;
    let q = 32;
    let ntru = Ntru::new(n, p, q);
    let (f, f_p, h) = loop {
        if let Some((f, f_p, h)) = ntru.generate_keys() {
            break (f, f_p, h);
        }
    };
    println!("f = {:?}", f);
    println!("f_p = {:?}", f_p);
    println!("h = {:?}", h);
    let m = Polynomial::new_from_ints(vec![-1, 0, 0, 0, -1, 1, 1]);
    println!("m = {:?}", m);
    let e = ntru.encrypt(m.clone(), h);
    println!("e = {:?}", e);
    let c = ntru.decrypt(e, f, f_p);
    println!("c = {:?}", c);
    assert_eq!(m, c);
}
