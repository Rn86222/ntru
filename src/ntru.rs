use crate::polynomial::*;
use num_rational::Rational64 as Rational;
use rand::Rng;

pub struct Ntru {
    n: usize,
    p: Int,
    q: Int,
}

impl Ntru {
    pub fn new(n: usize, p: Int, q: Int) -> Ntru {
        Ntru { n, p, q }
    }

    fn generate_random_numbers(n: usize, left: Int, right: Int) -> Vec<Int> {
        let mut rng = rand::thread_rng();
        let mut v = Vec::new();
        for _ in 0..n {
            v.push(rng.gen_range(left..=right));
        }
        v
    }

    pub fn generate_keys(&self) -> Option<(Polynomial, Polynomial, Polynomial)> {
        let f = Polynomial::new_from_ints(Self::generate_random_numbers(self.n, -1, 1));
        let g = Polynomial::new_from_ints(Self::generate_random_numbers(self.n, -1, 1));
        let f_p = f.clone().inv(self.p, self.n)?;
        let f_q = f.clone().inv(self.q, self.n)?;
        let pf_q = f_q.mul_const(Rational::from(self.p));
        let h = pf_q.cyclic_covolution(g, self.n).coef_modulo(self.q)?;
        Some((f, f_p, h))
    }

    pub fn encrypt(&self, m: Polynomial, h: Polynomial) -> Polynomial {
        let r = Polynomial::new_from_ints(Self::generate_random_numbers(self.n, -1, 1));
        let rh = r.cyclic_covolution(h, self.n);
        (rh + m).coef_modulo(self.q).unwrap()
    }

    pub fn decrypt(&self, e: Polynomial, f: Polynomial, f_p: Polynomial) -> Polynomial {
        let a = e
            .cyclic_covolution(f, self.n)
            .coef_modulo_into_interval(self.q, -self.q / 2, self.q / 2, true)
            .unwrap();
        let b = a
            .coef_modulo_into_interval(self.p, -self.p / 2, self.p / 2, false)
            .unwrap();
        b.cyclic_covolution(f_p, self.n)
            .coef_modulo_into_interval(self.p, -self.p / 2, self.p / 2, false)
            .unwrap()
    }
}
