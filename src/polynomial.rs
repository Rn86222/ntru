use std::{
    fmt::Debug,
    ops::{Add, Mul, Sub},
};

use num_rational::Rational64 as Rational;
pub type Int = i64;
type Coefficient = Rational;

fn extended_euclid(a: Int, b: Int) -> (Int, Int, Int) {
    if b == 0 {
        return (a, 1, 0);
    }
    let (d, x, y) = extended_euclid(b, a % b);
    (d, y, x - (a / b) * y)
}

fn mod_inv(a: Int, modulo: Int) -> Option<Int> {
    if modulo == 0 {
        panic!("modulo must not be zero");
    }
    let (d, x, _) = extended_euclid(a, modulo);
    if d == 1 {
        Some(x % modulo)
    } else {
        None
    }
}

fn fraction_mod(f: Rational, modulo: Int) -> Option<Rational> {
    let (d, _, _) = extended_euclid(*f.denom(), modulo);
    if d != 1 {
        return None;
    }
    Some(Rational::from(
        ((mod_inv(*f.denom(), modulo).unwrap() * f.numer() % modulo) + modulo) % modulo,
    ))
}

fn format_coefficient(coef: Coefficient) -> String {
    if coef.denom() == &1 {
        format!("{}", coef.numer())
    } else {
        format!("{}/{}", coef.numer(), coef.denom())
    }
}

#[derive(Clone)]
pub struct Polynomial {
    coefficients: Vec<Coefficient>,
}

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        self.coefficients == other.coefficients
    }
}

impl Eq for Polynomial {}

impl Polynomial {
    pub fn new(coefficients: Vec<Coefficient>) -> Self {
        Self { coefficients }
    }

    pub fn new_from_ints(coefficients: Vec<Int>) -> Self {
        Self::new(
            coefficients
                .iter()
                .map(|x| Rational::from(*x))
                .collect::<Vec<_>>(),
        )
    }

    fn trim(mut self) -> Self {
        if self.coefficients.is_empty() {
            return self;
        }
        for i in (0..self.coefficients.len()).rev() {
            if self.coefficients[i] != Rational::from(0) {
                self.coefficients.resize(i + 1, Rational::from(0));
                return self;
            }
        }
        Self::new(vec![Rational::from(0)])
    }

    fn resize(self, other: Self) -> (Self, Self) {
        let mut a = self.coefficients;
        let mut b = other.coefficients;
        let max = a.len().max(b.len());
        a.resize(max, Rational::from(0));
        b.resize(max, Rational::from(0));
        (Self::new(a), Self::new(b))
    }

    pub fn coef_modulo(self, modulo: Int) -> Option<Self> {
        let mut coefficients = vec![Rational::from(0); self.coefficients.len()];
        for (i, coef) in coefficients.iter_mut().enumerate() {
            *coef = fraction_mod(self.coefficients[i], modulo)?;
        }
        Some(Self::new(coefficients).trim())
    }

    fn div(self, other: Self) -> (Self, Self) {
        let other_coefficients = other.trim().coefficients;
        let other_deg = other_coefficients.len() - 1;
        let mut diviedend = self.trim();
        let mut self_deg = diviedend.coefficients.len() - 1;
        let (q, r) = if self_deg >= other_deg {
            let mut q = vec![Rational::from(0); self_deg - other_deg + 1];
            while self_deg >= other_deg && diviedend.coefficients != vec![Rational::from(0)] {
                let mut d = other_coefficients.clone();
                for _ in 0..self_deg - other_deg {
                    d.insert(0, Rational::from(0));
                }
                q[self_deg - other_deg] = diviedend.coefficients[self_deg] / d[d.len() - 1];
                d = d
                    .iter()
                    .map(|x| *x * q[self_deg - other_deg])
                    .collect::<Vec<_>>();
                diviedend = (diviedend - Self::new(d)).trim();
                if diviedend.coefficients.is_empty() {
                    break;
                }
                self_deg = diviedend.coefficients.len() - 1;
            }
            let r = diviedend;
            (Self::new(q), r)
        } else {
            (Self::new(vec![Rational::from(0)]), diviedend)
        };
        (q.trim(), r.trim())
    }

    fn extended_euclid(self, other: Self) -> (Self, Self, Self) {
        let mut switch = false;
        let a = self.trim();
        let b = other.trim();
        let (mut a1, mut b1) = if a.coefficients.len() <= b.coefficients.len() {
            (a, b)
        } else {
            switch = true;
            (b, a)
        };
        let mut qs = vec![];
        let mut rs = vec![];
        while b1.coefficients != vec![Rational::from(0)] {
            let (q, r) = a1.div(b1.clone());
            qs.push(q);
            rs.push(r.clone());
            a1 = b1;
            b1 = r;
        }
        let mut s = vec![Self::new(vec![Rational::from(0)]); qs.len() + 2];
        let mut t = vec![Self::new(vec![Rational::from(0)]); qs.len() + 2];
        s[0] = Self::new(vec![Rational::from(1)]);
        s[1] = Self::new(vec![Rational::from(0)]);
        t[0] = Self::new(vec![Rational::from(0)]);
        t[1] = Self::new(vec![Rational::from(1)]);
        for i in 2..s.len() {
            s[i] = s[i - 2].clone() - s[i - 1].clone() * qs[i - 2].clone();
            t[i] = t[i - 2].clone() - t[i - 1].clone() * qs[i - 2].clone();
        }
        let gcd_vals = rs[rs.len() - 2].coefficients.clone();
        let s_out = s[s.len() - 2].clone();
        let t_out = t[t.len() - 2].clone();

        let scale_factor = gcd_vals[gcd_vals.len() - 1];
        let gcd_vals = gcd_vals
            .iter()
            .map(|x| *x / scale_factor)
            .collect::<Vec<_>>();
        let s_out = s_out
            .coefficients
            .iter()
            .map(|x| *x / scale_factor)
            .collect::<Vec<_>>();
        let t_out = t_out
            .coefficients
            .iter()
            .map(|x| *x / scale_factor)
            .collect::<Vec<_>>();

        if switch {
            (Self::new(gcd_vals), Self::new(t_out), Self::new(s_out))
        } else {
            (Self::new(gcd_vals), Self::new(s_out), Self::new(t_out))
        }
    }

    #[allow(dead_code)]
    fn remain_modulo(self, other: Self, modulo: Int) -> Option<Self> {
        let (_, remain) = self.div(other);
        remain.coef_modulo(modulo)
    }

    pub fn inv(self, modulo: Int, n: usize) -> Option<Self> {
        let mut d = vec![Rational::from(0); n + 1];
        d[0] = Rational::from(-1);
        d[n] = Rational::from(1);
        let (gcd, s, _) = self.extended_euclid(Self::new(d));
        let s = s.coef_modulo(modulo);
        if gcd.coefficients != vec![Rational::from(1)] {
            None
        } else {
            s
        }
    }

    pub fn cyclic_covolution(self, other: Self, n: usize) -> Self {
        let mut coefficients = vec![Rational::from(0); n];
        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                coefficients[(i + j) % n] += self.coefficients[i] * other.coefficients[j];
            }
        }
        Self::new(coefficients).trim()
    }

    pub fn mul_const(self, c: Coefficient) -> Self {
        let mut coefficients = self.coefficients;
        for coef in coefficients.iter_mut() {
            *coef *= c;
        }
        Self::new(coefficients).trim()
    }

    pub fn coef_modulo_into_interval(
        self,
        modulo: Int,
        left: Int,
        right: Int,
        flag: bool,
    ) -> Option<Self> {
        let mut coefficients = vec![Rational::from(0); self.coefficients.len()];
        for (i, coef) in coefficients.iter_mut().enumerate() {
            *coef = fraction_mod(self.coefficients[i], modulo)?;
            if flag {
                if *coef < Rational::from(left) {
                    *coef += modulo;
                } else if *coef >= Rational::from(right) {
                    *coef -= modulo;
                }
            } else if *coef <= Rational::from(left) {
                *coef += modulo;
            } else if *coef > Rational::from(right) {
                *coef -= modulo;
            }
        }
        Some(Self::new(coefficients).trim())
    }
}

impl Add for Polynomial {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let (a, b) = self.resize(other);
        let mut coefficients = vec![Rational::from(0); a.coefficients.len()];
        for (i, coef) in coefficients.iter_mut().enumerate() {
            *coef = a.coefficients[i] + b.coefficients[i];
        }
        Self::new(coefficients).trim()
    }
}

impl Sub for Polynomial {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let (a, b) = self.resize(other);
        let mut coefficients = vec![Rational::from(0); a.coefficients.len()];
        for (i, coef) in coefficients.iter_mut().enumerate() {
            *coef = a.coefficients[i] - b.coefficients[i];
        }
        Self::new(coefficients).trim()
    }
}

impl Mul for Polynomial {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut coefficients =
            vec![Rational::from(0); self.coefficients.len() + other.coefficients.len() - 1];
        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                coefficients[i + j] += self.coefficients[i] * other.coefficients[j];
            }
        }
        Self::new(coefficients).trim()
    }
}

impl Debug for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for (i, coef) in self.coefficients.iter().enumerate() {
            if *coef == Rational::from(0) {
                if self.coefficients.len() == 1 {
                    s += "0";
                }
                continue;
            }
            if i == 0 {
                s += &format!("{} ", &format_coefficient(*coef));
            } else if i == 1 {
                if coef > &Rational::from(0) {
                    if coef == &Rational::from(1) {
                        s += "+ x ";
                    } else {
                        s += &format!("+ {}x ", format_coefficient(*coef));
                    }
                } else if coef == &Rational::from(-1) {
                    s += "- x ";
                } else {
                    s += &format!(" - {}x", format_coefficient(-*coef));
                }
            } else if coef > &Rational::from(0) {
                if coef == &Rational::from(1) {
                    s += &format!("+ x^{} ", i);
                } else {
                    s += &format!("+ {}x^{} ", format_coefficient(*coef), i);
                }
            } else if coef == &Rational::from(-1) {
                s += &format!("- x^{} ", i);
            } else {
                s += &format!("- {}x^{} ", format_coefficient(-*coef), i);
            }
        }
        write!(f, "{}", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let f = Polynomial::new_from_ints(vec![-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1]);
        let g = Polynomial::new_from_ints(vec![-1, 0, 1, 1, 0, 1, 0, 0, -1, 0, -1]);
        let p = 3;
        let q = 32;
        let n = 11;
        let f_p = f.clone().inv(p, n).unwrap();
        let f_q = f.clone().inv(q, n).unwrap();
        let pf_q = f_q.mul_const(Rational::from(p));
        let h = pf_q.cyclic_covolution(g, n).coef_modulo(q).unwrap();
        let m = Polynomial::new_from_ints(vec![-1, 0, 0, 1, -1, 0, 0, 0, -1, 1, 1]);
        let r = Polynomial::new_from_ints(vec![-1, 0, 1, 1, 1, -1, 0, -1]);
        let rh = r.cyclic_covolution(h, n).coef_modulo(q).unwrap();
        let e = (rh + m.clone()).coef_modulo(q).unwrap();
        let a = f
            .cyclic_covolution(e, n)
            .coef_modulo_into_interval(q, -q / 2, q / 2, true)
            .unwrap();
        let b = a
            .coef_modulo_into_interval(p, -p / 2, p / 2, false)
            .unwrap();
        let c = f_p
            .cyclic_covolution(b, n)
            .coef_modulo_into_interval(p, -p / 2, p / 2, false)
            .unwrap();
        assert_eq!(m, c);
    }
}
