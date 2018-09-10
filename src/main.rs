// This crate is licensed under Apache License 2.0 or CC BY 4.0 at your option.

extern crate num;
use num::{Zero, One, rational::Ratio, BigInt};

extern crate integer_partitions;
use integer_partitions::Partitions;

mod polynomial;
use polynomial::Polynomial;

type Int = BigInt;
type Rational = Ratio<Int>;
type RP = Polynomial<Rational>;

fn main() {
	let mut args = std::env::args().skip(1);
	match args.next() {
		None => {
			for p in 0.. {
				let qmax = if p == 0 { 0 } else { p*(p-1)/4 };
				for q in 0..=qmax {
					println!("p={}, q={}: {}", p, q, t(q, p))
				}
				println!("");
			}
		},
		Some(ref s) if s == "2" => {
			for q in 0.. {
				println!("{}", t(q, 2*q));
			}
		},
		_ => {
			println!("Run without command line arguments for A259976 or with argument '2' for A005368.");
		}
	}
}

fn rl(n: usize, a: &[usize]) -> Vec<usize> {
	let mut v = vec![0; n];
	for &x in a {
		v[x - 1] += 1;
	}
	v
}

fn xi(q: usize) -> RP {
// generates the polynomial which calculates the character of some symmetric group S_m
// corresponding to the partition [m-q, q] for a permutation with given set of cycles
// see F. D. Murnaghan, The theory of group representations, any edition, p. 144
	if q == 0 {
		return RP::one()
	}
	if q == 1 {
		return RP::from(vec![-Rational::one(), Rational::one()])
	}
	let mut s = RP::zero();
	let mut pp = Partitions::new(q - 1);
	while let Some(sigma) = pp.next() {
		let r = rl(q - 1, sigma);
		let mut term = RP::from(vec![Rational::from(-2 * Int::from(r[0]) - 1), Rational::one()]) * Rational::new(Int::one(), Int::from(r[0]+1));
		for (i, &v) in r.iter().enumerate() {
			let mut f = vec![0; q];
			f[i] = 1;
			for j in 1..=usize::from(v) {
				term = term * (RP::monomial(f.clone(), Rational::one()) - Rational::from(Int::from(j-1))) * Rational::new(Int::one(), Int::from(j));
			}
		}
		s = s + term;
	}
	let mut pp = Partitions::new(q);
	while let Some(beta) = pp.next() {
		let r = rl(q, beta);
		if !r[0].is_zero() {
			continue
		}
		let mut term = RP::one();
		for (i, &v) in r.iter().enumerate() {
			let mut f = vec![0; q];
			f[i] = 1;
			for j in 1..=v {
				term = term * (RP::monomial(f.clone(), Rational::one()) - Rational::from(Int::from(j-1))) * Rational::new(Int::one(), Int::from(j));
			}
		}
		s = s + term;
	}
	s
}

fn t(q: usize, p: usize) -> Int {
// https://oeis.org/A259976
// see the linked Merris & Watkins paper, pp. 539-541
	if q == 0 {
		return Int::one();
	}
	let xiq = xi(q); // cache?
	let mut r = Rational::zero();
	let mut ppp = Partitions::new(p);
	while let Some(pp) = ppp.next() {
		let t = // number of partitions with set of cycles pp divided by p!
			Rational::new(Int::one(),
			Int::from(pp.iter().product::<usize>()) *
			rl(p, pp).iter().map(|&x| { (1..=x).map(|n| { Int::from(n) }).product::<Int>() }).product::<Int>());
		let sigma = cycles(edges(from_cycles(pp)));
		r += xiq.eval(sigma) * t;
	}
	if !r.is_integer() { eprintln!("{}", r) };
	assert!(r.is_integer());
	r.to_integer()
}

fn from_cycles(pp: &[usize]) -> Vec<usize> {
// constructs a permutation from cycles
	let mut perm = Vec::new();
	let mut s = 0;
	for &cyc in pp {
		for i in 1..cyc {
			perm.push(s + i)
		}
		perm.push(s);
		s += cyc;
	}
	perm
}

fn cycles(p: Vec<usize>) -> Vec<Int> {
// counts cycles of each length in a permutation
	let mut r = vec![Int::zero(); p.len()];
	for i0 in 0..p.len() {
		let mut i = i0;
		let mut j = 0;
		while p[i] > i0 {
			i = p[i];
			j += 1;
		}
		if p[i] == i0 {
			r[j] += Int::one();
		}
	}
	r
}

fn edges(p: Vec<usize>) -> Vec<usize> {
// turns a permutations of vertices of K_p to a permutation of its edges
	let mut r = Vec::new();
	for (v2, w2) in p.iter().enumerate() {
		for (v1, w1) in p.iter().enumerate() {
			if v1 >= v2 {
				break
			}
			let (w1, w2) = if w1 < w2 { (w1, w2) } else { (w2, w1) };
			r.push(w2*(w2-1)/2+w1);
		}
	}
	r
}

#[cfg(test)]
mod tests {
	use super::*;
	
	#[test]
	fn a259976() {
		assert_eq!((0..=7).map(|q| { t(q, 6) }).collect::<Vec<_>>(), vec![1, 0, 1, 3, 4, 6, 6, 3].into_iter().map(|x| { Int::from(x) }).collect::<Vec<_>>());
		assert_eq!((0..=8).map(|q| { t(q, 7) }).collect::<Vec<_>>(), vec![1, 0, 1, 3, 5, 11, 20, 24, 32].into_iter().map(|x| { Int::from(x) }).collect::<Vec<_>>());
	}
}
