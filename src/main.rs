// This crate is licensed under Apache License 2.0 or CC BY 4.0 at your option.

extern crate num;
use num::{Zero, One, rational::Ratio};

extern crate integer_partitions;
use integer_partitions::Partitions;

extern crate permutohedron;
use permutohedron::{heap_recursive as permutations, factorial};

mod polynomial;
use polynomial::Polynomial;

type Int = i128;
type Rational = Ratio<Int>;
type RP = Polynomial<Rational>;

fn main() {
	for p in 3..=9 {
		let qmax = p*(p-1)/4;
		for q in 0..=qmax {
			println!("p={}, q={}: {}", p, q, t(q, p))
		}
		println!("");
	}
}

fn rl(n: usize, a: &[usize]) -> Vec<Int> {
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
		let mut term = RP::from(vec![(-2 * r[0] - 1).into(), 1.into()]) * Rational::new(1, r[0]+1);
		for (i, &v) in r.iter().enumerate() {
			let mut f = vec![0; q];
			f[i] = 1;
			for j in 1..=v {
				term = term * (RP::monomial(f.clone(), Rational::from(1)) - Rational::from(j-1)) * Rational::new(1, j);
			}
		}
		s = s + term;
	}
	let mut pp = Partitions::new(q);
	while let Some(beta) = pp.next() {
		let r = rl(q, beta);
		if r[0] != 0 {
			continue
		}
		let mut term = RP::one();
		for (i, &v) in r.iter().enumerate() {
			let mut f = vec![0; q];
			f[i] = 1;
			for j in 1..=v {
				term = term * (RP::monomial(f.clone(), Rational::from(1)) - Rational::from(j-1)) * Rational::new(1, j);
			}
		}
		s = s + term;
	}
	s
}

fn t(q: usize, p: usize) -> Int {
// https://oeis.org/A259976
// see the linked Merris & Watkins paper, pp. 539-541
	let xiq = xi(q); // cache?
	let mut r = Rational::zero();
	let mut range: Vec<_> = (0..p).collect();
	permutations(&mut range, |pp| {
		let sigma = cycles(edges(Vec::from(pp)));
		r += xiq.eval(sigma)
	});
	r /= factorial(p) as Int;
	if !r.is_integer() { println!("{}", r) };
	assert!(r.is_integer());
	*r.numer()
}

fn cycles(p: Vec<usize>) -> Vec<Int> {
// counts cycles of each length in a permutation
	let mut r = vec![0; p.len()];
	for i0 in 0..p.len() {
		let mut i = i0;
		let mut j = 0;
		while p[i] > i0 {
			i = p[i];
			j += 1;
		}
		if p[i] == i0 {
			r[j] += 1;
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