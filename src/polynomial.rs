use num::{Zero, One};
use std::{ops, collections::HashMap};

#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<T> {
	nvars: usize,
	factors: HashMap<Vec<usize>, T>,
}

impl<T> Polynomial<T> {
	pub fn _into_hashmap(self) -> HashMap<Vec<usize>, T> {
		self.factors
	}
	
	pub fn _as_hashmap(&self) -> &HashMap<Vec<usize>, T> {
		&self.factors
	}
	
	pub fn _nvars(&self) -> usize {
		self.nvars
	}
}

impl<T> Zero for Polynomial<T> where T: Zero + Clone + ops::AddAssign {
	fn zero() -> Self {
		Polynomial {
			nvars: 0,
			factors: HashMap::new(),
		}
	}
	
	fn is_zero(&self) -> bool {
		self.factors.iter().all(|(_, x)| { x.is_zero() })
	}
}

impl<T> One for Polynomial<T> where T: Zero + One + Clone + std::iter::Sum + PartialEq + ops::AddAssign {
	fn one() -> Self {
		Polynomial::monomial(Vec::new(), T::one())
	}
	
	fn is_one(&self) -> bool {
		self.factors.iter().all(|(v, x)| { x.is_zero() || x.is_one() && v.iter().all(|x| {x.is_zero()}) })
	}
}

impl<T> Polynomial<T> where T: Zero + Clone + ops::AddAssign {
	fn cleanup(mut self) -> Self {
		self.factors = self.factors.into_iter().filter(|(_, x)| { !x.is_zero() }).collect();
		self
	}

	pub fn monomial(v: Vec<usize>, factor: T) -> Self {
		if factor.is_zero() {
			return Self::zero()
		}
		let nvars = v.len();
		let mut factors = HashMap::new();
		factors.insert(v, factor);
		Polynomial {
			nvars,
			factors,
		}
	}
		
	fn extend(&mut self, l: usize) {
		let t = self.nvars;
		if l <= t {
			return
		}
		let tf = self.factors.drain().map(|(mut v, x)| { for _ in t..l { v.push(0); } (v, x) }).collect();
		self.factors = tf;
		self.nvars = l;
	}
}
	
impl<T> Polynomial<T> where T: Zero + Clone + One + ops::AddAssign {
	pub fn eval<X>(&self, xx: Vec<X>) -> T where X: Into<T> {
		let mut x: Vec<T> = xx.into_iter().map(|t| { t.into() }).collect();
		for _ in x.len()..self.nvars {
			x.push(T::zero());
		}
		let mut sum = T::zero();
		for (v, factor) in &self.factors {
			let mut m = factor.clone();
			for (t, &p) in x.iter().zip(v.iter()) {
				for _ in 0..p {
					m = m * t.clone();
				}
			}
			sum = sum + m;
		}
		sum
	}
}

impl<T> From<Vec<T>> for Polynomial<T> where T: Zero + Clone + ops::AddAssign {
	fn from(factors: Vec<T>) -> Self {
		Polynomial {
			nvars: 1,
			factors: factors.into_iter().enumerate().map(|(n, x)| { (vec!(n), x) }).collect()
		}
	}
}

impl<T> From<HashMap<Vec<usize>, T>> for Polynomial<T> where T: Zero + Clone + ops::AddAssign {
	fn from(mut factors: HashMap<Vec<usize>, T>) -> Self {
		if factors.is_empty() {
			return Polynomial::zero()
		}
		let nvars = factors.iter().map(|(v, _)| {v.len()}).max().unwrap();
		if nvars == factors.iter().map(|(v, _)| {v.len()}).min().unwrap() {
			return Polynomial {
				nvars,
				factors,
			}.cleanup()
		}
		let factors = factors.drain().map( |(mut v, f)| {
			let t = v.len();
			for _ in t..nvars {
				v.push(0)
			}
			(v, f)
		}).collect();
		Polynomial {
			nvars,
			factors,
		}.cleanup()
	}
}

impl<T> ops::Add for Polynomial<T> where T: Zero + Clone + ops::AddAssign {
	type Output = Polynomial<T>;
	fn add(mut self, mut rhs: Self) -> Self {
		self.extend(rhs.nvars);
		rhs.extend(self.nvars);
		assert!(self.nvars == rhs.nvars);
		for (k, v) in rhs.factors {
			if let Some(vv) = self.factors.get_mut(&k) {
				*vv += v;
				continue
			}
			self.factors.insert(k, v);
		}
		self.cleanup()
	}
}

impl<T> ops::Sub for Polynomial<T> where T:
	Zero
	+ ops::Sub<Output=T>
	+ ops::Neg<Output=T>
	+ Clone
	+ ops::AddAssign
	+ ops::SubAssign
{
	type Output = Polynomial<T>;
	fn sub(mut self, mut rhs: Self) -> Self {
		self.extend(rhs.nvars);
		rhs.extend(self.nvars);
		for (k, v) in rhs.factors {
			if let Some(vv) = self.factors.get_mut(&k) {
				*vv -= v;
				continue
			}
			self.factors.insert(k, -v);
		}
		self.cleanup()
	}
}

impl<T> ops::Neg for Polynomial<T> where T: ops::Neg<Output=T> + ops::AddAssign {
	type Output = Self;
	fn neg(mut self) -> Self {
		self.factors = self.factors.into_iter().map(|(x, v)| { (x, -v) }).collect();
		self
	}
}

impl<T> ops::Add<T> for Polynomial<T> where T: Zero + Clone + ops::AddAssign {
	type Output = Polynomial<T>;
	fn add(self, rhs: T) -> Self {
		self + Self::monomial(Vec::new(), rhs)
	}
}

impl<T> ops::Sub<T> for Polynomial<T> where T: Zero + ops::Neg<Output=T> + Clone + ops::AddAssign {
	type Output = Polynomial<T>;
	fn sub(self, rhs: T) -> Self {
		self + Self::monomial(Vec::new(), -rhs)
	}
}

impl<T> ops::Mul<T> for Polynomial<T> where T: Zero + ops::Mul<Output=T> + Clone + ops::AddAssign {
	type Output = Polynomial<T>;
	fn mul(mut self, rhs: T) -> Self {
		self.factors = self.factors.into_iter().map(|(x, v)| { (x, v * rhs.clone()) }).collect();
		self.cleanup()
	}
}

impl<T> ops::Mul for Polynomial<T> where T: Zero + ops::Mul<Output=T> + Clone + std::iter::Sum + ops::AddAssign {
	type Output = Polynomial<T>;
	fn mul(mut self, mut rhs: Self) -> Self {
		self.extend(rhs.nvars);
		rhs.extend(self.nvars);
		let mut res = Self::zero();
		for (k, v) in &self.factors {
			for (kk, vv) in &rhs.factors {
				let k2 = k.iter().zip(kk.iter()).map(|(x, y)| { x + y }).collect::<Vec<_>>();
				res = res + Self::monomial(k2, v.clone() * vv.clone());
			}
		}
		res.cleanup()
	}
}

#[cfg(test)]
mod tests {
	use super::*;
		
	#[test]
	fn eval() {
		let p = Polynomial::from(vec![4.0, 5.0, 6.0]);
		let x = 2;
		assert_eq!(p.eval(vec!(x)), 38.0);
		let x = 2.0;
		assert_eq!(p.eval(vec!(x)), 38.0);
	}
	
	#[test]
	fn extend() {
		let mut p = Polynomial::monomial(vec![], 2);
		p.extend(1);
		assert_eq!(p, Polynomial::monomial(vec![0], 2));
	}
	
	#[test]
	fn add() {
		let p1 = Polynomial::from(vec![3, -5, 2]);
		let p2 = Polynomial::from(vec![7, 3, -1, 3]);
		let s = Polynomial::from(vec![10, -2, 1, 3]);
		assert_eq!(p1.clone() + p2.clone(), s);
		assert_eq!(p2 + p1, s);
	}
	
	#[test]
	fn mul() {
		let p1 = Polynomial::from(vec![3, -5, 2]);
		let p2 = Polynomial::from(vec![7, 3, -1, 3]);
		let p = Polynomial::from(vec![21, -26, -4, 20, -17, 6]);
		assert_eq!(p1.clone() * p2.clone(), p);
		assert_eq!(p2 * p1, p);
	}
	
	#[test]
	fn pm_monom() {
		let p = Polynomial::from(vec![2, 3, -1]);
		assert_eq!(p.clone()+6, Polynomial::from(vec![8, 3, -1]));
		assert_eq!(p.clone()-3, Polynomial::from(vec![-1, 3, -1]));
	}
}
