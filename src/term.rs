use std::fmt;
use std::hash::Hash;
use std::ops::{Mul, Neg, Not};

use approx_collections::{ApproxEq, ApproxEqZero, ApproxHash, Precision};

use super::{Axes, Scalar};

/// Term in 2D CGA.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Term {
    /// Basis vectors.
    pub axes: Axes,
    /// Scalar coefficient.
    pub coef: Scalar,
}

impl fmt::Display for Term {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let &Self { axes, coef } = self;
        if axes == Axes::S {
            write!(f, "{coef}")
        } else {
            write!(f, "{coef}{axes}")
        }
    }
}

impl Neg for Term {
    type Output = Term;

    fn neg(self) -> Self::Output {
        let Self { axes, coef } = self;
        Term { axes, coef: -coef }
    }
}

impl Mul for Term {
    type Output = Term;

    #[inline(always)]
    fn mul(self, rhs: Self) -> Self::Output {
        let sign = self.axes * rhs.axes;
        Term {
            axes: self.axes ^ rhs.axes,
            coef: self.coef * rhs.coef * sign,
        }
    }
}

impl Mul<Scalar> for Term {
    type Output = Term;

    fn mul(self, rhs: Scalar) -> Self::Output {
        Term {
            axes: self.axes,
            coef: self.coef * rhs,
        }
    }
}

impl Not for Term {
    type Output = Term;

    fn not(self) -> Self::Output {
        self.dual()
    }
}

impl ApproxEq for Term {
    fn approx_eq(&self, other: &Self, prec: Precision) -> bool {
        self.axes == other.axes && self.coef.approx_eq(&other.coef, prec)
    }
}

impl ApproxEqZero for Term {
    fn approx_eq_zero(&self, prec: Precision) -> bool {
        self.coef.approx_eq_zero(prec)
    }
}

impl ApproxHash for Term {
    fn intern_floats<F: FnMut(&mut f64)>(&mut self, f: &mut F) {
        self.coef.intern_floats(f);
    }

    fn interned_eq(&self, other: &Self) -> bool {
        self.coef.interned_eq(&other.coef)
    }

    fn interned_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.axes.hash(state);
        self.coef.interned_hash(state);
    }
}

impl Term {
    /// Constructs a term.
    #[inline(always)]
    pub fn new(axes: Axes, coef: Scalar) -> Self {
        Self { axes, coef }
    }

    /// Constructs a scalar term.
    pub fn scalar(coef: Scalar) -> Self {
        Self::new(Axes::S, coef)
    }

    /// Constructs a pseudoscalar term.
    pub fn pseudoscalar(coef: Scalar) -> Self {
        Self::new(Axes::MPXY, coef)
    }

    /// Returns the dual of the term.
    #[must_use]
    pub fn dual(self) -> Self {
        self * Self::pseudoscalar(1.0)
    }

    /// Returns the antidual of the term.
    #[must_use]
    pub fn antidual(self) -> Self {
        self.dual() * Axes::PSEUDOSCALAR_SQUARED
    }

    /// Returns the reverse of the term.
    #[must_use]
    pub fn rev(self) -> Self {
        Self::new(self.axes, self.coef * self.axes.rev_sign())
    }

    /// Returns the conjugate of the term.
    #[must_use]
    pub fn conj(self) -> Self {
        Self::new(self.axes, self.coef * self.axes.conj_sign())
    }
}
