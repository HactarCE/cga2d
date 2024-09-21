use std::fmt;
use std::ops::{Mul, Neg, Not};

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

impl Not for Term {
    type Output = Term;

    fn not(self) -> Self::Output {
        self.dual()
    }
}

impl approx::AbsDiffEq for Term {
    type Epsilon = Scalar;

    fn default_epsilon() -> Self::Epsilon {
        Scalar::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.axes == other.axes && approx::AbsDiffEq::abs_diff_eq(&self.coef, &other.coef, epsilon)
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
        Self::pseudoscalar(1.0) * self
    }

    /// Returns the reverse of the term.
    #[must_use]
    pub fn reverse(self) -> Self {
        Self::new(self.axes, self.coef * self.axes.reverse())
    }

    /// Returns the reverse of the term.
    #[must_use]
    pub fn conjugate(self) -> Self {
        Self::new(self.axes, self.coef * self.axes.reverse())
    }
}
