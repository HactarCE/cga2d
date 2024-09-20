use std::ops::Mul;

use super::{Axes, Scalar};

/// Term in 2D CGA.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Term {
    /// Basis vectors.
    pub axes: Axes,
    /// Scalar coefficient.
    pub coef: Scalar,
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
