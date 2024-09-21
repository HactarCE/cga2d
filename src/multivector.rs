use std::{
    fmt,
    ops::{Add, Div, Mul, Sub},
};

use crate::ops::grade_project_and_sum_terms;

use super::{Axes, Scalar, Term};

/// Multivector supporting an arbitrary subset of terms.
pub trait Multivector:
    fmt::Debug
    + Copy
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Scalar, Output = Self>
    + Div<Scalar, Output = Self>
{
    /// Array type `[Term; N]` of terms in the blade.
    type Terms: Copy + AsRef<[Term]> + IntoIterator<Item = Term>;

    /// Dual blade type.
    type Dual: Multivector;

    /// Returns an array of the terms in the blade.
    fn terms(self) -> Self::Terms;

    /// Returns whether the terms of `self` and `other` have the same basis
    /// vectors.
    fn has_same_terms_as(self, other: Self) -> bool;

    /// Returns a coefficient of the blade.
    fn get(&self, axes: Axes) -> Option<&Scalar>;
    /// Returns a mutable reference to a coefficient of the blade.
    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar>;

    /// Returns a multivector equal to zero.
    fn zero() -> Self;

    /// Returns the dual of the blade, which results from right-multiplying the
    /// blade by the pseudoscalar.
    fn dual(self) -> Self::Dual {
        crate::ops::grade_project_and_sum_terms(self.terms().into_iter().map(|t| t.dual()))
    }

    /// Returns the antidual of the blade, which results from left-multiplying
    /// the blade by the pseudoscalar.
    fn antidual(self) -> Self::Dual {
        crate::ops::grade_project_and_sum_terms(self.terms().into_iter().map(|t| t.antidual()))
    }

    /// Returns the scalar dot product of two multivectors.
    fn dot(self, other: Self) -> Scalar {
        crate::ops::multiply_and_grade_project(self, other)
    }

    /// Returns the sum of the squares of all the components in the blade.
    fn mag2(self) -> Scalar {
        self.dot(self)
    }

    /// Returns the absolute value of the magnitude of the blade.
    fn abs_mag(self) -> Scalar {
        let mag2 = self.mag2();
        mag2.abs().sqrt()
    }

    /// Returns the magnitude of the blade. Real magnitudes are represented
    /// using positive numbers and imaginary magnitudes are represented using
    /// negative numbers.
    fn signed_mag(self) -> Scalar {
        let mag2 = self.mag2();
        mag2.abs().sqrt() * mag2.signum()
    }

    /// Normalizes a multivector with respect to its magnitude.
    fn normalize(self) -> Self {
        self / self.abs_mag()
    }

    /// Returns the termwise reverse of the multivector.
    fn rev(self) -> Self {
        crate::ops::grade_project_and_sum_terms(self.terms().into_iter().map(|t| t.rev()))
    }

    /// Returns the inverse of the multivector.
    fn inv(self) -> Self {
        let rev = self.rev();
        rev / self.dot(rev)
    }

    /// Returns the sandwich product of the multivector with `inner`.
    fn sandwich<M: Multivector>(self, inner: M) -> M {
        grade_project_and_sum_terms(
            itertools::iproduct!(
                self.terms().as_ref(),
                inner.terms().as_ref(),
                self.terms().as_ref()
            )
            .map(|(&l, &m, &r)| l * m * r.conj()),
        )
    }
}
