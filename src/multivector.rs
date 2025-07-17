use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

use approx_collections::{traits::*, Precision};

use super::{Axes, Scalar, Term};
use crate::ops::grade_project_and_sum_terms;

pub type Coefficients<M> =
    std::iter::Map<<<M as Multivector>::Terms as IntoIterator>::IntoIter, fn(Term) -> Scalar>;

/// Multivector supporting an arbitrary subset of terms.
pub trait Multivector:
    'static
    + fmt::Debug
    + Copy
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Scalar, Output = Self>
    + Div<Scalar, Output = Self>
    + ApproxEq
    + ApproxEqZero
    + ApproxHash
    + ForEachFloat
{
    /// Array type `[Term; N]` of terms in the blade.
    type Terms: Copy + AsRef<[Term]> + IntoIterator<Item = Term>;

    /// Dual blade type.
    type Dual: Multivector;

    /// Returns an array of the terms in the blade.
    fn terms(self) -> Self::Terms;

    /// Returns an iterator over the coefficients in the blade.
    fn coefs(self) -> Coefficients<Self> {
        self.terms().into_iter().map(|t| t.coef)
    }

    /// Returns whether the terms of `self` and `other` have the same basis
    /// vectors.
    ///
    /// This also determines whether [`Multivector::coefs()`] will represent the
    /// same terms for `self` and `other`.
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

    /// Rescales the multivector so that it represents the same object but at a
    /// canonical scale. Use this for comparison where different orientations
    /// should be considered distinct.
    ///
    /// This does **not** necessarily correspond to normalization.
    ///
    /// Uses [`Precision::DEFAULT`].
    #[must_use = "rescale_oriented() returns a mutated copy"]
    fn rescale_oriented(self) -> Self {
        self.rescale_oriented_with_prec(Precision::DEFAULT)
    }
    /// Rescales the multivector so that it represents the same object but at a
    /// canonical scale. Use this for comparison where different orientations
    /// should be considered distinct.
    ///
    /// This does **not** necessarily correspond to normalization.
    #[must_use = "rescale_oriented_with_prec() returns a mutated copy"]
    fn rescale_oriented_with_prec(self, prec: Precision) -> Self {
        match abs_max(self.coefs().filter(|x| !prec.eq_zero(x))) {
            Some(largest_coef) => self * largest_coef.abs().recip(),
            None => Self::zero(),
        }
    }

    /// Rescales the multivector so that it represents the same object but at a
    /// canonical scale and with a canonical orientation. Use this for
    /// comparison that ignores orientation.
    ///
    /// This does **not** necessarily correspond to normalization.
    ///
    /// Uses [`Precision::DEFAULT`].
    #[must_use = "rescale_unoriented() returns a mutated copy"]
    fn rescale_unoriented(self) -> Self {
        self.rescale_unoriented_with_prec(Precision::DEFAULT)
    }
    /// Rescales the multivector so that it represents the same object but at a
    /// canonical scale and with a canonical orientation. Use this for
    /// comparison that ignores orientation.
    ///
    /// This does **not** necessarily correspond to normalization.
    #[must_use = "rescale_unoriented_with_prec() returns a mutated copy"]
    fn rescale_unoriented_with_prec(self, prec: Precision) -> Self {
        let sign = self
            .coefs()
            .find(|x| !prec.eq_zero(x))
            .unwrap_or(0.0)
            .signum();
        match abs_max(self.coefs().filter(|x| !prec.eq_zero(x))) {
            Some(largest_coef) => self * largest_coef.recip() * sign,
            None => Self::zero(),
        }
    }
}

/// Returns the maximum absolute value of all the floats in the list, using the
/// same ordering as [`f64::total_cmp()`].
fn abs_max(floats: impl Iterator<Item = f64>) -> Option<f64> {
    // Positive floats compare like integers. See the comments in the source for
    // `f64::total_cmp()` for an explanation.
    floats.map(|x| x.abs().to_bits()).max().map(f64::from_bits)
}
