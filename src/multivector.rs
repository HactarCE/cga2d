use std::fmt;

use super::{Axes, Scalar, Term};

/// Multivector supporting an arbitrary subset of terms.
pub trait Multivector: fmt::Debug + Default + Copy + PartialEq {
    /// Array type `[Term; N]` of terms in the blade.
    type Terms: AsRef<[Term]> + IntoIterator<Item = Term>;

    /// Dual blade type.
    type Dual: Multivector;

    /// Returns an array of the terms in the blade.
    fn terms(self) -> Self::Terms;

    /// Returns a coefficient of the blade.
    fn get(&self, axes: Axes) -> Option<&Scalar>;
    /// Returns a mutable reference to a coefficient of the blade.
    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar>;

    /// Returns the dual of the blade, which results from right-multiplying the
    /// blade by the pseudoscalar.
    fn dual(self) -> Self::Dual {
        let mut ret = Self::Dual::default();
        for term in self.terms() {
            let new_term = term.dual();
            *ret.get_mut(new_term.axes).expect("bad dual") += new_term.coef;
        }
        ret
    }

    /// Returns the antidual of the blade, which results from left-multiplying
    /// the blade by the pseudoscalar.
    fn antidual(self) -> Self::Dual {
        let mut ret = Self::Dual::default();
        for term in self.terms() {
            let new_term = term.antidual();
            *ret.get_mut(new_term.axes).expect("bad dual") += new_term.coef;
        }
        ret
    }
}
