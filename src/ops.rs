#![allow(clippy::suspicious_arithmetic_impl)]

use std::fmt;
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{
    Add, AddAssign, BitAnd, BitXor, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Not, Shl,
    Sub, SubAssign,
};

use approx_collections::{ApproxEq, ApproxEqZero, ApproxHash, Precision};

use super::{
    Axes, Blade, Blade1, Blade2, Blade3, Flector, Multivector, Pseudoscalar, Rotoflector, Rotor,
    Scalar, Term,
};

/// Wedge operation between two blades.
pub trait Wedge<Rhs: Blade>: Blade {
    /// Output blade type.
    type Output: Blade;

    /// Returns the wedge of two blades.
    fn wedge(self, rhs: Rhs) -> <Self as Wedge<Rhs>>::Output {
        multiply_and_grade_project(self, rhs)
    }
}

macro_rules! impl_index_term {
    ($type:ty) => {
        impl Index<Axes> for $type {
            type Output = Scalar;

            fn index(&self, index: Axes) -> &Self::Output {
                self.get(index)
                    .expect(concat!("invalid index for ", stringify!($ty)))
            }
        }

        impl IndexMut<Axes> for $type {
            fn index_mut(&mut self, index: Axes) -> &mut Self::Output {
                self.get_mut(index)
                    .expect(concat!("invalid index for ", stringify!($ty)))
            }
        }
    };
}

impl_index_term!(Blade1);
impl_index_term!(Blade2);
impl_index_term!(Blade3);
impl_index_term!(Pseudoscalar);
impl_index_term!(Rotor);
impl_index_term!(Flector);
impl_index_term!(Rotoflector);

macro_rules! impl_display_terms {
    ($type:ty) => {
        impl fmt::Display for $type {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                let mut is_first = true;
                for term in self.terms() {
                    if !std::mem::take(&mut is_first) {
                        write!(f, " + ")?;
                    }
                    write!(f, "{term}")?;
                }
                Ok(())
            }
        }
    };
}

impl_display_terms!(Blade1);
impl_display_terms!(Blade2);
impl_display_terms!(Blade3);
impl_display_terms!(Pseudoscalar);
impl_display_terms!(Rotor);
impl_display_terms!(Flector);
impl_display_terms!(Rotoflector);

macro_rules! impl_neg {
    ($ty:ty) => {
        impl Neg for $ty {
            type Output = Self;

            fn neg(self) -> Self::Output {
                grade_project_and_sum_terms(self.terms().map(|t| -t))
            }
        }
    };
}

impl_neg!(Blade1);
impl_neg!(Blade2);
impl_neg!(Blade3);
impl_neg!(Pseudoscalar);
impl_neg!(Rotor);
impl_neg!(Flector);
impl_neg!(Rotoflector);

macro_rules! impl_add_sub_term_ops {
    ($type:ty) => {
        impl<T> Add<T> for $type
        where
            $type: AddAssign<T>,
        {
            type Output = $type;

            fn add(mut self, rhs: T) -> Self::Output {
                self += rhs;
                self
            }
        }
        impl<T> Sub<T> for $type
        where
            $type: SubAssign<T>,
        {
            type Output = $type;

            fn sub(mut self, rhs: T) -> Self::Output {
                self -= rhs;
                self
            }
        }

        impl AddAssign<Term> for $type {
            fn add_assign(&mut self, rhs: Term) {
                self[rhs.axes] += rhs.coef;
            }
        }

        impl SubAssign<Term> for $type {
            fn sub_assign(&mut self, rhs: Term) {
                self[rhs.axes] -= rhs.coef;
            }
        }

        impl AddAssign for $type {
            fn add_assign(&mut self, rhs: Self) {
                for term in rhs.terms() {
                    *self += term;
                }
            }
        }

        impl SubAssign for $type {
            fn sub_assign(&mut self, rhs: Self) {
                for term in rhs.terms() {
                    *self -= term;
                }
            }
        }
    };
}

impl_add_sub_term_ops!(Blade1);
impl_add_sub_term_ops!(Blade2);
impl_add_sub_term_ops!(Blade3);
impl_add_sub_term_ops!(Pseudoscalar);
impl_add_sub_term_ops!(Rotor);
impl_add_sub_term_ops!(Flector);
impl_add_sub_term_ops!(Rotoflector);

impl<T: Add<Term>> Add<T> for Term {
    type Output = T::Output;

    fn add(self, rhs: T) -> Self::Output {
        rhs + self
    }
}
impl<T: Neg<Output = T> + Add<Term>> Sub<T> for Term {
    type Output = <T as Add<Term>>::Output;

    fn sub(self, rhs: T) -> Self::Output {
        -rhs + self
    }
}

/// Implements `MulAssign` and `DivAssign` (which isn't possible for most of our
/// multiplication & division operations).
macro_rules! impl_mul_div_scalar_ops {
    ($type:ty) => {
        impl<T> Mul<T> for $type
        where
            $type: MulAssign<T>,
        {
            type Output = $type;

            fn mul(mut self, rhs: T) -> Self::Output {
                self *= rhs;
                self
            }
        }
        impl<T> Div<T> for $type
        where
            $type: DivAssign<T>,
        {
            type Output = $type;

            fn div(mut self, rhs: T) -> Self::Output {
                self /= rhs;
                self
            }
        }

        impl MulAssign<Scalar> for $type {
            fn mul_assign(&mut self, rhs: Scalar) {
                for term in self.terms() {
                    self[term.axes] *= rhs;
                }
            }
        }

        impl DivAssign<Scalar> for $type {
            fn div_assign(&mut self, rhs: Scalar) {
                for term in self.terms() {
                    self[term.axes] /= rhs;
                }
            }
        }
    };
}

impl_mul_div_scalar_ops!(Blade1);
impl_mul_div_scalar_ops!(Blade2);
impl_mul_div_scalar_ops!(Blade3);
impl_mul_div_scalar_ops!(Pseudoscalar);
impl_mul_div_scalar_ops!(Rotor);
impl_mul_div_scalar_ops!(Flector);
impl_mul_div_scalar_ops!(Rotoflector);

macro_rules! impl_blade_wedge {
    ($lhs:ty, $rhs:ty, $out:ty) => {
        impl Wedge<$rhs> for $lhs {
            type Output = $out;
        }
    };
}

impl_blade_wedge!(Scalar, Scalar, Scalar);
impl_blade_wedge!(Scalar, Blade1, Blade1);
impl_blade_wedge!(Scalar, Blade2, Blade2);
impl_blade_wedge!(Scalar, Blade3, Blade3);
impl_blade_wedge!(Scalar, Pseudoscalar, Pseudoscalar);

impl_blade_wedge!(Blade1, Scalar, Blade1);
impl_blade_wedge!(Blade1, Blade1, Blade2);
impl_blade_wedge!(Blade1, Blade2, Blade3);
impl_blade_wedge!(Blade1, Blade3, Pseudoscalar);

impl_blade_wedge!(Blade2, Scalar, Blade2);
impl_blade_wedge!(Blade2, Blade1, Blade3);
impl_blade_wedge!(Blade2, Blade2, Pseudoscalar);

impl_blade_wedge!(Blade3, Scalar, Blade3);
impl_blade_wedge!(Blade3, Blade1, Pseudoscalar);

impl_blade_wedge!(Pseudoscalar, Scalar, Pseudoscalar);

macro_rules! impl_multivector_binary_ops {
    ($lhs:ty, $rhs:ty, $out:ty, $trait:ident, $method:ident, $op:tt) => {
        impl $trait<$rhs> for $lhs {
            type Output = $out;

            fn $method(self, rhs: $rhs) -> Self::Output {
                multiply_and_grade_project(self, rhs)
            }
        }
    };

    (($lhs:ty) ^ ($rhs:ty) -> $out:ty) => {
        impl_multivector_binary_ops!($lhs, $rhs, $out, BitXor, bitxor, *);
    };
    (($lhs:ty) * ($rhs:ty) -> $out:ty) => {
        impl_multivector_binary_ops!($lhs, $rhs, $out, Mul, mul, *);
        impl_multivector_binary_ops!(($lhs) / ($rhs) -> $out); // also implement division
    };
    (($lhs:ty) << ($rhs:ty) -> $out:ty) => {
        impl_multivector_binary_ops!($lhs, $rhs, $out, Shl, shl, *);
    };
    (($lhs:ty) & ($rhs:ty) -> $out:ty) => {
        impl BitAnd<$rhs> for $lhs {
            type Output = $out;

            fn bitand(self, rhs: $rhs) -> Self::Output {
                (self.dual() ^ rhs.dual()).antidual()
            }
        }
    };
    (($lhs:ty) / ($rhs:ty) -> $out:ty) => {
        impl Div<$rhs> for $lhs {
            type Output = $out;

            fn div(self, rhs: $rhs) -> Self::Output {
                self * rhs.inv()
            }
        }
    };

    (($lhs:ty) $op:tt ($rhs:ty) -> $out:ty) => {
        compile_error!(concat!("unsupported operation: ", stringify!($op)));
    };

    ($(($lhs:ty) $op:tt ($rhs:ty) -> $out:ty);* $(;)?) => {
        $( impl_multivector_binary_ops!(($lhs) $op ($rhs) -> $out); )*
    };
}

impl_multivector_binary_ops!(
    (Scalar) * (Blade1) -> Blade1;
    (Scalar) * (Blade2) -> Blade2;
    (Scalar) * (Blade3) -> Blade3;
    (Scalar) * (Pseudoscalar) -> Pseudoscalar;

    (Blade1) * (Blade1) -> Rotor;
    (Blade1) * (Blade2) -> Flector;
    (Blade1) * (Blade3) -> Rotor;
    (Blade1) * (Pseudoscalar) -> Blade3;

    (Blade2) * (Blade1) -> Flector;
    (Blade2) * (Blade2) -> Rotor;
    (Blade2) * (Blade3) -> Flector;
    (Blade2) * (Pseudoscalar) -> Blade2;

    (Blade3) * (Blade1) -> Rotor;
    (Blade3) * (Blade2) -> Flector;
    (Blade3) * (Blade3) -> Rotor;
    (Blade3) * (Pseudoscalar) -> Blade1;

    (Pseudoscalar) * (Blade1) -> Flector;
    (Pseudoscalar) * (Blade2) -> Rotor;
    (Pseudoscalar) * (Blade3) -> Flector;
    (Pseudoscalar) * (Pseudoscalar) -> Scalar;

    (Scalar) * (Rotor) -> Rotor;
    (Blade1) * (Rotor) -> Flector;
    (Blade2) * (Rotor) -> Rotor;
    (Blade3) * (Rotor) -> Flector;

    (Scalar) * (Flector) -> Flector;
    (Blade1) * (Flector) -> Rotor;
    (Blade2) * (Flector) -> Flector;
    (Blade3) * (Flector) -> Rotor;

    (Pseudoscalar) * (Rotor) -> Rotor;
    (Pseudoscalar) * (Flector) -> Flector;

    (Rotor) * (Blade1) -> Flector;
    (Rotor) * (Blade2) -> Rotor;
    (Rotor) * (Blade3) -> Flector;
    (Rotor) * (Pseudoscalar) -> Rotor;
    (Rotor) * (Rotor) -> Rotor;
    (Rotor) * (Flector) -> Flector;

    (Flector) * (Blade1) -> Rotor;
    (Flector) * (Blade2) -> Flector;
    (Flector) * (Blade3) -> Rotor;
    (Flector) * (Pseudoscalar) -> Flector;
    (Flector) * (Rotor) -> Flector;
    (Flector) * (Flector) -> Rotor;

    (Scalar) * (Rotoflector) -> Rotoflector;
    (Blade1) * (Rotoflector) -> Rotoflector;
    (Blade2) * (Rotoflector) -> Rotoflector;
    (Blade3) * (Rotoflector) -> Rotoflector;
    (Pseudoscalar) * (Rotoflector) -> Rotoflector;
    (Rotor) * (Rotoflector) -> Rotoflector;
    (Flector) * (Rotoflector) -> Rotoflector;

    (Rotoflector) * (Blade1) -> Rotoflector;
    (Rotoflector) * (Blade2) -> Rotoflector;
    (Rotoflector) * (Blade3) -> Rotoflector;
    (Rotoflector) * (Pseudoscalar) -> Rotoflector;
    (Rotoflector) * (Rotor) -> Rotoflector;
    (Rotoflector) * (Flector) -> Rotoflector;

    (Rotoflector) * (Rotoflector) -> Rotoflector;

    (Scalar) ^ (Blade1) -> Blade1;
    (Blade1) ^ (Scalar) -> Blade1;

    (Scalar) ^ (Blade2) -> Blade2;
    (Blade1) ^ (Blade1) -> Blade2;
    (Blade2) ^ (Scalar) -> Blade2;

    (Scalar) ^ (Blade3) -> Blade3;
    (Blade1) ^ (Blade2) -> Blade3;
    (Blade2) ^ (Blade1) -> Blade3;
    (Blade3) ^ (Scalar) -> Blade3;

    (Scalar) ^ (Pseudoscalar) -> Pseudoscalar;
    (Blade1) ^ (Blade3) -> Pseudoscalar;
    (Blade2) ^ (Blade2) -> Pseudoscalar;
    (Blade3) ^ (Blade1) -> Pseudoscalar;
    (Pseudoscalar) ^ (Scalar) -> Pseudoscalar;

    (Pseudoscalar) & (Blade3) -> Blade3;
    (Blade3) & (Pseudoscalar) -> Blade3;

    (Pseudoscalar) & (Blade2) -> Blade2;
    (Blade3) & (Blade3) -> Blade2;
    (Blade2) & (Pseudoscalar) -> Blade2;

    (Pseudoscalar) & (Blade1) -> Blade1;
    (Blade3) & (Blade2) -> Blade1;
    (Blade2) & (Blade3) -> Blade1;
    (Blade1) & (Pseudoscalar) -> Blade1;

    (Pseudoscalar) & (Scalar) -> Scalar;
    (Blade3) & (Blade1) -> Scalar;
    (Blade2) & (Blade2) -> Scalar;
    (Blade1) & (Blade3) -> Scalar;
    (Scalar) & (Pseudoscalar) -> Scalar;

    (Blade1) << (Blade1) -> Scalar;
    (Blade1) << (Blade2) -> Blade1;
    (Blade1) << (Blade3) -> Blade2;
    (Blade1) << (Pseudoscalar) -> Blade3;

    (Blade2) << (Blade2) -> Scalar;
    (Blade2) << (Blade3) -> Blade1;
    (Blade2) << (Pseudoscalar) -> Blade2;

    (Blade3) << (Blade3) -> Scalar;
    (Blade3) << (Pseudoscalar) -> Blade1;

    (Pseudoscalar) << (Pseudoscalar) -> Blade1;
);

macro_rules! impl_multivector_dual {
    ($type:ty) => {
        impl Not for $type {
            type Output = <$type as Multivector>::Dual;

            fn not(self) -> Self::Output {
                self.dual()
            }
        }
    };
}

impl_multivector_dual!(Blade1);
impl_multivector_dual!(Blade2);
impl_multivector_dual!(Blade3);
impl_multivector_dual!(Pseudoscalar);
impl_multivector_dual!(Rotor);
impl_multivector_dual!(Flector);
impl_multivector_dual!(Rotoflector);

pub(crate) fn multiply_and_grade_project<L, R, O>(lhs: L, rhs: R) -> O
where
    L: Multivector,
    R: Multivector,
    O: Multivector,
{
    grade_project_and_sum_terms(
        itertools::iproduct!(lhs.terms().as_ref(), rhs.terms().as_ref()).map(|(&l, &r)| l * r),
    )
}

macro_rules! impl_sum_terms {
    ($type:ty) => {
        impl Sum<Term> for $type {
            fn sum<I: Iterator<Item = Term>>(iter: I) -> Self {
                grade_project_and_sum_terms(iter)
            }
        }
    };
}

impl_sum_terms!(Scalar);
impl_sum_terms!(Blade1);
impl_sum_terms!(Blade2);
impl_sum_terms!(Blade3);
impl_sum_terms!(Pseudoscalar);
impl_sum_terms!(Rotor);
impl_sum_terms!(Flector);
impl_sum_terms!(Rotoflector);

pub(crate) fn grade_project_and_sum_terms<M: Multivector>(
    terms: impl IntoIterator<Item = Term>,
) -> M {
    let mut ret = M::zero();
    for term in terms {
        if let Some(ret_coef) = ret.get_mut(term.axes) {
            *ret_coef += term.coef;
        }
    }
    ret
}

macro_rules! impl_approx_eq {
    ($type:ty) => {
        impl ApproxEq for $type {
            fn approx_eq(&self, other: &Self, prec: Precision) -> bool {
                self.has_same_terms_as(*other)
                    && std::iter::zip(self.terms(), other.terms())
                        .all(|(a, b)| prec.eq(a.coef, b.coef))
            }
        }

        impl ApproxEqZero for $type {
            fn approx_eq_zero(&self, prec: Precision) -> bool {
                self.terms().iter().all(|term| prec.eq_zero(term.coef))
            }
        }
    };
}

impl_approx_eq!(Blade1);
impl_approx_eq!(Blade2);
impl_approx_eq!(Blade3);
impl_approx_eq!(Pseudoscalar);
impl_approx_eq!(Rotor);
impl_approx_eq!(Flector);
impl_approx_eq!(Rotoflector);

macro_rules! impl_approx_hash {
    ($type:ty; $($field:ident),* $(,)?) => {
        impl ApproxHash for $type {
            fn intern_floats<F: FnMut(&mut f64)>(&mut self, f: &mut F) {
                // unpack to throw compile error if any fields are missing
                let Self { $($field),* } = self;
                $( f($field); )*
            }

            fn interned_eq(&self, other: &Self) -> bool {
                $( self.$field.to_bits() == other.$field.to_bits() )&&*
            }

            fn interned_hash<H: std::hash::Hasher>(&self, state: &mut H) {
                $( self.$field.to_bits().hash(state); )*
            }
        }
    };
}

impl_approx_hash!(Blade1; m, p, x, y);
impl_approx_hash!(Blade2; mp, mx, px, my, py, xy);
impl_approx_hash!(Blade3; mpx, mpy, mxy, pxy);
impl_approx_hash!(Pseudoscalar; mpxy);
impl_approx_hash!(Rotor; s, mp, mx, px, my, py, xy, mpxy);
impl_approx_hash!(Flector; m, p, x, y, mpx, mpy, mxy, pxy);
impl ApproxHash for Rotoflector {
    fn intern_floats<F: FnMut(&mut f64)>(&mut self, f: &mut F) {
        match self {
            Rotoflector::Zero => (),
            Rotoflector::Rotor(rotor) => {
                rotor.intern_floats(f);
                if (*rotor).interned_eq(&Rotor::zero()) {
                    *self = Rotoflector::Zero;
                }
            }
            Rotoflector::Flector(flector) => {
                flector.intern_floats(f);
                if (*flector).interned_eq(&Flector::zero()) {
                    *self = Rotoflector::Zero
                }
            }
        }
    }

    fn interned_eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Rotoflector::Zero, Rotoflector::Zero) => true,
            (Rotoflector::Rotor(r1), Rotoflector::Rotor(r2)) => r1.interned_eq(r2),
            (Rotoflector::Flector(f1), Rotoflector::Flector(f2)) => f1.interned_eq(f2),

            (Rotoflector::Zero, _) => false,
            (Rotoflector::Rotor(_), _) => false,
            (Rotoflector::Flector(_), _) => false,
        }
    }

    fn interned_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        std::mem::discriminant(self).hash(state);
        match self {
            Rotoflector::Zero => (),
            Rotoflector::Rotor(rotor) => rotor.interned_hash(state),
            Rotoflector::Flector(flector) => flector.interned_hash(state),
        }
    }
}
