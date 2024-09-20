use std::ops::{Add, AddAssign, BitAnd, BitXor, Index, IndexMut, Not, Shl, Sub, SubAssign};

use super::{Axes, Blade1, Blade2, Blade3, Multivector, Pseudoscalar, Scalar, Term};

macro_rules! impl_index_term {
    ($ty:ty) => {
        impl Index<Axes> for $ty {
            type Output = Scalar;

            fn index(&self, index: Axes) -> &Self::Output {
                self.get(index)
                    .expect(concat!("invalid index for ", stringify!($ty)))
            }
        }

        impl IndexMut<Axes> for $ty {
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

macro_rules! impl_multivector_binary_ops {
    ($lhs:ty, $rhs:ty, $out:ty, $trait:ident, $method:ident, $op:tt) => {
        impl $trait<$rhs> for $lhs {
            type Output = $out;

            fn $method(self, rhs: $rhs) -> Self::Output {
                let mut out = <$out>::default();
                for lhs_term in self.terms() {
                    for rhs_term in rhs.terms() {
                        let new_term = lhs_term $op rhs_term;
                        if let Some(out_coef) = out.get_mut(new_term.axes) {
                            *out_coef += new_term.coef
                        }
                    }
                }
                out
            }
        }
    };

    (($lhs:ty) ^ ($rhs:ty) -> $out:ty) => {
        impl_multivector_binary_ops!($lhs, $rhs, $out, BitXor, bitxor, *);
    };
    (($lhs:ty) * ($rhs:ty) -> $out:ty) => {
        impl_multivector_binary_ops!($lhs, $rhs, $out, Mul, mul, *);
    };
    (($lhs:ty) << ($rhs:ty) -> $out:ty) => {
        impl_multivector_binary_ops!($lhs, $rhs, $out, Shl, shl, *);
    };
    (($lhs:ty) & ($rhs:ty) -> $out:ty) => {
        impl BitAnd<$rhs> for $lhs {
            type Output = $out;

            fn bitand(self, rhs: $rhs) -> Self::Output {
                // In 2D, the sign is flipped compared to the proper antiwedge.
                (self.dual() ^ rhs.dual()).antidual()
            }
        }
    };

    ($(($lhs:ty) $op:tt ($rhs:ty) -> $out:ty);* $(;)?) => {
        $( impl_multivector_binary_ops!(($lhs) $op ($rhs) -> $out); )*
    };
}

impl_multivector_binary_ops!(
    (Blade1) ^ (Blade1) -> Blade2;

    (Blade1) ^ (Blade2) -> Blade3;
    (Blade2) ^ (Blade1) -> Blade3;

    (Blade1) ^ (Blade3) -> Pseudoscalar;
    (Blade2) ^ (Blade2) -> Pseudoscalar;
    (Blade3) ^ (Blade1) -> Pseudoscalar;

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

    (Blade3) & (Blade3) -> Blade2;

    (Blade3) & (Blade2) -> Blade1;
    (Blade2) & (Blade3) -> Blade1;

    (Blade3) & (Blade1) -> Scalar;
    (Blade2) & (Blade2) -> Scalar;
    (Blade1) & (Blade3) -> Scalar;
);

macro_rules! impl_multivector_dual {
    ($ty:ty) => {
        impl Not for $ty {
            type Output = <$ty as Multivector>::Dual;

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
