#![allow(missing_docs)]

use std::fmt;
use std::ops::Mul;

use super::Scalar;

bitflags::bitflags! {
    /// Subset of 2D CGA basis vectors.
    #[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
    #[cfg_attr(feature = "bytemuck", derive(bytemuck::Zeroable, bytemuck::NoUninit))]
    #[repr(C)]
    pub struct Axes: u8 {
        /// Scalar
        const S = 0b0000;

        /// e₋
        const M = 0b0001;

        /// e₊
        const P = 0b0010;

        /// x
        const X = 0b0100;

        /// y
        const Y = 0b1000;

        /// e₋e₊
        const MP = 0b0011;
        /// e₋x
        const MX = 0b0101;
        /// e₊x
        const PX = 0b0110;
        /// e₋e₊x
        const MPX = 0b0111;
        /// e₋y
        const MY = 0b1001;
        /// e₊y
        const PY = 0b1010;
        /// e₋e₊y
        const MPY = 0b1011;
        /// xy
        const XY = 0b1100;
        /// e₋xy
        const MXY = 0b1101;
        /// e₊xy
        const PXY = 0b1110;
        /// e₋e₊xy
        const MPXY = 0b1111;
    }
}

impl Axes {
    /// Returns the number of component basis vectors.
    pub fn grade(self) -> u8 {
        self.bits().count_ones() as u8
    }
}

impl fmt::Display for Axes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.contains(Self::M) {
            write!(f, "e₋")?;
        }
        if self.contains(Self::P) {
            write!(f, "e₊")?;
        }
        if self.contains(Self::X) {
            write!(f, "x")?;
        }
        if self.contains(Self::Y) {
            write!(f, "y")?;
        }
        Ok(())
    }
}

/// Sign of the geometric product of two terms, indexed by two `Axes`.
const GEOMETRIC_PRODUCT_SIGN_LUT: [u16; 16] = [
    0b0000000000000000,
    0b1010101010101010,
    0b1010101010101010,
    0b0000000000000000,
    0b0110011001100110,
    0b1100110011001100,
    0b1100110011001100,
    0b0110011001100110,
    0b1001011010010110,
    0b0011110000111100,
    0b0011110000111100,
    0b1001011010010110,
    0b1111000011110000,
    0b0101101001011010,
    0b0101101001011010,
    0b1111000011110000,
];

/// Sign of the reverse of a term, indexed by an `Axes`.
const REVERSE_SIGN_LUT: u16 = 0b0111111011101000;
const CONJUGATE_SIGN_LUT: u16 = 0b0001011101111110;

impl Mul for Axes {
    type Output = Scalar;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul_sign(rhs)
    }
}

impl Axes {
    /// Pseudoscalar axes.
    pub const PSEUDOSCALAR: Self = Self::MPXY;

    /// Sign of the pseudoscalar squared.
    pub const PSEUDOSCALAR_SQUARED: Scalar = Self::PSEUDOSCALAR.mul_sign(Self::PSEUDOSCALAR);

    /// Returns the sign of the reverse of the axes.
    pub const fn rev_sign(self) -> Scalar {
        get_bit_as_sign(REVERSE_SIGN_LUT, self.bits())
    }

    /// Returns the sign of the conjugate of the axes.
    pub const fn conj_sign(self) -> Scalar {
        get_bit_as_sign(CONJUGATE_SIGN_LUT, self.bits())
    }

    const fn mul_sign(self, rhs: Self) -> Scalar {
        get_bit_as_sign(GEOMETRIC_PRODUCT_SIGN_LUT[self.bits() as usize], rhs.bits())
    }
}

const fn get_bit_as_sign(bitmask: u16, index: u8) -> Scalar {
    match bitmask & (1 << index as u16) {
        0 => 1.0,
        _ => -1.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geometric_product_sign() {
        assert_eq!(Axes::X * Axes::X, 1.0);
        assert_eq!(Axes::X * Axes::Y, 1.0);
        assert_eq!(Axes::Y * Axes::X, -1.0);
        assert_eq!(Axes::X * Axes::M, -1.0);
        assert_eq!(Axes::X * Axes::P, -1.0);
        assert_eq!(Axes::P * Axes::P, 1.0);
        assert_eq!(Axes::M * Axes::M, -1.0);
    }

    #[test]
    fn test_axes_grade() {
        assert_eq!(Axes::S.grade(), 0);
        assert_eq!(Axes::X.grade(), 1);
        assert_eq!(Axes::Y.grade(), 1);
        assert_eq!(Axes::M.grade(), 1);
        assert_eq!(Axes::P.grade(), 1);
        assert_eq!(Axes::PY.grade(), 2);
        assert_eq!(Axes::MXY.grade(), 3);
        assert_eq!(Axes::MPXY.grade(), 4);
    }
}
