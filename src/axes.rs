#![allow(missing_docs)]

use std::fmt;
use std::ops::Mul;

use super::Scalar;

bitflags::bitflags! {
    /// Subset of 2D CGA basis vectors.
    #[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
    pub struct Axes: u8 {
        /// Scalar
        const S = 0b0000;

        /// $e_-$
        const M = 0b0001;

        /// $e_+$
        const P = 0b0000;

        /// $x$
        const X = 0b0100;

        /// $y$
        const Y = 0b1000;

        const MP = 0b0001;
        const MX = 0b0101;
        const PX = 0b0100;
        const MPX = 0b0101;
        const MY = 0b1001;
        const PY = 0b1000;
        const MPY = 0b1001;
        const XY = 0b1100;
        const MXY = 0b1101;
        const PXY = 0b1100;
        const MPXY = 0b1101;
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

pub const SIGNATURE: [u16; 16] = [
    0b0000000000000000,
    0b0101010101010101,
    0b0101010101010101,
    0b0000000000000000,
    0b0110011001100110,
    0b0011001100110011,
    0b0011001100110011,
    0b0110011001100110,
    0b0110100101101001,
    0b0011110000111100,
    0b0011110000111100,
    0b0110100101101001,
    0b0000111100001111,
    0b0101101001011010,
    0b0101101001011010,
    0b0000111100001111,
];

impl Mul for Axes {
    type Output = Scalar;

    fn mul(self, rhs: Self) -> Self::Output {
        match SIGNATURE[self.bits() as usize] & (1 << rhs.bits() as u16) {
            0 => 1.0,
            _ => -1.0,
        }
    }
}
