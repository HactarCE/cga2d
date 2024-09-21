use crate::Pseudoscalar;

use super::{Axes, Blade1, Blade2, Blade3, Multivector, Scalar, Term};

macro_rules! impl_from_terms {
    (($from:ty) -> $to:ty) => {
        impl From<$from> for $to {
            fn from(value: $from) -> Self {
                crate::ops::grade_project_and_sum_terms(value.terms())
            }
        }
    };
}

/// Multivector in either the even or odd subalgebra of 2D CGA, used to
/// represent an arbitrary conformal transformation.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Rotoflector {
    /// Degenerate transform.
    Zero,
    /// Orientation-preserving conformal transformation.
    Rotor(Rotor),
    /// Non-orientation-preserving conformal transformation.
    Flector(Flector),
}
impl Default for Rotoflector {
    fn default() -> Self {
        Self::ident()
    }
}
impl Multivector for Rotoflector {
    type Terms = [Term; 8];

    type Dual = Self;

    fn zero() -> Self {
        Self::Zero
    }

    fn has_same_terms_as(self, other: Self) -> bool {
        std::mem::discriminant(&self) == std::mem::discriminant(&other)
    }

    fn terms(self) -> Self::Terms {
        match self {
            Rotoflector::Zero => [Term::new(Axes::S, 0.0); 8],
            Rotoflector::Rotor(rotor) => rotor.terms(),
            Rotoflector::Flector(flector) => flector.terms(),
        }
    }

    fn get(&self, axes: Axes) -> Option<&Scalar> {
        match self {
            Rotoflector::Zero => None,
            Rotoflector::Rotor(rotor) => rotor.get(axes),
            Rotoflector::Flector(flector) => flector.get(axes),
        }
    }

    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar> {
        match self {
            Rotoflector::Zero => {
                *self = match axes.grade() % 2 {
                    0 => Rotoflector::Rotor(Rotor::default()),
                    _ => Rotoflector::Flector(Flector::default()),
                };
                self.get_mut(axes)
            }
            Rotoflector::Rotor(rotor) => rotor.get_mut(axes),
            Rotoflector::Flector(flector) => flector.get_mut(axes),
        }
    }
}
impl Rotoflector {
    /// Returns the identity transformation.
    pub fn ident() -> Self {
        Self::Rotor(Rotor::ident())
    }
}
impl From<Rotor> for Rotoflector {
    fn from(value: Rotor) -> Self {
        Self::Rotor(value)
    }
}
impl From<Flector> for Rotoflector {
    fn from(value: Flector) -> Self {
        Self::Flector(value)
    }
}
impl_from_terms!((Scalar) -> Rotoflector);
impl_from_terms!((Blade1) -> Rotoflector);
impl_from_terms!((Blade2) -> Rotoflector);
impl_from_terms!((Blade3) -> Rotoflector);
impl_from_terms!((Pseudoscalar) -> Rotoflector);

/// Multivector in the even subalgebra of 2D CGA, used to represent
/// orientation-preserving (i.e., non-reflecting) conformal transformations.
#[derive(Debug, Copy, Clone, PartialEq)]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Zeroable, bytemuck::NoUninit))]
#[repr(C)]
#[allow(missing_docs)]
pub struct Rotor {
    pub s: Scalar,

    pub mp: Scalar,
    pub mx: Scalar,
    pub px: Scalar,
    pub my: Scalar,
    pub py: Scalar,
    pub xy: Scalar,

    pub mpxy: Scalar,
}
impl Default for Rotor {
    fn default() -> Self {
        Self::ident()
    }
}
impl Multivector for Rotor {
    type Terms = [Term; 8];

    type Dual = Self;

    fn zero() -> Self {
        Self {
            s: 0.0,

            mp: 0.0,
            mx: 0.0,
            px: 0.0,
            my: 0.0,
            py: 0.0,
            xy: 0.0,

            mpxy: 0.0,
        }
    }

    fn terms(self) -> Self::Terms {
        [
            Term::new(Axes::S, self.s),
            Term::new(Axes::MP, self.mp),
            Term::new(Axes::MX, self.mx),
            Term::new(Axes::PX, self.px),
            Term::new(Axes::MY, self.my),
            Term::new(Axes::PY, self.py),
            Term::new(Axes::XY, self.xy),
            Term::new(Axes::MPXY, self.mpxy),
        ]
    }

    fn has_same_terms_as(self, _other: Self) -> bool {
        true
    }

    fn get(&self, axes: crate::Axes) -> Option<&Scalar> {
        match axes {
            Axes::S => Some(&self.s),
            Axes::MP => Some(&self.mp),
            Axes::MX => Some(&self.mx),
            Axes::PX => Some(&self.px),
            Axes::MY => Some(&self.my),
            Axes::PY => Some(&self.py),
            Axes::XY => Some(&self.xy),
            Axes::MPXY => Some(&self.mpxy),
            _ => None,
        }
    }

    fn get_mut(&mut self, axes: crate::Axes) -> Option<&mut Scalar> {
        match axes {
            Axes::S => Some(&mut self.s),
            Axes::MP => Some(&mut self.mp),
            Axes::MX => Some(&mut self.mx),
            Axes::PX => Some(&mut self.px),
            Axes::MY => Some(&mut self.my),
            Axes::PY => Some(&mut self.py),
            Axes::XY => Some(&mut self.xy),
            Axes::MPXY => Some(&mut self.mpxy),
            _ => None,
        }
    }
}
impl Rotor {
    /// Returns the identity transformation.
    pub fn ident() -> Self {
        Self {
            s: 1.0,
            ..Default::default()
        }
    }
}
impl_from_terms!((Scalar) -> Rotor);
impl_from_terms!((Blade2) -> Rotor);
impl_from_terms!((Pseudoscalar) -> Rotor);

/// Multivector in the odd subalgebra of 2D CGA, used to represent
/// non-orientation-preserving (i.e., reflecting) conformal transformations.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Zeroable, bytemuck::NoUninit))]
#[repr(C)]
#[allow(missing_docs)]
pub struct Flector {
    pub m: Scalar,
    pub p: Scalar,
    pub x: Scalar,
    pub y: Scalar,

    pub mpx: Scalar,
    pub mpy: Scalar,
    pub mxy: Scalar,
    pub pxy: Scalar,
}
impl Multivector for Flector {
    type Terms = [Term; 8];

    type Dual = Self;

    fn zero() -> Self {
        Self {
            m: 0.0,
            p: 0.0,
            x: 0.0,
            y: 0.0,

            mpx: 0.0,
            mpy: 0.0,
            mxy: 0.0,
            pxy: 0.0,
        }
    }

    fn terms(self) -> Self::Terms {
        [
            Term::new(Axes::M, self.m),
            Term::new(Axes::P, self.p),
            Term::new(Axes::X, self.x),
            Term::new(Axes::Y, self.y),
            Term::new(Axes::MPX, self.mpx),
            Term::new(Axes::MPY, self.mpy),
            Term::new(Axes::MXY, self.mxy),
            Term::new(Axes::PXY, self.pxy),
        ]
    }

    fn has_same_terms_as(self, _other: Self) -> bool {
        true
    }

    fn get(&self, axes: crate::Axes) -> Option<&Scalar> {
        match axes {
            Axes::M => Some(&self.m),
            Axes::P => Some(&self.p),
            Axes::X => Some(&self.x),
            Axes::Y => Some(&self.y),
            Axes::MPX => Some(&self.mpx),
            Axes::MPY => Some(&self.mpy),
            Axes::MXY => Some(&self.mxy),
            Axes::PXY => Some(&self.pxy),
            _ => None,
        }
    }

    fn get_mut(&mut self, axes: crate::Axes) -> Option<&mut Scalar> {
        match axes {
            Axes::M => Some(&mut self.m),
            Axes::P => Some(&mut self.p),
            Axes::X => Some(&mut self.x),
            Axes::Y => Some(&mut self.y),
            Axes::MPX => Some(&mut self.mpx),
            Axes::MPY => Some(&mut self.mpy),
            Axes::MXY => Some(&mut self.mxy),
            Axes::PXY => Some(&mut self.pxy),
            _ => None,
        }
    }
}
impl_from_terms!((Blade1) -> Flector);
impl_from_terms!((Blade3) -> Flector);
