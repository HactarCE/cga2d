#![allow(missing_docs)]

use super::{Axes, Multivector, Scalar, Term};

/// Multivector of a compile-time-known grade.
pub trait Blade: Multivector {
    /// Number of basis vectors comprising each term of the blade.
    const GRADE: u8;
}

impl Multivector for Scalar {
    type Terms = [Term; 1];

    type Dual = Pseudoscalar;

    fn terms(self) -> Self::Terms {
        [Term::new(Axes::S, self)]
    }

    fn get(&self, axes: Axes) -> Option<&Scalar> {
        match axes {
            Axes::S => Some(self),
            _ => None,
        }
    }

    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar> {
        match axes {
            Axes::S => Some(self),
            _ => None,
        }
    }
}
impl Blade for Scalar {
    const GRADE: u8 = 0;
}

/// 1-blade, used to represent points, vectors, and round points.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Blade1 {
    pub m: Scalar,
    pub p: Scalar,
    pub x: Scalar,
    pub y: Scalar,
}
impl Multivector for Blade1 {
    type Terms = [Term; 4];

    type Dual = Blade3;

    fn terms(self) -> [Term; 4] {
        [
            Term::new(Axes::M, self.m),
            Term::new(Axes::P, self.p),
            Term::new(Axes::X, self.x),
            Term::new(Axes::Y, self.y),
        ]
    }

    fn get(&self, axes: Axes) -> Option<&Scalar> {
        match axes {
            Axes::M => Some(&self.m),
            Axes::P => Some(&self.p),
            Axes::X => Some(&self.x),
            Axes::Y => Some(&self.y),
            _ => None,
        }
    }

    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar> {
        match axes {
            Axes::M => Some(&mut self.m),
            Axes::P => Some(&mut self.p),
            Axes::X => Some(&mut self.x),
            Axes::Y => Some(&mut self.y),
            _ => None,
        }
    }
}
impl Blade for Blade1 {
    const GRADE: u8 = 1;
}

/// 2-blade, used to represent point pairs (real and imaginary), tangent
/// vectors, and flat points.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Blade2 {
    pub mp: Scalar,
    pub mx: Scalar,
    pub px: Scalar,
    pub my: Scalar,
    pub py: Scalar,
    pub xy: Scalar,
}
impl Multivector for Blade2 {
    type Terms = [Term; 6];

    type Dual = Blade2;

    fn terms(self) -> [Term; 6] {
        [
            Term::new(Axes::MP, self.mp),
            Term::new(Axes::MX, self.mx),
            Term::new(Axes::PX, self.px),
            Term::new(Axes::MY, self.my),
            Term::new(Axes::PY, self.py),
            Term::new(Axes::XY, self.xy),
        ]
    }

    fn get(&self, axes: Axes) -> Option<&Scalar> {
        match axes {
            Axes::MP => Some(&self.mp),
            Axes::MX => Some(&self.mx),
            Axes::PX => Some(&self.px),
            Axes::MY => Some(&self.my),
            Axes::PY => Some(&self.py),
            Axes::XY => Some(&self.xy),
            _ => None,
        }
    }

    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar> {
        match axes {
            Axes::MP => Some(&mut self.mp),
            Axes::MX => Some(&mut self.mx),
            Axes::PX => Some(&mut self.px),
            Axes::MY => Some(&mut self.my),
            Axes::PY => Some(&mut self.py),
            Axes::XY => Some(&mut self.xy),
            _ => None,
        }
    }
}
impl Blade for Blade2 {
    const GRADE: u8 = 2;
}

/// 3-blade, used to represent circles (real and imaginary).
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Blade3 {
    pub mpx: Scalar,
    pub mpy: Scalar,
    pub mxy: Scalar,
    pub pxy: Scalar,
}
impl Multivector for Blade3 {
    type Terms = [Term; 4];

    type Dual = Blade1;

    fn terms(self) -> [Term; 4] {
        [
            Term::new(Axes::MPX, self.mpx),
            Term::new(Axes::MPY, self.mpy),
            Term::new(Axes::MXY, self.mxy),
            Term::new(Axes::PXY, self.pxy),
        ]
    }

    fn get(&self, axes: Axes) -> Option<&Scalar> {
        match axes {
            Axes::MPX => Some(&self.mpx),
            Axes::MPY => Some(&self.mpy),
            Axes::MXY => Some(&self.mxy),
            Axes::PXY => Some(&self.pxy),
            _ => None,
        }
    }

    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar> {
        match axes {
            Axes::MPX => Some(&mut self.mpx),
            Axes::MPY => Some(&mut self.mpy),
            Axes::MXY => Some(&mut self.mxy),
            Axes::PXY => Some(&mut self.pxy),
            _ => None,
        }
    }
}
impl Blade for Blade3 {
    const GRADE: u8 = 3;
}

/// 4-blade, used to represent pseduoscalar quantities.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Pseudoscalar {
    pub mpxy: Scalar,
}
impl Multivector for Pseudoscalar {
    type Terms = [Term; 1];

    type Dual = Scalar;

    fn terms(self) -> Self::Terms {
        [Term::new(Axes::MPXY, self.mpxy)]
    }

    fn get(&self, axes: Axes) -> Option<&Scalar> {
        match axes {
            Axes::MPXY => Some(&self.mpxy),
            _ => None,
        }
    }

    fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar> {
        match axes {
            Axes::MPXY => Some(&mut self.mpxy),
            _ => None,
        }
    }
}
impl Blade for Pseudoscalar {
    const GRADE: u8 = 4;
}
