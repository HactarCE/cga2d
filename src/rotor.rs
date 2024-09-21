use super::{Axes, Multivector, Scalar, Term};

/// Multivector in the even subalgebra of 2D CGA, used to represent
/// orientation-preserving (i.e., non-reflecting) conformal transformations.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
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
impl Multivector for Rotor {
    type Terms = [Term; 8];

    type Dual = Self;

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
