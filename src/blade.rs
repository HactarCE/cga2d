use super::{Axes, Multivector, Scalar, Term, Wedge};

/// ∞, also called nᵢ, representing the point at infinity.
pub const NI: Blade1 = Blade1 {
    m: 1.0,
    p: 1.0,
    x: 0.0,
    y: 0.0,
};

/// nₒ, representing the origin.
pub const NO: Blade1 = Blade1 {
    m: 0.5,
    p: -0.5,
    x: 0.0,
    y: 0.0,
};
/// Multivector of a compile-time-known grade.
pub trait Blade: Multivector {
    /// Number of basis vectors comprising each term of the blade.
    const GRADE: u8;

    /// Returns `!self ^ other`
    fn connect<R: Blade>(self, other: R) -> <Self::Dual as Wedge<R>>::Output
    where
        Self::Dual: Wedge<R>,
    {
        self.dual().wedge(other)
    }

    /// Returns whether the blade is "flat" (whether one of its simple component
    /// blades is ∞).
    ///
    /// `epsilon` is a small value used for comparison.
    fn is_flat(self, epsilon: Scalar) -> bool
    where
        Self: Wedge<Blade1>,
    {
        self.wedge(NI)
            .terms()
            .into_iter()
            .all(|t| t.coef.abs() <= epsilon)
    }
}

impl Multivector for Scalar {
    type Terms = [Term; 1];

    type Dual = Pseudoscalar;

    fn zero() -> Self {
        Self::default()
    }

    fn terms(self) -> Self::Terms {
        [Term::new(Axes::S, self)]
    }

    fn has_same_terms_as(self, _other: Self) -> bool {
        true
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
#[allow(missing_docs)]
pub struct Blade1 {
    pub m: Scalar,
    pub p: Scalar,
    pub x: Scalar,
    pub y: Scalar,
}
impl Multivector for Blade1 {
    type Terms = [Term; 4];

    type Dual = Blade3;

    fn zero() -> Self {
        Self::default()
    }

    fn has_same_terms_as(self, _other: Self) -> bool {
        true
    }

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
impl Blade1 {
    /// Returns the ∞ component of the blade (also called nᵢ).
    pub fn ni(self) -> Scalar {
        (self.m + self.p) / 2.0
    }
    /// Returns the nₒ component of the blade.
    pub fn no(self) -> Scalar {
        self.m - self.p
    }

    /// Returns the coordinates of a point in 2D Euclidean space.
    ///
    /// If the blade does not represent a conformal point, then the output may
    /// be meaningless.
    pub fn unpack_point(self) -> (Scalar, Scalar) {
        let no = self.no();
        (self.x / no, self.y / no)
    }

    /// Normalizes a point with respect to nₒ.
    #[must_use]
    pub fn normalize_point(self) -> Self {
        self / self.no()
    }
}

/// 2-blade, used to represent point pairs (real and imaginary), tangent
/// vectors, and flat points.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[allow(missing_docs)]
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

    fn zero() -> Self {
        Self::default()
    }

    fn has_same_terms_as(self, _other: Self) -> bool {
        true
    }

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
impl Blade2 {
    /// Returns the pair of points in a point pair, or `None` the point pair is
    /// imaginary.
    pub fn unpack_point_pair(self) -> Option<[Blade1; 2]> {
        let mag = self.mag2().sqrt();
        mag.is_finite().then(|| {
            [1.0, -1.0].map(|sign| {
                let multiplier = (NI << self).inv();
                (multiplier << self) + (sign * mag * multiplier)
            })
        })
    }
}

/// 3-blade, used to represent circles (real and imaginary).
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct Blade3 {
    pub mpx: Scalar,
    pub mpy: Scalar,
    pub mxy: Scalar,
    pub pxy: Scalar,
}
impl Multivector for Blade3 {
    type Terms = [Term; 4];

    type Dual = Blade1;

    fn zero() -> Self {
        Self::default()
    }

    fn has_same_terms_as(self, _other: Self) -> bool {
        true
    }

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
impl Blade3 {
    /// Returns the line or circle in Euclidean space.
    ///
    /// `epsilon` is a small value used for comparison.
    pub fn unpack(self, epsilon: Scalar) -> LineOrCircle {
        let dual = self.dual();
        if approx::abs_diff_eq!(dual.m, dual.p, epsilon = epsilon) {
            LineOrCircle::Line {
                a: dual.x,
                b: dual.y,
                c: dual.p,
            }
        } else {
            let no = dual.no();
            let cx = dual.x / no;
            let cy = dual.y / no;
            let r2 = dual.mag2() / (no * no);
            let r = r2.abs().sqrt() * r2.signum();
            LineOrCircle::Circle { cx, cy, r }
        }
    }
}
impl From<LineOrCircle> for Blade3 {
    fn from(value: LineOrCircle) -> Self {
        match value {
            LineOrCircle::Line { a, b, c } => Blade1 {
                m: c,
                p: c,
                x: a,
                y: b,
            }
            .dual(),
            LineOrCircle::Circle { cx, cy, r } => crate::circle(crate::point(cx, cy), r),
        }
    }
}

/// Euclidean line or circle.
#[derive(Debug, Copy, Clone, PartialEq)]
#[allow(missing_docs)]
pub enum LineOrCircle {
    /// Line described by the equation _ax+by=c_.
    Line { a: Scalar, b: Scalar, c: Scalar },
    /// Circle described by the equation _(x-cx)^2+(y-cy)^2=r^2_. If _r_ is
    /// negative, then the circle is imaginary.
    ///
    /// Imaginary circles are actually circles with an _imaginary_ radius (so
    /// their radius _squared_ is negative) but we represent them with a
    /// negative radius for convenience.
    ///
    /// The radius may be zero, such as in the case of the dual of a point.
    Circle { cx: Scalar, cy: Scalar, r: Scalar },
}

/// 4-blade, used to represent pseudoscalar quantities.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct Pseudoscalar {
    pub mpxy: Scalar,
}
impl Multivector for Pseudoscalar {
    type Terms = [Term; 1];

    type Dual = Scalar;

    fn zero() -> Self {
        Self::default()
    }

    fn has_same_terms_as(self, _other: Self) -> bool {
        true
    }

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
