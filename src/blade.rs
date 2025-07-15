use approx_collections::{ApproxSign, Precision};

use super::{Axes, Multivector, Scalar, Term, Wedge};

/// ∞, also called n ͚, representing the point at infinity.
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

macro_rules! impl_multivector {
    ($type:ty {
        dual: $dual:ty,
        len: $len:expr,
        terms: [$(($axes:path, $field:ident)),* $(,)?] $(,)?
    }) => {
        impl Multivector for $type {
            type Terms = [Term; $len];

            type Dual = $dual;

            fn zero() -> Self {
                Self::default()
            }

            fn has_same_terms_as(self, _other: Self) -> bool {
                true
            }

            fn terms(self) -> [Term; $len] {
                [$(Term::new($axes, self.$field)),*]
            }

            fn get(&self, axes: Axes) -> Option<&Scalar> {
                match axes {
                    $($axes => Some(&self.$field),)*
                    _ => None,
                }
            }

            fn get_mut(&mut self, axes: Axes) -> Option<&mut Scalar> {
                match axes {
                    $($axes => Some(&mut self.$field),)*
                    _ => None,
                }
            }
        }
    };
}

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
    /// Uses [`Precision::DEFAULT`].
    fn is_flat(self) -> bool
    where
        Self: Wedge<Blade1>,
    {
        self.is_flat_with_prec(Precision::DEFAULT)
    }

    /// Returns whether the blade is "flat" (whether one of its simple component
    /// blades is ∞).
    fn is_flat_with_prec(self, prec: Precision) -> bool
    where
        Self: Wedge<Blade1>,
    {
        self.wedge(NI).terms().into_iter().all(|t| prec.eq_zero(t))
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
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Zeroable, bytemuck::NoUninit))]
#[repr(C)]
#[allow(missing_docs)]
pub struct Blade1 {
    pub m: Scalar,
    pub p: Scalar,
    pub x: Scalar,
    pub y: Scalar,
}
impl_multivector!(Blade1 {
    dual: Blade3,
    len: 4,
    terms: [(Axes::M, m), (Axes::P, p), (Axes::X, x), (Axes::Y, y)],
});
impl Blade for Blade1 {
    const GRADE: u8 = 1;
}
impl Blade1 {
    /// Returns the ∞ component of the blade (also called n ͚).
    pub fn ni(self) -> Scalar {
        (self.m + self.p) / 2.0
    }
    /// Returns the nₒ component of the blade.
    pub fn no(self) -> Scalar {
        self.m - self.p
    }

    /// Returns the coordinates of a point in 2D Euclidean space.
    ///
    /// Uses [`Precision::DEFAULT`].
    pub fn unpack(self) -> Point {
        self.unpack_with_prec(Precision::DEFAULT)
    }
    /// Returns the coordinates of a point in 2D Euclidean space.
    pub fn unpack_with_prec(self, prec: Precision) -> Point {
        // TODO: handle other cases
        if prec.eq(self.m, self.p) && prec.eq_zero([self.x, self.y]) {
            Point::Infinity
        } else {
            let no = self.no();
            Point::Finite([self.x / no, self.y / no])
        }
    }

    /// Normalizes a point with respect to nₒ.
    #[must_use]
    pub fn normalize_point(self) -> Self {
        self / self.no()
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Point {
    Finite([Scalar; 2]),
    Infinity,
    RoundPointDual(LineOrCircle),
    CircleDual(LineOrCircle),
}
// TODO: `impl ApproxEq, ApproxHash for Point`
impl Default for Point {
    fn default() -> Self {
        Point::Finite([0.0; 2])
    }
}
impl Point {
    /// Returns the finite coordinates of the point, or `None` if it is not a
    /// real finite point.
    pub fn finite(self) -> Option<[Scalar; 2]> {
        match self {
            Point::Finite(xy) => Some(xy),
            _ => None,
        }
    }
}

/// 2-blade, used to represent point pairs (real and imaginary), tangent
/// vectors, and flat points.
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Zeroable, bytemuck::NoUninit))]
#[repr(C)]
#[allow(missing_docs)]
pub struct Blade2 {
    pub mp: Scalar,
    pub mx: Scalar,
    pub px: Scalar,
    pub my: Scalar,
    pub py: Scalar,
    pub xy: Scalar,
}
impl_multivector!(Blade2 {
    dual: Blade2,
    len: 6,
    terms: [
        (Axes::MP, mp),
        (Axes::MX, mx),
        (Axes::PX, px),
        (Axes::MY, my),
        (Axes::PY, py),
        (Axes::XY, xy),
    ],
});
impl Blade for Blade2 {
    const GRADE: u8 = 2;
}
impl Blade2 {
    /// Returns the real, imaginary, or tangent pair of points in a point pair.
    ///
    /// Uses [`Precision::DEFAULT`].
    pub fn unpack(self) -> PointPair {
        self.unpack_with_prec(Precision::DEFAULT)
    }

    /// Returns the real, imaginary, or tangent pair of points in a point pair.
    pub fn unpack_with_prec(self, prec: Precision) -> PointPair {
        let mag2 = self.mag2();

        match mag2.approx_sign(prec) {
            approx_collections::Sign::Negative => PointPair::Imaginary {},
            approx_collections::Sign::Zero => PointPair::Tangent {
                point: self.unpack_raw(0.0)[0],
            },
            approx_collections::Sign::Positive => PointPair::Real(self.unpack_raw(mag2.sqrt())),
        }
    }

    fn unpack_raw(self, mag: f64) -> [Blade1; 2] {
        // We need to use an arbitrary point that is not in the point
        // pair. Of these three points, there is guaranteed to be one
        // with a nonzero result.
        let candidates = [NI, NO, crate::point(1.0, 0.0)];
        let multiplier = candidates
            .iter()
            .map(|&arbitrary| arbitrary << self)
            .max_by(|a, b| f64::total_cmp(&a.mag2(), &b.mag2()))
            .expect("empty candidates list")
            .inv();

        [1.0, -1.0].map(|sign| (multiplier << self) + (sign * mag * multiplier))
    }

    /// Rotates a tangent vector counterclockwise by `angle` (in radians).
    pub fn rotate(self, angle: Scalar) -> Self {
        self * angle.cos() + self.dual() * angle.sin()
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum PointPair {
    Real([Blade1; 2]),
    #[non_exhaustive]
    Tangent {
        point: Blade1,
        // TODO: vector
    },
    #[non_exhaustive]
    Imaginary {
        // TODO: imaginary points
    },
}
// TODO: `impl ApproxEq, ApproxHash for PointPair`
impl PointPair {
    /// Returns the real point pair, or `None` if it is not real.
    pub fn real(self) -> Option<[Blade1; 2]> {
        match self {
            PointPair::Real(pair) => Some(pair),
            _ => None,
        }
    }
}

/// 3-blade, used to represent circles (real and imaginary).
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Zeroable, bytemuck::NoUninit))]
#[repr(C)]
#[allow(missing_docs)]
pub struct Blade3 {
    pub mpx: Scalar,
    pub mpy: Scalar,
    pub mxy: Scalar,
    pub pxy: Scalar,
}
impl_multivector!(Blade3 {
    dual: Blade1,
    len: 4,
    terms: [
        (Axes::MPX, mpx),
        (Axes::MPY, mpy),
        (Axes::MXY, mxy),
        (Axes::PXY, pxy),
    ],
});
impl Blade for Blade3 {
    const GRADE: u8 = 3;
}
impl Blade3 {
    /// Returns the line or circle in Euclidean space.
    ///
    /// Uses [`Precision::DEFAULT`].
    pub fn unpack(self) -> LineOrCircle {
        self.unpack_with_prec(Precision::DEFAULT)
    }

    /// Returns the line or circle in Euclidean space.
    pub fn unpack_with_prec(self, prec: Precision) -> LineOrCircle {
        let dual = self.dual();
        if prec.eq(dual.m, dual.p) {
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
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Zeroable, bytemuck::NoUninit))]
#[repr(C)]
#[allow(missing_docs)]
pub struct Pseudoscalar {
    pub mpxy: Scalar,
}
impl_multivector!(Pseudoscalar {
    dual: Scalar,
    len: 1,
    terms: [(Axes::MPXY, mpxy)],
});
impl Blade for Pseudoscalar {
    const GRADE: u8 = 4;
}
impl ApproxSign for Pseudoscalar {
    fn approx_sign(&self, prec: Precision) -> approx_collections::Sign {
        self.mpxy.approx_sign(prec)
    }
}
