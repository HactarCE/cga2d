use std::{
    hash::Hash,
    ops::{Mul, Neg},
};

use approx_collections::{
    ApproxEq, ApproxEqZero, ApproxHash, ApproxHasher, ApproxSign, ForEachFloat, Precision,
};

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
    fn no_eq_zero(self, prec: Precision) -> bool {
        prec.eq(self.m, self.p)
    }

    /// Returns the coordinates of a point in 2D Euclidean space, or `None` if
    /// the point is not on the conformal sphere.
    ///
    /// Uses [`Precision::DEFAULT`].
    pub fn unpack(self) -> Option<Point> {
        self.unpack_with_prec(Precision::DEFAULT)
    }
    /// Returns the coordinates of a point in 2D Euclidean space, or `None` if
    /// the point is not on the conformal sphere.
    pub fn unpack_with_prec(self, prec: Precision) -> Option<Point> {
        if self.no_eq_zero(prec) {
            if [self.x, self.y].approx_eq_zero(prec) {
                Some(Point::Infinity)
            } else {
                None
            }
        } else {
            if self.mag2().approx_eq_zero(prec) {
                let no = self.no();
                Some(Point::Finite([self.x / no, self.y / no]))
            } else {
                None
            }
        }
    }

    /// Returns the coordinates of a point in 2D Euclidean space.
    ///
    /// The result is undefined if the point is not on the conformal sphere.
    ///
    /// Uses [`Precision::DEFAULT`].
    pub fn unpack_unchecked(self) -> Point {
        self.unpack_unchecked_with_prec(Precision::DEFAULT)
    }
    /// Returns the coordinates of a point in 2D Euclidean space.
    ///
    /// The result is undefined if the point is not on the conformal sphere.
    pub fn unpack_unchecked_with_prec(self, prec: Precision) -> Point {
        if self.no_eq_zero(prec) {
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

/// Euclidean point.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Point {
    /// Finite point with X and Y coordinates.
    Finite([Scalar; 2]),
    /// Point at infinity.
    Infinity,
}
impl Default for Point {
    fn default() -> Self {
        Point::Finite([0.0; 2])
    }
}
impl From<Point> for Blade1 {
    fn from(value: Point) -> Self {
        value.to_blade()
    }
}
impl ApproxEq for Point {
    fn approx_eq(&self, other: &Self, prec: Precision) -> bool {
        match (self, other) {
            (Point::Finite(a), Point::Finite(b)) => prec.eq(a, b),
            (Point::Infinity, Point::Infinity) => true,

            (Point::Finite(_), _) | (_, Point::Finite(_)) => false,
        }
    }
}
impl ApproxHash for Point {
    fn approx_hash<H: ApproxHasher>(&self, state: &mut H) {
        std::mem::discriminant(self).hash(state);
        match self {
            Point::Finite(xy) => xy.approx_hash(state),
            Point::Infinity => (),
        }
    }
}
impl ForEachFloat for Point {
    fn for_each_float(&mut self, f: &mut impl FnMut(&mut f64)) {
        match self {
            Point::Finite(xy) => xy.for_each_float(f),
            Point::Infinity => (),
        }
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

    /// Returns a blade representing the point.
    pub fn to_blade(self) -> Blade1 {
        match self {
            Point::Finite([x, y]) => crate::point(x, y),
            Point::Infinity => NI,
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
    pub fn unpack(self) -> Dipole {
        self.unpack_with_prec(Precision::DEFAULT)
    }

    /// Returns the real, imaginary, or tangent pair of points in a point pair.
    pub fn unpack_with_prec(self, prec: Precision) -> Dipole {
        let mag2 = self.mag2();

        match mag2.approx_sign(prec) {
            approx_collections::Sign::Negative => {
                Dipole::Imaginary(self.dual().unpack_raw((-mag2).sqrt(), prec, [1.0, -1.0]))
            }
            approx_collections::Sign::Zero => {
                let point = self.unpack_raw(0.0, prec, [0.0])[0];
                let vector = match point {
                    Point::Finite(_) => [self.mx - self.px, self.my - self.py],
                    Point::Infinity => [self.mx, self.my], // (mx, my) = (px, py)
                };
                Dipole::Tangent(point, vector).normalize()
            }
            approx_collections::Sign::Positive => {
                Dipole::Real(self.unpack_raw(mag2.sqrt(), prec, [1.0, -1.0]))
            }
        }
    }

    /// Unpacks a real point pair or tangent vector.
    ///
    /// The result is undefined for an imaginary point pair.
    ///
    /// - To get both points of a real point, use `distances = [1.0, -1.0]`.
    /// - To get the point of a tangent vector, use `distances = [0.0]`.
    fn unpack_raw<const N: usize>(
        self,
        mag: f64,
        prec: Precision,
        distances: [f64; N],
    ) -> [Point; N] {
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

        distances.map(|sign| {
            ((multiplier << self) + (sign * mag * multiplier)).unpack_unchecked_with_prec(prec)
        })
    }

    /// Rotates a tangent vector counterclockwise by `angle` (in radians).
    pub fn rotate(self, angle: Scalar) -> Self {
        self * angle.cos() + self.dual() * angle.sin()
    }
}

/// Real or imaginary point pair in Euclidean space.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Dipole {
    /// Pair of real points.
    Real([Point; 2]),
    /// Tangent point, represented by a point and a Euclidean tangent vector.
    Tangent(Point, [Scalar; 2]),
    /// Imaginary point pair, represented by the points in its dual.
    Imaginary([Point; 2]),
}
impl From<Blade2> for Dipole {
    fn from(value: Blade2) -> Self {
        value.unpack()
    }
}
impl From<Dipole> for Blade2 {
    fn from(value: Dipole) -> Self {
        value.to_blade()
    }
}
impl ApproxEq for Dipole {
    fn approx_eq(&self, other: &Self, prec: Precision) -> bool {
        match (self, other) {
            (Dipole::Real(pp1), Dipole::Real(pp2)) => prec.eq(pp1, pp2),
            (Dipole::Tangent(p1, v1), Dipole::Tangent(p2, v2)) => {
                prec.eq(p1, p2) && prec.eq(v1, v2)
            }
            (Dipole::Imaginary(pp1), Dipole::Imaginary(pp2)) => prec.eq(pp1, pp2),

            (Dipole::Real(_), _) | (_, Dipole::Real(_)) => false,
            (Dipole::Tangent(_, _), _) | (_, Dipole::Tangent(_, _)) => false,
        }
    }
}
impl ApproxHash for Dipole {
    fn approx_hash<H: ApproxHasher>(&self, state: &mut H) {
        std::mem::discriminant(self).hash(state);
        match self {
            Dipole::Real(pp) => pp.approx_hash(state),
            Dipole::Tangent(p, v) => {
                p.approx_hash(state);
                v.approx_hash(state);
            }
            Dipole::Imaginary(pp) => pp.approx_hash(state),
        }
    }
}
impl ForEachFloat for Dipole {
    fn for_each_float(&mut self, f: &mut impl FnMut(&mut f64)) {
        match self {
            Dipole::Real(pp) => pp.for_each_float(f),
            Dipole::Tangent(p, v) => {
                p.for_each_float(f);
                v.for_each_float(f);
            }
            Dipole::Imaginary(pp) => pp.for_each_float(f),
        }
    }
}
impl Neg for Dipole {
    type Output = Dipole;

    fn neg(self) -> Self::Output {
        match self {
            Dipole::Real([p, q]) => Dipole::Real([q, p]),
            Dipole::Tangent(p, [vx, vy]) => Dipole::Tangent(p, [-vx, -vy]),
            Dipole::Imaginary([p, q]) => Dipole::Imaginary([q, p]),
        }
    }
}
impl Dipole {
    /// Normalizes tangent points so that the vector is unit-length.
    #[must_use = "normalize() returns a mutated copy"]
    pub fn normalize(self) -> Self {
        match self {
            Dipole::Tangent(p, [vx, vy]) => {
                let scale = vector_mag([vx, vy]).recip();
                Dipole::Tangent(p, [vx * scale, vy * scale])
            }
            _ => self,
        }
    }

    /// Returns the real point pair, or `None` if it is not real.
    pub fn real(self) -> Option<[Point; 2]> {
        match self {
            Dipole::Real(pair) => Some(pair),
            _ => None,
        }
    }

    /// Returns a blade representing the dipole.
    pub fn to_blade(self) -> Blade2 {
        match self {
            Dipole::Real([p1, p2]) => p1.to_blade() ^ p2.to_blade(),
            Dipole::Tangent(p, v) => crate::tangent_point(p, v),
            Dipole::Imaginary([p1, p2]) => (p1.to_blade() ^ p2.to_blade()).antidual(),
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
    pub fn unpack(self) -> Circle {
        self.unpack_with_prec(Precision::DEFAULT)
    }

    /// Returns the line or circle in Euclidean space.
    pub fn unpack_with_prec(self, prec: Precision) -> Circle {
        let dual = self.antidual();
        if dual.no_eq_zero(prec) {
            if [dual.x, dual.y].approx_eq_zero(prec) {
                Circle::Infinity(Orientation::from(dual.m)) // m = p
            } else {
                Circle::Line {
                    a: dual.x,
                    b: dual.y,
                    c: dual.p,
                }
                .normalize()
            }
        } else {
            let no = dual.no();
            let cx = dual.x / no;
            let cy = dual.y / no;
            let r2 = dual.mag2() / (no * no);
            let r = r2.abs().sqrt() * r2.signum();
            Circle::Circle {
                cx,
                cy,
                r,
                ori: Orientation::from(no),
            }
        }
    }
}

/// Euclidean line or circle.
#[derive(Debug, Copy, Clone, PartialEq)]
#[allow(missing_docs)]
pub enum Circle {
    /// Line described by the equation _ax+by=c_.
    ///
    /// The vector `(a, b)` represents the normal vector of the line. `a` and
    /// `b` cannot both be zero. When returned from [`Blade3::unpack()`] or
    /// [`Blade3::unpack_with_prec()`], the vector `(a, b)` is unit-length.
    Line { a: Scalar, b: Scalar, c: Scalar },
    /// Real or imaginary circle described by the equation
    /// _(x-cx)^2+(y-cy)^2=r^2_. If _r_ is positive, then the circle is real; if
    /// _r_ is negative, then the circle is imaginary.
    ///
    /// Imaginary circles are actually circles with an _imaginary_ radius (so
    /// their radius _squared_ is negative) but we represent them with a
    /// negative radius for convenience.
    ///
    /// The radius may be zero, in which case the circle is the dual of a real
    /// point.
    Circle {
        cx: Scalar,
        cy: Scalar,
        r: Scalar,
        ori: Orientation,
    },
    /// Circle with zero radius centered on the point at infinity.
    Infinity(Orientation),
}
impl From<Blade3> for Circle {
    fn from(value: Blade3) -> Self {
        value.unpack()
    }
}
impl From<Circle> for Blade3 {
    fn from(value: Circle) -> Self {
        value.to_blade()
    }
}
impl ApproxEq for Circle {
    fn approx_eq(&self, other: &Self, prec: Precision) -> bool {
        match (self, other) {
            (
                Circle::Line {
                    a: a1,
                    b: b1,
                    c: c1,
                },
                Circle::Line {
                    a: a2,
                    b: b2,
                    c: c2,
                },
            ) => prec.eq([a1, b1, c1], [a2, b2, c2]),
            (
                Circle::Circle {
                    cx: cx1,
                    cy: cy1,
                    r: r1,
                    ori: ori1,
                },
                Circle::Circle {
                    cx: cx2,
                    cy: cy2,
                    r: r2,
                    ori: ori2,
                },
            ) => prec.eq([cx1, cy1, r1], [cx2, cy2, r2]) && ori1 == ori2,
            (Circle::Infinity(ori1), Circle::Infinity(ori2)) => ori1 == ori2,

            (Circle::Line { .. }, _) | (_, Circle::Line { .. }) => false,
            (Circle::Circle { .. }, _) | (_, Circle::Circle { .. }) => false,
        }
    }
}
impl ApproxHash for Circle {
    fn approx_hash<H: ApproxHasher>(&self, state: &mut H) {
        std::mem::discriminant(self).hash(state);
        match *self {
            Circle::Line { a, b, c } => [a, b, c].approx_hash(state),
            Circle::Circle { cx, cy, r, ori } => {
                [cx, cy, r].approx_hash(state);
                ori.hash(state);
            }
            Circle::Infinity(ori) => ori.hash(state),
        }
    }
}
impl ForEachFloat for Circle {
    fn for_each_float(&mut self, f: &mut impl FnMut(&mut f64)) {
        match self {
            Circle::Line { a, b, c } => {
                f(a);
                f(b);
                f(c);
            }
            Circle::Circle { cx, cy, r, ori: _ } => {
                f(cx);
                f(cy);
                f(r);
            }
            Circle::Infinity(_) => (),
        }
    }
}
impl Neg for Circle {
    type Output = Circle;

    fn neg(self) -> Self::Output {
        match self {
            Circle::Line { a, b, c } => Circle::Line {
                a: -a,
                b: -b,
                c: -c,
            },
            Circle::Circle { cx, cy, r, ori } => Circle::Circle {
                cx,
                cy,
                r,
                ori: -ori,
            },
            Circle::Infinity(ori) => Circle::Infinity(-ori),
        }
    }
}
impl Circle {
    /// Normalizes lines so that the vector `(a, b)` is unit-length.
    #[must_use = "normalize() returns a mutated copy"]
    pub fn normalize(self) -> Self {
        match self {
            Circle::Line { a, b, c } => {
                let scale = vector_mag([a, b]).recip();
                Circle::Line {
                    a: a * scale,
                    b: b * scale,
                    c: c * scale,
                }
            }
            _ => self,
        }
    }

    /// Returns a blade representing the circle.
    pub fn to_blade(self) -> Blade3 {
        match self {
            Circle::Line { a, b, c } => crate::line(a, b, c),
            Circle::Circle { cx, cy, r, ori } => ori * crate::circle(crate::point(cx, cy), r),
            Circle::Infinity(ori) => ori * !crate::NI,
        }
    }
}

/// Orientation of an object.
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq, Hash)]
pub enum Orientation {
    /// Positive or counterclockwise.
    #[default]
    Pos = 1,
    /// Negative or clockwise.
    Neg = -1,
}
impl From<Scalar> for Orientation {
    fn from(value: Scalar) -> Self {
        if value.is_sign_positive() {
            Orientation::Pos
        } else {
            Orientation::Neg
        }
    }
}
impl<T: Neg<Output = T>> Mul<T> for Orientation {
    type Output = T;

    fn mul(self, rhs: T) -> Self::Output {
        match self {
            Orientation::Pos => rhs,
            Orientation::Neg => -rhs,
        }
    }
}
impl Neg for Orientation {
    type Output = Orientation;

    fn neg(self) -> Self::Output {
        match self {
            Orientation::Pos => Orientation::Neg,
            Orientation::Neg => Orientation::Pos,
        }
    }
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

fn vector_mag([x, y]: [Scalar; 2]) -> Scalar {
    ((x * x) + (y * y)).sqrt()
}
