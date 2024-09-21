//! [Conformal Geometric Algebra][cga-wiki] in 2D.
//!
//! [cga-wiki]: https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page

#![warn(missing_docs, rust_2018_idioms)]
#![warn(
    clippy::cargo,
    clippy::doc_markdown,
    clippy::if_then_some_else_none,
    clippy::manual_let_else,
    clippy::semicolon_if_nothing_returned,
    clippy::semicolon_inside_block,
    clippy::stable_sort_primitive,
    clippy::undocumented_unsafe_blocks,
    clippy::uninlined_format_args,
    clippy::unwrap_used
)]

mod axes;
mod blade;
mod multivector;
mod ops;
mod rotoflector;
mod term;

/// Traits and basic types (blades, NI, NO, rotor/flector/rotoflector).
pub mod prelude {
    pub use crate::traits::*;

    pub use crate::blade::{Blade1, Blade2, Blade3, Pseudoscalar, NI, NO};
    pub use crate::rotoflector::{Flector, Rotoflector, Rotor};
}

/// Traits.
pub mod traits {
    pub use crate::blade::Blade;
    pub use crate::multivector::Multivector;
    pub use crate::ops::Wedge;
}

pub use axes::Axes;
pub use blade::{Blade, Blade1, Blade2, Blade3, LineOrCircle, Pseudoscalar, NI, NO};
pub use multivector::Multivector;
pub use ops::Wedge;
pub use rotoflector::{Flector, Rotoflector, Rotor};
pub use term::Term;

/// 0-blade, used to represent scalar quantities.
pub type Scalar = f64;

/// Interpolates between two multivectors using an angle.
pub fn slerp<M: Multivector>(a: M, b: M, angle: Scalar) -> M {
    a.normalize() * angle.cos() + b.normalize() * angle.sin()
}

/// Lifts a Euclidean vector into 2D CGA.
pub fn vector(x: Scalar, y: Scalar) -> Blade1 {
    let m = 0.0;
    let p = 0.0;
    Blade1 { m, p, x, y }
}

/// Lifts a Euclidean point into 2D conformal space.
pub fn point(x: Scalar, y: Scalar) -> Blade1 {
    let mag2 = x * x + y * y;
    NO + vector(x, y) + 0.5 * mag2 * NI
}

/// Constructs a circle given a center point and a radius.
///
/// If `radius` is negative, constructs an imaginary circle.
pub fn circle(center: Blade1, radius: Scalar) -> Blade3 {
    !(center - 0.5 * radius * radius.abs() * NI)
}

/// Constructs a line from the equation _ax+by+c=0_.
pub fn line(a: Scalar, b: Scalar, c: Scalar) -> Blade3 {
    Blade1 {
        m: c,
        p: c,
        x: a,
        y: b,
    }
    .dual()
}

#[cfg(test)]
mod tests;
