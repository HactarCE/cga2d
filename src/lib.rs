//! [Conformal Geometric Algebra][cga-wiki] in 2D.
//!
//! [cga-wiki]:
//!     https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page
//!
//! ```rust
//! use cga2d::prelude::*;
//!
//! let p1 = cga2d::point(1.0, 3.0);
//! let p2 = cga2d::point(-3.0, 5.0);
//! let line = p1 ^ p2 ^ NI;
//!
//! let epsilon = 0.0001; // comparison threshold
//! assert!(line.is_flat(epsilon));
//!
//! assert_eq!(!(line ^ cga2d::point(-1.0, 4.0)), 0.0);
//!
//! let circ = cga2d::circle(cga2d::point(3.0, 1.5), 3.0);
//! assert_eq!(circ.sandwich(NI).unpack_point(), (3.0, 1.5));
//!
//! let rot90_ccw: Rotor = cga2d::line(1.0, 1.0, 0.0) * cga2d::line(1.0, 0.0, 0.0);
//! assert_eq!(rot90_ccw.sandwich(cga2d::point(3.0, 4.0)).unpack_point(), (-4.0, 3.0));
//! ```
//!
//! # Multivector types
//!
//! There is no unified multivector type. Instead, there is a [`Multivector`]
//! trait, implemented by several blades (which also implement the [`Blade`]
//! trait) and several rotoflectors.
//!
//! ## Blades
//!
//! | Blade type           | Grade | Used to represent                        |
//! |:-------------------- |:-----:|:---------------------------------------- |
//! | [`Scalar`] = [`f64`] |     0 | Scalar quantities                        |
//! | [`Blade1`]           |     1 | Points, vectors, round points            |
//! | [`Blade2`]           |     2 | Point pairs, tangent points, flat points |
//! | [`Blade3`]           |     3 | Circles (real & imaginary), lines        |
//! | [`Pseudoscalar`]     |     4 | Pseudoscalar quantities                  |
//!
//! ## Rotoflectors
//!
//! | Rotoflector type | Subalgebra  | Used to represent                                    |
//! |:---------------- |:-----------:|:---------------------------------------------------- |
//! | [`Rotor`]        | even        | orientation-preserving conformal transformations     |
//! | [`Flector`]      | odd         | non-orientation-preserving conformal transformations |
//! | [`Rotoflector`]  | even or odd | conformal transformations                            |
//!
//! Note that [`Rotoflector`] contains _either_ even terms _or_ odd terms. It is
//! not a general-purpose multivector.
//!
//! There is no general-purpose multivector type.
//!
//! # Construction
//!
//! The constants [`NI`] and [`NO`] contain 1-blades representing the point at
//! infinity and the origin respectively.
//!
//! ## From components
//!
//! All multivectors can be constructed directly via components.
//!
//! ```rust
//! # use cga2d::prelude::*;
//! let my_vector = Blade1 {
//!     m: 0.0,
//!     p: 0.0,
//!     x: 3.0,
//!     y: 4.0,
//! };
//! ```
//!
//! ## Helper functions
//!
//! There are also several convenience functions built-in for constructing
//! common multivectors, and these can be composed with operations.
//!
//! ```rust
//! # use cga2d::prelude::*;
//! // Type annotations are not required
//! let v: Blade1 = cga2d::vector(3.0, 4.0);
//! let center: Blade1 = cga2d::point(3.0, 4.0);
//! let real_circle: Blade3 = cga2d::circle(center, 7.0);
//! let imag_circle: Blade3 = cga2d::circle(center, -7.0);
//! let line: Blade3 = cga2d::line(3.0, 4.0, 2.0);
//!
//! let point_pair: Blade2 = cga2d::point(3.0, 4.0) ^ NO;
//! let flat_point: Blade2 = cga2d::point(3.0, 4.0) ^ NI;
//! ```
//!
//! Additionally, all multivectors can be constructed by summing terms. Terms
//! that cannot be represented by the multivector are discarded. I.e., the terms
//! are grade-projected.
//!
//! ## From terms
//!
//! ```rust
//! # use cga2d::prelude::*;
//! let vector: Blade1 = [
//!     cga2d::Term::new(cga2d::Axes::X, 3.0),
//!     cga2d::Term::new(cga2d::Axes::X, 4.0)
//! ]
//! .into_iter()
//! .sum();
//! ```
//!
//! ## From blades
//!
//! Rotoflectors can be constructed from [`Blade`]s of the appropriate grade.
//!
//! ```rust
//! # use cga2d::prelude::*;
//! let center = cga2d::point(3.0, 4.0);
//! let circle = cga2d::circle(center, 7.0);
//! let circle_inversion = Flector::from(circle);
//!
//! let central_inversion = Rotor::from(NI ^ NO);
//! let inverted_circle = central_inversion.sandwich(circle);
//! let epsilon = 0.0001; // comparison threshold
//! assert_eq!(inverted_circle.unpack(epsilon), cga2d::LineOrCircle::Circle {
//!     cx: -3.0,
//!     cy: -4.0,
//!     r: 7.0
//! });
//! ```
//!
//! # Operations
//!
//! ## Wedge product
//!
//! - [`Blade`] `^` [`Blade`] `->` [`Blade`]
//!
//! ## Antiwedge product
//!
//! - [`Blade`] `&` [`Blade`] `->` [`Blade`]
//!
//! ## Left contraction
//!
//! - [`Multivector`] `<<` [`Multivector`] `->` [`Multivector`]
//!
//! ## Negation
//!
//! - `-`[`Multivector`] `->` [`Multivector`]
//!
//! ## Dual
//!
//! - `!`[`Multivector`] `->` [`Multivector`]
//!
//! ## Scaling
//!
//! - [`Multivector`] `*` [`Scalar`] `->` [`Multivector`]
//! - [`Scalar`] `*` [`Multivector`] `->` [`Multivector`]
//! - [`Multivector`] `/` [`Scalar`] `->` [`Multivector`]
//!
//! ## Geometric product
//!
//! - [`Multivector`] `*` [`Multivector`] `->` [`Rotor`] (where sum of grades is
//!   even)
//! - [`Multivector`] `*` [`Multivector`] `->` [`Flector`] (where sum of grades
//!   is odd)
//!
//! ## Geometric product by inverse
//!
//! - [`Multivector`] `/` [`Multivector`] `->` [`Multivector`]
//!
//! ## Addition & subtraction
//!
//! - [`Multivector`] `+` [`Multivector`] `->` [`Multivector`] (must have same
//!   type)
//! - [`Multivector`] `+` [`Term`] `->` [`Multivector`] (panics if the
//!   multivector type doesn't support the term)
//! - [`Multivector`] `-` [`Term`] `->` [`Multivector`] (panics if the
//!   multivector type doesn't support the term)
//!
//! ## Indexing
//!
//! Indexing panics if the multivector type doesn't support the term. For a
//! non-panicking alternative, see [`Multivector::get()`] and
//! [`Multivector::get_mut()`].
//!
//! - [`Multivector`]`[`[`Axes`]`]` `->` [`Scalar`] (panics if the multivector
//!   type doesn't support the term)

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
