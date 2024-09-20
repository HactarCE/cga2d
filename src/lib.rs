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
mod term;

pub use axes::Axes;
pub use blade::{Blade, Blade1, Blade2, Blade3, Pseudoscalar};
pub use multivector::Multivector;
pub use term::Term;

/// 0-blade, used to represent scalar quantities.
pub type Scalar = f64;
