use std::fmt;

use crate::*;

const APPROX: Precision = Precision::DEFAULT;

mod objects;
mod transforms;

const SCALES: &[Scalar] = &[1.0, 10.0, -1.0, -10.0];
const NUMBERS: &[Scalar] = &[-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0];
const ORIENTATIONS: [Orientation; 2] = [Orientation::Pos, Orientation::Neg];

fn points() -> impl Clone + Iterator<Item = Point> {
    itertools::chain!(
        itertools::iproduct!(NUMBERS.iter(), NUMBERS.iter()).map(|(&x, &y)| Point::Finite([x, y])),
        [Point::Infinity],
    )
}
fn nonzero_vectors() -> impl Clone + Iterator<Item = [Scalar; 2]> {
    itertools::iproduct!(NUMBERS.iter(), NUMBERS.iter())
        .map(|(&x, &y)| [x, y])
        .filter(|xy| !APPROX.eq_zero(xy))
}
fn point_pairs() -> impl Clone + Iterator<Item = (Point, Point)> {
    itertools::iproduct!(points(), points()).filter(|(p, q)| !APPROX.eq(p, q))
}
fn dipoles() -> impl Clone + Iterator<Item = Dipole> {
    itertools::chain!(
        point_pairs().map(|(p, q)| Dipole::Real([p, q])),
        itertools::iproduct!(points(), nonzero_vectors()).map(|(p, v)| Dipole::Tangent(p, v)),
        point_pairs().map(|(p, q)| Dipole::Imaginary([p, q])),
    )
}
fn circles() -> impl Clone + Iterator<Item = Circle> {
    itertools::chain!(
        itertools::iproduct!(nonzero_vectors(), NUMBERS)
            .map(|([a, b], &c)| { Circle::Line { a, b, c } }),
        itertools::iproduct!(NUMBERS, NUMBERS, NUMBERS, ORIENTATIONS)
            .map(|(&cx, &cy, &r, ori)| { Circle::Circle { cx, cy, r, ori } }),
        ORIENTATIONS.map(Circle::Infinity),
    )
}

#[track_caller]
fn assert_mvec_eq<T: Multivector>(a: T, b: T) {
    let a = a.rescale_oriented();
    let b = b.rescale_oriented();
    assert!(APPROX.eq(a, b), "expected {a:?} and {b:?} to be equal")
}

#[track_caller]
fn assert_approx_eq<T: ApproxEq + fmt::Debug>(a: T, b: T) {
    assert!(APPROX.eq(&a, &b), "expected {a:?} and {b:?} to be equal",)
}
