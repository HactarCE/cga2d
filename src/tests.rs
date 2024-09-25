use std::fmt;

use super::*;

const EPS: Scalar = 0.0001;

#[track_caller]
fn assert_approx_eq<T: fmt::Debug + approx::AbsDiffEq<Epsilon = Scalar>>(a: T, b: T) {
    approx::assert_abs_diff_eq!(a, b, epsilon = EPS);
}

#[track_caller]
fn assert_eq_up_to_scale<M: Multivector>(a: M, b: M) {
    let max_axes = a
        .terms()
        .into_iter()
        .max_by(|t1, t2| f64::total_cmp(&t1.coef.abs(), &t2.coef.abs()))
        .unwrap()
        .axes;
    approx::assert_abs_diff_eq!(
        a / *a.get(max_axes).unwrap_or(&0.0),
        b / *b.get(max_axes).unwrap_or(&0.0),
        epsilon = EPS,
    );
}

const SCALARS: &[Scalar] = &[1.0, 10.0, -1.0, -10.0];

#[test]
fn test_point() {
    for &s in SCALARS {
        for (x, y) in [(0.0, 0.0), (3.0, -4.0)] {
            let p = s * point(x, y);
            assert!(!p.is_flat(EPS));

            let (x_out, y_out) = p.unpack_point();
            assert_approx_eq(x, x_out);
            assert_approx_eq(y, y_out);
        }
    }

    assert!(NI.is_flat(EPS));
}

#[test]
fn test_point_pair() {
    for &s in SCALARS {
        for (p1, p2) in [
            (crate::point(1.0, 2.0), crate::point(3.0, -4.0)),
            (NO, NI),
            (NI, NO),
            (point(1.0, 2.0), NI),
            (NI, point(1.0, 2.0)),
        ] {
            let pp = s * p1 ^ p2;

            assert_eq!(pp.is_flat(EPS), p1 == NI || p2 == NI);

            let [mut p1_out, mut p2_out] = pp.unpack_point_pair().unwrap();
            if s < 0.0 {
                std::mem::swap(&mut p1_out, &mut p2_out);
            }
            assert_eq_up_to_scale(p1, p1_out);
            assert_eq_up_to_scale(p2, p2_out);
        }
    }
}

#[test]
fn test_unpack_circle() {
    for &s in SCALARS {
        for ((cx, cy), r) in [((3.0, -4.0), 7.0), ((-1.0, 6.0), 0.0), ((3.0, -4.0), -9.0)] {
            let circ = s * circle(point(cx, cy), r);
            assert!(!circ.is_flat(EPS));

            match circ.unpack(EPS) {
                LineOrCircle::Line { .. } => panic!("expected circle"),
                LineOrCircle::Circle {
                    cx: cx_out,
                    cy: cy_out,
                    r: r_out,
                } => {
                    assert_approx_eq(cx, cx_out);
                    assert_approx_eq(cy, cy_out);
                    assert_approx_eq(r, r_out);
                }
            }
        }
    }
}

#[test]
fn test_unpack_line() {
    for &s in SCALARS {
        for (a, b, c) in [
            (-3.0, -5.0, 7.0),
            (-3.0, -5.0, -7.0),
            (-1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
        ] {
            let l = s * line(a, b, c);
            assert!(l.is_flat(EPS));

            if a != 0.0 {
                assert_approx_eq((point(c / a, 0.0) ^ l).dual(), 0.0);
            }
            if b != 0.0 {
                assert_approx_eq((point(0.0, c / b) ^ l).dual(), 0.0);
            }

            match l.unpack(EPS) {
                LineOrCircle::Line {
                    a: a_out,
                    b: b_out,
                    c: c_out,
                } => {
                    assert_approx_eq(a * c_out, a_out * c);
                    assert_approx_eq(b * c_out, b_out * c);
                }
                LineOrCircle::Circle { .. } => panic!("expected line"),
            }
        }
    }
}

#[test]
fn test_point_reflection() {
    let p = point(1.0, 4.0);

    let pp = p ^ point(6.0, -7.0);

    let center = point(3.0, 4.0);
    let circ = circle(center, 7.0);
    for &s in SCALARS {
        let central_inversion = s * Rotor::from(NI ^ NO);

        let (x, y) = central_inversion.sandwich(p).unpack_point();
        assert_approx_eq(x, -1.0);
        assert_approx_eq(y, -4.0);

        let [(x1, y1), (x2, y2)] = (-central_inversion.sandwich(pp))
            .unpack_point_pair()
            .unwrap()
            .map(|p| p.unpack_point());
        assert_approx_eq(x1, -1.0);
        assert_approx_eq(y1, -4.0);
        assert_approx_eq(x2, -6.0);
        assert_approx_eq(y2, 7.0);

        match central_inversion.sandwich(circ).unpack(EPS) {
            LineOrCircle::Line { .. } => {
                panic!("expected circle")
            }
            LineOrCircle::Circle { cx, cy, r } => {
                assert_approx_eq(cx, -3.0);
                assert_approx_eq(cy, -4.0);
                assert_approx_eq(r, 7.0);
            }
        }
    }
}

#[test]
fn test_tangent_vector_rotate() {
    let p = NO;

    let line_x = point(1.0, 0.0) ^ NI ^ NO;
    let tangent_vector_x = p << line_x;
    let tangent_vector_y = tangent_vector_x.rotate(std::f64::consts::FRAC_PI_2);

    let line_y = point(0.0, 1.0) ^ NI ^ NO;
    let expected_tangent_vector_y = p << line_y;

    assert_approx_eq(tangent_vector_y, expected_tangent_vector_y);
}

#[test]
fn test_rotate() {
    let r = rotate(std::f64::consts::FRAC_PI_2);
    assert_eq_up_to_scale(r.sandwich(point(1.0, -3.0)), point(3.0, 1.0));
}

#[test]
fn test_scale() {
    let r = scale(13.0);
    assert_eq_up_to_scale(r.sandwich(point(1.0, -3.0)), point(13.0, -39.0));
    assert_approx_eq(r.mag2(), 1.0);
}

#[test]
fn test_rotoflector() {
    let rot = Rotoflector::from(rotate(std::f64::consts::FRAC_PI_2)); // rotate 90 degree counterclockwise
    let refl = Rotoflector::from(line(1.0, 0.0, 1.0)); // reflect over x=1
    assert_eq_up_to_scale((rot * refl).sandwich(point(5.0, -2.0)), point(2.0, -3.0));
}
