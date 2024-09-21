use std::fmt;

use super::*;

const EPS: Scalar = 0.0001;

#[track_caller]
fn assert_approx_eq<T: fmt::Debug + approx::AbsDiffEq<Epsilon = Scalar>>(a: T, b: T) {
    approx::assert_abs_diff_eq!(a, b, epsilon = EPS);
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
        let ((x1, y1), (x2, y2)) = ((1.0, 2.0), (3.0, -4.0));

        let pp = s * point(x1, y1) ^ point(x2, y2);
        assert!(!pp.is_flat(EPS));

        let [mut out1, mut out2] = pp.unpack_point_pair().unwrap();
        if s < 0.0 {
            std::mem::swap(&mut out1, &mut out2);
        }

        let (x1_out, y1_out) = out1.unpack_point();
        let (x2_out, y2_out) = out2.unpack_point();
        assert_approx_eq(x1, x1_out);
        assert_approx_eq(y1, y1_out);
        assert_approx_eq(x2, x2_out);
        assert_approx_eq(y2, y2_out);
    }

    assert!((NO ^ NI).is_flat(EPS));
    assert!((point(1.0, 2.0) ^ NI).is_flat(EPS));
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

        let [(x1, y1), (x2, y2)] = central_inversion
            .sandwich(pp)
            .unpack_point_pair()
            .unwrap()
            .map(|p| p.unpack_point());
        assert_approx_eq(x1, 1.0);
        assert_approx_eq(y1, 4.0);
        assert_approx_eq(x2, 6.0);
        assert_approx_eq(y2, -7.0);

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
