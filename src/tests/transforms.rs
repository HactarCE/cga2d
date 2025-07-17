use super::*;

#[test]
fn test_point_reflection() {
    let p = point(1.0, 4.0);

    let pp = p ^ point(6.0, -7.0);

    let center = point(3.0, 4.0);
    let circ = circle(center, 7.0);
    for &s in SCALES {
        let central_inversion = s * Rotor::from(NI ^ NO);

        let [x, y] = central_inversion
            .sandwich(p)
            .unpack()
            .unwrap()
            .finite()
            .unwrap();
        assert_approx_eq(x, -1.0);
        assert_approx_eq(y, -4.0);

        let [[x1, y1], [x2, y2]] = (-central_inversion.sandwich(pp))
            .unpack()
            .real()
            .unwrap()
            .map(|p| p.finite().unwrap());
        assert_approx_eq(x1, -1.0);
        assert_approx_eq(y1, -4.0);
        assert_approx_eq(x2, -6.0);
        assert_approx_eq(y2, 7.0);

        match central_inversion.sandwich(circ).unpack_with_prec(APPROX) {
            Circle::Circle { cx, cy, r, ori } => {
                assert_approx_eq(cx, -3.0);
                assert_approx_eq(cy, -4.0);
                assert_approx_eq(r, 7.0);
                assert_eq!(ori, Orientation::Pos);
            }
            _ => panic!("expected circle"),
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
    assert_mvec_eq(r.sandwich(point(1.0, -3.0)), point(3.0, 1.0));
}

#[test]
fn test_scale() {
    let r = scale(13.0);
    assert_mvec_eq(r.sandwich(point(1.0, -3.0)), point(13.0, -39.0));
    assert_approx_eq(r.mag2(), 1.0);
}

#[test]
fn test_rotoflector() {
    let rot = Rotoflector::from(rotate(std::f64::consts::FRAC_PI_2)); // rotate 90 degree counterclockwise
    let refl = Rotoflector::from(line(1.0, 0.0, 1.0)); // reflect over x=1
    assert_mvec_eq((rot * refl).sandwich(point(5.0, -2.0)), point(2.0, -3.0));
}
