use super::*;

#[test]
fn test_blade1_unpack() {
    for p in points() {
        for &s in SCALES {
            let blade = p.to_blade();
            let unpacked = (blade * s).unpack().unwrap();

            // Point -> Blade1 -> Point
            assert_approx_eq(p, unpacked);

            // Blade1 -> Point -> Blade1
            assert_mvec_eq(blade, unpacked.to_blade());
        }
    }
}

#[test]
fn test_blade2_unpack() {
    for pp in dipoles() {
        let blade = pp.to_blade();

        let expected_blade = match pp {
            Dipole::Real([p, q]) => p.to_blade() ^ q.to_blade(),
            Dipole::Tangent(p, [vx, vy]) => match p {
                Point::Finite([px, py]) => {
                    point(px, py) << (point(px, py) ^ point(px - vx, py - vy) ^ NI)
                }
                Point::Infinity => NI << (NO ^ vector(vx, vy) ^ NI),
            },
            Dipole::Imaginary([p, q]) => (p.to_blade() ^ q.to_blade()).antidual(),
        };
        assert_mvec_eq(blade, expected_blade);

        for &s in SCALES {
            let unpacked = (blade * s).unpack();

            // PointPair -> Blade2 -> PointPair
            assert_approx_eq((Orientation::from(s) * pp).normalize(), unpacked);

            // Blade2 -> PointPair -> Blade2
            assert_mvec_eq(Orientation::from(s) * blade, unpacked.to_blade());
        }
    }
}

#[test]
fn test_blade3_unpack() {
    for c in circles() {
        let blade = c.to_blade();

        let expected_blade = match c {
            Circle::Line { a, b, c } => {
                let scale = c / (a * a + b * b);
                let px = a * scale;
                let py = b * scale;
                crate::point(px, py) ^ crate::point(px + b, py - a) ^ NI
            }
            Circle::Circle { cx, cy, r, ori } => {
                if r == 0.0 {
                    ori * !crate::point(cx, cy)
                } else {
                    let obj = crate::point(cx + r, cy)
                        ^ crate::point(cx, cy + ori * r)
                        ^ crate::point(cx - r, cy);
                    if r > 0.0 {
                        obj
                    } else {
                        ori * !(!obj ^ NI) ^ !obj
                    }
                }
            }
            Circle::Infinity(ori) => ori * !NI,
        };
        assert_mvec_eq(blade, expected_blade);

        for &s in SCALES {
            let unpacked = (blade * s).unpack();

            // Circle -> Blade3 -> Circle
            assert_approx_eq((Orientation::from(s) * c).normalize(), unpacked);

            // Blade3 -> Circle -> Blade3
            assert_mvec_eq(Orientation::from(s) * blade, unpacked.to_blade());
        }
    }
}

#[test]
fn test_rescale_unoriented() {
    let m = circle(point(2.0, 3.0), 4.0);

    // Ignores sign
    let m1 = m.rescale_unoriented();
    let m2 = (-m).rescale_unoriented();
    assert_eq!(m1, m2);

    // Idempotent
    assert!(APPROX.eq(m1, m1.rescale_unoriented()));

    // Equivalent to original (with optional sign flip)
    let pos_m = m.rescale_oriented();
    let neg_m = -m.rescale_oriented();
    assert!(APPROX.eq(pos_m, m1) || APPROX.eq(neg_m, m1));
    assert!(APPROX.eq(pos_m, m2) || APPROX.eq(neg_m, m2));
}

#[test]
fn test_rescale_zeros() {
    let zero = Scalar::zero();
    assert_eq!(zero, zero.rescale_oriented());
    assert_eq!(zero, zero.rescale_unoriented());

    let zero = Blade1::zero();
    assert_eq!(zero, zero.rescale_oriented());
    assert_eq!(zero, zero.rescale_unoriented());

    let zero = Blade2::zero();
    assert_eq!(zero, zero.rescale_oriented());
    assert_eq!(zero, zero.rescale_unoriented());

    let zero = Blade3::zero();
    assert_eq!(zero, zero.rescale_oriented());
    assert_eq!(zero, zero.rescale_unoriented());

    let zero = Pseudoscalar::zero();
    assert_eq!(zero, zero.rescale_oriented());
    assert_eq!(zero, zero.rescale_unoriented());
}

#[test]
fn test_is_flat() {
    for &s in SCALES {
        for p in points() {
            let expected_is_flat = p == Point::Infinity;
            assert_eq!((s * p.to_blade()).is_flat(), expected_is_flat);
        }

        for pp in dipoles() {
            let expected_is_flat = match pp {
                Dipole::Real(points) => points.contains(&Point::Infinity),
                Dipole::Tangent(point, _) => point == Point::Infinity,
                Dipole::Imaginary(_) => false,
            };
            assert_eq!((s * pp.to_blade()).is_flat(), expected_is_flat);
        }

        for c in circles() {
            let expected_is_flat = match c {
                Circle::Line { .. } | Circle::Infinity(_) => true,
                Circle::Circle { .. } => false,
            };
            assert_eq!((s * c.to_blade()).is_flat(), expected_is_flat);
        }
    }
}

#[test]
fn test_dipole_constructor() {
    for (p, q) in point_pairs() {
        let Point::Finite([px, py]) = p else {
            continue;
        };
        let Point::Finite([qx, qy]) = q else {
            continue;
        };
        let center = [(px + qx) / 2.0, (py + qy) / 2.0];
        let [dx, dy] = [qx - px, qy - py];
        let radius = (dx * dx + dy * dy).sqrt() / 2.0;
        let blade = crate::dipole(Point::Finite(center), [dx, dy], radius);
        assert_mvec_eq(p.to_blade() ^ q.to_blade(), blade);
    }
}
