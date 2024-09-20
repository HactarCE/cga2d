use super::*;

#[test]
fn test_point() {
    let p = point(3.0, -4.0);
    assert_eq!(p.unpack_point(), (3.0, -4.0));
}

#[test]
fn test_point_pair() {
    let (x1, y1) = (1.0, 2.0);
    let (x2, y2) = (3.0, -4.0);

    let p1 = point(x1, y1);
    let p2 = point(x2, y2);
    let [out1, out2] = (p1 ^ p2).unpack_point_pair();
    assert_eq!(out1.unpack_point(), (x1, y1));
    assert_eq!(out2.unpack_point(), (x2, y2));
}
