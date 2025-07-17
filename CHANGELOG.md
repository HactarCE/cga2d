# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Wedge products with `Scalar`
- Antiwedge products with `Pseudoscalar`
- `Blade::is_flat_with_prec()`
- `impl approx_collections::ApproxSign for Pseudoscalar`

### Changed

- Switched from `approx` to `approx_collections`
- Overhauled blade unpacking
  - There are now two methods on `Blade1`, `Blade2`, and `Blade3`: `unpack()` and `unpack_with_prec(prec: Precision)`
  - Unpacking a `Blade1` now returns a `Point` enum instead of `(Scalar, Scalar)`
  - Unpacking a `Blade2` now returns a `PointPair` enum instead of `Option<Blade1>`
  - Unpacking a `Blade3` now returns a `Circle` (renamed from `LineOrCircle`)
  - Added conversions between `Blade1`/`Blade2`/`Blade3` and `Point`/`PointPair`/`Circle` respectively
- Removed `epsilon` parameter from `Blade::is_flat()`
- Added `Blade::is_flat_with_prec(prec: Precision)`

## 0.4.0 - 2024-09-25

### Added

- `rotate(angle: Scalar) -> Rotor`
- `scale(factor: Scalar) -> Rotor`

### Fixed

- `Rotoflector` operations having wrong scalar component

## 0.3.0 - 2024-09-24

### Added

- Multiplication between rotors/flectors and blades

## 0.2.2 - 2024-09-22

### Added

- `approx::AbsDiffEq<Epsilon = Scalar>` bound on `Multivector`

### Fixed

- `.unpack_point_pair()` with flat points (one point is `NI`)

## 0.2.1 - 2024-09-22

### Fixed

- Stack overflow in `Rotor::ident()` and `Rotor::default()`

## 0.2.0 - 2024-09-21

### Added

- `bytemuck` support via feature flag

## 0.1.0 - 2024-09-21

- Initial release
