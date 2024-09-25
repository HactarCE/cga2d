# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- `rotate(angle: Scalar) -> Rotor`
- `scale(factor: Scalar) -> Rotor`

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
