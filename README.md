# cga2d

![Crates.io Version](https://img.shields.io/crates/v/cga2d) ![docs.rs](https://img.shields.io/docsrs/cga2d)

`cga2d` is a library for 2D Conformal Geometric Algebra with static types for various objects. It has traits for `Multivector` and `Blade` and types for blades of each grade, in addition to `Rotor`, `Flector`, and `Rotoflector`.

There is currently no general-purpose multivector type with all 16 components, but I am open to adding one in the future if there is use for one.

Read [the documentation][docs] for more details.

The scalar type is `f64`. I'm open to adding support for other scalar types, probably via feature flags. I'd prefer not to use generics, because that would greatly hinder ergonomics.

[docs]: https://docs.rs/cga2d/latest/cga2d/

## Known issues

- Not enough tests! I don't actually know if all the operations are implemented correctly.
