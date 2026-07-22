# Changelog

Notable changes to `insol`. The library is not versioned yet, so changes are
grouped by the refactor that introduced them.

## Unreleased

Substantial refactor. The API is not backwards compatible and results change —
see [Migrating](#migrating) below.

### Fixed

- **Daily-mean insolation was systematically wrong by up to 2.5 W/m².** The
  daily integral was evaluated by summing 24 hourly samples, which converges
  poorly because `max(cos z, 0)` has a kink at sunrise and sunset. It is now
  evaluated in closed form (Berger 1978), which is exact.
- **The array interfaces were wrong by a further 1.8–2.1 W/m².** They evaluated
  insolation at 61 fixed latitudes and cubic-splined to the requested ones, so
  they disagreed with the point interface for the same latitude.
- Together these biased global annual mean insolation low by 0.19 W/m²
  (341.105 against the required `S0/(4√(1-e²))` = 341.298).
- **`day_year` was accepted but ignored** below the interface: `ORBIT` fixed the
  tropical year at 360 days and the time of perihelion was derived in those
  units, so passing 365 shifted the seasons relative to the orbit and returned
  plausible but wrong values. 360, 365 and 366 are now supported properly and
  anything else is refused at initialization.
- Reading past the end of the orbital tables was possible for times beyond
  +1 Ma; the range is now checked.
- Error paths exited with status 0, so a caller checking the exit status saw a
  fatal error as success. They now use `error stop 1`.
- `calc_insol_day_internal` divided by a hardcoded `24.0` while the number of
  samples was a named parameter.
- The point and array interfaces disagreed on clipping negative values.

### Added

- `insol_class`, holding the orbital tables and optionally a latitude grid.
- `insol_init` / `insol_end`.
- Regression test suite (`make check`): 1440 reference points across both
  calendars, cross-validated against `scripts/make_reference.jl`, an
  independent Julia implementation sharing no code with the library; plus
  interface consistency, non-negativity and bounds, polar night, the equinox
  convention, and global energy balance.
- CI over three builds: default, `debug=1` with `-fcheck=all`, and
  `-DINSOL_SP`.
- `wp` compile-time precision switch for the interface (`-DINSOL_SP` selects
  single). All internal calculation stays double.

### Changed

- Build moved to [configme](https://github.com/fesmc/configme); sources to
  `src/`, tests to `tests/`, build products to `libinsol/`.
- The library no longer holds mutable module state, so it is safe to call from
  several threads.
- Bisection rather than a linear scan over the 6001-entry orbital table.
- Passing a fixed latitude grid to `insol_init` precomputes its geometry once.
  Against the previous implementation, per-day calls are ~21× faster at a
  single latitude, ~7× at 181, and ~2.8× at 2257 — while being exact.

### Removed

- The latitude spline and `interp1D.f90`, with no replacement: the analytic
  kernel is faster than the interpolation it existed to avoid.
- `sinsol_orbit.f`, a Fortran 77 duplicate of the same orbital solution.
- The vendored `ncio.f90` (3779 lines); the irradiance tests take it from
  fesm-utils.
- `test_data/` (6.9 MB), read only by the removed netCDF test program.
- `input/orbit_par_5000.txt`, unreferenced.

### Migrating

The old calls were functions taking latitudes per call:

```fortran
insol = calc_insol_day(day,lats,time_bp,fldr=fldr)
```

They are now subroutines taking a state object initialized once:

```fortran
type(insol_class) :: is
real(wp) :: insol(size(lats))

call insol_init(is,lats,fldr=fldr)      ! once
call calc_insol_day(is,insol,day,time_bp)
```

`S0` and `day_year` move from per-call arguments to `insol_init`. Point queries
keep an explicit latitude and need no grid:

```fortran
call insol_init(is)
call calc_insol_day(is,insol,day,lat,time_bp)
```

Any warm-up call made only to trigger the old lazy table loading (the state was
previously a module global initialized on first use) can be deleted —
`insol_init` does that explicitly.

Expect insolation values to shift by up to ~2.5 W/m² relative to previous
results. This is a correction, but it will move long paleo simulations.

## Earlier

Developed from 2013 onward as a component of paleoclimate simulations, later
extended with `irradiance.f90` for renewable energy work. See the git history.
