# Changelog

Notable changes to `insol`, newest first. Versions are annotated git tags.

## Unreleased

Additive: fixed orbits. Backwards compatible — `time_bp` calls are unchanged
and return bit-for-bit the same values.

Build: `libinsol/include/libinsol.a` is a file target of `make static` (was a phony
target that re-archived on every call, so dependants in a host build relinked each time).

### Added

- `orbit_pars_class`, the orbital elements (`ecc`, `obl`, `per`) at one moment,
  with two constructors: `orbit_pars_at(is,time_bp)` reads them from the tables
  (also takes an `orbit_class`), `orbit_pars_fixed(ecc,obl,per)` takes them as
  given, for a fixed orbit.
- Every routine taking `time_bp` — `calc_insol_day`, `calc_insol_days`,
  `calc_insol_ave`, `calc_insol_diurnal`, `calc_insol_diurnal_days`,
  `calc_daylength` — also takes an `orbit_pars_class` in its place. The
  `time_bp` form is `orbit_pars_at` followed by the elements form.

### Changed

- The step from the elements to the time and longitude of perihelion is
  factored out of `calc_orbital_par` into `calc_perihelion`, shared by table
  and fixed orbits. The calendar enters only there, so one set of elements is
  valid with any `day_year`.

## v1.1.0

Additive release: sub-daily and daylength routines, and a selectable mean-cosine
convention. Backwards compatible — existing calls are unchanged, and the default
`coszm_kind` reproduces previous results.

### Added

- `calc_insol_diurnal` / `calc_insol_diurnal_days`: the sub-daily insolation and
  cosine of the zenith angle sampled across the day, with the daily mean, min,
  max and mean cosine returned alongside. The mean and mean cosine come from the
  closed form; only the extremes depend on the sampling. Resolves on the rank of
  `solar` — point, 1D or 2D grid (plus 1D/2D for the multi-day form) — mirroring
  `calc_insol_day` / `calc_insol_days`. Arrays are laid out with the grid
  dimensions first and `n_hour` last.
- `calc_daylength` (point, 1D and 2D): hours the sun is above the horizon.
- A `coszm_kind` option on `insol_init` selecting how the mean cosine of the
  zenith angle is defined — `"flux"` (flux-weighted, the default and the dEBM
  convention) or `"daytime"` (the plain daytime mean, the CLIMBER-X convention).
  Resolved to an integer code once at init, so the hot path carries no string
  comparison.

## v1.0.0

Substantial refactor. The API is not backwards compatible and results change —
see [Migrating](#migrating) below.

### Fixed

- **Irradiance DNI carried a 30-minute timing bug.** The old `calc_hour_angle`
  applied a hardcoded −0.5 h offset, which skewed the zenith most at sunrise and
  sunset where `dni = fdir/cos(z)` is most sensitive. Repointing the geometry at
  `insol_orbit` removed it: DNI RMSE against NREL fell from 171 to 52 W/m².
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
- `insol_orbit`, the orbital solution as a module shared by the insolation core
  and the irradiance side, so both see the same Sun.
- Real-date solar geometry (`calc_solar_position_date`): the same orbit phased
  on a Gregorian date rather than the idealized paleo calendar, for present-day
  climate models and the irradiance decomposition. The equation of time falls
  out of the Kepler solve (validated −14 to +16 min against NREL).
- An elevation threshold on `calc_insol_day` (`elev_min`, `tday`, `coszm`
  optional arguments) for diurnal melt schemes such as dEBM: the daily mean
  restricted to when the sun is above a per-gridpoint elevation. `elev_min = 0`
  reproduces the plain daily mean exactly.
- Regression test suite (`make check`): 1440 reference points across both
  calendars, cross-validated against `scripts/make_reference.jl`, an
  independent Julia implementation sharing no code with the library; plus
  interface consistency, non-negativity and bounds, polar night, the equinox
  convention, and global energy balance.
- Tests for the irradiance side, which previously had none: DISC-estimated DNI
  scored against measured NREL NSRDB station data (~52 W/m² RMSE against a
  measured DNI RMS of ~558), plus a solar-position and Julian-Day suite.
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
- `irradiance` takes its geometry from `insol_orbit` and its `wp`/constants
  from `insol_defs` rather than carrying its own. Its DNI routines now take an
  `orbit_class` and a real date (year + UTC hour of year); the `time_zone` and
  `tisr` arguments are gone. The solar constant unifies on 1365 (from 1370).
- The library now has no external dependencies — no netCDF, no fesm-utils — and
  builds with only a Fortran compiler.

### Removed

- The latitude spline and `interp1D.f90`, with no replacement: the analytic
  kernel is faster than the interpolation it existed to avoid.
- `irradiance`'s own solar geometry: the Spencer (1971) declination and
  equation-of-time series, `calc_hour_angle`, and the duplicated
  top-of-atmosphere series, all superseded by `insol_orbit`.
- `sinsol_orbit.f`, a Fortran 77 duplicate of the same orbital solution.
- The vendored `ncio.f90` (3779 lines) and the `test_dni` ERA5 processing
  driver that used it, along with the fesm-utils dependency they required.
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
