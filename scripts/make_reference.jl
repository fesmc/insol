#!/usr/bin/env julia
#
# Generate reference insolation values for the insol regression tests.
#
# This is an *independent* reimplementation of the orbital solution and the
# daily-mean insolation integral. It deliberately shares no code with the
# Fortran library: agreement between the two is the actual test. The orbital
# element tables (Laskar et al. 2004) are the only shared input.
#
# Daily-mean insolation is evaluated in closed form (Berger 1978),
#
#     Q = S0 (a/r)^2 / pi * (H0 sin(phi) sin(dec) + cos(phi) cos(dec) sin(H0))
#
# with the half-day angle H0 = acos(-tan(phi) tan(dec)), rather than by hourly
# quadrature, so the reference carries no discretization error.
#
# Usage:  julia scripts/make_reference.jl [input_dir] [output_file]

using Printf

const S0_DEFAULT = 1365.0

# Day of the vernal equinox (March 21, 12:00 GMT) counted from the start of the
# year, per calendar. It depends only on the lengths of January and February,
# so `day_year` selects a labelling of the orbit, not a different orbit.
const TVEREX = Dict(360 => 80.5,    # 30 + 30 + 20.5
                    365 => 79.5,    # 31 + 28 + 20.5
                    366 => 80.5)    # 31 + 29 + 20.5

"Orbital element table from Laskar et al. (2004): time in ka, angles in degrees."
struct OrbitTable
    time::Vector{Float64}   # ka relative to 1950
    ecc::Vector{Float64}    # eccentricity
    obl::Vector{Float64}    # obliquity [deg]
    per::Vector{Float64}    # longitude of perihelion [deg]
end

"Parse one row of a LA2004 table (Fortran D-exponent notation)."
function parse_row(line)
    parts = split(line)
    length(parts) < 4 && error("malformed orbital table row: $line")
    return ntuple(i -> parse(Float64, replace(parts[i], 'D' => 'E', 'd' => 'e')), 4)
end

"""
    read_orbit_table(dir)

Read the past (INSOLN, -5 Ma..0) and future (INSOLP, 0..+1 Ma) tables and splice
them into a single record ordered from oldest to youngest. Both files include
t=0; the shared row is stored once.
"""
function read_orbit_table(dir)
    past   = parse_row.(readlines(joinpath(dir, "LA2004.INSOLN.5Ma.txt")))
    future = parse_row.(readlines(joinpath(dir, "LA2004.INSOLP.1Ma.txt")))

    # `past` runs 0 -> -5000 ka, so reverse it; drop its t=0 row in favour of
    # the identical first row of `future`.
    rows = vcat(reverse(past)[1:end-1], future)

    time = [r[1] for r in rows]
    ecc  = [r[2] for r in rows]
    obl  = [rad2deg(r[3]) for r in rows]
    per  = [rad2deg(r[4]) for r in rows]

    issorted(time) || error("spliced orbital table is not monotonic in time")
    return OrbitTable(time, ecc, obl, per)
end

"""
    orbital_par(tab, time_bp)

Interpolate the orbital elements to `time_bp` (years before 1950) and derive the
time of perihelion. Returns `(ecc, obl, per, tper, zavexpe)` with angles in
degrees except `zavexpe` (radians) and `tper` (days from start of year).
"""
function orbital_par(tab::OrbitTable, time_bp, day_year)
    t_ka = time_bp * 1e-3
    (t_ka < first(tab.time) || t_ka > last(tab.time)) &&
        error("time_bp = $time_bp ka outside orbital table range")

    n2 = searchsortedfirst(tab.time, t_ka)
    n2 = clamp(n2 == 1 ? 2 : n2, 2, length(tab.time))
    n1 = n2 - 1

    w = (t_ka - tab.time[n1]) / (tab.time[n2] - tab.time[n1])

    ecc = (1 - w) * tab.ecc[n1] + w * tab.ecc[n2]
    obl = (1 - w) * tab.obl[n1] + w * tab.obl[n2]

    # Perihelion is an angle: unwrap the branch cut before interpolating.
    p1, p2 = tab.per[n1], tab.per[n2]
    p2 < p1 && (p2 += 360.0)
    per = mod((1 - w) * p1 + w * p2, 360.0)

    # Time of perihelion, following Monin (1986): convert the perihelion angle
    # measured from the autumnal equinox into an eccentric then mean anomaly.
    ang     = mod(per - 180.0, 360.0)
    zavexpe = deg2rad(ang)
    ecc_anom = 2 * atan(sqrt((1 - ecc) / (1 + ecc)) * tan(zavexpe / 2))
    mean_anom = ecc_anom - ecc * sin(ecc_anom)          # Kepler's equation
    days = mod(mean_anom/(2pi)*day_year, day_year)
    tper = mod(days + TVEREX[day_year], day_year)

    return ecc, obl, per, tper, zavexpe
end

"""
    solar_geometry(day, ecc, obl, tper, zavexpe)

Solar declination [rad] and the inverse-square distance factor (a/r)^2 for a
given day of the year. Declination is held fixed across the day, matching the
convention of the library.
"""
function solar_geometry(day, ecc, obl, tper, zavexpe, day_year)
    ytime = day / day_year * 2pi

    # Solve Kepler's equation for the eccentric anomaly by Newton iteration.
    m = ytime - 2pi * tper / day_year
    e = m
    for _ in 1:100
        de = (m + ecc * (sin(e) - e * cos(e))) / (1 - ecc * cos(e)) - e
        e += de
        abs(de) < 1e-14 && break
    end

    dist2 = (1 / (1 - ecc * cos(e)))^2                  # (a/r)^2

    true_anom = 2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(e / 2))
    lambda    = true_anom + zavexpe                      # true longitude
    decl      = asin(sin(deg2rad(obl)) * sin(lambda))

    return decl, dist2
end

"""
    insol_day(lat, decl, dist2; S0)

Daily-mean insolation [W/m2] at latitude `lat` [deg], in closed form.
"""
function insol_day(lat, decl, dist2; S0 = S0_DEFAULT)
    sinlat, coslat = sincos(deg2rad(lat))
    sindec, cosdec = sincos(decl)

    denom = coslat * cosdec
    cosH0 = abs(denom) < 1e-14 ? -sign(sinlat * sindec) : -sinlat * sindec / denom

    H0 = cosH0 >= 1 ? 0.0 :          # polar night
         cosH0 <= -1 ? pi :          # polar day
         acos(cosH0)

    return S0 * dist2 / pi * (H0 * sinlat * sindec + coslat * cosdec * sin(H0))
end

"""
    coszm_day(lat, decl, dist2)

Insolation-weighted daily-mean cosine of the solar zenith angle, i.e.
`<cos^2 z> / <cos z>` over the sunlit part of the day. Returns 0 in polar night.
"""
function coszm_day(lat, decl, dist2)
    sinlat, coslat = sincos(deg2rad(lat))
    sindec, cosdec = sincos(decl)

    a = sinlat * sindec
    b = coslat * cosdec

    cosH0 = abs(b) < 1e-14 ? -sign(a) : -a / b
    H0 = cosH0 >= 1 ? 0.0 : cosH0 <= -1 ? pi : acos(cosH0)
    H0 == 0 && return 0.0

    # integral of cos(z) and cos(z)^2 over [-H0, H0], with cos(z) = a + b cos(w)
    i1 = 2 * (a * H0 + b * sin(H0))
    i2 = 2 * ((a^2 + b^2 / 2) * H0 + 2 * a * b * sin(H0) + b^2 / 4 * sin(2 * H0))

    return i1 <= 0 ? 0.0 : i2 / i1
end

function main()
    input_dir = length(ARGS) >= 1 ? ARGS[1] : "input"
    outfile   = length(ARGS) >= 2 ? ARGS[2] : "tests/reference/insol_reference.txt"

    tab = read_orbit_table(input_dir)

    # Reference points chosen to exercise: both hemispheres, the equator, the
    # poles (polar day and polar night), solstices and equinoxes, and times
    # spanning the full range of the orbital tables.
    times = [0.0, -6000.0, -21000.0, -115000.0, -1.0e6, -4.9e6, 5.0e5, 1.0e6]
    lats  = [-90.0, -65.0, -45.0, -23.5, 0.0, 23.5, 45.0, 65.0, 90.0]

    # Both supported calendar lengths, so the reference also pins down the
    # vernal equinox convention for each.
    daysets = Dict(360 => [1, 46, 80, 105, 172, 200, 264, 300, 355, 360],
                   365 => [1, 46, 79, 105, 172, 200, 264, 300, 355, 365])

    open(outfile, "w") do io
        println(io, "# Reference daily-mean insolation and mean cosine zenith angle.")
        println(io, "# Generated by scripts/make_reference.jl -- an independent")
        println(io, "# implementation of the Berger (1978) analytic daily integral.")
        println(io, "# Orbital elements: Laskar et al. (2004), tables in input/.")
        println(io, "#")
        @printf(io, "# S0 = %.4f W/m2\n", S0_DEFAULT)
        println(io, "#")
        @printf(io, "#%11s %9s %10s %6s %22s %22s\n",
                "time_bp", "day_year", "lat", "day", "insol", "coszm")

        for day_year in sort(collect(keys(daysets)))
            for time_bp in times
                ecc, obl, per, tper, zavexpe = orbital_par(tab, time_bp, day_year)
                for day in daysets[day_year]
                    decl, dist2 = solar_geometry(day, ecc, obl, tper, zavexpe, day_year)
                    for lat in lats
                        q  = insol_day(lat, decl, dist2)
                        cz = coszm_day(lat, decl, dist2)
                        @printf(io, "%12.1f %9d %10.4f %6d %22.14e %22.14e\n",
                                time_bp, day_year, lat, day, q, cz)
                    end
                end
            end
        end
    end

    n = length(times) * sum(length(d) for d in values(daysets)) * length(lats)
    println("Wrote $n reference values to $outfile")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
