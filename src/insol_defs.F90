module insol_defs
    ! Precision kinds and physical constants shared by the insol core.
    !
    ! All internal calculation is done in double precision (dp). The orbital
    ! solution solves Kepler's equation by Newton iteration and differences
    ! angles accumulated over millions of years, so single precision is not
    ! adequate there and is not offered.
    !
    ! Only the *interface* is generic: `wp` sets the kind of the insolation
    ! values exchanged with the caller, so a single-precision host model need
    ! not convert at every call site. Build with -DINSOL_SP to select single
    ! precision; the default is double.
    !
    ! Note: irradiance.f90 deliberately does not use this module. It carries its
    ! own kinds and constants (including a different solar constant, 1370 rather
    ! than 1365 W/m2) inherited from the renewable-energy work it was written
    ! for. Unifying the two would silently change its results, so the two halves
    ! of the library are kept separate until that side is revisited.

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: sp = kind(1.0)

#ifdef INSOL_SP
    integer, parameter :: wp = sp
#else
    integer, parameter :: wp = dp
#endif

    ! Mathematical constants
    real(dp), parameter :: pi = 2.0_dp*acos(0.0_dp)
    real(dp), parameter :: degrees_to_radians = pi/180.0_dp
    real(dp), parameter :: radians_to_degrees = 180.0_dp/pi

    ! Default solar constant [W/m2]
    real(dp), parameter :: S0_default = 1365.0_dp

    ! Length of the year in days. The orbital solution fixes the vernal equinox
    ! at day 80.5 of a 360-day tropical year (see calc_orbital_par), so this is
    ! the only value the current solution supports.
    integer,  parameter :: day_year_default = 360

    private
    public :: sp, dp, wp
    public :: pi, degrees_to_radians, radians_to_degrees
    public :: S0_default, day_year_default

end module insol_defs
