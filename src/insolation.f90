! The MIT License (MIT)
!
! Copyright (c) 2014 Alexander Robinson and Mahé Perrette
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module insolation

    use insol_defs, only : sp, dp, wp, pi, deg_to_rad => degrees_to_radians, &
                           S0_default, day_year_default
    use insol_orbit, only : orbit_class, orbit_init, orbit_end, &
                            calc_orbital_par, calc_solar_geometry

    implicit none

    ! How the mean cosine of the solar zenith angle is defined. The two differ
    ! only in weighting: "flux" weights each instant by the insolation it carries
    ! (the dEBM convention, exact via I2/I1); "daytime" is the unweighted mean
    ! over the part of the day the sun is up (a + b sinH0/H0, the convention
    ! CLIMBER-X's albedo schemes expect). Fixed once per state at initialization.
    integer, parameter :: COSZM_KIND_FLUX    = 0
    integer, parameter :: COSZM_KIND_DAYTIME = 1

    type insol_class
        ! The orbital solution, shared with the rest of the library.
        type(orbit_class) :: orb

        ! Configuration, fixed at initialization
        real(dp) :: S0                       ! Solar constant [W/m2]
        integer  :: coszm_kind               ! COSZM_KIND_FLUX or COSZM_KIND_DAYTIME

        ! Optional latitude grid. When one is given, sin and cos of every
        ! latitude are precomputed here rather than recomputed on every call --
        ! for a fixed model grid this is the dominant cost of the daily kernel.
        ! Stored flat, with the original shape kept so output can be reshaped.
        logical  :: grid
        integer  :: n_lat                    ! Total number of grid points
        integer  :: n1, n2                   ! Grid shape (n2 = 0 for a 1D grid)
        real(dp), allocatable :: sinlat(:)
        real(dp), allocatable :: coslat(:)
    end type

    ! A state is initialized either bare -- valid only for point queries -- or
    ! with a 1D or 2D latitude grid, which additionally enables the grid
    ! routines below.
    interface insol_init
        module procedure insol_init_nogrid
        module procedure insol_init_1D
        module procedure insol_init_2D
    end interface

    ! Daily-mean insolation. Resolved on the rank of `insol`: rank 1 or 2 with
    ! no latitude argument uses the grid held in `is`; rank 0 with an explicit
    ! latitude is a one-off point query needing no grid.
    interface calc_insol_day
        module procedure calc_insol_day_1D
        module procedure calc_insol_day_2D
        module procedure calc_insol_day_pt
    end interface

    ! Daily-mean insolation for several days at once. `insol` carries one
    ! extra trailing dimension for the days.
    interface calc_insol_days
        module procedure calc_insol_days_1D
        module procedure calc_insol_days_2D
    end interface

    ! Insolation averaged over a range of days (eg. a season or the year).
    interface calc_insol_ave
        module procedure calc_insol_ave_1D
        module procedure calc_insol_ave_2D
    end interface

    ! Sub-daily (diurnal) insolation and zenith angle over the 1D grid, sampled
    ! at n_hour points across the day, with the daily mean/min/max and the mean
    ! cosine returned alongside. `calc_insol_diurnal` does one day; the `_days`
    ! form fills a whole set of days in one call.
    interface calc_insol_diurnal
        module procedure calc_insol_diurnal_1D
        module procedure calc_insol_diurnal_2D
        module procedure calc_insol_diurnal_pt
    end interface

    interface calc_insol_diurnal_days
        module procedure calc_insol_diurnal_days_1D
        module procedure calc_insol_diurnal_days_2D
    end interface

    ! Length of the daylit part of the day, in hours.
    interface calc_daylength
        module procedure calc_daylength_1D
        module procedure calc_daylength_2D
        module procedure calc_daylength_pt
    end interface

    private
    ! Re-exported so a caller needs only `use insolation` to declare variables
    ! of the kind the interface exchanges.
    public :: sp, dp, wp
    public :: insol_class
    public :: insol_init
    public :: insol_end
    public :: calc_insol_day
    public :: calc_insol_days
    public :: calc_insol_ave
    public :: calc_insol_diurnal
    public :: calc_insol_diurnal_days
    public :: calc_daylength

contains

! =====================================================================================================
! Initialization

    subroutine insol_init_nogrid(is,fldr,S0,day_year,coszm_kind)
        ! Load the orbital tables without a latitude grid. The resulting state
        ! supports point queries only.

        implicit none

        type(insol_class), intent(INOUT) :: is
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year
        character(len=*),  intent(IN), optional :: coszm_kind

        call insol_config(is,fldr,S0,day_year,coszm_kind)

        is%grid  = .FALSE.
        is%n_lat = 0
        is%n1    = 0
        is%n2    = 0

        return

    end subroutine insol_init_nogrid

    subroutine insol_init_1D(is,lats,fldr,S0,day_year,coszm_kind)
        ! Load the orbital tables and precompute the geometry of a 1D latitude grid.

        implicit none

        type(insol_class), intent(INOUT) :: is
        real(dp),          intent(IN)    :: lats(:)
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year
        character(len=*),  intent(IN), optional :: coszm_kind

        call insol_config(is,fldr,S0,day_year,coszm_kind)
        call insol_set_grid(is,reshape(lats,[size(lats)]),size(lats),0)

        return

    end subroutine insol_init_1D

    subroutine insol_init_2D(is,lats,fldr,S0,day_year,coszm_kind)
        ! Load the orbital tables and precompute the geometry of a 2D latitude grid.

        implicit none

        type(insol_class), intent(INOUT) :: is
        real(dp),          intent(IN)    :: lats(:,:)
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year
        character(len=*),  intent(IN), optional :: coszm_kind

        call insol_config(is,fldr,S0,day_year,coszm_kind)
        call insol_set_grid(is,reshape(lats,[size(lats)]),size(lats,1),size(lats,2))

        return

    end subroutine insol_init_2D

    subroutine insol_end(is)
        ! Release the orbital tables and any cached grid geometry.

        implicit none

        type(insol_class), intent(INOUT) :: is

        call orbit_end(is%orb)

        if (allocated(is%sinlat)) deallocate(is%sinlat)
        if (allocated(is%coslat)) deallocate(is%coslat)

        is%grid  = .FALSE.
        is%n_lat = 0

        return

    end subroutine insol_end

! =====================================================================================================
! Daily-mean insolation

    subroutine calc_insol_day_pt(is,insol,day,lat,time_bp,elev_min,tday,coszm)
        ! Daily-mean insolation at one latitude. Needs no latitude grid, so it
        ! works with a state initialized by insol_init(is) alone.
        ! time_bp = years before present (present == 1950)
        !
        ! `insol` is the mean over the whole 24-hour day. With `elev_min` given,
        ! only the part of the day when the sun stands higher than `elev_min`
        ! degrees above the horizon contributes -- what diurnal melt schemes such
        ! as dEBM (Krebs-Kanzow et al., 2018) need. `elev_min = 0` reproduces the
        ! plain daily mean exactly. The optional outputs describe that window:
        ! `tday` is the fraction of the day the sun is above `elev_min`, and
        ! `coszm` the flux-weighted mean cosine of the zenith angle over it. The
        ! mean flux *during* the window (dEBM's SW_Phi) is then insol / tday.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: lat
        real(dp),          intent(IN)  :: time_bp
        real(dp),          intent(IN),  optional :: elev_min   ! [degrees]
        real(wp),          intent(OUT), optional :: tday       ! fraction of day above elev_min
        real(wp),          intent(OUT), optional :: coszm      ! flux-weighted mean cos(zenith)

        ! Local variables. One point is treated as a grid of length one; the
        ! named arrays (rather than array-constructor actuals) keep the compiler
        ! from making a temporary for the contiguous kernel dummy.
        real(dp) :: decl, dist2, latrad
        real(dp) :: sinlat1(1), coslat1(1), elev1(1)
        real(wp) :: insol1(1), tday1(1), coszm1(1)

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        latrad     = lat*deg_to_rad
        sinlat1(1) = sin(latrad)
        coslat1(1) = cos(latrad)
        elev1(1)   = 0.0_dp
        if (present(elev_min)) elev1(1) = elev_min

        call calc_insol_kernel(insol1,sinlat1,coslat1,decl,dist2,is%S0, &
                               elev_min=elev1,tday=tday1,coszm=coszm1, &
                               coszm_kind=is%coszm_kind)

        insol = insol1(1)
        if (present(tday))  tday  = tday1(1)
        if (present(coszm)) coszm = coszm1(1)

        return

    end subroutine calc_insol_day_pt

    subroutine calc_insol_day_1D(is,insol,day,time_bp,elev_min,tday,coszm)
        ! Daily-mean insolation over the 1D latitude grid held in `is`. See
        ! calc_insol_day_pt for the meaning of the optional arguments; here
        ! `elev_min` is one threshold per grid point, as dEBM estimates it
        ! locally.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:)
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp
        real(dp),          intent(IN),  optional :: elev_min(:)
        real(wp),          intent(OUT), optional :: tday(:)
        real(wp),          intent(OUT), optional :: coszm(:)

        ! Local variables
        real(dp) :: decl, dist2

        call check_grid(is,size(insol),1,"calc_insol_day")
        if (present(elev_min)) call check_elev_min(size(elev_min),is%n_lat,"calc_insol_day")

        ! The solar geometry does not depend on latitude, so it is evaluated
        ! once and reused for every grid point.
        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        ! Branch rather than forward the optionals: an absent optional passed on
        ! to another optional dummy is a gfortran hazard, and it keeps the common
        ! call on the fast, argument-free path.
        if (present(elev_min)) then
            call calc_insol_kernel(insol,is%sinlat,is%coslat,decl,dist2,is%S0, &
                                   elev_min=elev_min,tday=tday,coszm=coszm, &
                                   coszm_kind=is%coszm_kind)
        else
            call calc_insol_kernel(insol,is%sinlat,is%coslat,decl,dist2,is%S0, &
                                   tday=tday,coszm=coszm,coszm_kind=is%coszm_kind)
        end if

        return

    end subroutine calc_insol_day_1D

    subroutine calc_insol_day_2D(is,insol,day,time_bp,elev_min,tday,coszm)
        ! Daily-mean insolation over the 2D latitude grid held in `is`. See
        ! calc_insol_day_pt for the optional arguments.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:,:)
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp
        real(dp),          intent(IN),  optional :: elev_min(:,:)
        real(wp),          intent(OUT), optional :: tday(:,:)
        real(wp),          intent(OUT), optional :: coszm(:,:)

        ! Local variables
        real(dp) :: decl, dist2
        real(dp) :: elev1D(is%n_lat)
        real(wp) :: insol1D(is%n_lat), tday1D(is%n_lat), coszm1D(is%n_lat)

        call check_grid(is,size(insol),2,"calc_insol_day")
        if (present(elev_min)) call check_elev_min(size(elev_min),is%n_lat,"calc_insol_day")

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        ! Flatten the 2D optionals to the grid vector, call once, reshape back.
        ! Passing an absent optional straight through keeps the fast path.
        if (present(elev_min)) then
            elev1D = reshape(elev_min,[is%n_lat])
            call calc_insol_kernel(insol1D,is%sinlat,is%coslat,decl,dist2,is%S0, &
                                   elev_min=elev1D,tday=tday1D,coszm=coszm1D, &
                                   coszm_kind=is%coszm_kind)
        else
            call calc_insol_kernel(insol1D,is%sinlat,is%coslat,decl,dist2,is%S0, &
                                   tday=tday1D,coszm=coszm1D,coszm_kind=is%coszm_kind)
        end if

        insol = reshape(insol1D,[is%n1,is%n2])
        if (present(tday))  tday  = reshape(tday1D, [is%n1,is%n2])
        if (present(coszm)) coszm = reshape(coszm1D,[is%n1,is%n2])

        return

    end subroutine calc_insol_day_2D

! =====================================================================================================
! Daily-mean insolation over several days

    subroutine calc_insol_days_1D(is,insol,days,time_bp)
        ! Daily-mean insolation over the 1D latitude grid, for each day in `days`.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:,:)     ! [lat,day]
        integer,           intent(IN)  :: days(:)
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        integer  :: q
        real(dp) :: decl, dist2

        call check_grid(is,size(insol,1),1,"calc_insol_days")
        call check_ndays(size(insol,2),size(days),"calc_insol_days")

        do q = 1, size(days)
            call calc_solar_geometry_day(is,decl,dist2,days(q),time_bp)
            call calc_insol_grid(is,insol(:,q),decl,dist2)
        end do

        return

    end subroutine calc_insol_days_1D

    subroutine calc_insol_days_2D(is,insol,days,time_bp)
        ! Daily-mean insolation over the 2D latitude grid, for each day in `days`.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:,:,:)   ! [x,y,day]
        integer,           intent(IN)  :: days(:)
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        integer  :: q
        real(dp) :: decl, dist2
        real(wp) :: insol1D(is%n_lat)

        call check_grid(is,size(insol,1)*size(insol,2),2,"calc_insol_days")
        call check_ndays(size(insol,3),size(days),"calc_insol_days")

        do q = 1, size(days)
            call calc_solar_geometry_day(is,decl,dist2,days(q),time_bp)
            call calc_insol_grid(is,insol1D,decl,dist2)
            insol(:,:,q) = reshape(insol1D,[is%n1,is%n2])
        end do

        return

    end subroutine calc_insol_days_2D

! =====================================================================================================
! Insolation averaged over a range of days

    subroutine calc_insol_ave_1D(is,insol,day_range,time_bp,step)
        ! Mean insolation over an inclusive range of days, on the 1D grid.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:)
        integer,           intent(IN)  :: day_range(2)
        real(dp),          intent(IN)  :: time_bp
        integer,           intent(IN), optional :: step

        ! Local variables
        integer  :: day, ndays, day_step
        real(dp) :: decl, dist2
        real(wp) :: insol_day(is%n_lat)

        call check_grid(is,size(insol),1,"calc_insol_ave")
        call check_day_range(day_range,"calc_insol_ave")

        day_step = 1
        if (present(step)) day_step = step

        insol = 0.0_wp
        ndays = 0

        do day = day_range(1), day_range(2), day_step
            call calc_solar_geometry_day(is,decl,dist2,day,time_bp)
            call calc_insol_grid(is,insol_day,decl,dist2)
            insol = insol + insol_day
            ndays = ndays + 1
        end do

        insol = insol / real(ndays,wp)

        return

    end subroutine calc_insol_ave_1D

    subroutine calc_insol_ave_2D(is,insol,day_range,time_bp,step)
        ! Mean insolation over an inclusive range of days, on the 2D grid.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:,:)
        integer,           intent(IN)  :: day_range(2)
        real(dp),          intent(IN)  :: time_bp
        integer,           intent(IN), optional :: step

        ! Local variables
        integer  :: day, ndays, day_step
        real(dp) :: decl, dist2
        real(wp) :: insol1D(is%n_lat), insol_day(is%n_lat)

        call check_grid(is,size(insol),2,"calc_insol_ave")
        call check_day_range(day_range,"calc_insol_ave")

        day_step = 1
        if (present(step)) day_step = step

        insol1D = 0.0_wp
        ndays   = 0

        do day = day_range(1), day_range(2), day_step
            call calc_solar_geometry_day(is,decl,dist2,day,time_bp)
            call calc_insol_grid(is,insol_day,decl,dist2)
            insol1D = insol1D + insol_day
            ndays   = ndays + 1
        end do

        insol = reshape(insol1D,[is%n1,is%n2]) / real(ndays,wp)

        return

    end subroutine calc_insol_ave_2D

! =====================================================================================================
! Sub-daily (diurnal) insolation

    subroutine calc_insol_diurnal_core(is,day,time_bp,sinlat,coslat, &
                                       solar,cosz,solarm,solarmin,solarmax,coszm)
        ! Shared worker: fills every diurnal and aggregate output for one day over
        ! a flat latitude vector. The public wrappers reshape their point or 2D
        ! arguments onto this vector, so all the geometry lives here once. All
        ! arguments are mandatory, so the wrappers forward slices without tripping
        ! the absent-optional hazard.
        !
        ! The diurnal cycle is sampled at n_hour = size(solar,2) equal slots over
        ! the day. The instantaneous flux onto a horizontal surface is
        ! S0 (a/r)^2 max(cos z, 0) with cos z = a + b cos(h) -- the same geometry
        ! the closed-form daily mean integrates -- so the samples and the mean are
        ! mutually consistent. The daily mean and mean cosine are taken from that
        ! exact closed form rather than the sampling, and so do not depend on
        ! n_hour; only min and max, which have no closed form, come from the
        ! samples.
        !
        ! Layout is [lat, hour] -- latitude leading (contiguous, the inner loop)
        ! and hour trailing -- so the wrappers can reshape the latitude extent to
        ! a 2D grid and add a day dimension without a transpose.

        implicit none

        type(insol_class), intent(IN)  :: is
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp
        real(dp),          intent(IN),  contiguous :: sinlat(:)     ! [lat]
        real(dp),          intent(IN),  contiguous :: coslat(:)     ! [lat]
        real(wp),          intent(OUT) :: solar(:,:)     ! [lat,hour] insolation [W/m2]
        real(wp),          intent(OUT) :: cosz(:,:)      ! [lat,hour] cos(zenith) []
        real(wp),          intent(OUT) :: solarm(:)      ! [lat] daily-mean insolation
        real(wp),          intent(OUT) :: solarmin(:)    ! [lat] daily-min insolation
        real(wp),          intent(OUT) :: solarmax(:)    ! [lat] daily-max insolation
        real(wp),          intent(OUT) :: coszm(:)       ! [lat] mean cos(zenith)

        ! Local variables
        integer  :: j, k, nh, nlat
        real(dp) :: decl, dist2, sindec, cosdec, s0d, cosH, hourang, cz
        real(dp) :: a(size(sinlat)), b(size(sinlat))

        nlat = size(sinlat)
        nh   = size(solar,2)

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        sindec = sin(decl)
        cosdec = cos(decl)
        s0d    = is%S0*dist2

        do j = 1, nlat
            a(j) = sinlat(j)*sindec
            b(j) = coslat(j)*cosdec
            solarmin(j) =  huge(1.0_wp)
            solarmax(j) = -huge(1.0_wp)
        end do

        do k = 1, nh
            ! Hour angle of sample k: nh equal slots across the day, with midnight
            ! at -pi and local solar noon at 0, sampled at slot centres so no
            ! sample lands exactly on the discontinuity. Independent of latitude,
            ! so it is evaluated once per hour, outside the inner loop.
            hourang = 2.0_dp*pi*(real(k,dp)-0.5_dp)/real(nh,dp) - pi
            cosH    = cos(hourang)
            do j = 1, nlat
                cz = max(a(j) + b(j)*cosH, 0.0_dp)
                cosz(j,k)  = real(cz,wp)
                solar(j,k) = real(s0d*cz,wp)
                if (solar(j,k) < solarmin(j)) solarmin(j) = solar(j,k)
                if (solar(j,k) > solarmax(j)) solarmax(j) = solar(j,k)
            end do
        end do

        call calc_insol_kernel(solarm,sinlat,coslat,decl,dist2,is%S0, &
                               coszm=coszm,coszm_kind=is%coszm_kind)

        return

    end subroutine calc_insol_diurnal_core

    subroutine calc_insol_diurnal_pt(is,day,lat,time_bp,solar,cosz, &
                                     solarm,solarmin,solarmax,coszm)
        ! Diurnal insolation and cosine of the zenith angle at one latitude. Needs
        ! no grid, so it works on a bare state. `solar` and `cosz` are [n_hour];
        ! the aggregates are scalars. The mean cosine follows the state's
        ! coszm_kind.

        implicit none

        type(insol_class), intent(IN)  :: is
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: lat
        real(dp),          intent(IN)  :: time_bp
        real(wp),          intent(OUT) :: solar(:)       ! [hour]
        real(wp),          intent(OUT) :: cosz(:)        ! [hour]
        real(wp),          intent(OUT), optional :: solarm
        real(wp),          intent(OUT), optional :: solarmin
        real(wp),          intent(OUT), optional :: solarmax
        real(wp),          intent(OUT), optional :: coszm

        ! Local variables. One point is a grid of length one.
        integer  :: nh
        real(dp) :: latrad
        real(dp) :: sinlat1(1), coslat1(1)
        real(wp) :: solarm_l(1), solarmin_l(1), solarmax_l(1), coszm_l(1)
        real(wp), allocatable :: solar2(:,:), cosz2(:,:)

        call check_diurnal(size(solar),size(cosz),1,1,"calc_insol_diurnal")

        nh = size(solar)
        allocate(solar2(1,nh),cosz2(1,nh))

        latrad     = lat*deg_to_rad
        sinlat1(1) = sin(latrad)
        coslat1(1) = cos(latrad)

        call calc_insol_diurnal_core(is,day,time_bp,sinlat1,coslat1,solar2,cosz2, &
                                     solarm_l,solarmin_l,solarmax_l,coszm_l)

        solar = solar2(1,:)
        cosz  = cosz2(1,:)
        if (present(solarm))   solarm   = solarm_l(1)
        if (present(solarmin)) solarmin = solarmin_l(1)
        if (present(solarmax)) solarmax = solarmax_l(1)
        if (present(coszm))    coszm    = coszm_l(1)

        deallocate(solar2,cosz2)

        return

    end subroutine calc_insol_diurnal_pt

    subroutine calc_insol_diurnal_1D(is,day,time_bp,solar,cosz, &
                                     solarm,solarmin,solarmax,coszm)
        ! Diurnal insolation and cosine of the zenith angle over the 1D grid, for
        ! one day. `solar` and `cosz` are [n_lat, n_hour]; the number of sub-daily
        ! samples is taken from size(solar,2). The daily mean, min, max and mean
        ! cosine are optional. The mean cosine follows the state's coszm_kind.

        implicit none

        type(insol_class), intent(IN)  :: is
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp
        real(wp),          intent(OUT) :: solar(:,:)
        real(wp),          intent(OUT) :: cosz(:,:)
        real(wp),          intent(OUT), optional :: solarm(:)
        real(wp),          intent(OUT), optional :: solarmin(:)
        real(wp),          intent(OUT), optional :: solarmax(:)
        real(wp),          intent(OUT), optional :: coszm(:)

        ! Local variables
        real(wp) :: solarm_l(is%n_lat), solarmin_l(is%n_lat)
        real(wp) :: solarmax_l(is%n_lat), coszm_l(is%n_lat)

        call check_grid(is,size(solar,1),1,"calc_insol_diurnal")
        call check_diurnal(size(solar,2),size(cosz,2),size(solar,1),size(cosz,1), &
                           "calc_insol_diurnal")

        call calc_insol_diurnal_core(is,day,time_bp,is%sinlat,is%coslat,solar,cosz, &
                                     solarm_l,solarmin_l,solarmax_l,coszm_l)

        if (present(solarm))   solarm   = solarm_l
        if (present(solarmin)) solarmin = solarmin_l
        if (present(solarmax)) solarmax = solarmax_l
        if (present(coszm))    coszm    = coszm_l

        return

    end subroutine calc_insol_diurnal_1D

    subroutine calc_insol_diurnal_2D(is,day,time_bp,solar,cosz, &
                                     solarm,solarmin,solarmax,coszm)
        ! Diurnal insolation and cosine of the zenith angle over the 2D grid, for
        ! one day. `solar` and `cosz` are [n1, n2, n_hour]; the aggregates are
        ! [n1, n2]. The grid is flattened to the internal latitude vector and the
        ! results reshaped back, so the physics is identical to the 1D path.

        implicit none

        type(insol_class), intent(IN)  :: is
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp
        real(wp),          intent(OUT) :: solar(:,:,:)
        real(wp),          intent(OUT) :: cosz(:,:,:)
        real(wp),          intent(OUT), optional :: solarm(:,:)
        real(wp),          intent(OUT), optional :: solarmin(:,:)
        real(wp),          intent(OUT), optional :: solarmax(:,:)
        real(wp),          intent(OUT), optional :: coszm(:,:)

        ! Local variables
        integer  :: nh
        real(wp) :: solarm_l(is%n_lat), solarmin_l(is%n_lat)
        real(wp) :: solarmax_l(is%n_lat), coszm_l(is%n_lat)
        real(wp), allocatable :: solar2(:,:), cosz2(:,:)

        call check_grid(is,size(solar,1)*size(solar,2),2,"calc_insol_diurnal")
        call check_diurnal(size(solar,3),size(cosz,3), &
                           size(solar,1)*size(solar,2),size(cosz,1)*size(cosz,2), &
                           "calc_insol_diurnal")

        nh = size(solar,3)
        allocate(solar2(is%n_lat,nh),cosz2(is%n_lat,nh))

        call calc_insol_diurnal_core(is,day,time_bp,is%sinlat,is%coslat,solar2,cosz2, &
                                     solarm_l,solarmin_l,solarmax_l,coszm_l)

        solar = reshape(solar2,[is%n1,is%n2,nh])
        cosz  = reshape(cosz2, [is%n1,is%n2,nh])
        if (present(solarm))   solarm   = reshape(solarm_l,  [is%n1,is%n2])
        if (present(solarmin)) solarmin = reshape(solarmin_l,[is%n1,is%n2])
        if (present(solarmax)) solarmax = reshape(solarmax_l,[is%n1,is%n2])
        if (present(coszm))    coszm    = reshape(coszm_l,   [is%n1,is%n2])

        deallocate(solar2,cosz2)

        return

    end subroutine calc_insol_diurnal_2D

    subroutine calc_insol_diurnal_days_1D(is,days,time_bp,solar,cosz, &
                                          solarm,solarmin,solarmax,coszm)
        ! As calc_insol_diurnal, for each day in `days` at once. `solar` and
        ! `cosz` are [n_lat, n_days, n_hour]; the optional aggregates are
        ! [n_lat, n_days].

        implicit none

        type(insol_class), intent(IN)  :: is
        integer,           intent(IN)  :: days(:)
        real(dp),          intent(IN)  :: time_bp
        real(wp),          intent(OUT) :: solar(:,:,:)
        real(wp),          intent(OUT) :: cosz(:,:,:)
        real(wp),          intent(OUT), optional :: solarm(:,:)
        real(wp),          intent(OUT), optional :: solarmin(:,:)
        real(wp),          intent(OUT), optional :: solarmax(:,:)
        real(wp),          intent(OUT), optional :: coszm(:,:)

        ! Local variables
        integer  :: q, nh
        real(wp) :: solarm_l(is%n_lat), solarmin_l(is%n_lat)
        real(wp) :: solarmax_l(is%n_lat), coszm_l(is%n_lat)
        real(wp), allocatable :: solar2(:,:), cosz2(:,:)

        call check_ndays(size(solar,2),size(days),"calc_insol_diurnal_days")
        call check_grid(is,size(solar,1),1,"calc_insol_diurnal_days")
        call check_diurnal(size(solar,3),size(cosz,3),size(solar,1),size(cosz,1), &
                           "calc_insol_diurnal_days")

        nh = size(solar,3)
        allocate(solar2(is%n_lat,nh),cosz2(is%n_lat,nh))

        do q = 1, size(days)
            call calc_insol_diurnal_core(is,days(q),time_bp,is%sinlat,is%coslat, &
                                         solar2,cosz2, &
                                         solarm_l,solarmin_l,solarmax_l,coszm_l)
            solar(:,q,:) = solar2
            cosz(:,q,:)  = cosz2
            if (present(solarm))   solarm(:,q)   = solarm_l
            if (present(solarmin)) solarmin(:,q) = solarmin_l
            if (present(solarmax)) solarmax(:,q) = solarmax_l
            if (present(coszm))    coszm(:,q)    = coszm_l
        end do

        deallocate(solar2,cosz2)

        return

    end subroutine calc_insol_diurnal_days_1D

    subroutine calc_insol_diurnal_days_2D(is,days,time_bp,solar,cosz, &
                                          solarm,solarmin,solarmax,coszm)
        ! As calc_insol_diurnal on the 2D grid, for each day in `days`. `solar`
        ! and `cosz` are [n1, n2, n_days, n_hour]; the aggregates are
        ! [n1, n2, n_days].

        implicit none

        type(insol_class), intent(IN)  :: is
        integer,           intent(IN)  :: days(:)
        real(dp),          intent(IN)  :: time_bp
        real(wp),          intent(OUT) :: solar(:,:,:,:)
        real(wp),          intent(OUT) :: cosz(:,:,:,:)
        real(wp),          intent(OUT), optional :: solarm(:,:,:)
        real(wp),          intent(OUT), optional :: solarmin(:,:,:)
        real(wp),          intent(OUT), optional :: solarmax(:,:,:)
        real(wp),          intent(OUT), optional :: coszm(:,:,:)

        ! Local variables
        integer  :: q, nh
        real(wp) :: solarm_l(is%n_lat), solarmin_l(is%n_lat)
        real(wp) :: solarmax_l(is%n_lat), coszm_l(is%n_lat)
        real(wp), allocatable :: solar2(:,:), cosz2(:,:)

        call check_ndays(size(solar,3),size(days),"calc_insol_diurnal_days")
        call check_grid(is,size(solar,1)*size(solar,2),2,"calc_insol_diurnal_days")
        call check_diurnal(size(solar,4),size(cosz,4), &
                           size(solar,1)*size(solar,2),size(cosz,1)*size(cosz,2), &
                           "calc_insol_diurnal_days")

        nh = size(solar,4)
        allocate(solar2(is%n_lat,nh),cosz2(is%n_lat,nh))

        do q = 1, size(days)
            call calc_insol_diurnal_core(is,days(q),time_bp,is%sinlat,is%coslat, &
                                         solar2,cosz2, &
                                         solarm_l,solarmin_l,solarmax_l,coszm_l)
            solar(:,:,q,:) = reshape(solar2,[is%n1,is%n2,nh])
            cosz(:,:,q,:)  = reshape(cosz2, [is%n1,is%n2,nh])
            if (present(solarm))   solarm(:,:,q)   = reshape(solarm_l,  [is%n1,is%n2])
            if (present(solarmin)) solarmin(:,:,q) = reshape(solarmin_l,[is%n1,is%n2])
            if (present(solarmax)) solarmax(:,:,q) = reshape(solarmax_l,[is%n1,is%n2])
            if (present(coszm))    coszm(:,:,q)    = reshape(coszm_l,   [is%n1,is%n2])
        end do

        deallocate(solar2,cosz2)

        return

    end subroutine calc_insol_diurnal_days_2D

! =====================================================================================================
! Daylength

    subroutine calc_daylength_pt(is,daylength,day,lat,time_bp)
        ! Hours the sun is above the horizon at one latitude. This is the plain
        ! daily window (elevation floor zero), so daylength = 24 * tday.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: daylength   ! [hours]
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: lat
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        real(dp) :: decl, dist2, latrad
        real(dp) :: sinlat1(1), coslat1(1)
        real(wp) :: insol1(1), tday1(1)

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        latrad     = lat*deg_to_rad
        sinlat1(1) = sin(latrad)
        coslat1(1) = cos(latrad)

        call calc_insol_kernel(insol1,sinlat1,coslat1,decl,dist2,is%S0,tday=tday1)

        daylength = 24.0_wp*tday1(1)

        return

    end subroutine calc_daylength_pt

    subroutine calc_daylength_1D(is,daylength,day,time_bp)
        ! Daylength over the 1D grid.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: daylength(:)   ! [hours]
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        real(dp) :: decl, dist2
        real(wp) :: insol_l(is%n_lat)

        call check_grid(is,size(daylength),1,"calc_daylength")

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)
        call calc_insol_kernel(insol_l,is%sinlat,is%coslat,decl,dist2,is%S0, &
                               tday=daylength)

        daylength = 24.0_wp*daylength

        return

    end subroutine calc_daylength_1D

    subroutine calc_daylength_2D(is,daylength,day,time_bp)
        ! Daylength over the 2D grid.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: daylength(:,:)   ! [hours]
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        real(dp) :: decl, dist2
        real(wp) :: insol_l(is%n_lat), tday1D(is%n_lat)

        call check_grid(is,size(daylength),2,"calc_daylength")

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)
        call calc_insol_kernel(insol_l,is%sinlat,is%coslat,decl,dist2,is%S0, &
                               tday=tday1D)

        daylength = reshape(24.0_wp*tday1D,[is%n1,is%n2])

        return

    end subroutine calc_daylength_2D

    subroutine check_diurnal(nh_solar,nh_cosz,nlat_solar,nlat_cosz,routine)
        ! Validate the sub-daily arrays of a diurnal call: at least one hour
        ! sample, and cosz matching solar in both the hour and latitude extents.
        ! Grid membership is checked separately by the caller via check_grid.

        implicit none

        integer,          intent(IN) :: nh_solar, nh_cosz, nlat_solar, nlat_cosz
        character(len=*), intent(IN) :: routine

        if (nh_solar < 1) then
            write(*,*) trim(routine)//":: Error: need at least one hour sample."
            error stop 1
        end if

        if (nh_cosz /= nh_solar .or. nlat_cosz /= nlat_solar) then
            write(*,*) trim(routine)//":: Error: cosz and solar shapes differ."
            error stop 1
        end if

        return

    end subroutine check_diurnal

    subroutine check_elev_min(n_given,n_lat,routine)

        implicit none

        integer,          intent(IN) :: n_given, n_lat
        character(len=*), intent(IN) :: routine

        if (n_given /= n_lat) then
            write(*,*) trim(routine)//":: Error: elev_min does not match the grid."
            write(*,*) "grid points:    ", n_lat
            write(*,*) "elev_min size:  ", n_given
            error stop 1
        end if

        return

    end subroutine check_elev_min

! =====================================================================================================
! Internal helpers

    subroutine insol_config(is,fldr,S0,day_year,coszm_kind)
        ! Apply configuration defaults and load the orbital tables. The mean
        ! cosine convention is resolved from its string once, here, so the hot
        ! path stores only the integer code.

        implicit none

        type(insol_class), intent(INOUT) :: is
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year
        character(len=*),  intent(IN), optional :: coszm_kind

        is%S0 = S0_default
        if (present(S0)) is%S0 = S0

        is%coszm_kind = COSZM_KIND_FLUX
        if (present(coszm_kind)) is%coszm_kind = resolve_coszm_kind(coszm_kind)

        call orbit_init(is%orb,fldr,day_year)

        return

    end subroutine insol_config

    function resolve_coszm_kind(name) result(kind)
        ! Map the public string onto the internal code, case-insensitively.
        ! An unrecognised value is a programming error and stops the run.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: kind

        ! Local variables
        integer :: i, c
        character(len=len(name)) :: lname

        lname = name
        do i = 1, len(lname)
            c = iachar(lname(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) lname(i:i) = achar(c+32)
        end do

        select case (trim(adjustl(lname)))
            case ("flux")
                kind = COSZM_KIND_FLUX
            case ("daytime")
                kind = COSZM_KIND_DAYTIME
            case default
                write(*,*) "insol_init:: Error: unknown coszm_kind '"//trim(name)//"'."
                write(*,*) 'Expected "flux" (flux-weighted) or "daytime" (daytime mean).'
                error stop 1
        end select

    end function resolve_coszm_kind

    subroutine insol_set_grid(is,lats,n1,n2)
        ! Precompute and store the geometry of a latitude grid. `lats` arrives
        ! flattened; n2 = 0 marks a 1D grid.

        implicit none

        type(insol_class), intent(INOUT) :: is
        real(dp),          intent(IN)    :: lats(:)
        integer,           intent(IN)    :: n1, n2

        ! Local variables
        integer  :: j
        real(dp) :: latrad

        if (allocated(is%sinlat)) deallocate(is%sinlat)
        if (allocated(is%coslat)) deallocate(is%coslat)

        is%n_lat = size(lats)
        is%n1    = n1
        is%n2    = n2
        is%grid  = .TRUE.

        allocate(is%sinlat(is%n_lat),is%coslat(is%n_lat))

        do j = 1, is%n_lat
            latrad       = lats(j)*deg_to_rad
            is%sinlat(j) = sin(latrad)
            is%coslat(j) = cos(latrad)
        end do

        return

    end subroutine insol_set_grid

    subroutine calc_insol_grid(is,insol,decl,dist2)
        ! Plain daily mean at every point of the cached grid, used by the
        ! multi-day routines. The elevation-threshold path is handled directly
        ! in calc_insol_day_1D/2D: forwarding absent optional arrays through this
        ! extra layer is a gfortran hazard best avoided.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:)
        real(dp),          intent(IN)  :: decl, dist2

        call calc_insol_kernel(insol,is%sinlat,is%coslat,decl,dist2,is%S0)

        return

    end subroutine calc_insol_grid

    subroutine calc_insol_kernel(insol,sinlat,coslat,decl,dist2,S0,elev_min,tday,coszm,coszm_kind)
        ! The grid loop, taking its arrays directly rather than reaching through
        ! the state object. Declaring them contiguous tells the compiler the
        ! iterations are independent and unit-strided, which it cannot assume for
        ! allocatable components accessed through a derived type.
        !
        ! The daily mean is the closed-form integral of max(cos z, 0) over the
        ! day (Berger 1978). With elev_min present the lower limit becomes the
        ! hour angle at which the sun reaches that elevation, h(Phi):
        !
        !     cos(z) = a + b cos(h),   a = sin(lat)sin(dec), b = cos(lat)cos(dec)
        !     cos h(Phi) = (sin(Phi) - a)/b
        !
        ! and, integrating over |h| < h(Phi),
        !
        !     insol = S0 (a/r)^2 / pi * (a h + b sin h)      (24-hour mean)
        !     tday  = h/pi                                    (fraction of day)
        !     coszm = I2/I1  with  I1 = a h + b sin h,        (flux-weighted
        !             I2 = (a^2 + b^2/2) h + 2ab sin h         mean cos zenith)
        !                  + (b^2/4) sin 2h
        !
        ! With coszm_kind = COSZM_KIND_DAYTIME the mean is instead unweighted over
        ! the same window, coszm = I1/h = a + b sinH0/H0, which is what CLIMBER-X's
        ! albedo schemes expect. Absent, the weighting is flux (I2/I1).
        !
        ! elev_min absent is elev_min = 0: the plain sunrise-to-sunset daily mean.

        implicit none

        real(wp), intent(OUT) :: insol(:)
        real(dp), intent(IN), contiguous :: sinlat(:)
        real(dp), intent(IN), contiguous :: coslat(:)
        real(dp), intent(IN)  :: decl, dist2, S0
        real(dp), intent(IN),  optional, contiguous :: elev_min(:)
        real(wp), intent(OUT), optional :: tday(:)
        real(wp), intent(OUT), optional :: coszm(:)
        integer,  intent(IN),  optional :: coszm_kind

        ! Local variables
        integer  :: j, kind
        real(dp) :: sindec, cosdec, fac
        real(dp) :: a, b, sin_phi, cosH0, H0, sinH0, i1, i2
        logical  :: have_elev, want_tday, want_coszm

        ! Hoist everything that does not vary over the grid out of the loop.
        sindec = sin(decl)
        cosdec = cos(decl)
        fac    = S0*dist2/pi

        have_elev  = present(elev_min)
        want_tday  = present(tday)
        want_coszm = present(coszm)

        kind = COSZM_KIND_FLUX
        if (present(coszm_kind)) kind = coszm_kind

        do j = 1, size(sinlat)

            a = sinlat(j)*sindec
            b = coslat(j)*cosdec

            if (have_elev) then
                sin_phi = sin(elev_min(j)*deg_to_rad)
            else
                sin_phi = 0.0_dp
            end if

            ! Hour angle at which the elevation reaches the threshold. At the
            ! poles b vanishes and the elevation is constant through the day, so
            ! the threshold is either met all day or not at all.
            if (abs(b) .lt. 1.0e-14_dp) then
                cosH0 = sign(1.0_dp,sin_phi-a)
            else
                cosH0 = (sin_phi-a)/b
            end if

            if (cosH0 .ge. 1.0_dp) then
                H0    = 0.0_dp
                sinH0 = 0.0_dp
            else if (cosH0 .le. -1.0_dp) then
                H0    = pi
                sinH0 = 0.0_dp
            else
                H0    = acos(cosH0)
                sinH0 = sqrt(1.0_dp-cosH0*cosH0)
            end if

            i1 = a*H0 + b*sinH0
            insol(j) = real(fac*i1,wp)

            if (want_tday) tday(j) = real(H0/pi,wp)

            if (want_coszm) then
                if (i1 .gt. 0.0_dp) then
                    if (kind == COSZM_KIND_DAYTIME) then
                        ! Unweighted mean over the daylight window: I1/H0.
                        coszm(j) = real(i1/H0,wp)
                    else
                        ! Flux-weighted mean: I2/I1.
                        i2 = (a*a + 0.5_dp*b*b)*H0 + 2.0_dp*a*b*sinH0 &
                             + 0.25_dp*b*b*sin(2.0_dp*H0)
                        coszm(j) = real(i2/i1,wp)
                    end if
                else
                    coszm(j) = 0.0_wp
                end if
            end if

        end do

        return

    end subroutine calc_insol_kernel

    subroutine calc_solar_geometry_day(is,decl,dist2,day,time_bp)
        ! Orbital elements for `time_bp` and the solar geometry for `day`.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(dp),          intent(OUT) :: decl, dist2
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        real(dp) :: PER, ECC, XOBCH, TPERI, ZAVEXPE

        call calc_orbital_par(is%orb,time_bp,PER,ECC,XOBCH,TPERI,ZAVEXPE)
        call calc_solar_geometry(decl,dist2,day,is%orb%day_year,ECC,XOBCH,TPERI,ZAVEXPE)

        return

    end subroutine calc_solar_geometry_day

    subroutine check_grid(is,n_given,rank_given,routine)
        ! Fail loudly when a grid routine is called on a state that has no grid,
        ! or when the output array does not match the grid it was initialized with.

        implicit none

        type(insol_class), intent(IN) :: is
        integer,           intent(IN) :: n_given, rank_given
        character(len=*),  intent(IN) :: routine

        if (.not. is%grid) then
            write(*,*) trim(routine)//":: Error: no latitude grid available."
            write(*,*) "This state was initialized without one. Either call"
            write(*,*) "insol_init(is,lats,...) with the grid, or use the point"
            write(*,*) "interface calc_insol_day(is,insol,day,lat,time_bp)."
            error stop 1
        end if

        if (rank_given == 2 .and. is%n2 == 0) then
            write(*,*) trim(routine)//":: Error: state holds a 1D latitude grid,"
            write(*,*) "but a 2D output array was given."
            error stop 1
        end if

        if (rank_given == 1 .and. is%n2 /= 0) then
            write(*,*) trim(routine)//":: Error: state holds a 2D latitude grid,"
            write(*,*) "but a 1D output array was given."
            error stop 1
        end if

        if (n_given /= is%n_lat) then
            write(*,*) trim(routine)//":: Error: output array does not match the grid."
            write(*,*) "grid points: ", is%n_lat
            write(*,*) "output size: ", n_given
            error stop 1
        end if

        return

    end subroutine check_grid

    subroutine check_ndays(n_given,n_days,routine)

        implicit none

        integer,          intent(IN) :: n_given, n_days
        character(len=*), intent(IN) :: routine

        if (n_given /= n_days) then
            write(*,*) trim(routine)//":: Error: day dimension of the output array"
            write(*,*) "does not match the number of days requested."
            write(*,*) "days requested:   ", n_days
            write(*,*) "output dimension: ", n_given
            error stop 1
        end if

        return

    end subroutine check_ndays

    subroutine check_day_range(day_range,routine)

        implicit none

        integer,          intent(IN) :: day_range(2)
        character(len=*), intent(IN) :: routine

        if (day_range(2) .lt. day_range(1)) then
            write(*,*) trim(routine)//":: Error: day_range must give a positive"
            write(*,*) "total of days to calculate."
            write(*,*) "day_range: ", day_range
            error stop 1
        end if

        return

    end subroutine check_day_range


! =====================================================================================================

end module

