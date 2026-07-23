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

    type insol_class
        ! The orbital solution, shared with the rest of the library.
        type(orbit_class) :: orb

        ! Configuration, fixed at initialization
        real(dp) :: S0                       ! Solar constant [W/m2]

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

contains

! =====================================================================================================
! Initialization

    subroutine insol_init_nogrid(is,fldr,S0,day_year)
        ! Load the orbital tables without a latitude grid. The resulting state
        ! supports point queries only.

        implicit none

        type(insol_class), intent(INOUT) :: is
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year

        call insol_config(is,fldr,S0,day_year)

        is%grid  = .FALSE.
        is%n_lat = 0
        is%n1    = 0
        is%n2    = 0

        return

    end subroutine insol_init_nogrid

    subroutine insol_init_1D(is,lats,fldr,S0,day_year)
        ! Load the orbital tables and precompute the geometry of a 1D latitude grid.

        implicit none

        type(insol_class), intent(INOUT) :: is
        real(dp),          intent(IN)    :: lats(:)
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year

        call insol_config(is,fldr,S0,day_year)
        call insol_set_grid(is,reshape(lats,[size(lats)]),size(lats),0)

        return

    end subroutine insol_init_1D

    subroutine insol_init_2D(is,lats,fldr,S0,day_year)
        ! Load the orbital tables and precompute the geometry of a 2D latitude grid.

        implicit none

        type(insol_class), intent(INOUT) :: is
        real(dp),          intent(IN)    :: lats(:,:)
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year

        call insol_config(is,fldr,S0,day_year)
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
                               elev_min=elev1,tday=tday1,coszm=coszm1)

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
                                   elev_min=elev_min,tday=tday,coszm=coszm)
        else
            call calc_insol_kernel(insol,is%sinlat,is%coslat,decl,dist2,is%S0, &
                                   tday=tday,coszm=coszm)
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
                                   elev_min=elev1D,tday=tday1D,coszm=coszm1D)
        else
            call calc_insol_kernel(insol1D,is%sinlat,is%coslat,decl,dist2,is%S0, &
                                   tday=tday1D,coszm=coszm1D)
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

    subroutine insol_config(is,fldr,S0,day_year)
        ! Apply configuration defaults and load the orbital tables.

        implicit none

        type(insol_class), intent(INOUT) :: is
        character(len=*),  intent(IN), optional :: fldr
        real(dp),          intent(IN), optional :: S0
        integer,           intent(IN), optional :: day_year

        is%S0 = S0_default
        if (present(S0)) is%S0 = S0

        call orbit_init(is%orb,fldr,day_year)

        return

    end subroutine insol_config

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

    subroutine calc_insol_kernel(insol,sinlat,coslat,decl,dist2,S0,elev_min,tday,coszm)
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
        ! elev_min absent is elev_min = 0: the plain sunrise-to-sunset daily mean.

        implicit none

        real(wp), intent(OUT) :: insol(:)
        real(dp), intent(IN), contiguous :: sinlat(:)
        real(dp), intent(IN), contiguous :: coslat(:)
        real(dp), intent(IN)  :: decl, dist2, S0
        real(dp), intent(IN),  optional, contiguous :: elev_min(:)
        real(wp), intent(OUT), optional :: tday(:)
        real(wp), intent(OUT), optional :: coszm(:)

        ! Local variables
        integer  :: j
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
                    i2 = (a*a + 0.5_dp*b*b)*H0 + 2.0_dp*a*b*sinH0 &
                         + 0.25_dp*b*b*sin(2.0_dp*H0)
                    coszm(j) = real(i2/i1,wp)
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

