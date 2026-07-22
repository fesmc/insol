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

    use insol_defs, only : dp, wp, pi, deg_to_rad => degrees_to_radians, &
                           S0_default, day_year_default

    implicit none

    type insol_class
        ! Orbital element table (Laskar et al. 2004), loaded once by insol_init
        ! and read-only thereafter.
        integer :: n_orbit
        real(dp), allocatable :: time(:)     ! [ka] relative to 1950
        real(dp), allocatable :: eccm(:)     ! Eccentricity
        real(dp), allocatable :: xobm(:)     ! Obliquity [degrees]
        real(dp), allocatable :: perm(:)     ! Longitude of perihelion [degrees]

        ! Configuration, fixed at initialization
        real(dp) :: S0                       ! Solar constant [W/m2]
        integer  :: day_year                 ! Days in the year

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

        if (allocated(is%time))   deallocate(is%time)
        if (allocated(is%eccm))   deallocate(is%eccm)
        if (allocated(is%xobm))   deallocate(is%xobm)
        if (allocated(is%perm))   deallocate(is%perm)
        if (allocated(is%sinlat)) deallocate(is%sinlat)
        if (allocated(is%coslat)) deallocate(is%coslat)

        is%grid  = .FALSE.
        is%n_lat = 0

        return

    end subroutine insol_end

! =====================================================================================================
! Daily-mean insolation

    subroutine calc_insol_day_pt(is,insol,day,lat,time_bp)
        ! Daily-mean insolation at one latitude. Needs no latitude grid, so it
        ! works with a state initialized by insol_init(is) alone.
        ! time_bp = years before present (present == 1950)

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: lat
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        real(dp) :: decl, dist2

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        insol = real(calc_insol_day_internal(lat,decl,dist2,is%S0),wp)

        return

    end subroutine calc_insol_day_pt

    subroutine calc_insol_day_1D(is,insol,day,time_bp)
        ! Daily-mean insolation over the 1D latitude grid held in `is`.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:)
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        real(dp) :: decl, dist2

        call check_grid(is,size(insol),1,"calc_insol_day")

        ! The solar geometry does not depend on latitude, so it is evaluated
        ! once and reused for every grid point.
        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        call calc_insol_grid(is,insol,decl,dist2)

        return

    end subroutine calc_insol_day_1D

    subroutine calc_insol_day_2D(is,insol,day,time_bp)
        ! Daily-mean insolation over the 2D latitude grid held in `is`.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:,:)
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp

        ! Local variables
        real(dp) :: decl, dist2
        real(wp) :: insol1D(is%n_lat)

        call check_grid(is,size(insol),2,"calc_insol_day")

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        call calc_insol_grid(is,insol1D,decl,dist2)

        insol = reshape(insol1D,[is%n1,is%n2])

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

        is%day_year = day_year_default
        if (present(day_year)) is%day_year = day_year

        call load_orbital_table(is,fldr)

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
        ! Evaluate the daily kernel at every point of the cached grid, reusing
        ! the precomputed sin/cos of each latitude.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:)
        real(dp),          intent(IN)  :: decl, dist2

        call calc_insol_kernel(insol,is%sinlat,is%coslat,decl,dist2,is%S0)

        return

    end subroutine calc_insol_grid

    subroutine calc_insol_kernel(insol,sinlat,coslat,decl,dist2,S0)
        ! The grid loop, taking its arrays directly rather than reaching through
        ! the state object. Declaring them contiguous tells the compiler the
        ! iterations are independent and unit-strided, which it cannot assume for
        ! allocatable components accessed through a derived type.

        implicit none

        real(wp), intent(OUT) :: insol(:)
        real(dp), intent(IN), contiguous :: sinlat(:)
        real(dp), intent(IN), contiguous :: coslat(:)
        real(dp), intent(IN)  :: decl, dist2, S0

        ! Local variables
        integer  :: j
        real(dp) :: sindec, cosdec, fac
        real(dp) :: denom, cosH0, H0, sinH0

        ! Hoist everything that does not vary over the grid out of the loop.
        sindec = sin(decl)
        cosdec = cos(decl)
        fac    = S0*dist2/pi

        do j = 1, size(sinlat)

            denom = coslat(j)*cosdec
            if (abs(denom) .lt. 1.0e-14_dp) then
                cosH0 = sign(1.0_dp,-sinlat(j)*sindec)
            else
                cosH0 = -sinlat(j)*sindec/denom
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

            insol(j) = real(fac*(H0*sinlat(j)*sindec + coslat(j)*cosdec*sinH0),wp)

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

        call calc_orbital_par(is,time_bp,PER,ECC,XOBCH,TPERI,ZAVEXPE)
        call calc_solar_geometry(decl,dist2,day,is%day_year,ECC,XOBCH,TPERI,ZAVEXPE)

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
            stop
        end if

        if (rank_given == 2 .and. is%n2 == 0) then
            write(*,*) trim(routine)//":: Error: state holds a 1D latitude grid,"
            write(*,*) "but a 2D output array was given."
            stop
        end if

        if (rank_given == 1 .and. is%n2 /= 0) then
            write(*,*) trim(routine)//":: Error: state holds a 2D latitude grid,"
            write(*,*) "but a 1D output array was given."
            stop
        end if

        if (n_given /= is%n_lat) then
            write(*,*) trim(routine)//":: Error: output array does not match the grid."
            write(*,*) "grid points: ", is%n_lat
            write(*,*) "output size: ", n_given
            stop
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
            stop
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
            stop
        end if

        return

    end subroutine check_day_range


! =====================================================================================================

    subroutine calc_solar_geometry(decl,dist2,day,day_max,ECC,XOBCH,TPERI,ZAVEXPE)
        ! Solar declination and the inverse-square distance factor (a/r)^2 for a
        ! given day of the year.
        !
        ! Neither quantity depends on the time of day: ORBIT derives the
        ! declination from the time of year alone, and the clock angle enters
        ! only the zenith parameters PZEN2/PZEN3. It is therefore evaluated once
        ! per day rather than once per hour, and reused for every latitude.

        implicit none

        real(dp), intent(OUT) :: decl        ! Solar declination [radians]
        real(dp), intent(OUT) :: dist2       ! (a/r)^2, ratio of solar constant to its annual mean
        integer,  intent(IN)  :: day
        integer,  intent(IN)  :: day_max
        real(dp), intent(IN)  :: ECC, XOBCH, TPERI, ZAVEXPE

        ! Local variables
        real(dp) :: PCLOCK, PYTIME
        real(dp) :: PDISSE, PZEN1, PZEN2, PZEN3, PRAE

        PYTIME = dble(day)/dble(day_max)*2.0_dp*pi
        PCLOCK = 0.0_dp

        call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
                   PCLOCK,PYTIME,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)

        ! PZEN1 = sin(declination); PDISSE = (a/r)^2
        decl  = asin(max(-1.0_dp,min(1.0_dp,PZEN1)))
        dist2 = PDISSE

        return

    end subroutine calc_solar_geometry

    function calc_insol_day_internal(lat,decl,dist2,S0) result(solarm)
        ! Daily-mean insolation on a horizontal surface [W/m2].
        !
        ! The daily mean of the solar zenith angle has a closed form (Berger,
        ! 1978), so the integral is evaluated exactly rather than by summing over
        ! the hours of the day:
        !
        !   solarm = S0 (a/r)^2 / pi * (H0 sin(phi) sin(dec) + cos(phi) cos(dec) sin(H0))
        !
        ! with the half-day angle H0 = acos(-tan(phi) tan(dec)), i.e. the hour
        ! angle at sunrise and sunset. H0 = 0 in polar night and H0 = pi in polar
        ! day, so both are handled without a special case.
        !
        ! Hourly summation was used previously (adapted from CLIMBER2.4 routine
        ! 'sinsol' by A. Ganopolski). Because the integrand is clipped at zero at
        ! sunrise and sunset, it is not smooth there, and 24 samples per day left
        ! an error of up to 2.5 W/m2 -- large enough to bias the global annual
        ! mean below S0/(4 sqrt(1-e^2)). The closed form removes it entirely.

        implicit none

        real(dp), intent(IN) :: lat         ! Latitude [degrees North]
        real(dp), intent(IN) :: decl        ! Solar declination [radians]
        real(dp), intent(IN) :: dist2       ! (a/r)^2
        real(dp), intent(IN) :: S0          ! Solar constant [W/m2]
        real(dp) :: solarm

        ! Local variables
        real(dp) :: latrad

        latrad = lat*deg_to_rad

        solarm = calc_insol_internal(sin(latrad),cos(latrad),decl,dist2,S0)

        return

    end function calc_insol_day_internal

    pure function calc_insol_internal(sinlat,coslat,decl,dist2,S0) result(solarm)
        ! As calc_insol_day_internal, but taking the sine and cosine of the
        ! latitude directly. For a fixed grid these are precomputed once at
        ! initialization, which removes two transcendental evaluations per grid
        ! point per day -- the dominant cost of this kernel.
        !
        implicit none

        real(dp), intent(IN) :: sinlat      ! sin(latitude)
        real(dp), intent(IN) :: coslat      ! cos(latitude)
        real(dp), intent(IN) :: decl        ! Solar declination [radians]
        real(dp), intent(IN) :: dist2       ! (a/r)^2
        real(dp), intent(IN) :: S0          ! Solar constant [W/m2]
        real(dp) :: solarm

        ! Local variables
        real(dp) :: sindec, cosdec
        real(dp) :: denom, cosH0, H0, sinH0

        sindec = sin(decl)
        cosdec = cos(decl)

        ! Hour angle at sunrise/sunset. At the poles cos(phi) = 0 and the
        ! quotient is undefined; the sign of sin(phi)sin(dec) then decides
        ! between polar day and polar night.
        denom = coslat*cosdec
        if (abs(denom) .lt. 1.0e-14_dp) then
            cosH0 = sign(1.0_dp,-sinlat*sindec)
        else
            cosH0 = -sinlat*sindec/denom
        end if

        if (cosH0 .ge. 1.0_dp) then
            H0    = 0.0_dp                  ! Polar night: sun never rises
            sinH0 = 0.0_dp
        else if (cosH0 .le. -1.0_dp) then
            H0    = pi                      ! Polar day: sun never sets
            sinH0 = 0.0_dp
        else
            H0 = acos(cosH0)
            ! sin(H0) >= 0 for H0 in [0,pi], so the positive root is correct.
            ! Cheaper than evaluating sin() a second time.
            sinH0 = sqrt(1.0_dp-cosH0*cosH0)
        end if

        solarm = S0*dist2/pi * (H0*sinlat*sindec + coslat*cosdec*sinH0)

        return

    end function calc_insol_internal

    subroutine calc_orbital_par(is,time,per,ecc,obl,tper,zan)
        !
        ! Calculates Earth's orbital parameters for the current time before 1950 (time)
        !
        ! Adapted from Fortran77 subroutine BERGOR written by  
        ! S. J. Lorenz, University of Bremen, 08-03-94.
        !
        ! Orbital parameters are loaded from a precalculated table of orbital values
        !
        implicit none

        type(insol_class), intent(IN)  :: is
        real(dp),          intent(IN)  :: time
        real(dp),          intent(OUT) :: per, ecc, obl, tper, zan

        integer :: N1, N2, NMID
        real (dp) :: time_ka, ef
        real (dp) :: perm1, perm2
        real (dp) :: ang, tgnu, sqecc, tian
        real (dp) :: days, e, epc

        ! Get current time in ka
        time_ka = time*1d-3

        ! The tables span -5 Ma to +1 Ma; extrapolating beyond them is
        ! meaningless, and reading past the end used to run off the array.
        if (time_ka .lt. is%time(1) .or. time_ka .gt. is%time(is%n_orbit)) then
            write(*,*) "calc_orbital_par:: Error: time_bp outside the range of the"
            write(*,*) "orbital parameter tables."
            write(*,*) "time_bp [yr]:   ", time
            write(*,*) "table range [yr]:", is%time(1)*1d3, is%time(is%n_orbit)*1d3
            stop
        end if

        ! Determine the indices bracketing the current time by bisection. The
        ! table holds several thousand entries, so a linear scan from the start
        ! costs more than the interpolation it feeds.
        N1 = 1
        N2 = is%n_orbit
        do while (N2-N1 .gt. 1)
            NMID = (N1+N2)/2
            if (is%time(NMID) .gt. time_ka) then
                N2 = NMID
            else
                N1 = NMID
            end if
        end do

        ! Get weighting of each index
        ef = (time_ka-is%time(N1)) / (is%time(N2)-is%time(N1))

        ! Get current eccentricity, obliquity and perihelion
        ecc = (1.0_dp-ef)*is%eccm(N1) + ef*is%eccm(N2)
        obl = (1.0_dp-ef)*is%xobm(N1) + ef*is%xobm(N2)

        perm1 = is%perm(N1)
        perm2 = is%perm(N2)
        if (perm2 .lt. perm1) perm2 = perm2 + 360.0_dp
        
        per = (1.0_dp-ef)*perm1 + ef*perm2
        if (per .gt. 360.0_dp) per = per - 360.0_dp

        !  calculate time of perihelion in year:
        !  MONIN, A.S.: An Introduction to the Theory of Climate (Reidel Publ.)
        !   - PER: perihelion from aut. equinox in degrees (BERGER)
        !   - ECC: eccentricity of Earth's orbit (BERGER)

        ang = per - 180.0_dp
        if (ang .lt. 0.0_dp) ang = ang+360.0_dp
        zan   = ang / 180.0_dp*pi                !   =zavexpe
        tgnu  = tan(zan/2.0_dp)
        sqecc = dsqrt((1.0_dp-ecc)/(1.0_dp+ecc))
        e     = 2.0_dp*atan(sqecc*tgnu)          !   Eccentri!   Anomaly

        !   time angle in radians of perihelion from vernal equinox

        tian = e - ecc*sin(e)
        days = tian*180.0_dp/pi                    !   360 day year only
        if (days .lt. 0.0_dp) days = days+360.0_dp !   days from ver.eq. to perh.

        !   time in days from begin of year: vernal equinox fixed at 3/21, 12 GMT
        !    = 80.5 days in 360 day year

        tper = days+80.5_dp
        if (tper .gt. 360.0_dp) tper = tper-360.0_dp

        epc = ecc*100.0_dp
      
        return 

    end subroutine calc_orbital_par 

    subroutine load_orbital_table(is,fldr)
        ! Load the orbital parameter tables (Laskar et al. 2004) covering
        ! -5 Ma to +1 Ma into the state.
        !
        ! The two files are spliced into one record running from oldest to
        ! youngest. Both begin at t = 0, so the past table is read in reverse
        ! and its t = 0 row is overwritten by the identical first row of the
        ! future table.
        !
        ! Files are read from subdirectory 'input' unless 'fldr' is given.

        implicit none

        type(insol_class), intent(INOUT) :: is
        character(len=*),  intent(IN), optional :: fldr

        ! Local variables
        integer, parameter :: n_past   = 5001      ! rows in LA2004.INSOLN.5Ma.txt
        integer, parameter :: n_future = 1001      ! rows in LA2004.INSOLP.1Ma.txt

        character(len=512) :: input_dir
        character(len=512) :: filen, filep
        integer :: n, unit, io

        ! Define the input directory where orbital parameter files are located
        input_dir = "input"
        if (present(fldr)) input_dir = trim(fldr)

        ! First make sure all vectors are deallocated
        if (allocated(is%time)) deallocate(is%time)
        if (allocated(is%eccm)) deallocate(is%eccm)
        if (allocated(is%perm)) deallocate(is%perm)
        if (allocated(is%xobm)) deallocate(is%xobm)

        ! Allocate vectors to length of input parameter files
        is%n_orbit = n_past + n_future - 1
        allocate(is%time(is%n_orbit),is%eccm(is%n_orbit))
        allocate(is%perm(is%n_orbit),is%xobm(is%n_orbit))

        ! Specify file names of negative and positive time orbital parameters
        filen = trim(input_dir)//'/LA2004.INSOLN.5Ma.txt'
        filep = trim(input_dir)//'/LA2004.INSOLP.1Ma.txt'

        ! ==== PAST (NEGATIVE) TIMES ====
        open(newunit=unit,file=trim(filen),status="old",action="read",iostat=io)
        if (io /= 0) then
            write(*,*) "load_orbital_table:: Error: cannot open "//trim(filen)
            write(*,*) "Pass the location of the orbital tables as insol_init(...,fldr=)."
            stop
        end if

        do n = 1, n_past
            read(unit,*) is%time(n_past-n+1), is%eccm(n_past-n+1), &
                         is%xobm(n_past-n+1), is%perm(n_past-n+1)
        end do

        close(unit)

        ! ==== FUTURE (POSITIVE) TIMES ====
        open(newunit=unit,file=trim(filep),status="old",action="read",iostat=io)
        if (io /= 0) then
            write(*,*) "load_orbital_table:: Error: cannot open "//trim(filep)
            write(*,*) "Pass the location of the orbital tables as insol_init(...,fldr=)."
            stop
        end if

        do n = 1, n_future
            read(unit,*) is%time(n_past-1+n), is%eccm(n_past-1+n), &
                         is%xobm(n_past-1+n), is%perm(n_past-1+n)
        end do

        close(unit)

        ! Convert from radians to degrees
        is%xobm = is%xobm * 180.0_dp/pi
        is%perm = is%perm * 180.0_dp/pi

        return

    end subroutine load_orbital_table

    subroutine ORBIT(ECC,XOBCH,TPERI,ZAVEXPE,  &
                     PCLOCK,PYTIME,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)
    !   *ORBIT* - FOR SOLAR ORBITAL PARAMETERS.
!
!       J.F.GELEYN     E.C.M.W.F.     03/06/82.
!
!       CHANGED TO 18000YBP-CONDITIONS BY LAUTENSCHLAGER
!                       MPI F. METEOR. HAMBURG  4.87
!       CHANGED TO OPTIONAL ORBITAL PARAMETERS BY GRAF, MPI 1/91
!       CHANGED AGAIN TO A MORE EXACT SOLUTION BY HOFFMANN, MPI 7/93
!       REVISED BY LORENZ, UNIVERSITY OF BREMEN, 10/93
!       -- last change: 94-08-03
!       -- Converted to FORTRAN90 by A. Robinson, 12/2013
!
!       PURPOSE.
!       --------
!
!            THIS ROUTINE COMPUTES THE SOLAR CONSTANT NORMALISED BY ITS
!       ANNUAL MEAN AND THREE ORBITAL PARAMETERS DEPENDING ON THE TIME
!       OF THE DAY AS WELL AS OF THE YEAR (BOTH IN RADIANS). TIME SCALE
!       ORIGIN IS 01/01 00 GMT OF THE DESIGNED TIME DISC. THE DIFFERENT
!       TIME DISCS ARE SCALED ON 03/21 12 GMT AS THE TIME OF THE VERNAL
!       EQUINOX. ALSO RETURNED IS A CONSTANT FOR THE EFFECT OF THE EARTH'S
!       CURVATURE ON THE COSINE OF THE SOLAR ZENITH ANGLE.
!
!       THE ORBITAL PARAMETERS CAN BE TAKEN FROM PROGRAM BERGOR, WRITTEN
!       BY LORENZ, UNIVERSITY OF BREMEN, 8/94
!
!       THERE ARE SEVEN DUMMY ARGUMENTS: 
!       
!       *PCLOCK* IS THE TIME OF THE DAY
!       *PYTIME* IS THE TIME OF THE YEAR (BOTH IN RADIANS).
!       *PDISSE* IS THE RATIO OF THE SOLAR CONSTANT TO ITS ANNUAL MEAN
!       *PZEN1*, *PZEN2* AND *PZEN3* ARE ZENITHAL PARAMETERS.
!       *PRAE* IS THE RATIO OF THE HEIGHT OF THE EQUIVALENT ATMOSPHERE 
!       TO THE RADIUS OF THE EARTH.

!       *ZCDIS*, *ZCDEC* AND *ZCEQT* ARE ARRAYS FOR A SECOND ORDER
!       *FOURIER DEVELOPMENT FOR RESPECTIVELY: SOLAR CONSTANT, SOLAR
!       DECLINATION AND EQUATION OF TIME.
!       THEY ARE VALID FOR TODAY'S ORBITAL PARAMETERS ONLY (LORENZ, 11/93).
!
!       DIMENSION ZCDIS(5),ZCDEC(5),ZCEQT(5)
!       DATA ZCDIS/+1.000110,+0.034221,+0.001280,+0.000719,+0.000077/
!       DATA ZCDEC/+0.006918,-0.399912,+0.070257,-0.006758,+0.000907/
!       DATA ZCEQT/+0.000075,+0.001868,-0.032077,-0.014615,-0.040849/
!
!       *ZRAE* IS THE VALUE FOR *PRAE*.
!
        real (dp), parameter :: ZRAE = 0.1277e-02_dp   

        real(dp) :: TTROP 

        real (dp) :: BTIME, ECC, XOBCH, TPERI, ZAVEXPE
        real (dp) :: PCLOCK, PYTIME, PDISSE, PZEN1, PZEN2, PZEN3, PRAE

        real (dp) :: ZCLOCK, ZYTIME

        real (dp) :: time, eold, enew, eps, epc, zeps, cose, E
        real (dp) :: ZDISSE, zsqecc, ztgean, znu, zlambda, xobche
        real (dp) :: zsinde, zdecli, ZZEN1, ZZEN2, ZZEN3 
        integer   :: niter 

        ! == PRELIMINARY SETTINGS ==

        ZCLOCK = PCLOCK
        ZYTIME = PYTIME
        TTROP  = 360.0_dp     !       TROPICAL YEAR

!       DAY OF VERNAL EQUINOX: DEFINED ON 21. OF MARCH, 12.00 GMT
!       (BEGIN OF YEAR: 01. OF JANUARY, 0.00 GMT: ZYTIME=0.0, THEREAFTER
!       E. G. 2. OF JANUARY 12.00 GMT: TIME IN DAYS FROM BEGIN OF YEAR = 1.5)
!       - NOT USED IN THE EXACT FORMULATION, *ZAVEXPE* USED INSTEAD
!       TVEREX=80.5              ! TROPICAL YEAR: 360 DAYS
!       TVEREX=79.5              ! TROPICAL YEAR: 365 DAYS

!        ***  ORBITAL PARAMETERS : VALUES FOR DIFFERENT TIME DISCS
!             THESE VALUES CAN BE TAKEN FROM PROGRAM BERGOR (S. LORENZ)
!
!       2. COMPUTATIONS
!          ------------

!       ZC1YT=COS(ZYTIME)
!       ZS1YT=SIN(ZYTIME)
!       ZC2YT=ZC1YT**2-ZS1YT**2
!       ZS2YT=2.*ZS1YT*ZC1YT
!
!       CALCULATE ECCENTRI!   ANOMALY:
!
!       LESS EXACT METHOD OF H. GRAF: SIN(E)=E
!       E=(ZYTIME - 2*PI*TPERI/TTROP) / (1.-ECC)
!
!       EXACT CALCULATION WITH KEPLER'S LAW (ELLIPSE)
!       (MONIN, 1986 : AN INTRODUCTION TO THE THEORY OF CLIMATE)
!       USE FIRST DERIVATIVE (NEWTON) TO CALCULATE SOLUTION OF
!       EQUATION OF ECCENTRI!   ANOMALY *E*
!
        time  = ZYTIME-2.0_dp*pi*TPERI/TTROP 
        eold  = time/(1.0_dp-ECC)
        enew  = time
        eps   = 1.e-6_dp 
        niter = 0
        zeps  = eold - enew

        do while ( dabs(zeps) .gt. eps ) 
            zeps = eold - enew
            if (niter .ge. 30) then 
                write(*,*) ' SUBROUTINE *ORBIT* - eccentri!   anomaly not found!'
                write(*,*) ' ERROR IN   *ORBIT* -- STOP'
                stop
            end if 

            niter = niter+1
            eold  = enew
            cose  = cos(enew)
            enew  = (time+ECC*(sin(enew)-enew*cose))/(1.0_dp-ECC*cose)
        end do 

        E      = enew
        ZDISSE = (1.0_dp/(1.0_dp-ECC*cos(E)))**2

!       Change in the calculation of the declination.
!       Used are not the formulas from Paltridge/Platt as in the ECHAM model
!       with fixed constants for contemporanious conditions
!       but the exact equations from Monin (s.a. - Keplers law)
!       are solved. Day of vernal equinox is fixed for a 360 day year on the
!       21. Maerz noon, with start of year at 01/01/00 GMT (s.a.).

        zsqecc = sqrt((1_dp+ECC)/(1_dp-ECC))
        ztgean = tan(E/2_dp)

!       znu: true anomaly (actual angle of Earth's position from perihelion)
!       zlambda: true longitude of the Earth (actual angle from vernal equinox)

        znu     = 2.0_dp*atan(zsqecc*ztgean)
        zlambda = znu + zavexpe
        xobche  = xobch/180.0_dp*pi
        zsinde  = sin(xobche)*sin(zlambda)
        zdecli  = asin(zsinde)

        ZZEN1 = sin(ZDECLI)
        ZZEN2 = cos(ZDECLI)*cos(ZCLOCK)
        ZZEN3 = cos(ZDECLI)*sin(ZCLOCK)

        ! 3. RETURN
        !    ------

        PDISSE = ZDISSE
        PZEN1  = ZZEN1
        PZEN2  = ZZEN2
        PZEN3  = ZZEN3
        PRAE   = ZRAE
 
        return

    end subroutine ORBIT 

end module

