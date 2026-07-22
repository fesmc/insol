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

    ! Insolation while the sun is above a minimum elevation angle. Resolved on
    ! the rank of `insol` in the same way as calc_insol_day.
    interface calc_insol_above
        module procedure calc_insol_above_1D
        module procedure calc_insol_above_2D
        module procedure calc_insol_above_pt
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
    public :: calc_insol_above

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
! Insolation while the sun is above a minimum elevation

    subroutine calc_insol_above_pt(is,insol,tday,coszm,day,lat,time_bp,elev_min)
        ! Insolation at one latitude, restricted to the part of the day when the
        ! sun stands higher than `elev_min` above the horizon.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol     ! Mean flux while above elev_min [W/m2]
        real(wp),          intent(OUT) :: tday      ! Fraction of the day above elev_min
        real(wp),          intent(OUT) :: coszm     ! Flux-weighted mean cos(zenith)
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: lat
        real(dp),          intent(IN)  :: time_bp
        real(dp),          intent(IN)  :: elev_min  ! [degrees above the horizon]

        ! Local variables
        real(dp) :: decl, dist2, latrad
        real(wp) :: insol1(1), tday1(1), coszm1(1)

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        latrad = lat*deg_to_rad

        call calc_insol_above_kernel(insol1,tday1,coszm1,[sin(latrad)],[cos(latrad)], &
                                     [elev_min],decl,dist2,is%S0)

        insol = insol1(1)
        tday  = tday1(1)
        coszm = coszm1(1)

        return

    end subroutine calc_insol_above_pt

    subroutine calc_insol_above_1D(is,insol,tday,coszm,day,time_bp,elev_min)
        ! As calc_insol_above_pt, over the 1D latitude grid held in `is`.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:)
        real(wp),          intent(OUT) :: tday(:)
        real(wp),          intent(OUT) :: coszm(:)
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp
        real(dp),          intent(IN)  :: elev_min(:)   ! One threshold per grid point

        ! Local variables
        real(dp) :: decl, dist2

        call check_grid(is,size(insol),1,"calc_insol_above")
        call check_elev_min(size(elev_min),is%n_lat,"calc_insol_above")

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        call calc_insol_above_kernel(insol,tday,coszm,is%sinlat,is%coslat, &
                                     elev_min,decl,dist2,is%S0)

        return

    end subroutine calc_insol_above_1D

    subroutine calc_insol_above_2D(is,insol,tday,coszm,day,time_bp,elev_min)
        ! As calc_insol_above_pt, over the 2D latitude grid held in `is`.

        implicit none

        type(insol_class), intent(IN)  :: is
        real(wp),          intent(OUT) :: insol(:,:)
        real(wp),          intent(OUT) :: tday(:,:)
        real(wp),          intent(OUT) :: coszm(:,:)
        integer,           intent(IN)  :: day
        real(dp),          intent(IN)  :: time_bp
        real(dp),          intent(IN)  :: elev_min(:,:)

        ! Local variables
        real(dp) :: decl, dist2
        real(wp) :: insol1D(is%n_lat), tday1D(is%n_lat), coszm1D(is%n_lat)

        call check_grid(is,size(insol),2,"calc_insol_above")
        call check_elev_min(size(elev_min),is%n_lat,"calc_insol_above")

        call calc_solar_geometry_day(is,decl,dist2,day,time_bp)

        call calc_insol_above_kernel(insol1D,tday1D,coszm1D,is%sinlat,is%coslat, &
                                     reshape(elev_min,[is%n_lat]),decl,dist2,is%S0)

        insol = reshape(insol1D,[is%n1,is%n2])
        tday  = reshape(tday1D, [is%n1,is%n2])
        coszm = reshape(coszm1D,[is%n1,is%n2])

        return

    end subroutine calc_insol_above_2D

    subroutine calc_insol_above_kernel(insol,tday,coszm,sinlat,coslat,elev_min, &
                                       decl,dist2,S0)
        ! Insolation over the part of the day when the solar elevation exceeds
        ! elev_min, in closed form.
        !
        ! Written for melt schemes that resolve the diurnal cycle: dEBM
        ! (Krebs-Kanzow et al., 2018) melts only while the sun is above a
        ! critical elevation Phi, estimated locally, which is why the threshold
        ! is taken per grid point rather than as a scalar.
        !
        ! With cos(z) = a + b cos(h), a = sin(lat)sin(dec), b = cos(lat)cos(dec),
        ! the sun is above elevation Phi for |h| < h_phi where
        !
        !     cos(h_phi) = (sin(Phi) - a)/b
        !
        ! and the integrals over that window are
        !
        !     I1 = int cos(z)   dh = 2( a h_phi + b sin(h_phi) )
        !     I2 = int cos(z)^2 dh = 2( (a^2 + b^2/2) h_phi + 2ab sin(h_phi)
        !                               + (b^2/4) sin(2 h_phi) )
        !
        ! giving the mean flux over the window as S0 (a/r)^2 I1/(2 h_phi), the
        ! fraction of the day as h_phi/pi, and the flux-weighted mean cosine of
        ! the zenith angle -- the quantity that sets the effective air mass for
        ! the energy actually received -- as I2/I1.
        !
        ! elev_min = 0 recovers sunrise-to-sunset, for which insol*tday equals
        ! the daily mean returned by calc_insol_day.

        implicit none

        real(wp), intent(OUT) :: insol(:)
        real(wp), intent(OUT) :: tday(:)
        real(wp), intent(OUT) :: coszm(:)
        real(dp), intent(IN), contiguous :: sinlat(:)
        real(dp), intent(IN), contiguous :: coslat(:)
        real(dp), intent(IN) :: elev_min(:)
        real(dp), intent(IN) :: decl, dist2, S0

        ! Local variables
        integer  :: j
        real(dp) :: sindec, cosdec, flux
        real(dp) :: a, b, sin_phi, cos_h, h, sin_h, i1, i2

        sindec = sin(decl)
        cosdec = cos(decl)
        flux   = S0*dist2

        do j = 1, size(sinlat)

            a = sinlat(j)*sindec
            b = coslat(j)*cosdec

            sin_phi = sin(elev_min(j)*deg_to_rad)

            ! At the poles b vanishes and the elevation does not vary through the
            ! day: the sun is either above the threshold for the whole day or not
            ! at all.
            if (abs(b) .lt. 1.0e-14_dp) then
                if (a .gt. sin_phi) then
                    cos_h = -1.0_dp
                else
                    cos_h =  1.0_dp
                end if
            else
                cos_h = (sin_phi-a)/b
            end if

            if (cos_h .ge. 1.0_dp) then
                h     = 0.0_dp              ! Never rises above elev_min
                sin_h = 0.0_dp
            else if (cos_h .le. -1.0_dp) then
                h     = pi                  ! Above elev_min all day
                sin_h = 0.0_dp
            else
                h     = acos(cos_h)
                sin_h = sqrt(1.0_dp-cos_h*cos_h)
            end if

            tday(j) = real(h/pi,wp)

            if (h .le. 0.0_dp) then
                insol(j) = 0.0_wp
                coszm(j) = 0.0_wp
                cycle
            end if

            i1 = 2.0_dp*(a*h + b*sin_h)
            i2 = 2.0_dp*((a*a + 0.5_dp*b*b)*h + 2.0_dp*a*b*sin_h &
                         + 0.25_dp*b*b*sin(2.0_dp*h))

            if (i1 .le. 0.0_dp) then
                insol(j) = 0.0_wp
                coszm(j) = 0.0_wp
            else
                insol(j) = real(flux*i1/(2.0_dp*h),wp)
                coszm(j) = real(i2/i1,wp)
            end if

        end do

        return

    end subroutine calc_insol_above_kernel

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


end module

