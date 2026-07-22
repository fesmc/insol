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

module insol_orbit
    ! Earth's orbital solution: the shared geometric foundation of the library.
    !
    ! Holds the orbital element tables of Laskar et al. (2004) and derives from
    ! them, for any time between -5 Ma and +1 Ma, the quantities that both the
    ! daily-mean insolation (`insolation`) and the surface irradiance
    ! decomposition (`irradiance`) are built on: solar declination and the
    ! inverse-square distance factor.
    !
    ! Keeping this in one place means both consumers see the same orbit. The
    ! irradiance side previously carried its own declination, a Fourier series
    ! (Spencer 1971) valid only for the present day and differing from this
    ! solution by up to 0.33 degrees.

    use insol_defs, only : dp, pi, deg_to_rad => degrees_to_radians

    implicit none

    type orbit_class
        ! Orbital element table, loaded once by orbit_init and read-only after.
        integer :: n_orbit
        real(dp), allocatable :: time(:)     ! [ka] relative to 1950
        real(dp), allocatable :: eccm(:)     ! Eccentricity
        real(dp), allocatable :: xobm(:)     ! Obliquity [degrees]
        real(dp), allocatable :: perm(:)     ! Longitude of perihelion [degrees]

        integer :: day_year                  ! Days in the year (360, 365, 366)
    end type

    private
    public :: orbit_class
    public :: orbit_init
    public :: orbit_end
    public :: calc_orbital_par
    public :: calc_solar_geometry
    public :: calc_solar_geometry_at
    public :: calc_solar_position

contains

    subroutine orbit_init(orb,fldr,day_year)
        ! Load the orbital tables and fix the calendar.

        implicit none

        type(orbit_class), intent(INOUT) :: orb
        character(len=*),  intent(IN), optional :: fldr
        integer,           intent(IN), optional :: day_year

        orb%day_year = 360
        if (present(day_year)) orb%day_year = day_year

        call check_day_year(orb%day_year)

        call load_orbital_table(orb,fldr)

        return

    end subroutine orbit_init

    subroutine orbit_end(orb)
        ! Release the orbital tables.

        implicit none

        type(orbit_class), intent(INOUT) :: orb

        if (allocated(orb%time)) deallocate(orb%time)
        if (allocated(orb%eccm)) deallocate(orb%eccm)
        if (allocated(orb%xobm)) deallocate(orb%xobm)
        if (allocated(orb%perm)) deallocate(orb%perm)

        return

    end subroutine orbit_end

    subroutine calc_solar_geometry(decl,dist2,day,day_max,ECC,XOBCH,TPERI,ZAVEXPE)
        ! Solar declination and the inverse-square distance factor (a/r)^2 for a
        ! given whole day of the year.
        !
        ! Neither quantity depends on the time of day, so the daily-mean
        ! insolation evaluates this once per day and reuses it for every
        ! latitude.

        implicit none

        real(dp), intent(OUT) :: decl        ! Solar declination [radians]
        real(dp), intent(OUT) :: dist2       ! (a/r)^2, ratio of solar constant to its annual mean
        integer,  intent(IN)  :: day
        integer,  intent(IN)  :: day_max
        real(dp), intent(IN)  :: ECC, XOBCH, TPERI, ZAVEXPE

        ! Local variables
        real(dp) :: eqtime

        call calc_solar_geometry_at(decl,dist2,eqtime,real(day,dp),day_max, &
                                    ECC,XOBCH,TPERI,ZAVEXPE)

        return

    end subroutine calc_solar_geometry

    subroutine calc_solar_geometry_at(decl,dist2,eqtime,day_time,day_max, &
                                      ECC,XOBCH,TPERI,ZAVEXPE)
        ! Solar declination, distance factor and equation of time at an arbitrary
        ! instant, given as elapsed days since 00:00 on 1 January.
        !
        ! Follows Monin (1986): solve Kepler's equation for the eccentric anomaly,
        ! from which the true anomaly and hence the true solar longitude follow
        ! exactly. Declination is then the projection of that longitude onto the
        ! equator, tilted by the obliquity.
        !
        ! The equation of time -- apparent minus mean solar time -- comes from the
        ! same solve: it is the difference between the mean sun's longitude and
        ! the true sun's right ascension, so it captures both the eccentricity
        ! and the obliquity contribution without a separate empirical series.
        !
        ! Replaces the ECMWF subroutine ORBIT (Geleyn 1982, and successive
        ! revisions through Lorenz 1994), which computed the same declination and
        ! distance but exposed neither the anomalies nor the equation of time.

        implicit none

        real(dp), intent(OUT) :: decl        ! Solar declination [radians]
        real(dp), intent(OUT) :: dist2       ! (a/r)^2
        real(dp), intent(OUT) :: eqtime      ! Equation of time [radians of hour angle]
        real(dp), intent(IN)  :: day_time    ! Elapsed days since 1 Jan 00:00
        integer,  intent(IN)  :: day_max     ! Days in the year
        real(dp), intent(IN)  :: ECC         ! Eccentricity
        real(dp), intent(IN)  :: XOBCH       ! Obliquity [degrees]
        real(dp), intent(IN)  :: TPERI       ! Time of perihelion [days from year start]
        real(dp), intent(IN)  :: ZAVEXPE     ! Longitude of perihelion from vernal equinox [radians]

        ! Local variables
        integer, parameter :: niter_max = 30
        real(dp), parameter :: tol = 1.0e-12_dp

        integer  :: niter
        real(dp) :: ytime, mean_anom, ecc_anom, ecc_anom_new, de
        real(dp) :: true_anom, lambda, obliq, mean_lon, right_asc

        ytime = day_time/real(day_max,dp)*2.0_dp*pi

        ! Mean anomaly: angle of a fictitious sun moving at constant rate,
        ! measured from perihelion.
        mean_anom = ytime - 2.0_dp*pi*TPERI/real(day_max,dp)

        ! Solve Kepler's equation E - e sin(E) = M by Newton iteration.
        ecc_anom = mean_anom
        do niter = 1, niter_max
            ecc_anom_new = (mean_anom + ECC*(sin(ecc_anom)-ecc_anom*cos(ecc_anom))) &
                           / (1.0_dp - ECC*cos(ecc_anom))
            de = ecc_anom_new - ecc_anom
            ecc_anom = ecc_anom_new
            if (abs(de) .lt. tol) exit
        end do

        if (niter .gt. niter_max) then
            write(*,*) "calc_solar_geometry_at:: Error: eccentric anomaly did not converge."
            write(*,*) "eccentricity: ", ECC
            write(*,*) "mean anomaly: ", mean_anom
            error stop 1
        end if

        ! Distance factor: the solar constant scales as (a/r)^2.
        dist2 = (1.0_dp/(1.0_dp - ECC*cos(ecc_anom)))**2

        ! True anomaly, then true longitude measured from the vernal equinox.
        true_anom = 2.0_dp*atan(sqrt((1.0_dp+ECC)/(1.0_dp-ECC))*tan(ecc_anom/2.0_dp))
        lambda    = true_anom + ZAVEXPE

        obliq = XOBCH*deg_to_rad
        decl  = asin(max(-1.0_dp,min(1.0_dp,sin(obliq)*sin(lambda))))

        ! Equation of time. The mean sun advances uniformly along the equator;
        ! the true sun's position on the equator is its right ascension, which
        ! lags or leads because the orbit is eccentric and the ecliptic tilted.
        mean_lon  = mean_anom + ZAVEXPE
        right_asc = atan2(cos(obliq)*sin(lambda),cos(lambda))

        eqtime = mean_lon - right_asc

        ! Wrap into (-pi,pi]: the difference is at most a few degrees, so any
        ! larger value is an accumulated multiple of a full turn.
        eqtime = eqtime - 2.0_dp*pi*anint(eqtime/(2.0_dp*pi))

        return

    end subroutine calc_solar_geometry_at

    subroutine calc_solar_position(orb,cosz,time_bp,day_time,lat,lon,decl,eqtime)
        ! Cosine of the solar zenith angle at an instant and a location.
        !
        ! `day_time` is elapsed days since 00:00 UTC on 1 January, so 1.5 is
        ! noon on 2 January. `lon` is degrees east. Negative values of cosz mean
        ! the sun is below the horizon and are returned as given, so a caller can
        ! distinguish night from a sun exactly on the horizon.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(dp),          intent(OUT) :: cosz
        real(dp),          intent(IN)  :: time_bp     ! [yr] relative to 1950
        real(dp),          intent(IN)  :: day_time    ! Elapsed days since 1 Jan 00:00 UTC
        real(dp),          intent(IN)  :: lat         ! [degrees north]
        real(dp),          intent(IN)  :: lon         ! [degrees east]
        real(dp),          intent(OUT), optional :: decl
        real(dp),          intent(OUT), optional :: eqtime

        ! Local variables
        real(dp) :: PER, ECC, XOBCH, TPERI, ZAVEXPE
        real(dp) :: decl_now, dist2, eqtime_now
        real(dp) :: hour_utc, solar_time, hour_angle
        real(dp) :: latrad

        call calc_orbital_par(orb,time_bp,PER,ECC,XOBCH,TPERI,ZAVEXPE)
        call calc_solar_geometry_at(decl_now,dist2,eqtime_now,day_time,orb%day_year, &
                                    ECC,XOBCH,TPERI,ZAVEXPE)

        ! Hour angle: zero at local solar noon, advancing 2*pi per day. Mean
        ! solar time at the meridian, corrected to apparent solar time by the
        ! equation of time.
        hour_utc   = (day_time - floor(day_time))*2.0_dp*pi
        solar_time = hour_utc + lon*deg_to_rad + eqtime_now
        hour_angle = solar_time - pi

        latrad = lat*deg_to_rad

        cosz = sin(latrad)*sin(decl_now) + cos(latrad)*cos(decl_now)*cos(hour_angle)

        if (present(decl))   decl   = decl_now
        if (present(eqtime)) eqtime = eqtime_now

        return

    end subroutine calc_solar_position

    subroutine calc_orbital_par(orb,time,per,ecc,obl,tper,zan)
        !
        ! Calculates Earth's orbital parameters for the current time before 1950 (time)
        !
        ! Adapted from Fortran77 subroutine BERGOR written by  
        ! S. J. Lorenz, University of Bremen, 08-03-94.
        !
        ! Orbital parameters are loaded from a precalculated table of orbital values
        !
        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(dp),          intent(IN)  :: time
        real(dp),          intent(OUT) :: per, ecc, obl, tper, zan

        integer :: N1, N2, NMID
        real (dp) :: time_ka, ef
        real (dp) :: perm1, perm2
        real (dp) :: ang, tgnu, sqecc, tian
        real (dp) :: days, e, epc
        real (dp) :: day_year

        ! Get current time in ka
        time_ka = time*1d-3

        ! The tables span -5 Ma to +1 Ma; extrapolating beyond them is
        ! meaningless, and reading past the end used to run off the array.
        if (time_ka .lt. orb%time(1) .or. time_ka .gt. orb%time(orb%n_orbit)) then
            write(*,*) "calc_orbital_par:: Error: time_bp outside the range of the"
            write(*,*) "orbital parameter tables."
            write(*,*) "time_bp [yr]:   ", time
            write(*,*) "table range [yr]:", orb%time(1)*1d3, orb%time(orb%n_orbit)*1d3
            error stop 1
        end if

        ! Determine the indices bracketing the current time by bisection. The
        ! table holds several thousand entries, so a linear scan from the start
        ! costs more than the interpolation it feeds.
        N1 = 1
        N2 = orb%n_orbit
        do while (N2-N1 .gt. 1)
            NMID = (N1+N2)/2
            if (orb%time(NMID) .gt. time_ka) then
                N2 = NMID
            else
                N1 = NMID
            end if
        end do

        ! Get weighting of each index
        ef = (time_ka-orb%time(N1)) / (orb%time(N2)-orb%time(N1))

        ! Get current eccentricity, obliquity and perihelion
        ecc = (1.0_dp-ef)*orb%eccm(N1) + ef*orb%eccm(N2)
        obl = (1.0_dp-ef)*orb%xobm(N1) + ef*orb%xobm(N2)

        perm1 = orb%perm(N1)
        perm2 = orb%perm(N2)
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
        day_year = real(orb%day_year,dp)

        days = tian/(2.0_dp*pi)*day_year           !   days from ver.eq. to perh.
        if (days .lt. 0.0_dp) days = days+day_year

        !   time in days from begin of year: vernal equinox fixed at 3/21, 12 GMT

        tper = days+calc_day_vernal_equinox(orb%day_year)
        if (tper .gt. day_year) tper = tper-day_year

        epc = ecc*100.0_dp

        return

    end subroutine calc_orbital_par 

    subroutine load_orbital_table(orb,fldr)
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

        type(orbit_class), intent(INOUT) :: orb
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
        if (allocated(orb%time)) deallocate(orb%time)
        if (allocated(orb%eccm)) deallocate(orb%eccm)
        if (allocated(orb%perm)) deallocate(orb%perm)
        if (allocated(orb%xobm)) deallocate(orb%xobm)

        ! Allocate vectors to length of input parameter files
        orb%n_orbit = n_past + n_future - 1
        allocate(orb%time(orb%n_orbit),orb%eccm(orb%n_orbit))
        allocate(orb%perm(orb%n_orbit),orb%xobm(orb%n_orbit))

        ! Specify file names of negative and positive time orbital parameters
        filen = trim(input_dir)//'/LA2004.INSOLN.5Ma.txt'
        filep = trim(input_dir)//'/LA2004.INSOLP.1Ma.txt'

        ! ==== PAST (NEGATIVE) TIMES ====
        open(newunit=unit,file=trim(filen),status="old",action="read",iostat=io)
        if (io /= 0) then
            write(*,*) "load_orbital_table:: Error: cannot open "//trim(filen)
            write(*,*) "Pass the location of the orbital tables as insol_init(...,fldr=)."
            error stop 1
        end if

        do n = 1, n_past
            read(unit,*) orb%time(n_past-n+1), orb%eccm(n_past-n+1), &
                         orb%xobm(n_past-n+1), orb%perm(n_past-n+1)
        end do

        close(unit)

        ! ==== FUTURE (POSITIVE) TIMES ====
        open(newunit=unit,file=trim(filep),status="old",action="read",iostat=io)
        if (io /= 0) then
            write(*,*) "load_orbital_table:: Error: cannot open "//trim(filep)
            write(*,*) "Pass the location of the orbital tables as insol_init(...,fldr=)."
            error stop 1
        end if

        do n = 1, n_future
            read(unit,*) orb%time(n_past-1+n), orb%eccm(n_past-1+n), &
                         orb%xobm(n_past-1+n), orb%perm(n_past-1+n)
        end do

        close(unit)

        ! Convert from radians to degrees
        orb%xobm = orb%xobm * 180.0_dp/pi
        orb%perm = orb%perm * 180.0_dp/pi

        return

    end subroutine load_orbital_table

    function calc_day_vernal_equinox(day_year) result(tverex)
        ! Day of the vernal equinox, counted from the start of the year, for
        ! each supported calendar.
        !
        ! The orbital solution pins the vernal equinox to March 21, 12:00 GMT.
        ! Which day number that is depends only on the lengths of January and
        ! February in the calendar being used, so `day_year` selects a labelling
        ! of the orbit rather than changing the orbit itself.

        implicit none

        integer, intent(IN) :: day_year
        real(dp) :: tverex

        select case (day_year)

            case (360)
                tverex = 80.5_dp        ! 30 + 30 + 20.5

            case (365)
                tverex = 79.5_dp        ! 31 + 28 + 20.5

            case (366)
                tverex = 80.5_dp        ! 31 + 29 + 20.5

            case default
                ! Unreachable: insol_config rejects other values at init.
                tverex = 0.0_dp

        end select

        return

    end function calc_day_vernal_equinox

    subroutine check_day_year(day_year)
        ! Only calendars whose vernal equinox day is defined are supported.
        ! Silently accepting others would shift the seasons relative to the
        ! orbit, which is how this argument used to behave.

        implicit none

        integer, intent(IN) :: day_year

        select case (day_year)

            case (360, 365, 366)
                ! Supported

            case default
                write(*,*) "insol_init:: Error: unsupported day_year."
                write(*,*) "day_year: ", day_year
                write(*,*) "Supported values are 360, 365 and 366."
                error stop 1

        end select

        return

    end subroutine check_day_year

 

end module insol_orbit
