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
                   PCLOCK,PYTIME,dble(day_max),PDISSE,PZEN1,PZEN2,PZEN3,PRAE)

        ! PZEN1 = sin(declination); PDISSE = (a/r)^2
        decl  = asin(max(-1.0_dp,min(1.0_dp,PZEN1)))
        dist2 = PDISSE

        return

    end subroutine calc_solar_geometry

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

    subroutine ORBIT(ECC,XOBCH,TPERI,ZAVEXPE,  &
                     PCLOCK,PYTIME,TTROP,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)
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


        real (dp) :: BTIME, ECC, XOBCH, TPERI, ZAVEXPE
        real (dp) :: PCLOCK, PYTIME, PDISSE, PZEN1, PZEN2, PZEN3, PRAE

        ! Length of the tropical year, in the same day units as TPERI and
        ! PYTIME. Supplied by the caller rather than fixed at 360 so that the
        ! solution can be labelled on a 365- or 366-day calendar.
        real (dp) :: TTROP

        real (dp) :: ZCLOCK, ZYTIME

        real (dp) :: time, eold, enew, eps, epc, zeps, cose, E
        real (dp) :: ZDISSE, zsqecc, ztgean, znu, zlambda, xobche
        real (dp) :: zsinde, zdecli, ZZEN1, ZZEN2, ZZEN3 
        integer   :: niter 

        ! == PRELIMINARY SETTINGS ==

        ZCLOCK = PCLOCK
        ZYTIME = PYTIME

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
                error stop 1
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

end module insol_orbit
