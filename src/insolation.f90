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

    use interp1D 

    implicit none 

    integer,  parameter :: dp  = kind(1.0d0)

    real (dp), parameter :: pi = 2.d0*acos(0.d0)
    real (dp), parameter :: deg_to_rad = pi/180.0_dp
    
    type orbit_class
        integer :: n_orbit, n_lat
        real (dp), dimension(:), allocatable :: time, eccm, perm, xobm
        real (dp), dimension(:), allocatable :: lats, solarm, coszm 
    end type 

    type(orbit_class) :: OPAR 

    interface calc_insol_day
        module procedure calc_insol_day_pt  
        module procedure calc_insol_day_1D
        module procedure calc_insol_day_2D 
    end interface 

    interface calc_insol_days
        module procedure calc_insol_days_pt
        module procedure calc_insol_days_1D
        module procedure calc_insol_days_2D
    end interface

    private
    public :: calc_insol_ave_pt
    public :: calc_insol_day 
    public :: calc_insol_days

contains 

! =====================================================================================================
! Interface functions for routines for annual mean insolation

    function calc_insol_ave_pt(day_range,lat,time_bp,S0,day_year,fldr,step) result(insol)
        ! Given a range of days of the year, latitude and time before present,
        ! Calculate the average insolation at each latitude for the day_range
        ! time_bp = years before present (present==1950)

        implicit none 

        integer,  intent(IN) :: day_range(2)
        real(dp), intent(IN) :: lat
        real(dp), intent(IN) :: time_bp
        real(dp), intent(IN), optional :: S0
        integer,  intent(IN), optional :: day_year
        character(len=*), optional :: fldr 
        integer,  intent(IN), optional :: step 
        real(dp) :: insol
        
        ! Local variables
        integer, parameter :: nh = 24 
        integer :: j, h, day, ndays 
        real(dp) :: PER, ECC, XOBCH, TPERI, ZAVEXPE
        real(dp) :: PRAE, PCLOCK, PYTIME 
        real(dp) :: PDISSE(nh), PZEN1(nh), PZEN2(nh), PZEN3(nh)
        real(dp) :: S0_value 
        integer  :: day_max 
        integer  :: day_step 

        
        if (day_range(2) .le. day_range(1)) then 
            write(*,*) "calc_insol_ave_pt:: Error: day_range must give a positive total of days to calculate."
            write(*,*) "day_range: ", day_range 
            stop 
        end if 

        ! Determine the step to use for the integration
        day_step = 1
        if (present(step)) day_step = step 

        ! Define the solar constant
        S0_value = 1365.0_dp
        if (present(S0)) S0_value = S0  

        ! Define year length in days 
        day_max = 360
        if (present(day_year)) day_max = day_year

        ! Calculate orbital parameters for current year 
        call calc_orbital_par(time_bp,PER,ECC,XOBCH,TPERI,ZAVEXPE,fldr)

        insol = 0.0 
        ndays = 0

        do day = day_range(1), day_range(2), day_step

            ! Get hourly zenith angle values for input into daily insol function
            PYTIME = dble(day)/dble(day_max)*2.0_dp*pi
            do h = 1, nh
                PCLOCK = dble(h)/dble(nh)*2.0_dp*pi
                call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
                        PCLOCK,PYTIME,PDISSE(h),PZEN1(h),PZEN2(h),PZEN3(h),PRAE)
            end do 

            ! Calculate daily insolation at latitude of interest
            insol = insol + calc_insol_day_internal(lat,PDISSE,PZEN1,PZEN2,S0_value)

            ndays = ndays + 1 
        end do 

        ! Get the average 
        insol = insol / real(ndays,dp)

        return 

    end function calc_insol_ave_pt


! =====================================================================================================
! Interface functions for routines for a single day
    
    function calc_insol_day_pt(day,lat,time_bp,S0,day_year,fldr) result(insol)
        ! Given day of year, latitudes and time before present,
        ! Calculate the daily insolation values at each latitude 
        ! time_bp = years before present (present==1950)

        implicit none 

        integer,  intent(IN) :: day
        real(dp), intent(IN) :: lat
        real(dp), intent(IN) :: time_bp
        real(dp), intent(IN), optional :: S0
        integer,  intent(IN), optional :: day_year
        character(len=*), optional :: fldr 
        real(dp) :: insol
        
        ! Local variables
        integer :: j, h 
        integer, parameter :: nh = 24 
        real(dp) :: PER, ECC, XOBCH, TPERI, ZAVEXPE
        real(dp) :: PRAE, PCLOCK, PYTIME 
        real(dp) :: PDISSE(nh), PZEN1(nh), PZEN2(nh), PZEN3(nh)
        real(dp) :: S0_value 
        integer  :: day_max 
        
        ! Define the solar constant
        S0_value = 1365.0_dp
        if (present(S0)) S0_value = S0  

        ! Define year length in days 
        day_max = 360
        if (present(day_year)) day_max = day_year

        ! Calculate orbital parameters for current year 
        call calc_orbital_par(time_bp,PER,ECC,XOBCH,TPERI,ZAVEXPE,fldr)

        ! Get hourly zenith angle values for input into daily insol function
        PYTIME = dble(day)/dble(day_max)*2.0_dp*pi
        do h = 1, nh
            PCLOCK = dble(h)/dble(nh)*2.0_dp*pi
            call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
                       PCLOCK,PYTIME,PDISSE(h),PZEN1(h),PZEN2(h),PZEN3(h),PRAE)
        end do 

        ! Calculate daily insolation at latitude of interest
        insol = calc_insol_day_internal(lat,PDISSE,PZEN1,PZEN2,S0_value)

        return 

    end function calc_insol_day_pt

    function calc_insol_day_1D(day,lats,time_bp,S0,day_year,fldr) result(insol)
        ! Given day of year, latitudes and time before present,
        ! Calculate the daily insolation values at each latitude 
        ! time_bp = years before present (present==1950)

        implicit none 

        integer   :: day
        real (dp) :: lats(:) 
        real (dp) :: time_bp
        real (dp) :: insol(size(lats))

        integer, parameter :: nh = 24 

        real (dp) :: PER, ECC, XOBCH, TPERI, ZAVEXPE
        real (dp) :: PRAE, PCLOCK, PYTIME 
        real (dp) :: PDISSE(nh), PZEN1(nh), PZEN2(nh), PZEN3(nh)

        integer :: j, h 

        real (dp), optional :: S0
        real (dp)           :: S0_value 

        integer, optional   :: day_year
        integer             :: day_max 

        character(len=*), optional :: fldr 

        ! Define the solar constant
        S0_value = 1365.0_dp
        if (present(S0)) S0_value = S0  

        ! Define year length in days 
        day_max = 360
        if (present(day_year)) day_max = day_year

        ! Calculate orbital parameters for current year 
        call calc_orbital_par(time_bp,PER,ECC,XOBCH,TPERI,ZAVEXPE,fldr)

        ! Get hourly zenith angle values for input into daily insol function
        PYTIME = dble(day)/dble(day_max)*2.0_dp*pi
        do h = 1, nh
            PCLOCK = dble(h)/dble(nh)*2.0_dp*pi
            call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
                       PCLOCK,PYTIME,PDISSE(h),PZEN1(h),PZEN2(h),PZEN3(h),PRAE)
        end do 

if (.TRUE.) then 
        ! ================================================
        ! Using spline interpolation 
        ! ================================================
        ! Calculate daily insolation at predefined latitude values,
        ! then interpolate via spline to get insolation at desired latitudes
        do j = 1, OPAR%n_lat
            OPAR%solarm(j) = calc_insol_day_internal(OPAR%lats(j),PDISSE,PZEN1,PZEN2,S0_value)
        end do 
        insol(:) = interp_spline(OPAR%lats,OPAR%solarm,lats)

else
        ! ================================================
        ! Direct calculation at each latitude (no spline)
        ! ================================================
        ! Calculate daily insolation at each latitude
        do j = 1, size(lats)
            insol(j) = calc_insol_day_internal(lats(j),PDISSE,PZEN1,PZEN2,S0_value)
        end do

end if 

        where (insol .lt. 0.0_dp) insol = 0.0_dp 
        
        return 

    end function calc_insol_day_1D

    function calc_insol_day_2D(day,lats,time_bp,S0,day_year,fldr) result(insol)
        ! Wrapper function to calculate 2D array of insolation values
        ! for a given day, given a 2D array of latitude values
        ! time_bp = years before present (present==1950)

        implicit none 

        integer,  intent(IN) :: day
        real(dp), intent(IN) :: lats(:,:) 
        real(dp), intent(IN) :: time_bp
        real(dp), intent(IN), optional :: S0
        integer,  intent(IN), optional :: day_year
        character(len=*), intent(IN), optional :: fldr 
        real (dp) :: insol(size(lats,1),size(lats,2))
        
        ! Local variables 
        integer :: n 
        real(dp), allocatable :: lats1D(:)
        real(dp), allocatable :: insol1D(:)
        
        n = size(lats,1)*size(lats,2)
        allocate(lats1D(n),insol1D(n))

        lats1D = reshape(lats,[n])
        insol1D = calc_insol_day_1D(day,lats1D,time_bp,S0,day_year,fldr)
        insol = reshape(insol1D,[size(lats,1),size(lats,2)])

        return 

    end function calc_insol_day_2D

! =====================================================================================================
! Interface functions for routines with a vector of days

    function calc_insol_days_pt(days,lat,time_bp,S0,day_year,fldr) result(insol)
        ! Given day of year, latitudes and time before present,
        ! Calculate the daily insolation values at each latitude 
        ! time_bp = years before present (present==1950)

        implicit none 

        integer,  intent(IN) :: days(:)
        real(dp), intent(IN)  :: lat
        real(dp), intent(IN)  :: time_bp
        real(dp), intent(IN), optional :: S0
        integer,  intent(IN), optional :: day_year
        character(len=*), intent(IN), optional :: fldr 
        real(dp) :: insol(size(days))
        
        ! Local variables
        integer :: q 

        do q = 1, size(days,1)
            insol(q) = calc_insol_day_pt(days(q),lat,time_bp,S0,day_year,fldr)
        end do 

        return 

    end function calc_insol_days_pt
    
    function calc_insol_days_1D(days,lats,time_bp,S0,day_year,fldr) result(insol)
        ! Given day of year, latitudes and time before present,
        ! Calculate the daily insolation values at each latitude 
        ! time_bp = years before present (present==1950)

        implicit none 

        integer,  intent(IN) :: days(:)
        real(dp), intent(IN) :: lats(:) 
        real(dp), intent(IN) :: time_bp
        real(dp), intent(IN), optional :: S0
        integer, intent(IN),  optional :: day_year
        character(len=*), intent(IN), optional :: fldr 
        real(dp) :: insol(size(lats,1),size(days,1))
        
        
        ! Local variables
        integer :: q 

        do q = 1, size(days,1)
            insol(:,q) = calc_insol_day_1D(days(q),lats,time_bp,S0,day_year,fldr)
        end do 

        return 

    end function calc_insol_days_1D

    function calc_insol_days_2D(days,lats,time_bp,S0,day_year,fldr) result(insol)
        ! Wrapper function to calculate 2D array of insolation values
        ! for a given day, given a 2D array of latitude values
        ! time_bp = years before present (present==1950)

        implicit none 

        integer,  intent(IN) :: days(:)
        real(dp), intent(IN) :: lats(:,:) 
        real(dp), intent(IN) :: time_bp
        real(dp), intent(IN), optional :: S0
        integer,  intent(IN), optional :: day_year
        character(len=*), intent(IN), optional :: fldr 
        real (dp) :: insol(size(lats,1),size(lats,2),size(days,1))
        
        ! Local variables 
        integer :: n, q 
        real(dp), allocatable :: lats1D(:)
        real(dp), allocatable :: insol1D(:)
        
        n = size(lats,1)*size(lats,2)
        allocate(lats1D(n),insol1D(n))

        lats1D = reshape(lats,[n])

        do q = 1, size(days,1)
            insol1D = calc_insol_day_1D(days(q),lats1D,time_bp,S0,day_year,fldr)
            insol(:,:,q) = reshape(insol1D,[size(lats,1),size(lats,2)])
        end do 

        return 

    end function calc_insol_days_2D

! =====================================================================================================

    function calc_insol_day_analytical_internal(day,time_bp,lat,S0) result(solarm)
        ! Calculate the daily mean insolation value,
        ! given day of year, time before present, latitude and solar constant,

        implicit none 

        real(dp), intent(IN) :: time_bp
        real(dp), intent(IN) :: lat 
        real(dp), intent(IN) :: S0
        real (dp) :: solarm

        ! Local variables
        integer   :: day, h
        
        integer, parameter :: nh = 24
        integer, parameter :: day_max = 360 

        
        return 

    end function calc_insol_day_analytical_internal

    function calc_insol_day_internal(lat,PDISSE,PZEN1,PZEN2,S0) result(solarm)
        ! Given latitude and sun hourly angle,
        ! Calculate the daily insolation value 
        ! Adapted from CLIMBER2.4 subroutine 'sinsol' by A. Ganopolski

        implicit none 

        integer   :: day, h, nh   
        real (dp), intent(IN) :: lat, S0  
        real (dp), intent(IN) :: PDISSE(:), PZEN1(:), PZEN2(:)
        real (dp) :: latrad, coszm, cosp, cosn, S
        real (dp) :: solarm

        ! How many hours in the day?
        nh = size(PDISSE)

        ! Daily insolation is calculated by hourly integration for each day
        solarm = 0.0_dp
        coszm  = 0.0_dp 
        latrad = lat*deg_to_rad 

        do h = 1, nh 

            cosp = PZEN1(h)*dsin(latrad) + PZEN2(h)*dcos(latrad)
            cosn = max(cosp,0.0_dp)

            S = S0*cosn*PDISSE(h)
            solarm = solarm + S
            coszm  = coszm + S*cosn

        end do

        ! Get daily average insolation and zenith angle 
        solarm = solarm/24.0_dp
        if (solarm .gt. 0.0_dp) then
            coszm = coszm/(solarm*24.0_dp)
        else
            coszm = 0.0_dp 
        endif

        return 

    end function calc_insol_day_internal

    subroutine calc_orbital_par(time,per,ecc,obl,tper,zan,fldr)
        !
        ! Calculates Earth's orbital parameters for the current time before 1950 (time)
        !
        ! Adapted from Fortran77 subroutine BERGOR written by  
        ! S. J. Lorenz, University of Bremen, 08-03-94.
        !
        ! Orbital parameters are loaded from a precalculated table of orbital values
        !
        implicit none 

        real (dp) :: time, per, ecc, obl, tper, zan 
        character(len=*), optional :: fldr 

        integer :: N1, N2 
        real (dp) :: time_ka, ef 
        real (dp) :: perm1, perm2 
        real (dp) :: ang, tgnu, sqecc, tian
        real (dp) :: days, e, epc

        ! First, initialize global orbital parameter data if necessary
        if ( .not. allocated(OPAR%time) ) call orbital_par_init(fldr)

        ! Get current time in ka 
        time_ka = time*1d-3 

        ! Determine indices bracketing current time
        do N1 = 1, OPAR%n_orbit
            if (OPAR%time(N1) .gt. time_ka) exit 
        end do 
        N1 = N1 - 1
        N2 = N1 + 1 

        ! Get weighting of each index
        ef = (time_ka-OPAR%time(N1)) / (OPAR%time(N2)-OPAR%time(N1))
        
        ! Get current eccentricity, obliquity and perihelion
        ecc = (1.0_dp-ef)*OPAR%eccm(N1) + ef*OPAR%eccm(N2)
        obl = (1.0_dp-ef)*OPAR%xobm(N1) + ef*OPAR%xobm(N2)

        perm1 = OPAR%perm(N1)
        perm2 = OPAR%perm(N2)
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

    subroutine orbital_par_init(fldr)
        ! Initialize the global module variable OPAR,
        ! which holds the orbital parameters from -5 Ma => 1 Ma
        ! Files should be located in subdirectory 'input' or 
        ! the path 'fldr' should be given as an argument

        implicit none 

        character (len=*), optional :: fldr 
        character (len=512) :: input_dir
        character (len=512) :: filen, filep
        integer :: n 

        ! Define the input directory where orbital parameter files are located
        input_dir = "input"
        if (present(fldr)) input_dir = trim(fldr)

        ! First make sure all vectors are deallocated
        if (allocated(OPAR%time)) deallocate(OPAR%time)
        if (allocated(OPAR%eccm)) deallocate(OPAR%eccm)
        if (allocated(OPAR%perm)) deallocate(OPAR%perm)
        if (allocated(OPAR%xobm)) deallocate(OPAR%xobm)
        
        ! Allocate vectors to length of input parameter files
        OPAR%n_orbit = 6001 
        allocate(OPAR%time(OPAR%n_orbit),OPAR%eccm(OPAR%n_orbit))
        allocate(OPAR%perm(OPAR%n_orbit),OPAR%xobm(OPAR%n_orbit))

        ! Specify file names of negative and positive time orbital parameters
        filen = trim(input_dir)//'/LA2004.INSOLN.5Ma.txt'
        filep = trim(input_dir)//'/LA2004.INSOLP.1Ma.txt'

        ! ==== PAST (NEGATIVE) TIMES ====
        write(*,*) "Loading orbital parameters from file: "//trim(filen)
        open (1,file=trim(filen))

        do n=1,5001
            read (1,*) OPAR%time(5001-n+1), OPAR%eccm(5001-n+1), OPAR%xobm(5001-n+1), OPAR%perm(5001-n+1)
        enddo

        close (1)

        ! ==== FUTURE (POSITIVE) TIMES ====
        write(*,*) "Loading orbital parameters from file: "//trim(filep)
        open (1,file=trim(filep))

        do n=1,1001
            read (1,*) OPAR%time(5000+n), OPAR%eccm(5000+n), OPAR%xobm(5000+n), OPAR%perm(5000+n)
        enddo

        close (1)

        ! Convert from radians to degrees      
        OPAR%xobm = OPAR%xobm * 180.0/pi
        OPAR%perm = OPAR%perm * 180.0/pi

        ! Also initialize the default latitude vector for insol calculations
        ! (actual lat values are interpolated from these generic calcs)
        OPAR%n_lat = 61 
        if (allocated(OPAR%lats))   deallocate(OPAR%lats)
        if (allocated(OPAR%solarm)) deallocate(OPAR%solarm)
        if (allocated(OPAR%coszm))  deallocate(OPAR%coszm)
        allocate(OPAR%lats(OPAR%n_lat),OPAR%solarm(OPAR%n_lat),OPAR%coszm(OPAR%n_lat))

        do n = 1, OPAR%n_lat
            OPAR%lats(n) = -90.0_dp + (n-1)*180.0_dp/dble(OPAR%n_lat-1)
        end do 

        
        return 

    end subroutine orbital_par_init 

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

