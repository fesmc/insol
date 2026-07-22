! Tests for the solar geometry layer (insol_orbit).
!
! The daily-mean insolation tests in test_insol.f90 exercise the idealized
! calendar used for paleoclimate. These cover the real-date path instead, which
! is what a present-day climate model or the irradiance decomposition uses, and
! they are scored against NREL's own solar position calculation.
!
! Exits with status 1 if any check fails.

program test_orbit

    use insol_defs, only : dp, pi, radians_to_degrees
    use insol_orbit

    implicit none

    character(len=*), parameter :: fldr     = "input"
    character(len=*), parameter :: ref_file = "tests/reference/nrel_stations.txt"

    ! Bounds on agreement with the NSRDB solar zenith angle. NREL's algorithm is
    ! accurate to ~0.0003 degrees, so in principle this could be far tighter.
    ! In practice the comparison floors out near 0.19 deg RMS and no correction
    ! moves it: the best-fitting time offset is exactly zero, a refraction
    ! correction changes it by 0.005 deg, and the residual is scatter rather than
    ! bias. The likely cause is that a half-hourly NSRDB record represents its
    ! interval rather than an instant, which no instantaneous calculation can
    ! reproduce. The bounds are set to catch regressions, not to certify accuracy.
    real(dp), parameter :: tol_rms = 0.25_dp
    real(dp), parameter :: tol_max = 0.50_dp

    integer :: n_pass = 0
    integer :: n_fail = 0

    write(*,*)
    write(*,*) "=============================================="
    write(*,*) " insol_orbit tests"
    write(*,*) "=============================================="

    call test_julian_day()
    call test_solar_position_nrel()
    call test_equation_of_time()

    write(*,*)
    write(*,"(a,i0,a,i0,a)") "  ", n_pass, " passed, ", n_fail, " failed"
    write(*,*)

    if (n_fail > 0) then
        write(*,*) "FAILED"
        stop 1
    end if

    write(*,*) "PASSED"

contains

    subroutine test_julian_day()
        ! Worked examples from Meeus (1998), chapter 7.

        implicit none

        call check_value("Julian Day 2000-01-01", calc_julian_day(2000,1,1), 2451544.5_dp)
        call check_value("Julian Day 1987-01-27", calc_julian_day(1987,1,27), 2446822.5_dp)
        call check_value("Julian Day 1957-10-04", calc_julian_day(1957,10,4), 2436115.5_dp)
        call check_value("Julian Day 1900-01-01", calc_julian_day(1900,1,1),  2415020.5_dp)

        ! The Gregorian leap rule: 1900 was not a leap year, 2000 was.
        call check_value("1900 is not a leap year", &
                         calc_julian_day(1900,3,1)-calc_julian_day(1900,2,28), 1.0_dp)
        call check_value("2000 is a leap year", &
                         calc_julian_day(2000,3,1)-calc_julian_day(2000,2,28), 2.0_dp)

    end subroutine test_julian_day

    subroutine test_solar_position_nrel()
        ! Compare the real-date solar position against the NSRDB station records.

        implicit none

        type(orbit_class) :: orb
        character(len=512) :: line
        character(len=32)  :: station
        real(dp) :: lat, lon, tz, doy, hr, mn, hour_utc
        real(dp) :: zen_ref, ghi, dni, dhi, temp, pres
        real(dp) :: cosz, zen, err, sum2, err_max
        integer  :: year, io, unit, n

        call orbit_init(orb,fldr=fldr,day_year=365)

        n = 0
        sum2 = 0.0_dp
        err_max = 0.0_dp

        open(newunit=unit,file=ref_file,status="old",iostat=io)
        if (io /= 0) then
            call fail("solar position vs NREL", "cannot open "//ref_file)
            return
        end if

        do
            read(unit,'(a)',iostat=io) line
            if (io /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            read(line,*) station, lat, lon, tz, year, doy, hr, mn, hour_utc, &
                         zen_ref, ghi, dni, dhi, temp, pres

            ! hour_utc counts from 1 January 00:00 UTC, so it can be handed to
            ! the date routine directly as an hour offset from that date.
            call calc_solar_position_date(orb,cosz,year,1,1,hour_utc,lat,lon)

            zen = acos(max(-1.0_dp,min(1.0_dp,cosz)))*radians_to_degrees

            ! Score only where the sun is well up. Near and below the horizon the
            ! comparison is dominated by refraction and by how the record
            ! represents its interval, neither of which this routine models.
            if (zen_ref .lt. 85.0_dp) then
                err  = abs(zen-zen_ref)
                sum2 = sum2 + err*err
                if (err > err_max) err_max = err
                n = n + 1
            end if
        end do

        close(unit)

        call orbit_end(orb)

        if (n == 0) then
            call fail("solar position vs NREL", "no daytime records read")
            return
        end if

        write(*,"(a,i0,a,f8.4,a,f8.4)") "     NREL daytime records = ", n, &
                                        ",  RMS = ", sqrt(sum2/real(n,dp)), &
                                        ",  max = ", err_max

        call check("solar position vs NREL (RMS)", sqrt(sum2/real(n,dp)), tol_rms, &
                   "RMS zenith error [deg]")
        call check("solar position vs NREL (max)", err_max, tol_max, &
                   "max zenith error [deg]")

    end subroutine test_solar_position_nrel

    subroutine test_equation_of_time()
        ! The equation of time swings between about -14 and +16 minutes, with
        ! four zero crossings a year. Getting the sign or the obliquity term
        ! wrong changes the shape, so check the extremes rather than a mean.

        implicit none

        type(orbit_class) :: orb
        real(dp) :: cosz, eqt, eqt_min, eqt_max
        integer  :: day

        call orbit_init(orb,fldr=fldr,day_year=365)

        eqt_min =  1.0e30_dp
        eqt_max = -1.0e30_dp

        do day = 1, 365
            call calc_solar_position_date(orb,cosz,2019,1,1,real(day-1,dp)*24.0_dp, &
                                          0.0_dp,0.0_dp,eqtime=eqt)
            ! Radians of hour angle -> minutes of time.
            eqt = eqt*radians_to_degrees*4.0_dp
            eqt_min = min(eqt_min,eqt)
            eqt_max = max(eqt_max,eqt)
        end do

        call orbit_end(orb)

        write(*,"(a,f8.3,a,f8.3,a)") "     equation of time range = ", eqt_min, &
                                     " to ", eqt_max, " min"

        call check_range("equation of time minimum", eqt_min, -15.0_dp, -13.0_dp)
        call check_range("equation of time maximum", eqt_max,  15.5_dp,  17.5_dp)

    end subroutine test_equation_of_time

    ! ------------------------------------------------------------------
    ! Check framework
    ! ------------------------------------------------------------------

    subroutine check(name, value, tol, what)

        implicit none

        character(len=*), intent(IN) :: name, what
        real(dp),         intent(IN) :: value, tol

        if (.not. (value <= tol)) then
            call fail(name, trim(what)//" = "//trim(rtoa(value))// &
                            " exceeds tolerance "//trim(rtoa(tol)))
        else
            call pass(name)
        end if

    end subroutine check

    subroutine check_value(name, value, expected)

        implicit none

        character(len=*), intent(IN) :: name
        real(dp),         intent(IN) :: value, expected

        if (.not. (abs(value-expected) <= 1.0e-9_dp)) then
            call fail(name, "got "//trim(rtoa(value))// &
                            ", expected "//trim(rtoa(expected)))
        else
            call pass(name)
        end if

    end subroutine check_value

    subroutine check_range(name, value, lo, hi)

        implicit none

        character(len=*), intent(IN) :: name
        real(dp),         intent(IN) :: value, lo, hi

        if (.not. (value >= lo .and. value <= hi)) then
            call fail(name, trim(rtoa(value))//" outside ["// &
                            trim(rtoa(lo))//","//trim(rtoa(hi))//"]")
        else
            call pass(name)
        end if

    end subroutine check_range

    subroutine pass(name)

        implicit none

        character(len=*), intent(IN) :: name

        n_pass = n_pass + 1
        write(*,"(a,a)") "  [ ok ] ", name

    end subroutine pass

    subroutine fail(name, reason)

        implicit none

        character(len=*), intent(IN) :: name, reason

        n_fail = n_fail + 1
        write(*,"(a,a)") "  [FAIL] ", name
        write(*,"(a,a)") "         ", reason

    end subroutine fail

    function rtoa(x) result(str)

        implicit none

        real(dp), intent(IN) :: x
        character(len=32)    :: str

        if (abs(x) >= 1.0e5_dp) then
            write(str,'(f14.1)') x
        else
            write(str,'(es12.4)') x
        end if
        str = adjustl(str)

    end function rtoa

end program test_orbit
