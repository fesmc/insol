! Tests for the irradiance module: estimating direct normal irradiance (DNI)
! from global horizontal irradiance (GHI).
!
! Unlike the insolation core, the decomposition models (DISC and the rest) are
! empirical estimators, not closed forms. There is no exact answer to check
! against -- only measurements. So these tests score the estimate against the
! NREL NSRDB station DNI and assert a bound on the error, and they exist mainly
! to pin the behaviour so a refactor cannot silently change it.
!
! Exits with status 1 if any check fails.

program test_irradiance

    use irradiance

    implicit none

    character(len=*), parameter :: ref_file = "tests/reference/nrel_stations.txt"

    ! Bound on the RMSE of estimated vs measured DNI, over daytime records.
    ! DISC is an empirical estimator; against station data its error is a large
    ! fraction of DNI, so this bounds the estimate rather than certifying it.
    ! Set from the measured baseline with headroom, to catch regressions.
    real(wp), parameter :: tol_rmse = 185.0_wp

    integer :: n_pass = 0
    integer :: n_fail = 0

    write(*,*)
    write(*,*) "=============================================="
    write(*,*) " irradiance tests"
    write(*,*) "=============================================="

    call test_dni_vs_nrel()

    write(*,*)
    write(*,"(a,i0,a,i0,a)") "  ", n_pass, " passed, ", n_fail, " failed"
    write(*,*)

    if (n_fail > 0) then
        write(*,*) "FAILED"
        stop 1
    end if

    write(*,*) "PASSED"

contains

    subroutine test_dni_vs_nrel()
        ! Estimate DNI from the station GHI with DISC and compare to measured DNI.

        implicit none

        character(len=512) :: line
        character(len=32)  :: station
        real(wp) :: lat, lon, tz, hour_utc, zen_ref
        real(wp) :: ghi, dni_ref, dhi, temp, pres_mbar
        real(wp) :: dni, pres_pa, err, sum2, sum_ref2, bias
        integer  :: year, doy, hr, mn
        integer  :: io, unit, n

        n = 0
        sum2 = 0.0_wp
        sum_ref2 = 0.0_wp
        bias = 0.0_wp

        open(newunit=unit,file=ref_file,status="old",iostat=io)
        if (io /= 0) then
            call fail("DNI vs NREL", "cannot open "//ref_file)
            return
        end if

        do
            read(unit,'(a)',iostat=io) line
            if (io /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            read(line,*) station, lat, lon, tz, year, doy, hr, mn, hour_utc, &
                         zen_ref, ghi, dni_ref, dhi, temp, pres_mbar

            ! Score only sunlit records with real irradiance. DISC is undefined
            ! at night and meaningless where the sun is near the horizon.
            if (zen_ref >= 85.0_wp) cycle
            if (ghi <= 0.0_wp) cycle

            ! NSRDB pressure is in millibar; the module expects pascals.
            pres_pa = pres_mbar*100.0_wp

            ! hour_utc counts from 00:00 UTC on 1 January, so pass it as the hour
            ! of the year with a zero time zone. tisr < 0 tells DISC to compute
            ! the top-of-atmosphere irradiance itself.
            call calc_dni_disc(dni,hour_utc,lat,lon,0.0_wp,ghi,pres_pa,-1.0_wp, &
                               print_diagnostics=.FALSE.)

            err      = dni-dni_ref
            sum2     = sum2 + err*err
            sum_ref2 = sum_ref2 + dni_ref*dni_ref
            bias     = bias + err
            n = n + 1
        end do

        close(unit)

        if (n == 0) then
            call fail("DNI vs NREL", "no scorable records read")
            return
        end if

        write(*,"(a,i0)")            "     scored records      = ", n
        write(*,"(a,f8.2,a)")        "     RMSE(DNI)           = ", sqrt(sum2/real(n,wp)), " W/m2"
        write(*,"(a,f8.2,a)")        "     mean bias           = ", bias/real(n,wp), " W/m2"
        write(*,"(a,f8.2,a)")        "     RMS of measured DNI = ", sqrt(sum_ref2/real(n,wp)), " W/m2"

        call check("DNI vs NREL (RMSE)", sqrt(sum2/real(n,wp)), tol_rmse, &
                   "RMSE of estimated DNI [W/m2]")

    end subroutine test_dni_vs_nrel

    subroutine check(name, value, tol, what)

        implicit none

        character(len=*), intent(IN) :: name, what
        real(wp),         intent(IN) :: value, tol

        if (.not. (value <= tol)) then
            call fail(name, trim(what)//" exceeds tolerance")
        else
            call pass(name)
        end if

    end subroutine check

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

end program test_irradiance
