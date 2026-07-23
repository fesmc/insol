module irradiance
    ! Surface irradiance decomposition for renewable-energy applications:
    ! estimating direct normal irradiance (DNI) from global horizontal
    ! irradiance (GHI), chiefly via the DISC model (Maxwell 1987).
    !
    ! The solar geometry -- declination, equation of time, zenith angle, the
    ! Sun-Earth distance factor -- comes from insol_orbit, the same orbital
    ! solution the insolation core uses. Earlier this module carried its own
    ! present-day Fourier series (Spencer 1971) with a fixed 365-day calendar
    ! and a different solar constant; sharing insol_orbit removes that
    ! duplication and lets the estimate use real dates and orbital elements.
    !
    ! The decomposition models themselves (DISC and the rest) are empirical and
    ! kept in the working precision wp; their coefficients are fits to data, so
    ! double precision would add no meaningful accuracy. The geometry they rely
    ! on is computed in double precision inside insol_orbit.

    use insol_defs, only : sp, dp, wp, pi, &
                           degrees_to_radians, radians_to_degrees, S0_default
    use insol_orbit, only : orbit_class, calc_solar_position_date

    implicit none

    ! Missing value and aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,wp)
    real(wp), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT

    ! Solar constant [W/m2]. 1365 is the classic value used across the
    ! insolation core (insol_defs). The modern best estimate is closer to
    ! 1361 W/m2 (Kopp & Lean 2011); this module previously used 1370, the value
    ! assumed when the DISC coefficients were originally fitted, so DNI
    ! estimates shift slightly with the change to 1365.
    real(wp), parameter :: S0 = real(S0_default,wp)

    private
    public :: sp, dp, wp, MISSING_VALUE, MV, pi, degrees_to_radians, radians_to_degrees
    public :: calc_dni_from_fdir
    public :: calc_dni_from_fdir_60min
    public :: calc_dni_60min
    public :: calc_dni_disc
    public :: calc_zenith_angle

contains
    
    subroutine calc_dni_from_fdir(orb,dni,year,hour,lat,lon,fdir,print_diagnostics)
        ! DNI from the direct (beam) horizontal flux fdir, as dni = fdir/cos(z),
        ! using the Blanc & Wald (2016) effective zenith angle over the hour.
        ! `hour` is hours since 00:00 UTC on 1 January of `year`.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(wp),          intent(OUT) :: dni        ! [W/m2]
        integer,           intent(IN)  :: year
        real(wp),          intent(IN)  :: hour       ! Hours since 1 Jan 00:00 UTC
        real(wp),          intent(IN)  :: lat        ! [Degrees North]
        real(wp),          intent(IN)  :: lon        ! [Degrees East]
        real(wp),          intent(IN)  :: fdir       ! [W/m2]
        logical,           intent(IN)  :: print_diagnostics

        ! Local variables
        real(wp) :: zenith_angle_eff
        real(wp) :: cosz

        ! Calculate effective zenith angle for this hour
        call calc_effective_zenith_angle(orb,zenith_angle_eff,year,hour,lat,lon)

        ! Calculate DNI 
        if (zenith_angle_eff .lt. 89.0_wp) then 

            cosz = cos(zenith_angle_eff*pi/180.0_wp)
            dni = fdir / cosz 

        else

            dni = 0.0_wp 

        end if

        return

    end subroutine calc_dni_from_fdir

    subroutine calc_dni_from_fdir_60min(orb,dni,year,hour,lat,lon,fdir,fdir_m1,print_diagnostics)
        ! Hourly-mean DNI from the direct horizontal flux, integrating the solar
        ! position over the preceding hour. `hour` is hours since 00:00 UTC on
        ! 1 January of `year`; fdir_m1 is the flux one hour earlier.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(wp),          intent(OUT) :: dni        ! [W/m2]
        integer,           intent(IN)  :: year
        real(wp),          intent(IN)  :: hour       ! Hours since 1 Jan 00:00 UTC
        real(wp),          intent(IN)  :: lat        ! [Degrees North]
        real(wp),          intent(IN)  :: lon        ! [Degrees East]
        real(wp),          intent(IN)  :: fdir       ! [W/m2]
        real(wp),          intent(IN)  :: fdir_m1    ! [W/m2]
        logical,           intent(IN)  :: print_diagnostics

        ! Local variables 
        integer  :: k, nmin
        real(wp) :: hour_val
        real(wp) :: fdir_now
        real(wp) :: zenith_angle_now
        real(wp) :: cosz 

        real(wp), allocatable :: dni_vals(:)


        ! Do 60 minutes of calculations by default
        nmin = 60

        allocate(dni_vals(nmin))

        dni_vals = 0.0_wp 

        do k = 1, nmin 

if (.TRUE.) then 
    ! This procedure seems more consistent with defintion of fdir in the file,
    ! which is the cumulative fdir over the whole hour (which we then divide by sec_per_hr). 
    ! So we can assume fdir is linearly changing over the past hour to calculate a more
    ! accurate value of dni for each minute of the last hour. 
    ! This approach doesn't work so well for NREL data using GHI data, which 
    ! perhaps is not cumulative but direct measurements. 

            ! Get the current hour value
            hour_val = hour - real(nmin,wp)/60_wp + real(k,wp)/60_wp

            ! Linear interpolation to get fdir now
            fdir_now = fdir_m1 + (fdir-fdir_m1) * real(k,wp) / real(nmin,wp)

else 
            hour_val = hour + (-(nmin/2.0_wp) + real(k,wp)) /60_wp
            fdir_now = fdir 
end if 

            ! Get the current zenith angle
            call calc_zenith_angle(orb,zenith_angle_now,year,hour_val,lat,lon)

            if (zenith_angle_now .lt. 89.0_wp) then

                ! Get cos(zenith_angle)
                cosz = cos(zenith_angle_now*pi/180.d0)

                dni_vals(k) = fdir_now / cosz 

            else 

                dni_vals(k) = 0.0_wp 
            
            end if 

            !write(*,*) "calc_dni_from_fdir_60min: ", hour_vals(k), dni_vals(k)

        end do

        ! Return the mean DNI value for this hour
        dni = sum(dni_vals)/real(nmin,wp)

        return

    end subroutine calc_dni_from_fdir_60min

    subroutine calc_effective_zenith_angle(orb,zenith_angle_eff,year,hour,lat,lon)
        ! Calculate the effective zenith angle for a given hour
        ! following the algorithm A3 of:
        ! Blanc, P. and Wald, L.: On the effective solar zenith and
        ! azimuth angles to use with measurements of hourly
        ! irradiation, Adv. Sci. Res., 13, 1-6,
        ! https://doi.org/10.5194/asr-13-1-2016, 2016.

        ! Get the average zenith angle of the previous hour, only
        ! using values during daylight hours. `hour` is hours since 00:00 UTC
        ! on 1 January of `year`.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(wp),          intent(OUT) :: zenith_angle_eff
        integer,           intent(IN)  :: year
        real(wp),          intent(IN)  :: hour        ! Hours since 1 Jan 00:00 UTC
        real(wp),          intent(IN)  :: lat
        real(wp),          intent(IN)  :: lon

        ! Local variables
        integer  :: k, nmin
        real(wp) :: hour_now
        real(wp) :: wt_tot

        real(wp), allocatable :: za(:)
        real(wp), allocatable :: wt(:)

        ! Do 60 minutes of calculations by default
        nmin = 60

        allocate(za(nmin))
        allocate(wt(nmin))

        za = 0.0_wp
        wt = 0.0_wp

        do k = 1, nmin

            ! Get the current hour value (from t-1 to t)
            hour_now = hour - real(nmin,wp)/60_wp + real(k,wp)/60_wp

            ! Calculate the zenith angle at this time
            call calc_zenith_angle(orb,za(k),year,hour_now,lat,lon)

            if (za(k) .lt. 89.0_wp) then
                ! Angle is less than 90deg, include this minute in the weighting
                wt(k) = 1.0_wp 
            end if 

        end do 

        ! Get total weight 
        wt_tot = sum(wt) 


        if (wt_tot .gt. 0.0_wp) then 
            ! Some daylight minutes were found, calculate the 
            ! effective zenith angle as the weighted mean 
            ! over the hour. 

            zenith_angle_eff = sum(wt*za) / wt_tot 

        else
            ! Simply fill in angle with the hour value

            call calc_zenith_angle(orb,zenith_angle_eff,year,hour,lat,lon)

        end if

        return

    end subroutine calc_effective_zenith_angle

    subroutine calc_dni_60min(orb,dni,year,hour,lat,lon,ghi,ghi_m1,pres,print_diagnostics)
        ! Hourly-mean DNI from GHI: DISC applied minute by minute over the hour
        ! with the solar position advancing, then averaged. `hour` is hours since
        ! 00:00 UTC on 1 January of `year`.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(wp),          intent(OUT) :: dni        ! [W/m2]
        integer,           intent(IN)  :: year
        real(wp),          intent(IN)  :: hour       ! Hours since 1 Jan 00:00 UTC
        real(wp),          intent(IN)  :: lat        ! [Degrees North]
        real(wp),          intent(IN)  :: lon        ! [Degrees East]
        real(wp),          intent(IN)  :: ghi        ! [W/m2] Current hour GHI
        real(wp),          intent(IN)  :: ghi_m1     ! [W/m2] Previous hour GHI
        real(wp),          intent(IN)  :: pres       ! [Pa]
        logical,           intent(IN)  :: print_diagnostics

        ! Local variables 
        integer  :: k, nmin
        real(wp) :: hour_val
        real(wp) :: ghi_now

        real(wp), allocatable :: dni_vals(:)


        ! Do 60 minutes of calculations by default
        nmin = 60

        allocate(dni_vals(nmin))

        dni_vals = 0.0_wp 

        do k = 1, nmin 
            

if (.FALSE.) then 
    ! Assume ghi is for hour, and ghi_m1 is for hour-1
    ! This method does worse for NREL station data! 

            ! Get the current hour value
            hour_val = hour - real(nmin,wp)/60_wp + real(k,wp)/60_wp

            ! Linear interpolation to get fdir now
            ghi_now = ghi_m1 + (ghi-ghi_m1) * real(k,wp) / real(nmin,wp)

            call calc_dni_disc(orb,dni_vals(k),year,hour_val,lat,lon,ghi_now,pres,print_diagnostics)

else
    ! Assume current ghi is valid for the middle of the hour
    ! This method is better for NREL station data, continuing with this method for now,
    ! as it also doesn't change the ERA5 results compared to the FDIR method

            hour_val = hour + (-(nmin/2.0_wp) + real(k,wp)) /60_wp
            call calc_dni_disc(orb,dni_vals(k),year,hour_val,lat,lon,ghi,pres,print_diagnostics)

end if 

            !write(*,*) "calc_dni_60min: ", hour_val, dni_vals(k)

        end do

        ! Return the mean DNI value for this hour
        dni = sum(dni_vals)/real(nmin,wp)

        return

    end subroutine calc_dni_60min

    subroutine calc_dni_disc(orb,dni,year,hour,lat,lon,ghi,pres,print_diagnostics)
        ! Direct Insolation Simulation Code (DISC) model
        ! Maxwell (1987), SERI Technical Report
        ! https://www.nrel.gov/grid/solar-resource/disc.html
        !
        ! `hour` is hours since 00:00 UTC on 1 January of `year`. Solar geometry
        ! and the top-of-atmosphere irradiance come from insol_orbit; only the
        ! empirical clear-sky transmittance below is DISC's own.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(wp),          intent(OUT) :: dni        ! [W/m2]
        integer,           intent(IN)  :: year
        real(wp),          intent(IN)  :: hour       ! Hours since 1 Jan 00:00 UTC
        real(wp),          intent(IN)  :: lat        ! [Degrees North]
        real(wp),          intent(IN)  :: lon        ! [Degrees East]
        real(wp),          intent(IN)  :: ghi        ! [W/m2]
        real(wp),          intent(IN)  :: pres       ! [Pa]
        logical,           intent(IN)  :: print_diagnostics

        ! Local variables
        real(wp) :: pressure

        real(wp) :: I0                      ! TOA normal irradiance
        real(wp) :: zenith_angle
        real(wp) :: cos_zenith_angle
        real(wp) :: m_air_sl, m_air, kt, kn, knc, del_kn
        real(wp) :: a, b, c
        real(wp) :: kt2, kt3

        real(wp) :: dhi

        real(wp), parameter :: standard_pressure = 101325_wp

        ! Limit kt to a reasonable value, usually between 0 and 1 (see comments below).
        real(wp), parameter :: kt_max  = 1.0_wp
        real(wp), parameter :: ghi_tol = 1e-2

        if (pres .gt. 0.0_wp) then
            ! Use pressure provided as an argument
            pressure = pres
        else
            ! Assume standard pressure for now [Pa]
            pressure = standard_pressure
        end if

        ! Solar geometry and TOA normal irradiance from the shared orbit
        call calc_geometry(orb,year,hour,lat,lon,zenith_angle,cos_zenith_angle,I0)

        ! Calculate air mass
        if (zenith_angle < 87) then

            ! Use Bird & Hulstrom coefficient values (DEFAULT)
            ! Bird, R., Hulstrom, R., A Simplified Clear Sky Model for Direct and Diffuse Insolation on Horizontal Surfaces, 1981
            ! NREL SERI/TP-642-761
            m_air_sl = 1.0 / (cos_zenith_angle + 0.15 * (93.885-zenith_angle)**(-1.253))

            ! Use Kasten / Sandia Version
            ! Kasten, F. Discussions on the Relative Air Mass, Light Research Technology 25, 129.
            !m_air_sl = 1.0 / (cos_zenith_angle + 0.5057 * (96.080-zenith_angle)**(-1.634))

        else
            m_air_sl = 0.0
        end if

        ! Scale by pressure too to get actual estimated air mass
        m_air = m_air_sl * (pressure / standard_pressure)

        ! Clearness index: GHI over the horizontal TOA irradiance.
        call calc_clearness_index(kt,cos_zenith_angle,I0,ghi,kt_max)


        ! Calculate powers of kt
        kt2 = kt*kt
        kt3 = kt2*kt

        ! Calculate coefficients A, B and C 
        if (kt .gt. 0.6) then 

            a = -5.743 +  21.77*kt -  27.49*kt2 + 11.56*kt3
            b =  41.40 - 118.50*kt +  66.05*kt2 + 31.90*kt3
            c = -47.01 + 184.20*kt - 222.00*kt2 + 73.81*kt3

        else if (kt .gt. 0.0) then

            a =  0.512 - 1.560*kt + 2.286*kt2 - 2.222*kt3
            b =  0.370 + 0.962*kt
            c = -0.280 + 0.932*kt - 2.048*kt2

        else 

            a = 0.0
            b = 0.0 
            c = 0.0 

        end if 

        ! Calculate the clear-sky atmospheric transmissivity 

        if (kt .gt. 0.0) then 
            knc    = 0.886 - 0.122*m_air + 0.0121*m_air**2.0 &
                            - 0.000653*m_air**3.0 + 0.000014*m_air**4.0
            
            del_kn = a + b * exp(c * m_air)
            
            kn     = max(knc - del_kn,0.0_wp)

            ! Finally calculate DNI 
            dni = I0*kn 

        else 

            knc    = 0.0 
            del_kn = 0.0 
            kn     = 0.0 

            ! Set DNI to zero
            dni = 0.0 

        end if

        ! Further corrections
        if (ghi .lt. 0.0 .or. dni .lt. 0.0) then

            dni = 0.0

        end if

        if (print_diagnostics) then

            write(*,*) "=== calc_dni_disc ==="
            write(*,*)

            write(*,*) "year             = ", year
            write(*,*) "hour             = ", hour
            write(*,*) "lat              = ", lat
            write(*,*) "lon              = ", lon
            write(*,*) "zenith_angle     = ", zenith_angle
            write(*,*) "cos_zenith_angle = ", cos_zenith_angle
            write(*,*) "pressure         = ", pressure
            write(*,*) "pressure ratio   = ", (pressure / standard_pressure)

            write(*,*) "m_air            = ", m_air
            write(*,*) "ghi              = ", ghi
            write(*,*) "I0               = ", I0
            write(*,*) "kt               = ", kt
            write(*,*) "a                = ", a 
            write(*,*) "b                = ", b 
            write(*,*) "c                = ", c 
            write(*,*) "del_kn           = ", del_kn
            write(*,*) "kn               = ", kn 
            write(*,*) "knc              = ", knc 
            write(*,*) "dni              = ", dni
            write(*,*) 

        end if 

        return 

    end subroutine calc_dni_disc

    subroutine calc_dni_boland2013(dni,ghi)
        ! Logistic model for DNI calculation
        ! Boland et al. (2013, Renewable and Sustainable Energy Reviews)
        ! https://www.sciencedirect.com/science/article/pii/S1364032113005637

        implicit none

        real(wp), intent(OUT) :: dni 
        real(wp), intent(IN)  :: ghi 

        ! Local variables 
        real(wp) :: p 
        real(wp) :: day_angle 
        real(wp) :: cos_zenith_angle
        real(wp) :: I0 
        real(wp) :: I0_horizontal
        real(wp) :: kt 

        real(wp), parameter :: N = 0.006
        real(wp), parameter :: M = 4.38 

        real(wp), parameter :: b1 = 7.75
        real(wp), parameter :: b2 = 1.185
        real(wp), parameter :: b3 = 1.05
        real(wp), parameter :: b4 = 0.004
        real(wp), parameter :: b5 = -0.003



        
        !p = -b1*kt - b2*ast - b3*alpha - b4*kt_day + b5*psi 

        dni = (N*M) / (N+(M-N)*exp(p))

        return

    end subroutine calc_dni_boland2013

    subroutine calc_dni_from_diffuse(dni,ghi,dhi,zenith_angle)

        implicit none

        real(wp), intent(OUT) :: dni 
        real(wp), intent(IN)  :: ghi
        real(wp), intent(IN)  :: dhi 
        real(wp), intent(IN)  :: zenith_angle

        ! Boland et al (2013), Eq. 7 
        dni = (ghi - dhi) / cos(zenith_angle*pi/180.0_wp)

        return

    end subroutine calc_dni_from_diffuse

    subroutine calc_dhi_brl(dhi,kt,ghi)
        ! Logistic model for DHI calculation
        ! Ridley et al. (2010, Renewable Energy)
        ! https://www.sciencedirect.com/science/article/pii/S0960148109003012

        implicit none

        real(wp), intent(OUT) :: dhi 
        real(wp), intent(IN)  :: kt 
        real(wp), intent(IN)  :: ghi

        ! Local variables 
        real(wp) :: p 
        real(wp) :: diffuse_frac 

        ! Simple logistic function that depends only on kt 
        ! Ridley et al (2010), Eq. 3

        p   = -5.0033 + 8.6025*kt 
        diffuse_frac = 1.0_wp / (1.0_wp + exp(p))

        dhi = diffuse_frac*ghi

        return

    end subroutine calc_dhi_brl

    subroutine calc_clearness_index(kt,cos_zenith_angle,I0,ghi,kt_max)
        ! Clearness index: measured GHI over the horizontal top-of-atmosphere
        ! irradiance, cos(zenith) * I0. The TOA irradiance is supplied by the
        ! caller from the shared orbit rather than recomputed here.

        implicit none

        real(wp), intent(OUT) :: kt
        real(wp), intent(IN)  :: cos_zenith_angle
        real(wp), intent(IN)  :: I0                 ! TOA normal irradiance [W/m2]
        real(wp), intent(IN)  :: ghi
        real(wp), intent(IN)  :: kt_max

        ! Local variables
        real(wp) :: I0_horizontal

        I0_horizontal = cos_zenith_angle*I0

        ! As the zenith angle approaches 90 degrees (sunrise/sunset) the
        ! horizontal TOA irradiance vanishes and the ratio blows up, so guard
        ! the division and cap kt at a physical maximum (~1).
        if (I0_horizontal .gt. 0.0_wp) then
            kt = min(ghi/I0_horizontal,kt_max)
        else
            kt = 0.0_wp
        end if

        return

    end subroutine calc_clearness_index

    subroutine calc_zenith_angle(orb,zenith_angle,year,hour,lat,lon)
        ! Solar zenith angle [degrees] from the shared orbital solution.
        !
        ! `hour` is hours elapsed since 00:00 UTC on 1 January of `year`, the
        ! convention used throughout this module. Longitude is degrees east.
        ! Declination, equation of time and the hour angle are all handled by
        ! insol_orbit; there is no time-zone argument, as the timestamp is UTC.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        real(wp),          intent(OUT) :: zenith_angle
        integer,           intent(IN)  :: year
        real(wp),          intent(IN)  :: hour       ! Hours since 1 Jan 00:00 UTC
        real(wp),          intent(IN)  :: lat        ! Degrees north
        real(wp),          intent(IN)  :: lon        ! Degrees east

        ! Local variables
        real(wp) :: cos_zenith, I0

        call calc_geometry(orb,year,hour,lat,lon,zenith_angle,cos_zenith,I0)

        return

    end subroutine calc_zenith_angle

    subroutine calc_geometry(orb,year,hour,lat,lon,zenith_angle,cos_zenith,I0)
        ! Solar geometry for a location and time: the zenith angle, its cosine,
        ! and the top-of-atmosphere normal irradiance I0 = S0 (a/r)^2. All three
        ! come from one call to the shared orbital solution.

        implicit none

        type(orbit_class), intent(IN)  :: orb
        integer,           intent(IN)  :: year
        real(wp),          intent(IN)  :: hour, lat, lon
        real(wp),          intent(OUT) :: zenith_angle    ! [degrees]
        real(wp),          intent(OUT) :: cos_zenith
        real(wp),          intent(OUT) :: I0              ! [W/m2]

        ! Local variables
        real(dp) :: cosz, dist2

        ! Anchor to 1 January and let the elapsed hours carry the day of year.
        call calc_solar_position_date(orb,cosz,year,1,1,real(hour,dp), &
                                      real(lat,dp),real(lon,dp),dist2=dist2)

        cos_zenith   = real(cosz,wp)
        zenith_angle = real(acos(max(-1.0_dp,min(1.0_dp,cosz)))*radians_to_degrees,wp)
        I0           = real(S0*dist2,wp)

        return

    end subroutine calc_geometry

end module irradiance
