
! To compile:
! gfortran -o test_insol.x insolation.f90 test_insol_pt.f90

program test_dni 

    use irradiance  
    use ncio 

    character(len=56) :: station
    integer :: year0, year1
    integer :: year 

    integer :: narg 
    character(len=56) :: loc_now 
    character(len=56) :: year_now 
    real(wp) :: time_zone


    character(len=56) :: run_case 


    run_case = "era5"   ! "point", "nrel", "era5" 


    ! ====================================================================

    select case(trim(run_case))


    case("point")

        ! Run a diagnostic test with specific parameter values for testing
        ! (hard-coded inside routine below) 
        call calc_point_dni_test()

    case("nrel") 

        ! Run a calculation of NREL station data for a given set of years
        !station = "667526"
        !station = "677404"
        station = "848945"
        year0 = 1998
        year1 = 2020 

        do year = year0, year1

            call calc_nrel_station(fldr="data/nrel/NRELNSRDBViewerData", &
                                        stn=station,year=year,print_diagnostics=.FALSE.)

        end do 

    case("era5") 

        ! ===== Load location and year from command line arguments =====
        narg = command_argument_count()

        if (narg .ne. 2) then 
            write(*,*) "test_dni:: Error: location and year must be provided as an argument."
            stop 
        end if 

        call get_command_argument(1,loc_now)
        call get_command_argument(2,year_now)

        ! No time zone needed for ERA5 data, all times are UTC and consistent
        ! with a UTC time zone, it seems (ajr, 2022-02-09). This was checked by 
        ! using the solar zenith angle calculator here:
        ! https://midcdmz.nrel.gov/solpos/solpos.html
        ! To get zenith angles < 90 in the time period when GHI coming 
        ! out of the ERA5 data was > 0, a time_zone=0 needs to be specified. 
        ! If the "right" time zone is specified, then the zenith_angle is 
        ! offset by that amount and becomes inconsistent with the times at 
        ! which GHI > 0. 
        ! Note: it appears that the NREL station data does need to include
        ! this time_zone offset, despite apparently having time in UTC. 
if (.TRUE.) then 

        time_zone = 0.0_wp 

else 

        ! Define time zone depending on location 
        select case(trim(loc_now))

        case("fuentes")
            time_zone = 1.0_wp 
        case("dubai") 
            time_zone = 4.0_wp 

        case DEFAULT 

            write(5,*) "Error: time_zone not yet defined in code for: ", trim(loc_now)
            stop 

        end select

end if 

        call calc_era5_dni(loc_now,year_now,time_zone)

    end select


contains 

    subroutine calc_point_dni_test()
        ! Diagnostic testing, e.g., to compare with values in NREL Excel file smalldisc.xls
        ! https://www.nrel.gov/grid/solar-resource/disc.html

        implicit none

        ! Local variables 
        real(wp) :: lat, lon, time_zone
        real(wp) :: hour
        real(wp) :: zenith_angle 
        real(wp) :: pres 
        real(wp) :: ghi
        real(wp) :: dni

        real(wp) :: ghi_now 
        real(wp) :: ghi_m1 

        ! === Specify input variable values ===
        
        lon       = -75.0
        lat       =  40.0
        time_zone = -5.0 
        hour      =  9.0
        ghi       = 100.0
        pres      = 840.0 * 100.0  ! [mb] => [Pa]

        ghi_now   = ghi
        ghi_m1    = ghi

        ! === Calculate Zenith Angle and DNI ===
        ! (with print_diagnostics=.TRUE., zenith_angle and dni calculation info will be printed to screen)

        !call calc_zenith_angle(zenith_angle,hour,lat,lon,time_zone,print_diagnostics=.TRUE.)
        !call calc_dni_disc(dni,hour,lat,lon,time_zone,ghi,pres,-1.0_wp,print_diagnostics=.TRUE.)
        
        ! Interpolate time over one hour
        call calc_dni_60min(dni,hour,lat,lon,time_zone,ghi_now,ghi_m1,pres,-1.0_wp,print_diagnostics=.TRUE.)

        stop " === Done. === "

        return

    end subroutine calc_point_dni_test

    subroutine calc_nrel_station(fldr,stn,year,print_diagnostics)

        implicit none

        character(len=*), intent(IN) :: fldr
        character(len=*), intent(IN) :: stn
        integer, intent(IN) :: year
        logical, intent(IN) :: print_diagnostics 

        ! Local variables 
        character(len=512) :: filepath_info
        character(len=512) :: filepath_in
        character(len=512) :: filepath_out
        character(len=4)   :: year_str
        integer  :: id 
        real(wp) :: lat, lon, time_zone, elev
        real(wp) :: hour
        real(wp) :: zenith_angle 
        real(wp) :: temp 
        real(wp) :: pres 
        real(wp) :: ghi
        real(wp) :: dhi 
        real(wp) :: dni
        real(wp) :: ghi_now 
        real(wp) :: ghi_m1 
        real(wp) :: zenith_angle_pred
        real(wp) :: dni_pred 

        integer  :: io_unit, n, stat
        integer  :: io_unit_out   
        character(len=56) :: tmp(10) 

        ! Define filenames 
        write(year_str,"(i4)") year 
        filepath_info = trim(fldr)//"/"//trim(stn)//"-info.txt"
        filepath_in   = trim(fldr)//"/"//trim(stn)//"-"//year_str//".txt"
        filepath_out  = trim(fldr)//"/"//trim(stn)//"-dni-"//year_str//".txt"
        
        if (print_diagnostics) then 
            write(*,*) "filepath_info: ", trim(filepath_info)
            write(*,*) "filepath_in:   ", trim(filepath_in)
            write(*,*) "filepath_out:  ", trim(filepath_out)
        end if

        ! === Read in the station information first ===

        io_unit = 99 
        open(io_unit,status="old",file=trim(filepath_info))
        read(io_unit,*) tmp(1:5)            ! Read header
        read(io_unit,*) id, lat, lon, time_zone, elev 
        close(io_unit) 

        if (print_diagnostics) then
            write(*,*) "id:        ", id 
            write(*,*) "lat:       ", lat 
            write(*,*) "lon:       ", lon 
            write(*,*) "time_zone: ", time_zone 
            write(*,*) "elev:      ", elev 
            write(*,*) 
        else 
            write(*,"(a,i10.0)") trim(stn), year
        end if

        ! === Read in the input data, calculate dni and go to next line ===

        ! Define output file 
        io_unit_out = 100 
        open(io_unit_out,status="unknown",file=trim(filepath_out))
        write(io_unit_out,"(a,2x,a,2x,a)") "hour", "zenith_angle", "dni" 

        io_unit = 99 
        open(io_unit,status="old",file=trim(filepath_in))
        read(io_unit,*) tmp(1:7)            ! Read header

        if (print_diagnostics) then
            write(*,*) "hour, zenith_angle (in/pred), dni (in/pred): "
        end if 

        do n = 1, 1000000

            ! Read in this line 
            read(io_unit,*,iostat=stat) hour, ghi, dni, dhi, zenith_angle, temp, pres
            
            if (stat .ne. 0) then 
                ! End of file, exit 
                exit
            end if 

            ! Store ghi 
            if (n .eq. 1) then 
                ghi_m1 = ghi
            else 
                ghi_m1 = ghi_now 
            end if 

            ghi_now = ghi

            ! Perform calculations
            call calc_zenith_angle(zenith_angle_pred,hour,lat,lon,time_zone,print_diagnostics=.FALSE.) 
            !call calc_dni_disc(dni_pred,hour,lat,lon,time_zone,ghi,pres,-1.0_wp,print_diagnostics=.FALSE.)
            
            ! Interpolate time over one hour
            call calc_dni_60min(dni_pred,hour,lat,lon,time_zone,ghi_now,ghi,pres,-1.0_wp,print_diagnostics=.FALSE.)

            if (print_diagnostics) then
                ! Print summary 
                write(*,"(f10.0,4f14.2)") hour, zenith_angle, zenith_angle_pred, dni, dni_pred
            end if 

            ! Write to file too 
            write(io_unit_out,"(f10.0,3f14.2)") hour, zenith_angle_pred, dni_pred 

        end do 

        close(io_unit) 
        close(io_unit_out) 

        return

    end subroutine calc_nrel_station

    subroutine calc_era5_dni(loc_now,year_now,time_zone)

        implicit none 

        character(len=*), intent(IN) :: loc_now 
        character(len=*), intent(IN) :: year_now 
        real(wp),         intent(IN) :: time_zone 

        ! Local variables 
        character(len=512) :: filename_in 
        character(len=512) :: filename_out  

        integer  :: i, j, k, nx, ny, nt
        real(wp) :: fdir_now
        real(wp) :: fdir_m1
        real(wp) :: ghi_now
        real(wp) :: ghi_m1
        
        real(wp), allocatable :: lon(:) 
        real(wp), allocatable :: lat(:) 
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: hour(:) 
        real(wp), allocatable :: tisr(:,:,:) 
        real(wp), allocatable :: ghi(:,:,:) 
        real(wp), allocatable :: fdir(:,:,:) 
        real(wp), allocatable :: pres(:,:,:) 
        real(wp), allocatable :: zenith_angle(:,:,:) 
        real(wp), allocatable :: dni(:,:,:) 

        real(wp), allocatable :: dni_mod(:,:,:) 
        real(wp), allocatable :: fdir_mod(:,:,:) 

        real(wp), parameter :: radiation_tol = 1e-2 

        ! === Define filenames === 

        !filename_in  = "data/era5/southspain/era5_southspain_2019.nc" 
        !filename_out = "data/era5/southspain/era5_southspain_2019_dni.nc"
        
        filename_in  = "data/era5/"//trim(loc_now)//"/era5_"//trim(loc_now)//"_"//trim(year_now)//".nc" 
        filename_out = "data/era5/"//trim(loc_now)//"/era5-dni_"//trim(loc_now)//"_"//trim(year_now)//".nc" 
        
        ! === Get data dimensions === 

        nx = nc_size(filename_in,"longitude")
        ny = nc_size(filename_in,"latitude")
        nt = nc_size(filename_in,"time") 

        allocate(lon(nx))
        allocate(lat(ny))
        allocate(time(nt))
        allocate(hour(nt))
        allocate(tisr(nx,ny,nt))
        allocate(ghi(nx,ny,nt)) 
        allocate(fdir(nx,ny,nt)) 
        allocate(pres(nx,ny,nt))
        allocate(zenith_angle(nx,ny,nt)) 
        allocate(dni(nx,ny,nt))

        allocate(dni_mod(nx,ny,nt))
        allocate(fdir_mod(nx,ny,nt))

        ! === Read in dataset === 

        call nc_read(filename_in,"longitude",lon)
        call nc_read(filename_in,"latitude",lat)
        call nc_read(filename_in,"time",time)
        call nc_read(filename_in,"tisr",tisr)       ! TOA Incident solar radiation
        call nc_read(filename_in,"ssrd",ghi)        ! Surface global horizontal irradiance == [GHI]
        call nc_read(filename_in,"fdir",fdir)       ! Surface direct horizontal irradiance == [DNI*cos(zenith_angle)]
        call nc_read(filename_in,"sp",pres)         ! Surface pressure [Pa]

        ! Define the hour_of_year, by subtracting minimum value loaded 
        ! from NetCDF file. 
        ! Note: the first hour of the year should be 00:00:00 UTC.
        ! Therefore, hour of the year (hour) should start with hour=0
        hour = time - minval(time)                  ! First hour == 0.0

        tisr = tisr / 3600.0_wp                     ! [J/m2] => [W/m2]
        ghi  = ghi  / 3600.0_wp                     ! [J/m2] => [W/m2]
        fdir = fdir / 3600.0_wp                     ! [J/m2] => [W/m2] 

        where(ghi  .lt. radiation_tol) ghi  = 0.0_wp 
        where(tisr .lt. radiation_tol) tisr = 0.0_wp 
        where(fdir .lt. radiation_tol) fdir = 0.0_wp 

        ! ajr: testing!!!

        write(*,*) "range(lon):  ", minval(lon),  maxval(lon)
        write(*,*) "range(lat):  ", minval(lat),  maxval(lat)
        write(*,*) "range(time): ", minval(time), maxval(time)
        write(*,*) "range(tisr): ", minval(tisr), maxval(tisr)
        write(*,*) "range(ghi):  ", minval(ghi),  maxval(ghi)
        write(*,*) "range(fdir): ", minval(fdir), maxval(fdir)
        write(*,*) "range(sp):   ", minval(pres), maxval(pres)


        ! === Calculate Zenith Angle and DNI ===

        zenith_angle = 0.0_wp 
        dni          = 0.0_wp 
        dni_mod      = 0.0_wp 
        fdir_mod     = 0.0_wp 

        do i = 1, nx 
        do j = 1, ny 

            do k = 1, nt
                call calc_zenith_angle(zenith_angle(i,j,k),hour(k),lat(j),lon(i),time_zone,print_diagnostics=.FALSE.) 
                !call calc_dni_disc(dni(i,j,k),hour(k),lat(j),lon(i),time_zone,ghi(i,j,k),pres(i,j,k),-1.0_wp,print_diagnostics=.FALSE.)
                
                ! Get current and previous values of ghi and fdir 

                ghi_now = ghi(i,j,k) 
                if (k .eq. 1) then 
                    ghi_m1 = ghi(i,j,k)
                else 
                    ghi_m1 = ghi(i,j,k-1)
                end if 

                fdir_now = fdir(i,j,k) 
                if (k .eq. 1) then 
                    fdir_m1 = fdir(i,j,k)
                else 
                    fdir_m1 = fdir(i,j,k-1)
                end if 

                ! Calculate DNI directly from FDIR and the effective zenith angle
                call calc_dni_from_fdir(dni(i,j,k),hour(k),lat(j),lon(i),time_zone,fdir_now,print_diagnostics=.FALSE.)
                
                ! Calculate DNI directly from FDIR and the zenith angle over 60min
                !call calc_dni_from_fdir_60min(dni(i,j,k),hour(k),lat(j),lon(i),time_zone,fdir_now,fdir_m1,print_diagnostics=.FALSE.)
                
                ! Calculate DNI from DISC model interpolated over 60min
                call calc_dni_60min(dni_mod(i,j,k),hour(k),lat(j),lon(i),time_zone,ghi_now,ghi_m1,pres(i,j,k),-1.0_wp,print_diagnostics=.FALSE.)

                if (zenith_angle(i,j,k) .lt. 89.0) then 
                    ! Calculate fdir=dni*cos(z) too, for comparison with ERA5 FDIR variable
                    fdir_mod(i,j,k) = dni_mod(i,j,k)*cos(zenith_angle(i,j,k)*pi/180.0_wp)
                end if 

            end do

        end do 
        end do 


        ! === Write to output file ===

        call nc_create(filename_out)
        call nc_write_dim(filename_out,"longitude",x=lon)
        call nc_write_dim(filename_out,"latitude",x=lat)
        call nc_write_dim(filename_out,"time",x=hour)
        
        call nc_write(filename_out,"zenith_angle",  zenith_angle,   dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"tisr",          tisr,           dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"ghi",           ghi,            dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"fdir",          fdir,           dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"sp",            pres,           dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"dni",           dni,            dim1="longitude",dim2="latitude",dim3="time")

        call nc_write(filename_out,"dni_mod",       dni_mod,        dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"fdir_mod",      fdir_mod,       dim1="longitude",dim2="latitude",dim3="time")
        
        return

    end subroutine calc_era5_dni

end program test_dni