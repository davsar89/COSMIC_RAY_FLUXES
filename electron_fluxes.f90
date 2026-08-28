program electron_fluxes
   use, intrinsic :: iso_fortran_env, only: real64, error_unit, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use integration, only: log_simpson_integration
   use utility, only: create_log_grid, ID_to_string_mapping
   implicit none

   integer, parameter :: max_particle_id = 33
   integer, parameter :: npoints_ener = 1025 ! odd: 1024 intervals for composite Simpson integration
   real(real64), parameter :: default_emax = 1.0e14_real64
   real(real64), parameter :: annihilation_energy = 0.511_real64
   integer, parameter :: atmosphere_standard = 0
   integer, parameter :: atmosphere_msis = 1

   integer :: ip, iyear, imonth, iday
   integer :: solar_status, depth_status, atmosphere_model
   integer :: nargs, ii
   real(real64) :: emin, emax, alti, glat, glong, geometry
   real(real64) :: solar_w, rigidity, depth, continuum_flux, line_511_flux, total_flux
   real(real64), dimension(npoints_ener) :: energy_grid, flux_cm2_s_MeV
   logical :: solar_override_given
   character(len=32) :: particle_name, atmosphere_name, solar_source

   interface
      function getHP(year, month, day, status_code) result(value)
         import real64
         integer, intent(in) :: year, month, day
         integer, intent(out) :: status_code
         real(real64) :: value
      end function getHP

      function getr(latitude, longitude) result(value)
         import real64
         real(real64), intent(in) :: latitude, longitude
         real(real64) :: value
      end function getr

      function getd_model(altitude, latitude, model, status_code) result(value)
         import real64
         real(real64), intent(in) :: altitude, latitude
         integer, intent(in) :: model
         integer, intent(out) :: status_code
         real(real64) :: value
      end function getd_model

      function getSpec(particle_id, solar_index, cutoff_rigidity, atmospheric_depth, energy, g) result(value)
         import real64
         integer, intent(in) :: particle_id
         real(real64), intent(in) :: solar_index, cutoff_rigidity, atmospheric_depth, energy, g
         real(real64) :: value
      end function getSpec

      function get511flux(solar_index, cutoff_rigidity, atmospheric_depth) result(value)
         import real64
         real(real64), intent(in) :: solar_index, cutoff_rigidity, atmospheric_depth
         real(real64) :: value
      end function get511flux
   end interface

   ! The most recent valid daily value in the bundled FFPtable.day is 2019-05-27.
   iyear = 2019
   imonth = 5
   iday = 27
   geometry = 0.15_real64
   atmosphere_model = atmosphere_standard
   solar_override_given = .false.
   solar_w = 0.0_real64
   emax = default_emax

   nargs = command_argument_count()
   if (nargs < 5 .or. nargs > 11) then
      call print_usage()
      call fail('Expected between 5 and 11 command-line arguments.')
   end if

   call read_integer_argument(1, 'particle ID', ip)
   call read_real_argument(2, 'minimum energy', emin)
   call read_real_argument(3, 'altitude', alti)
   call read_real_argument(4, 'latitude', glat)
   call read_real_argument(5, 'longitude', glong)

   if (nargs >= 6) call read_integer_argument(6, 'year', iyear)
   if (nargs >= 7) call read_integer_argument(7, 'month', imonth)
   if (nargs >= 8) call read_integer_argument(8, 'day', iday)
   if (nargs >= 9) call read_real_argument(9, 'geometry parameter', geometry)
   if (nargs >= 10) call read_atmosphere_argument(10, atmosphere_model)
   if (nargs >= 11) then
      call read_real_argument(11, 'solar W-index override', solar_w)
      solar_override_given = .true.
   end if

   call validate_inputs(ip, emin, emax, alti, glat, glong, iyear, imonth, iday, &
      geometry, atmosphere_model, solar_override_given, solar_w)

   call ID_to_string_mapping(ip, particle_name)
   if (atmosphere_model == atmosphere_standard) then
      atmosphere_name = 'US Standard 1976'
   else
      atmosphere_name = 'NRLMSISE-00 table'
   end if

   if (solar_override_given) then
      solar_status = 0
      solar_source = 'user override'
   else
      solar_w = getHP(iyear, imonth, iday, solar_status)
      select case (solar_status)
       case (1)
         solar_source = 'daily neutron monitor'
       case (2)
         solar_source = 'annual Usoskin value'
       case (3)
         solar_source = 'daily GLE-adjusted value'
       case (4)
         call fail('No bundled solar-activity data are available for the requested date. ' // &
            'Pass a supported date or provide argument 11 as a solar W-index override.')
       case (5)
         call fail('The requested calendar date is not valid.')
       case default
         call fail('The solar-activity tables could not be read.')
      end select
   end if

   rigidity = getr(glat, glong)
   if (.not. ieee_is_finite(rigidity) .or. rigidity < 0.0_real64) then
      call fail('Cut-off rigidity calculation failed for the requested location.')
   end if

   depth = getd_model(alti, glat, atmosphere_model, depth_status)
   if (depth_status /= 0 .or. .not. ieee_is_finite(depth) .or. depth <= 0.0_real64) then
      call fail('Atmospheric-depth calculation failed for the requested altitude/model.')
   end if

   call create_log_grid(emin, emax, energy_grid, npoints_ener)
   do ii = 1, npoints_ener
      flux_cm2_s_MeV(ii) = getSpec(ip, solar_w, rigidity, depth, energy_grid(ii), geometry)
      if (.not. ieee_is_finite(flux_cm2_s_MeV(ii))) then
         call fail('The spectral model returned a non-finite value.')
      end if
      if (flux_cm2_s_MeV(ii) < 0.0_real64) then
         call fail('The spectral model returned a negative differential flux.')
      end if
   end do

   continuum_flux = log_simpson_integration(energy_grid, flux_cm2_s_MeV, npoints_ener)
   line_511_flux = 0.0_real64
   if (ip == 33 .and. emin <= annihilation_energy .and. annihilation_energy <= emax) then
      line_511_flux = get511flux(solar_w, rigidity, depth)
      if (.not. ieee_is_finite(line_511_flux) .or. line_511_flux < 0.0_real64) then
         call fail('The 511 keV photon-line model returned an invalid flux.')
      end if
   end if

   total_flux = continuum_flux + line_511_flux
   if (.not. ieee_is_finite(total_flux) .or. total_flux < 0.0_real64) then
      call fail('The integrated flux is invalid.')
   end if

   write(output_unit, '(A,I0,2A)') 'Particle ID: ', ip, ' (', trim(particle_name) // ')'
   if (ip >= 1 .and. ip <= 28) then
      write(output_unit, '(A,ES16.8,A)') 'Minimum energy: ', emin, ' MeV/nucleon'
   else
      write(output_unit, '(A,ES16.8,A)') 'Minimum energy: ', emin, ' MeV'
   end if
   if (ip >= 1 .and. ip <= 28) then
      write(output_unit, '(A,ES16.8,A)') 'Maximum energy: ', emax, ' MeV/nucleon'
   else
      write(output_unit, '(A,ES16.8,A)') 'Maximum energy: ', emax, ' MeV'
   end if
   write(output_unit, '(A,ES16.8,A)') 'Altitude: ', alti, ' km'
   write(output_unit, '(A,ES16.8,A)') 'Latitude: ', glat, ' degrees'
   write(output_unit, '(A,ES16.8,A)') 'Longitude: ', glong, ' degrees'
   write(output_unit, '(A,I4.4,A,I2.2,A,I2.2)') 'Solar date: ', iyear, '-', imonth, '-', iday
   write(output_unit, '(A,ES16.8,3A,I0,A)') 'Solar W index: ', solar_w, ' (', trim(solar_source), &
      ', status ', solar_status, ')'
   write(output_unit, '(A,ES16.8)') 'Vertical cut-off rigidity (GV): ', rigidity
   write(output_unit, '(A,ES16.8)') 'Atmospheric depth (g/cm2): ', depth
   write(output_unit, '(A,A)') 'Atmosphere model: ', trim(atmosphere_name)
   write(output_unit, '(A,ES16.8)') 'Geometry parameter: ', geometry
   write(output_unit, '(A,ES24.16,A)') 'Continuum energy-integrated flux: ', continuum_flux, ' cm^-2 s^-1'
   if (ip == 33) then
      write(output_unit, '(A,ES24.16,A)') '511 keV line contribution: ', line_511_flux, ' cm^-2 s^-1'
   end if
   write(output_unit, '(A,ES24.16,A)') 'Energy integrated flux: ', total_flux, ' cm^-2 s^-1'

   ! Stable, machine-readable record consumed by run_on_grid.py.
   write(output_unit, '(A,I0,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3,' // &
      'A,I0,A,I0,A,I0,A,ES24.16E3,A,I0,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3,A,I0,' // &
      'A,ES24.16E3,A,ES24.16E3,A,ES24.16E3)') &
      'RESULT_CSV,', ip, ',', emin, ',', emax, ',', alti, ',', glat, ',', glong, ',', &
      iyear, ',', imonth, ',', iday, ',', solar_w, ',', solar_status, ',', rigidity, ',', &
      depth, ',', geometry, ',', atmosphere_model, ',', continuum_flux, ',', line_511_flux, ',', total_flux

contains

   subroutine print_usage()
      write(error_unit, '(A)') 'Usage:'
      write(error_unit, '(A)') '  ./electron_fluxes PARTICLE_ID E_MIN ALT_KM LAT_DEG LON_DEG' // &
         ' [YEAR [MONTH [DAY [GEOMETRY [ATMOSPHERE [SOLAR_W]]]]]]'
      write(error_unit, '(A)') '  ATMOSPHERE: standard|us76|0 or msis|nrlmsise|1'
      write(error_unit, '(A)') '  Defaults: date=2019-05-27, geometry=0.15, atmosphere=standard.'
      write(error_unit, '(A)') '  SOLAR_W overrides the date lookup when supplied.'
   end subroutine print_usage

   subroutine fail(message)
      character(len=*), intent(in) :: message
      write(error_unit, '(A)') 'ERROR: ' // trim(message)
      flush(error_unit)
      error stop 1
   end subroutine fail

   subroutine read_integer_argument(index, description, value)
      integer, intent(in) :: index
      character(len=*), intent(in) :: description
      integer, intent(out) :: value
      character(len=256) :: argument
      integer :: status, ios

      call get_command_argument(index, value=argument, status=status)
      if (status /= 0 .or. len_trim(argument) == 0) then
         call fail('Missing ' // trim(description) // ' argument.')
      end if
      read(argument, *, iostat=ios) value
      if (ios /= 0) call fail('Invalid integer for ' // trim(description) // ': ' // trim(argument))
   end subroutine read_integer_argument

   subroutine read_real_argument(index, description, value)
      integer, intent(in) :: index
      character(len=*), intent(in) :: description
      real(real64), intent(out) :: value
      character(len=256) :: argument
      integer :: status, ios

      call get_command_argument(index, value=argument, status=status)
      if (status /= 0 .or. len_trim(argument) == 0) then
         call fail('Missing ' // trim(description) // ' argument.')
      end if
      read(argument, *, iostat=ios) value
      if (ios /= 0 .or. .not. ieee_is_finite(value)) then
         call fail('Invalid finite number for ' // trim(description) // ': ' // trim(argument))
      end if
   end subroutine read_real_argument

   subroutine read_atmosphere_argument(index, model)
      integer, intent(in) :: index
      integer, intent(out) :: model
      character(len=256) :: argument
      integer :: status

      call get_command_argument(index, value=argument, status=status)
      if (status /= 0 .or. len_trim(argument) == 0) call fail('Missing atmosphere argument.')
      call lowercase(argument)
      select case (trim(adjustl(argument)))
       case ('0', 'standard', 'us', 'us76', 'us-standard')
         model = atmosphere_standard
       case ('1', 'msis', 'nrlmsise', 'nrlmsise-00')
         model = atmosphere_msis
       case default
         call fail('Unknown atmosphere model: ' // trim(argument))
      end select
   end subroutine read_atmosphere_argument

   subroutine lowercase(text)
      character(len=*), intent(inout) :: text
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) text(i:i) = achar(code + 32)
      end do
   end subroutine lowercase

   subroutine validate_inputs(particle_id, min_energy, max_energy, altitude, latitude, longitude, &
      year, month, day, g, model, has_solar_override, solar_override)
      integer, intent(in) :: particle_id, year, month, day, model
      real(real64), intent(in) :: min_energy, max_energy, altitude, latitude, longitude, g, solar_override
      logical, intent(in) :: has_solar_override

      if (particle_id < 0 .or. particle_id > max_particle_id) then
         call fail('Particle ID must be between 0 and 33.')
      end if
      if (.not. ieee_is_finite(min_energy) .or. min_energy <= 0.0_real64) then
         call fail('Minimum energy must be finite and greater than zero.')
      end if
      if (min_energy >= max_energy) then
         call fail('Minimum energy must be smaller than the fixed maximum energy.')
      end if
      if (.not. ieee_is_finite(latitude) .or. latitude < -90.0_real64 .or. latitude > 90.0_real64) then
         call fail('Latitude must be between -90 and 90 degrees.')
      end if
      if (.not. ieee_is_finite(longitude) .or. longitude < -180.0_real64 .or. longitude > 180.0_real64) then
         call fail('Longitude must be between -180 and 180 degrees.')
      end if
      if (.not. valid_date(year, month, day)) call fail('The requested calendar date is invalid.')
      if (.not. valid_geometry(g)) then
         call fail('Invalid geometry parameter. Use 0..1 for ground, a negative aircraft value ' // &
            '(except -10), or a value >=10 for ideal/no-Earth modes.')
      end if
      if (model == atmosphere_standard) then
         if (.not. ieee_is_finite(altitude) .or. altitude < -0.5_real64 .or. altitude > 86.0_real64) then
            call fail('US Standard Atmosphere altitude must be within -0.5 to 86 km.')
         end if
      else if (model == atmosphere_msis) then
         if (.not. ieee_is_finite(altitude) .or. altitude < -0.6_real64 .or. altitude > 99.0_real64) then
            call fail('NRLMSISE table altitude must be within -0.6 to 99 km.')
         end if
      else
         call fail('Internal error: invalid atmosphere model.')
      end if
      if (has_solar_override) then
         if (.not. ieee_is_finite(solar_override) .or. solar_override < 0.0_real64) then
            call fail('Solar W-index override must be finite and nonnegative.')
         end if
      end if
   end subroutine validate_inputs

   logical function valid_geometry(g)
      real(real64), intent(in) :: g
      real(real64), parameter :: tolerance = 1.0e-12_real64

      valid_geometry = .false.
      if (.not. ieee_is_finite(g)) return
      if (g >= 0.0_real64 .and. g <= 1.0_real64) valid_geometry = .true.
      if (g < 0.0_real64 .and. abs(g + 10.0_real64) > tolerance) valid_geometry = .true.
      if (g >= 10.0_real64) valid_geometry = .true.
   end function valid_geometry

   logical function valid_date(year, month, day)
      integer, intent(in) :: year, month, day
      integer, dimension(12) :: days_in_month
      logical :: leap_year

      valid_date = .false.
      if (year < 1 .or. month < 1 .or. month > 12) return
      leap_year = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. mod(year, 400) == 0
      days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      if (leap_year) days_in_month(2) = 29
      if (day < 1 .or. day > days_in_month(month)) return
      valid_date = .true.
   end function valid_date

end program electron_fluxes
