program electron_fluxes
!  Calculate cosmic-ray fluxes in the atmosphere based on PARMA model
   parameter(npart=33) ! number of applicable particle
   implicit real*8 (a-h, o-z)
   dimension IangPart(0:npart)
   integer, parameter :: npoints_ener = 1024
   integer, parameter :: npoints_alt = 5
   real(kind=8) :: emin, emax
   character(len=100) :: arg         ! Buffer for reading command line arguments
   integer :: iostat                 ! I/O status indicator for error handling
   real(kind=8), dimension(npoints_ener) :: energy_grid
   real(kind=8), dimension(npoints_ener) :: flux_cm2_s_MeV
   real(kind=8), dimension(npoints_alt) :: altitude_grid
   character(len=10) :: name

   data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID

   ip=31 ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
   iyear=2024  ! Year
   imonth=7    ! Month
   iday=29     ! Date
   glat=21.0   ! Latitude (deg), -90 =< glat =< 90
   glong=-81.0 ! Longitude (deg), -180 =< glong =< 180
   g=0.15      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
   ! g does not make a difference for electrons between 10 and 20 km altitude

   alti=16.0   ! Altitude (km)
   ang=1.0     ! cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction).

   ee=1e9    ! Energy (MeV/n)

! Set the minimum and maximum energy values
   emin = 0.3 ! MeV, can be changed by input argument
   emax = 1e14

! Check command-line arguments
   call get_command_argument(1, arg, iostat)
   if (len_trim(arg) > 0) then
      read(arg, *) ip
   else
      print *, 'Please enter a particle ID number (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)'
      read *, ip
   endif

   call get_command_argument(2, arg, iostat)
   if (len_trim(arg) > 0) then
      read(arg, *) emin
   else
      print *, 'Please enter an energy threshold in MeV:'
      read *, emin
   endif

   call get_command_argument(3, arg, iostat)
   if (len_trim(arg) > 0) then
      read(arg, *) alti
   else
      print *, 'Please enter an altitude in kilometers:'
      read *, alti
   endif

   call intToStringMapping(ip, name)

   ! Output the entered values to confirm
   print *, 'Particle ID: ', ip ,'(',trim(name),')'
   print '(A, F10.2, A)', ' Energy Threshold: ', emin, ' MeV'
   print '(A, F17.1, A)', ' Altitude: ', alti, '  kilometers'

   s=getHP(iyear,imonth,iday,idummy) ! Solar activity (W index) of the day
   r=getr(glat,glong)                ! Vertical cut-off rigidity (GV)
   d=getd(alti,glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976

! Call the subroutine to fill the grid
   call create_log_grid(emin, emax, energy_grid, npoints_ener)

   call create_lin_grid(10, 20, altitude_grid, npoints_alt)

! Output the results
   !call output_array(energy_grid, npoints_ener)

   Flux = getSpec(ip,s,r,d,ee,g)
   print '(A, ES8.2, A, ES10.2)', ' Angular Integrated Flux (/cm2/s/(MeV/n)) for ', ee, ' MeV =>', Flux
   !if(IangPart(ip).ne.0) then ! Angular distribution is available for the particle
   !   DifFlux=Flux*getSpecAngFinal(iangpart(ip),s,r,d,e,g,ang)
   !   write(6,*) 'Angular Differential Flux(/cm2/s/(MeV/n)/sr)=',DifFlux
   !endif

   do ii = 1, npoints_ener
      flux_cm2_s_MeV(ii) = getSpec(ip,s,r,d,energy_grid(ii),g)
   end do

   !call output_array(flux_cm2_s_MeV, npoints_ener)

   !Calculate the total flux using the trapezoidal rule
   !total_flux = trapezoidal_integration(energy_grid, flux_cm2_s_MeV, npoints_ener)
   total_flux = simpson_rule(energy_grid, flux_cm2_s_MeV, npoints_ener)

   ! Output the result
   print *, 'Energy integrated flux : ', total_flux, ' cm^-2 s^-1'


end program electron_fluxes


