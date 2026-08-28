program test_integration
   use integration, only: log_simpson_integration, simpson_integration
   use utility, only: create_log_grid
   implicit none

   integer, parameter :: n = 1025
   real(kind=8), dimension(n) :: energy, flux
   real(kind=8), dimension(3) :: uneven_x, uneven_y
   real(kind=8), dimension(4) :: even_x, even_y
   real(kind=8) :: value, expected

   call create_log_grid(1.0d0, 100.0d0, energy, n)

   flux = 1.0d0 / energy
   value = log_simpson_integration(energy, flux, n)
   expected = log(100.0d0)
   call assert_close('log-grid integral of 1/E', value, expected, 1.0d-12)

   flux = 1.0d0
   value = log_simpson_integration(energy, flux, n)
   expected = 99.0d0
   call assert_close('log-grid integral of one', value, expected, 1.0d-8)

   uneven_x = [0.0d0, 1.0d0, 3.0d0]
   uneven_y = uneven_x**2
   value = simpson_integration(uneven_x, uneven_y, size(uneven_x))
   expected = 9.0d0
   call assert_close('unequal-grid quadratic', value, expected, 1.0d-12)

   even_x = [0.0d0, 1.0d0, 2.0d0, 3.0d0]
   even_y = 1.0d0
   value = simpson_integration(even_x, even_y, size(even_x))
   expected = 3.0d0
   call assert_close('even-point final interval', value, expected, 1.0d-12)

   print '(A)', 'Integration tests passed.'

contains

   subroutine assert_close(label, actual, target, tolerance)
      character(len=*), intent(in) :: label
      real(kind=8), intent(in) :: actual, target, tolerance
      if (abs(actual - target) > tolerance) then
         write(*, '(A,2ES24.16)') trim(label) // ' failed; actual/expected: ', actual, target
         error stop 1
      end if
   end subroutine assert_close

end program test_integration
