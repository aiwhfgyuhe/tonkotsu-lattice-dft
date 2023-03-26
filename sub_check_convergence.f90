  subroutine check_convergence (density_prev, density, x_size,  y_size,  z_size, eps_converge, logi_converge)
!$ USE OMP_LIB
  implicit none

  integer                          :: x     ,  y     ,  z
  integer         , intent(in)     :: x_size,  y_size,  z_size

  real (kind = 8 ), intent(in)     :: density_prev ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )
  real (kind = 8 ), intent(in)     :: density      ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )
  real (kind = 8 ), intent(in)     :: eps_converge
  logical         , intent(out)    :: logi_converge

  real (kind = 8 )                 :: diff, err

  err = 0.0d0

!$OMP PARALLEL DO PRIVATE(x, y, diff) REDUCTION(+:err)
  do z = 1, z_size
      do y = 1, y_size
          do x = 1, x_size
              diff = density(x, y, z) - density_prev(x, y, z)
              err = err + diff * diff
          end do
      end do
  end do
!$OMP END PARALLEL DO

  err = err / dble ( x_size * y_size * z_size )

  if ( err < eps_converge * eps_converge) then
     logi_converge = .true.
  else
     logi_converge = .false.
  end if

  end subroutine check_convergence
