  subroutine update_density ( structure, density_prev, density, x_size,  y_size,  z_size, temperature, ediff, alpha, activity, n_scf_loop_max)
!$ USE OMP_LIB
  implicit none

  integer                          :: x     ,  y     ,  z
  integer         , intent(in)     :: x_size,  y_size,  z_size
  real (kind = 8 ), intent(in)     :: temperature, ediff, alpha, activity
  integer         , intent(in)     :: n_scf_loop_max

  integer         , intent(in)     :: structure    ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )
  real (kind = 8 ), intent(in)     :: density_prev ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )
  real (kind = 8 ), intent(out)    :: density      ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )
  real (kind = 8 )                 :: c, s


!$OMP PARALLEL DO PRIVATE(x, y, s, c)
  do z = 1, z_size
      do y = 1,  y_size;  do x = 1, x_size
          if( structure (x, y, z) .eq. 0) then

              s =    density_prev(x - 1, y    , z    ) + alpha * structure (x - 1, y    , z    ) &
  &                + density_prev(x + 1, y    , z    ) + alpha * structure (x + 1, y    , z    ) &
  &                + density_prev(x    , y - 1, z    ) + alpha * structure (x    , y - 1, z    ) &
  &                + density_prev(x    , y + 1, z    ) + alpha * structure (x    , y + 1, z    ) &
  &                + density_prev(x    , y    , z - 1) + alpha * structure (x    , y    , z - 1) &
  &                + density_prev(x    , y    , z + 1) + alpha * structure (x    , y    , z + 1)

              c = activity * dexp(s / temperature)

              density(x, y, z) = c / (1.0d0 + c)

          end if

      end do;  end do
  end do
!$OMP END PARALLEL DO


  return
  end subroutine update_density

