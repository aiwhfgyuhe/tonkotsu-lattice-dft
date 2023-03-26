  subroutine calc_ave_density                                &
  & (          x_size,              y_size,          z_size, &
  &    x_min_internal,      y_min_internal,  z_min_internal, &
  &    x_max_internal,      y_max_internal,  z_max_internal, &
  &           density,   num_pore_internal,        activity, &
  &                                                file_num )
!$ USE OMP_LIB
  implicit none

  integer                       :: x             ,  y             ,  z
  integer         , intent(in)  :: x_size        ,  y_size        ,  z_size
  integer         , intent(in)  :: x_min_internal,  y_min_internal,  z_min_internal
  integer         , intent(in)  :: x_max_internal,  y_max_internal,  z_max_internal
  integer         , intent(in)  :: num_pore_internal

  real (kind = 8 ), intent(in)  :: density  ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )
  real (kind = 8 ), intent(in)  :: activity
  integer         , intent(in)  :: file_num
  real (kind = 8 )              :: sum

  sum = 0.0d0

!$OMP PARALLEL DO PRIVATE(x, y) REDUCTION(+:sum)
  do z = z_min_internal, z_max_internal
      do y = y_min_internal, y_max_internal
          do x = x_min_internal, x_max_internal
              sum = sum + density(x, y, z)
          end do
      end do
  end do
!$OMP END PARALLEL DO

  write( file_num, * )  activity, sum, sum / dble (num_pore_internal)

  return
  end subroutine calc_ave_density
