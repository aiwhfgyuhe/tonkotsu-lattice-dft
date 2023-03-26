  subroutine calc_omega                                      &
  & (          x_size,              y_size,          z_size, &
  &    x_min_internal,      y_min_internal,  z_min_internal, &
  &    x_max_internal,      y_max_internal,  z_max_internal, &
  &         structure,             density,     temperature, &
  &             alpha,            activity,         entropy, &
  &         energy_ff,           energy_sf,         chemipo )
!$ USE OMP_LIB
  implicit none

  integer                          :: x     ,  y     ,  z
  integer         , intent(in)     :: x_size,  y_size,  z_size
  integer         , intent(in)     :: x_min_internal,  y_min_internal,  z_min_internal
  integer         , intent(in)     :: x_max_internal,  y_max_internal,  z_max_internal

  real (kind = 8 ), intent(in)     :: temperature,  alpha, activity

  integer         , intent(in)     :: structure    ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )
  real (kind = 8 ), intent(in)     :: density      ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )

  real (kind = 8 ), intent(out)    :: entropy, energy_ff, energy_sf, chemipo

  entropy = 0.0d0
  energy_ff = 0.0d0
  energy_sf = 0.0d0
  chemipo = 0.0d0

!$OMP PARALLEL DO PRIVATE(x, y) REDUCTION(+:entropy, energy_ff, energy_sf)
  do z = z_min_internal, z_max_internal
      do y = y_min_internal,  y_max_internal;  do x = x_min_internal,  x_max_internal

          if ( ( density (x, y, z) .gt. 0.0d0 ) .and. ( density (x, y, z) .lt. 1.0d0 ) ) then
              entropy = entropy +          density (x, y, z)   * dlog (         density (x, y, z) ) &
  &                             + (1.0d0 - density (x, y, z) ) * dlog ( 1.0d0 - density (x, y, z) )
          end if

          if( structure (x, y, z) .eq. 0) then

              energy_ff = energy_ff - 0.5d0 * density (x, y, z  ) &
  &                      * (    density (x - 1, y    , z    )     &
  &                           + density (x + 1, y    , z    )     &
  &                           + density (x    , y - 1, z    )     &
  &                           + density (x    , y + 1, z    )     &
  &                           + density (x    , y    , z - 1)     &
  &                           + density (x    , y    , z + 1)     &
  &                        )

              energy_sf = energy_sf - alpha *  density (x, y, z  ) &
  &                       * (    structure (x - 1, y    , z    )   &
  &                            + structure (x + 1, y    , z    )   &
  &                            + structure (x    , y - 1, z    )   &
  &                            + structure (x    , y + 1, z    )   &
  &                            + structure (x    , y    , z - 1)   &
  &                            + structure (x    , y    , z + 1)   &
  &                         )

              chemipo = chemipo - density (x, y, z) * temperature * dlog (activity)

          end if

      end do;  end do
  end do
!$OMP END PARALLEL DO

  end subroutine calc_omega
