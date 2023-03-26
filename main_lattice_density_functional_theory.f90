  program lattice_density_functional_theory
!$ USE OMP_LIB
  implicit none

  integer                      :: x             ,  y             ,  z
  integer                      :: x_size        ,  y_size        ,  z_size
  integer                      :: x_min_internal,  y_min_internal,  z_min_internal
  integer                      :: x_max_internal,  y_max_internal,  z_max_internal
  integer                      :: bc_x          ,  bc_y          ,  bc_z
  integer                      :: n_scf_loop_max
  integer                      :: ierr

  integer                      :: num_pore_internal

  integer                      :: i

  real (kind = 8)              :: alpha
  real (kind = 8)              :: temperature
  real (kind = 8)              :: eps_converge
  real (kind = 8)              :: activity

  real (kind = 8)              :: grand_potential, entropy, energy_ff, energy_sf, chemipo

  logical                      :: logi_input_density
  logical                      :: logi_output_density
  logical                      :: logi_converge

  integer        , allocatable :: structure       (:,:,:)
  real (kind = 8), allocatable :: structure_real  (:,:,:)
  real (kind = 8), allocatable :: density         (:,:,:)
  real (kind = 8), allocatable :: density_prev    (:,:,:)

  character(len = 5)           :: filename

!!$  write(*,*) dlog   (2.718281828459d0)
!!$  write(*,*) dlog10 (2.718281828459d0)
!!$
!!$  write(*,*) dlog   (10.0d0)
!!$  write(*,*) dlog10 (10.0d0)

  open ( 9, file='./param.txt')
  open (10, file='./input_structure.txt')
  open (11, file='./input_density.txt')

  open (31, file='./output_structure_final.txt')
  open (41, file='./output_density_final.txt')
  open (51, file='./output_ave_dens_internal.txt')
  open (61, file='./output_omega.txt')

#ifdef debug
  open (410, file='./debug_structure_all.txt')
  open (420, file='./debug_output_omega_all.txt')
  write(420,'(a71)') '# i, activity,  entropy, energy_ff, energy_sf, chemipo, grand_potential'
#endif

  write( 61,'(a71)') '# i, activity,  entropy, energy_ff, energy_sf, chemipo, grand_potential'


  read ( 9, * ) x_size,  y_size,  z_size
  read ( 9, * ) x_min_internal,  y_min_internal,  z_min_internal
  read ( 9, * ) x_max_internal,  y_max_internal,  z_max_internal
  read ( 9, * ) alpha
  read ( 9, * ) temperature
  read ( 9, * ) eps_converge
  read ( 9, * ) bc_x  ,  bc_y  ,  bc_z
  read ( 9, * ) logi_input_density
  read ( 9, * ) n_scf_loop_max
  read ( 9, * )

  write( *, '(a39,3i10)'  ) ' # x_size ,  y_size ,  z_size        = ', x_size,  y_size,  z_size
  write( *, '(a39,3i10)'  ) ' # x_min,  y_min,  z_min (internal)  = ', x_min_internal,  y_min_internal,  z_min_internal
  write( *, '(a39,3i10)'  ) ' # x_max,  y_max,  z_max (internal)  = ', x_max_internal,  y_max_internal,  z_max_internal
  write( *, '(a39,f20.10)') ' # alpha       [eps_ff]              = ', alpha
  write( *, '(a39,f20.10)') ' # temperature [1/k_B]               = ', temperature
  write( *, '(a39,e20.10)') ' # eps_converge                      = ', eps_converge
  write( *, '(a85)'       ) ' ## if ave ( (dens - dens_prev ) ** 2 ) < eps_converge ** 2, exit scf calculation. ##'
  write( *, '(a39,3i10)'  ) ' # bc_x  ,  bc_y  ,  bc_z            = ', bc_x  ,  bc_y  ,  bc_z
  write( *, '(a37)'       ) ' ## bc = 1 : mirror,  bc = 2 : pbc ##'
  write( *, '(a39,l5)'    ) ' # logi_input_density                = ', logi_input_density
  write( *, *             )

  allocate ( structure       ( 0:x_size + 1, 0:y_size + 1, 0:z_size + 1) )
  allocate ( structure_real  ( 0:x_size + 1, 0:y_size + 1, 0:z_size + 1) )
  allocate ( density         ( 0:x_size + 1, 0:y_size + 1, 0:z_size + 1) )
  allocate ( density_prev    ( 0:x_size + 1, 0:y_size + 1, 0:z_size + 1) )

  do z = 0, z_size + 1  ; do y = 0, y_size + 1 ; do x = 0 , x_size + 1
      structure  (x, y, z) = 0
      density    (x, y, z) = 0.0d0
  end do; end do; end do

  do z = 1, z_size;  do y = 1, y_size
      read(10, *) ( structure (x, y, z), x = 1, x_size)
  end do;  end do


  num_pore_internal =   ( x_max_internal - x_min_internal + 1 ) &
  &                   * ( y_max_internal - y_min_internal + 1 ) &
  &                   * ( z_max_internal - z_min_internal + 1 )

  do z = z_min_internal, z_max_internal
      do y = y_min_internal, y_max_internal
          do x = x_min_internal, x_max_internal
              num_pore_internal = num_pore_internal - structure (x, y, z)
          end do
      end do
  end do

  write( *, * ) '# num_pore_internal = ', num_pore_internal
  write( *, * ) '# Reading input_density.txt ... '

  if ( logi_input_density ) then
      do z = 1, z_size
          do y = 1, y_size
              read(11, *) ( density   (x, y, z), x = 1, x_size)
          end do
!         write(*,*) z, '/', z_size
      end do
  end if
  write(*,*) '# Done.'


  structure_real (:, :, :) = dble ( structure (:, :, :) )
  call boundary_condition (  structure_real,  x_size,  y_size,  z_size,  bc_x,  bc_y,  bc_z )
  structure (:, :, :) = nint ( structure_real (:, :, :) )


  do while (.true.)

      read(9, *,iostat=ierr) activity, logi_output_density
      if ( ierr .ne. 0 ) exit
      write(*,*) ierr, activity, logi_output_density

      i = 0
      do while (.true.)

          i = i + 1
          if ( i .ge. n_scf_loop_max ) then
              write(*,*) '# i = ', i
              write(*,*) '# L-DFT not converged. '
              stop
          end if

          call boundary_condition (  density,  x_size,  y_size,  z_size,  bc_x,  bc_y,  bc_z )
!$OMP WORKSHARE
          density_prev (:, :, :) = density (:, :, :)
!$OMP END WORKSHARE

          call update_density ( structure, density_prev, density,  x_size,  y_size,  z_size, temperature, eps_converge, alpha, activity, n_scf_loop_max)
          write(*,*) i

#ifdef debug
          call calc_omega                                    &
  & (          x_size,              y_size,          z_size, &
  &    x_min_internal,      y_min_internal,  z_min_internal, &
  &    x_max_internal,      y_max_internal,  z_max_internal, &
  &         structure,             density,     temperature, &
  &             alpha,            activity,         entropy, &
  &         energy_ff,           energy_sf,         chemipo )

          grand_potential = temperature * entropy + energy_ff +  energy_sf +  chemipo
          write(420,'(i10, 6e25.15)') i, activity, entropy, energy_ff, energy_sf, chemipo, grand_potential
#endif

          call check_convergence (density_prev, density, x_size,  y_size,  z_size, eps_converge, logi_converge)
          if ( logi_converge ) exit

      end do
      write(*,*) '@', i, 'scf loop converged.'


      call calc_ave_density                                  &
  & (          x_size,              y_size,          z_size, &
  &    x_min_internal,      y_min_internal,  z_min_internal, &
  &    x_max_internal,      y_max_internal,  z_max_internal, &
  &           density,   num_pore_internal,        activity, &
  &                                                      51 )

      call calc_omega                                        &
  & (          x_size,              y_size,          z_size, &
  &    x_min_internal,      y_min_internal,  z_min_internal, &
  &    x_max_internal,      y_max_internal,  z_max_internal, &
  &         structure,             density,     temperature, &
  &             alpha,            activity,         entropy, &
  &         energy_ff,           energy_sf,         chemipo )
!     call calc_omega  ( structure,  density, x_size,  y_size,  z_size, temperature, alpha, activity, entropy, energy_ff, energy_sf, chemipo)
      grand_potential = temperature * entropy + energy_ff +  energy_sf +  chemipo
      write(61,'(i10, 6e25.15)') i, activity,  entropy, energy_ff, energy_sf, chemipo, grand_potential


      if ( logi_output_density ) then

          write(filename,'(i5.5)') nint (activity * 100000.d0)
          open(100,file = './density'//filename//'.txt')

          do z = 1, z_size;  do y = 1, y_size
              do x = 1, x_size
                  write(100, '(f12.7,$)') density (x, y, z)
              end do
              write(100, '(a)') ''
          end do;  end do

          close (100)

      end if

  end do

  write( *, * ) '# L-DFT Done.'
  write( *, * ) '# Writing output files ... '


  do z = 1, z_size;  do y = 1, y_size
      do x = 1, x_size
          write(31, '(i2,$)'   ) structure (x, y, z)
          write(41, '(f12.7,$)') density   (x, y, z)
      end do
      write(31, '(a)') ''
      write(41, '(a)') ''
  end do;  end do


#ifdef debug
  do z = 0, z_size + 1
      do y = 0, y_size + 1
          do x = 0, x_size + 1
              write(410, '(i2,$)') structure (x, y, z)
          end do
          write(410, '(a)') ''
      end do
  end do
#endif

  deallocate ( structure )
  deallocate ( density   )

  write( *, * ) '# Calculation finished successfully. '

  stop
  end program lattice_density_functional_theory
