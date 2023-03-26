  subroutine boundary_condition (  voxel_data,  x_size,  y_size,  z_size,  bc_x,  bc_y,  bc_z )
  implicit none

  integer                          :: x     ,  y     ,  z
  integer         , intent(in)     :: x_size,  y_size,  z_size
  integer         , intent(in)     :: bc_x  ,  bc_y  ,  bc_z
  real (kind = 8 ), intent(inout)  :: voxel_data ( 0:x_size + 1,  0:y_size + 1,  0:z_size + 1 )

!  bc = 1 : mirror,  bc = 2 : pbc
!
!**********  face  ****************************************

  do z = 1, z_size;  do y = 1, y_size
      if (bc_x .eq. 1 ) then
          voxel_data(         0,          y,          z) = voxel_data(     1,      y,      z)
          voxel_data(x_size + 1,          y,          z) = voxel_data(x_size,      y,      z)
      else if  (bc_x .eq. 2 ) then
          voxel_data(         0,          y,          z) = voxel_data(x_size,      y,      z)
          voxel_data(x_size + 1,          y,          z) = voxel_data(     1,      y,      z)
      else
          write(*,*) "error in bc_x"
          write(*,*) bc_x
          stop
      end if
  end do;  end do

  do z = 1, z_size;  do x = 1, x_size
      if (bc_y .eq. 1 ) then
          voxel_data(         x,          0,          z) = voxel_data(     x,      1,      z)
          voxel_data(         x, y_size + 1,          z) = voxel_data(     x, y_size,      z)
      else if  (bc_y .eq. 2 ) then
          voxel_data(         x,          0,          z) = voxel_data(     x, y_size,      z)
          voxel_data(         x, y_size + 1,          z) = voxel_data(     x,      1,      z)
      else
          write(*,*) "error in bc_y"
          write(*,*) bc_y
          stop
      end if
  end do;  end do

  do y = 1, y_size;  do x = 1, x_size
      if (bc_z .eq. 1 ) then
          voxel_data(         x,          y,          0) = voxel_data(     x,      y,      1)
          voxel_data(         x,          y, z_size + 1) = voxel_data(     x,      y, z_size)
      else if  (bc_z .eq. 2 ) then
          voxel_data(         x,          y,          0) = voxel_data(     x,      y, z_size)
          voxel_data(         x,          y, z_size + 1) = voxel_data(     x,      y,      1)
      else
          write(*,*) "error in bc_z"
          write(*,*) bc_z
          stop
      end if
  end do;  end do

!**********  edge  ****************************************

  if ( ( bc_x .eq. 2 ) .and. ( bc_y .eq. 2 ) ) then
      do z = 1, z_size
          voxel_data(         0,          0,          z) = voxel_data(x_size, y_size,      z)
          voxel_data(         0, y_size + 1,          z) = voxel_data(x_size,      1,      z)
          voxel_data(x_size + 1,          0,          z) = voxel_data(     1, y_size,      z)
          voxel_data(x_size + 1, y_size + 1,          z) = voxel_data(     1,      1,      z)
      end do
  end if


  if ( ( bc_y .eq. 2 ) .and. ( bc_z .eq. 2 ) ) then
      do x = 1, x_size
          voxel_data(         x,          0,          0) = voxel_data(     x, y_size, z_size)
          voxel_data(         x,          0, z_size + 1) = voxel_data(     x, y_size,      1)
          voxel_data(         x, y_size + 1,          0) = voxel_data(     x,      1, z_size)
          voxel_data(         x, y_size + 1, z_size + 1) = voxel_data(     x,      1,      1)
      end do
  end if


  if ( ( bc_z .eq. 2 ) .and. ( bc_x .eq. 2 ) ) then
      do y = 1, y_size
          voxel_data(         0,          y,          0) = voxel_data(x_size,      y, z_size)
          voxel_data(x_size + 1,          y,          0) = voxel_data(     1,      y, z_size)
          voxel_data(         0,          y, z_size + 1) = voxel_data(x_size,      y,      1)
          voxel_data(x_size + 1,          y, z_size + 1) = voxel_data(     1,      y,      1)
      end do
  end if

!**********  vertex  ****************************************

  if ( ( bc_x .eq. 2 ) .and. ( bc_y .eq. 2 ) .and. ( bc_z .eq. 2 ) ) then
      voxel_data(         0,          0,          0) = voxel_data(x_size, y_size, z_size)
      voxel_data(x_size + 1,          0,          0) = voxel_data(     1, y_size, z_size)
      voxel_data(         0, y_size + 1,          0) = voxel_data(x_size,      1, z_size)
      voxel_data(x_size + 1, y_size + 1,          0) = voxel_data(     1,      1, z_size)
      voxel_data(         0,          0, z_size + 1) = voxel_data(x_size, y_size,      1)
      voxel_data(x_size + 1,          0, z_size + 1) = voxel_data(     1, y_size,      1)
      voxel_data(         0, y_size + 1, z_size + 1) = voxel_data(x_size,      1,      1)
      voxel_data(x_size + 1, y_size + 1, z_size + 1) = voxel_data(     1,      1,      1)
  end if

  return
  end subroutine boundary_condition
