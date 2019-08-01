!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c FFTS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine fftwrk(zr,zk)

  use scratch_mod
  use mpi_mod
  use mpi
  use fftw_mod
  use iso_c_binding
  !use ftw2d_mod
  !use ftw1d_mod
  implicit none
  include 'param.inc'
  include 'fftw3.f03'

  complex, dimension(iktx,ikty,iktzp) :: zk
  real,    dimension(n1d,n3d,n2dp)    :: zr
  real :: norm
  integer :: h,i,j,k,l,ikstart
  integer, parameter :: n1n3=n1d*n3d
  integer, parameter :: iktxz=iktx*iktz,iktxy=iktx*ikty,iktyp=ikty/npe
  integer, parameter :: iktxyzp=iktx*ikty*iktzp

!  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
!  common/FTW2D/plan2_rk,plan2_kr
!  common/FTW1D/plan1_rk,plan1_kr

  integer :: ierror
!  integer :: mype
!  common/mpi/mype

!  if( .not. scratch_alloue ) then
!        allocate( zkt(iktx,iktz,iktyp), zk1(iktx,ikty,iktzp), zk2(iktx,ikty,iktzp) )
!  endif

  norm = float(n1*n2*n3) !#grids

  ! Do 2D (x,z) transforms at n2 levels all at once 
  ! Note, output data has dimensions iktx,iktz,iktyp 
  ! zk1 is scratch
  ! real to k-space

  do h = 1, n2pe, n1n3
   ! call dfftw_plan_dft_r2c(plan, m, n, in array, out array, flags)
     call dfftw_plan_dft_r2c_2d(plan2_rk, n1, n3, zr, zk1, fftw_measure)
     call dfftw_execute(plan2_rk)
  enddo
  ! call rfftwnd_f77_real_to_complex(plan, howmany, in array, istride, idist, out array, ostride, odist)
  ! call rfftwnd_f77_real_to_complex(plan2_rk,n2pe,zr,1,n1n3,zk1,1,iktxz)

  ! 2D transforms are in-place so output is in zk
  ! but has dimensions iktx,iktz,iktyp
  ! copy to zkt in prep for transpose
!!  zkt=zk ! this doesn't work on sharcnet!  Do this instead:  (sharcnet: shared hierarchical academic research computing network)
   do i=1,iktxyzp
     zkt(i,1,1)=zk(i,1,1)
   enddo

  ! Transpose zkt(iktx,iktz,iktyp) -> zk(iktx,ikty,iktzp)
  ! Transpose output from previous step to have dimensions iktx,ikty,iktzp
!if (MPI == 1 ) then
  call mpitranspose(zkt,iktx,iktz,iktyp,zk,ikty,iktzp,npe,zk1,zk2) ! k-space to real
!else
!  call serialtranspose(zkt,zk,iktx,iktz,ikty)
!endif

  ! Do remaining 1D (y) transforms at iktx*iktz rows
  do l=1,iktx, 1
     do k=1,iktzp
       ikstart=1+(k-1)*iktxy
     ! call dfftw_plan_dft_1d(plan, size, in array, out array, direction, flags)
       call dfftw_plan_dft_1d(plan1_rk, iktx, zk(ikstart,1,1),zk1,fftw_forward,fftw_measure)
       call dfftw_execute(plan1_rk)
       ! call fftw_f77(plan, howmany, in array, )
       ! call fftw_f77(plan1_rk,iktx,zk(ikstart,1,1),iktx,1,zk1,iktx,1)
     enddo
  enddo
  
  ! Normalize
  zk=zk/norm
end subroutine fftwrk


subroutine fftwkr(zr,zk)

  use scratch_mod
  use mpi_mod
  use mpi
  !use ftw1d_mod
  !use ftw2d_mod
  use fftw_mod
  use iso_c_binding
  implicit none
  include 'param.inc'
  include 'fftw3.f03'
!  include 'mpif.h'

  complex, dimension(iktx,ikty,iktzp) :: zk
  real,    dimension(n1d,n3d,n2dp)    :: zr
  integer :: h,i,j,k,ikstart
  integer, parameter :: n1n3=n1d*n3d
  integer, parameter :: iktxz=iktx*iktz,iktxy=iktx*ikty,iktyp=ikty/npe
  integer, parameter :: iktxyzp=iktx*ikty*iktzp

!  integer*8 :: plan2_rk,plan2_kr,plan1_rk,plan1_kr
!  common/FTW2D/plan2_rk,plan2_kr
!  common/FTW1D/plan1_rk,plan1_kr

  integer :: ierror
!  integer :: mype, ierror
!  common/mpi/mype


!  call mpi_barrier(MPI_COMM_WORLD,ierror)
!  if (mype .eq. 0) print*, 'In fftwkr before allocate'
!  call mpi_barrier(MPI_COMM_WORLD,ierror)


!  if( .not. scratch_alloue ) then
!        allocate( zkt(iktx,iktz,iktyp), zk1(iktx,ikty,iktzp), zk2(iktx,ikty,iktzp) )
!  endif

!  call mpi_barrier(MPI_COMM_WORLD,ierror)
!  if (mype .eq. 0) print*, 'In fftwkr after allocate' 
!  call mpi_barrier(MPI_COMM_WORLD,ierror)


  call realit(zk)

  ! Do 1D (y) transforms at iktx*iktz rows
  do j=1,iktx,1
     do k=1,iktzp
       ikstart=1+(k-1)*iktxy
       ! call dfftw_plan_dft_1d(plan, size, in array, out array, direction, flag)
       call dfftw_plan_dft_1d(plan1_kr, iktx, zk(ikstart,1,1), zk1, fftw_backward, fftw_measure)
       call dfftw_execute(plan1_kr)
       ! call fftw_f77(plan1_kr,iktx,zk(ikstart,1,1),iktx,1,zk1,iktx,1)
     enddo
  enddo

  ! Transpose zk(iktx,ikty,iktzp) -> zkt(iktx,iktz,iktyp)
  ! Note, output of transpose has dimensions iktx,iktz,iktyp but store it in zk
!if (MPI == 1) then
  call mpitranspose(zk,iktx,ikty,iktzp,zkt,iktz,iktyp,npe,zk1,zk2)
!else
!  call serialtranspose(zk,zkt,iktx,ikty,iktz)
!endif

  ! Copy transposed array back into zk in preparation for in-place 2D transforms
!! zk=zkt ! this doesnt work on sharcnet
  do i=1,iktxyzp
     zk(i,1,1)=zkt(i,1,1)
  enddo

  ! Finally do 2D (x,z) transforms at n2 levels all at once
  ! Note, input data in zk has dimensions iktx,iktz,iktyp
  ! zk1 is scrtatch
  ! call dfftw_plan_dft_c2r(plan, rank, size array, in array, out array, flags)

  !do h=1,n2pe,1
     call dfftw_plan_dft_c2r_2d(plan2_kr,n1,n3,zk,zk1,fftw_measure)
     call dfftw_execute(plan2_kr)
     ! call rfftwnd_f77_complex_to_real(plan2_kr,n2pe,zk,1,iktxz,zk1,1,n1n3)
  !enddo

end subroutine fftwkr
