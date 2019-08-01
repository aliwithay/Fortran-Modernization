module mic_mod
    ! for mic block
  implicit none
      integer :: ihydro,isolu,colcal,colgr
      integer :: sendfol
      real(8) :: ka1,ka2,diffvnd1,diffvnd2
      real(8) :: delt,ks
      integer :: nstop,micdtout,nout,dsdout
      real(8) :: pp,sp,rhoa,thetapp,qvpp,seq
      real(8) :: thetap,qvp,deltaqp
      real    :: cellw,slabw,oneoverh,h
      real(8) :: vol
      real    :: rmth,rmqv,radmin
      integer,parameter :: maxinfocol=18 !number of variables to be transered between processors.
      integer, allocatable, dimension(:) :: ixp,iyp,izp
      real, allocatable, dimension(:)    :: dtovertaup
      real, allocatable, dimension(:)    :: uudrop,vvdrop,wwdrop,flowu,flowv,floww,fv
      real, allocatable, dimension(:)    :: yrel !ndroppe
      integer,allocatable, dimension(:)  :: idpfol
      real, allocatable, dimension(:) :: xfols,yfols,zfols,rfols, ssfols
      integer, allocatable, dimension(:) :: idpfolr
      real, allocatable, dimension(:)    :: xfol,yfol,zfol,rfol,ssfol !ndrop/everynd
      real, allocatable, dimension(:,:)  :: infocols, infocolr!(ndroppe,maxinfocol)
      real, allocatable, dimension(:) :: xfolr,yfolr,zfolr,rfolr,ssfolr 
      integer, allocatable, dimension(:) :: head !ncell+m*m
      real, allocatable, dimension(:)    :: maxutemp !npe
      integer, allocatable, dimension(:,:) :: dcids !(5,2)
      integer, allocatable, dimension(:) :: dsd,dsdtot !nbins dsdtot only in mype0
      integer, allocatable, dimension(:) :: dsd_log2,ccn_log2,dsdtot_log2,ccntot_log2 !nbins dsdtot only in mype0
      integer, allocatable, dimension(:) :: idp
      integer, allocatable, dimension(:) :: idpfols
      real, allocatable, dimension(:)    :: x,y,z,r,r_ccn,kappa,udrop,vdrop,wdrop,dr3,supersatd
      contains
      real function log2(rad)
      real :: rad
        log2 = log(rad)/log(2.d0)
      end function log2

end module mic_mod


