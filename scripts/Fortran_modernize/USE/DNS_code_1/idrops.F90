  SUBROUTINE idrops(seed,rm)
! This subroutine determines the initial position and size of all droplets & its dry size
!! for droplets locations & ID# & random # generator
  use mic_mod
  implicit none
   include 'param.inc'
  include 'mpif.h'

  ! --- local --
  integer :: i,entr
  real :: ran1
  real,external :: rangen
  integer :: seed,id,ii,nsets
  integer, allocatable, dimension(:) :: nrad
  real :: delm
  real(8) :: rm,bin_factor
  integer :: jx,jy,jz
  
  integer :: mype,ierror
  common/mpi/ mype

 111 format(a20,i3)

   
!//////////Initialize everything to zero. 
  x = 0.0
  y = 0.0
  yrel = 0.0
  z = 0.0
  r = 0.0
  dr3 = 0.0
  delm = h*real(n)/real(m) !2deltx
!/////////define the 3D-coords for each droplet in one core///////////
  do id=1,ndropp
     idp(id) = mype*ndropp+id
     x(id) = rangen(seed)*h*N1*1.0d0
     y(id) = mype*h*n2pe + rangen(seed)*h*n2pe*1.0d0
     yrel(id) = y(id) - slabw*h
     if (yrel(id) .le. 0.0) then ! HAS TOO SMALL Of A YREL at ini
        yrel(id) = 0.0 - yrel(id) + h/100.0
        y(id) = yrel(id) + slabw*h
     elseif (yrel(id) .ge. real(n2pe)*h) then ! HAS TOO BIG OF A YREL at ini
        yrel(id) =  2.0d0*h*real(n2pe) - yrel(id) - h/100.0
        y(id) = yrel(id) + slabw*h
     endif
     z(id) = rangen(seed)*h*N3*1.0d0
  enddo


!/////////size dispersion (DSD)/////////////
   !currently most disp>1000 only works for N=128 eps=500 npe=64
   !for those disp, same shape of natural CCN as 6094 case but added seedings
   
  if (disp .eq. 1) then
      allocate(nrad(1)) 
         r_ccn=1.d-7
         kappa = 0.7053 !ammonium sulfate
      if (isolu .ge. 1) then
         if (mype .eq.0 ) print*,'retrieve initial droplet size'
         call wetradius() 
      else
         if (mype .eq. 0) print*,'isolu is off, prescribed droplet size'
         r=20.d-6
      endif
  elseif (disp .eq. 2) then
      allocate(nrad(2))
      if(isolu .ge. 1)then
      do id=1,ndropp,2
         r_ccn(id)=1.d-8
         r_ccn(id+1)=2.d-8
         kappa(id) = 1.33 !NaCl
         kappa(id+1) = 0.72 !Ammonium sulfate
      enddo
      else
      do id=1,ndropp,2
         r(id)=5.d-6
         r(id+1)=10.d-6
      enddo
      endif !isolu
  elseif (disp .eq. 6) then !r=5,10,15,20,25,30
      allocate(nrad(6))
      r_ccn=1.d-7
      kappa=0.7
      if (isolu .ge. 1) then
         call wetradius()
      else
      	do id = 1,ndropp,6
	 do ii=1,6
	    r(id+ii-1)=5.d-6+real(ii-1)*5.d-6
	 enddo!ii
      	enddo!id
      endif
  elseif (disp .eq. 6094) then
      !only work for N128 E500 
 	   !same shape as 12188 case but with number concentration of 86 cm^-3
      ! allocate droplet size (wet & dry)
      allocate(nrad(25))
      nrad(1) = 341
      nrad(2) = 884
      nrad(3) = 5285
      nrad(4) = 6035
      nrad(5) = 6070!69	  !6
      nrad(6) = 6079!18   !7
      nrad(7) = 6088!20	  !8
      nrad(8) = 6091!5	  !10
      nrad(9) = 6092!3    !11
      nrad(10)= 6093
      nrad(11) =6094!1   !12
     r(1     : nrad(1))     = 2.d-6
     r(1+nrad(1):nrad(2))   = 3.d-6
     r(1+nrad(2):nrad(3))   =4.d-6
     r(1+nrad(3):nrad(4))   =5.d-6
     r(1+nrad(4):nrad(5))   =6.d-6
     r(1+nrad(5):nrad(6))   =7.d-6
     r(1+nrad(6):nrad(7))   =8.d-6
     r(1+nrad(7):nrad(8))   =10.d-6
     r(1+nrad(8):nrad(9))   =12.d-6
     r(1+nrad(9):nrad(10))  =14.d-6
     r(1+nrad(10):nrad(11)) =16.d-6
     nrad(1)=341
     r_ccn(1        :nrad(1))   = 1.5119d-8
     nrad(2)=635
     nrad(3)=884
     r_ccn(1+nrad(1):nrad(2))   = 1.9049d-8
     r_ccn(1+nrad(2):nrad(3))   = 2.4000d-8
     nrad(4)=1096
     nrad(5)=1289
     nrad(6)=1507
     nrad(7)=1821
     nrad(8)=2312
     nrad(9)=3010
     nrad(10)=3843
     nrad(11)=4651
     nrad(12)=5285
     r_ccn(1+nrad(3):nrad(4))   = 3.0238d-8
     r_ccn(1+nrad(4):nrad(5))   = 3.8098d-8
     r_ccn(1+nrad(5):nrad(6))   = 4.80d-8
     r_ccn(1+nrad(6):nrad(7))   = 6.05d-8
     r_ccn(1+nrad(7):nrad(8))   = 7.62d-8
     r_ccn(1+nrad(8):nrad(9))   = 9.6d-8
     r_ccn(1+nrad(9):nrad(10))  = 1.21d-7
     r_ccn(1+nrad(10):nrad(11)) = 1.52d-7
     r_ccn(1+nrad(11):nrad(12)) = 1.92d-7
     nrad(13)=5687
     nrad(14)=5897
     nrad(15)=5993
     nrad(16)=6035
     r_ccn(1+nrad(12):nrad(13)) = 2.419d-7
     r_ccn(1+nrad(13):nrad(14)) = 3.048d-7
     r_ccn(1+nrad(14):nrad(15)) = 3.84d-7
     r_ccn(1+nrad(15):nrad(16)) = 4.84d-7
     nrad(17)=6056
     nrad(18)=6070
     r_ccn(1+nrad(16):nrad(17)) = 6.096d-7
     r_ccn(1+nrad(17):nrad(18)) = 7.68d-7
     nrad(19)=6078
     r_ccn(1+nrad(18):nrad(19)) = 9.676d-7
     nrad(20)=6084
     nrad(21)=6088
     r_ccn(1+nrad(19):nrad(20)) = 1.2191d-6
     r_ccn(1+nrad(20):nrad(21)) = 1.536d-6
     nrad(22)=6091
     r_ccn(1+nrad(21):nrad(22)) = 1.9352d-6
     nrad(23)=6092
     r_ccn(1+nrad(22):nrad(23)) = 2.4382d-6
     nrad(24)=6093
     nrad(25)=6094
     r_ccn(1+nrad(23):nrad(24)) = 3.072d-6
     r_ccn(1+nrad(24):nrad(25)) = 4.8765d-6

  elseif (disp .eq. 6796) then
      !only work for N=128 eps=500 npe=64
      !same shape as 6094 case but added 10 cm^-3 seeding
      allocate(nrad(26))
   !! allocate droplet size
      !natural aerosol
      nrad(1) = 341
      nrad(2) = 884
      nrad(3) = 5285
      nrad(4) = 6035
      nrad(5) = 6070!69	  !6
      nrad(6) = 6079!18   !7
      nrad(7) = 6088!20	  !8
      nrad(8) = 6091!5	  !10
      nrad(9) = 6092!3    !11
      nrad(10)= 6093
      nrad(11) =6094!1   !12
      r(1     : nrad(1))     = 2.d-6
      r(1+nrad(1):nrad(2))   = 3.d-6
      r(1+nrad(2):nrad(3))   =4.d-6
      r(1+nrad(3):nrad(4))   =5.d-6
      r(1+nrad(4):nrad(5))   =6.d-6
      r(1+nrad(5):nrad(6))   =7.d-6
      r(1+nrad(6):nrad(7))   =8.d-6
      r(1+nrad(7):nrad(8))   =10.d-6
      r(1+nrad(8):nrad(9))   =12.d-6
      r(1+nrad(9):nrad(10))  =14.d-6
      r(1+nrad(10):nrad(11)) =16.d-6
      kappa(1:nrad(11)) = 1.33 ! hygroscopicity (NaCl=1.33, (NH4)2SO4=0.72)
      !seeded aerosol
      nrad(12) =6796
      r(1+nrad(11):nrad(12)) =4.d-6
      kappa(1+nrad(11):nrad(12)) = 0.72
   !! allocate dry CCN size
      ! natural 
      nrad(1)=341
      r_ccn(1        :nrad(1))   = 1.5119d-8
      nrad(2)=635
      nrad(3)=884
      r_ccn(1+nrad(1):nrad(2))   = 1.9049d-8
      r_ccn(1+nrad(2):nrad(3))   = 2.4000d-8
      nrad(4)=1096
      nrad(5)=1289
      nrad(6)=1507
      nrad(7)=1821
      nrad(8)=2312
      nrad(9)=3010
      nrad(10)=3843
      nrad(11)=4651
      nrad(12)=5285
      r_ccn(1+nrad(3):nrad(4))   = 3.0238d-8
      r_ccn(1+nrad(4):nrad(5))   = 3.8098d-8
      r_ccn(1+nrad(5):nrad(6))   = 4.80d-8
      r_ccn(1+nrad(6):nrad(7))   = 6.05d-8
      r_ccn(1+nrad(7):nrad(8))   = 7.62d-8
      r_ccn(1+nrad(8):nrad(9))   = 9.6d-8
      r_ccn(1+nrad(9):nrad(10))  = 1.21d-7
      r_ccn(1+nrad(10):nrad(11)) = 1.52d-7
      r_ccn(1+nrad(11):nrad(12)) = 1.92d-7
      nrad(13)=5687
      nrad(14)=5897
      nrad(15)=5993
      nrad(16)=6035
      r_ccn(1+nrad(12):nrad(13)) = 2.419d-7
      r_ccn(1+nrad(13):nrad(14)) = 3.048d-7
      r_ccn(1+nrad(14):nrad(15)) = 3.84d-7
      r_ccn(1+nrad(15):nrad(16)) = 4.84d-7
      nrad(17)=6056
      nrad(18)=6070
      r_ccn(1+nrad(16):nrad(17)) = 6.096d-7
      r_ccn(1+nrad(17):nrad(18)) = 7.68d-7
      nrad(19)=6078
      r_ccn(1+nrad(18):nrad(19)) = 9.676d-7
      nrad(20)=6084
      nrad(21)=6088
      r_ccn(1+nrad(19):nrad(20)) = 1.2191d-6
      r_ccn(1+nrad(20):nrad(21)) = 1.536d-6
      nrad(22)=6091
      r_ccn(1+nrad(21):nrad(22)) = 1.9352d-6
      nrad(23)=6092
      r_ccn(1+nrad(22):nrad(23)) = 2.4382d-6
      nrad(24)=6093
      nrad(25)=6094
      r_ccn(1+nrad(23):nrad(24)) = 3.072d-6
      r_ccn(1+nrad(24):nrad(25)) = 4.8765d-6
      !seeded
      nrad(26)=6796
      r_ccn(1+nrad(25):nrad(26)) = 9.6d-8

  elseif (disp .eq. 6800) then
          !only work for N=128 eps=500 npe=64
          !GCCN
     allocate(nrad(26))
     nrad(1) = 341
     nrad(2) = 884
     nrad(3) = 5285
     nrad(4) = 6035
     nrad(5) = 6070!69	  !6
     nrad(6) = 6079!18   !7
     nrad(7) = 6088!20	  !8
     nrad(8) = 6091!5	  !10
     nrad(9) = 6092!3    !11
     nrad(10)= 6093
     nrad(11) =6094!1   !12
     nrad(12) =6796!seeding 
     r(1     : nrad(1))     = 2.d-6
     r(1+nrad(1):nrad(2))   = 3.d-6
     r(1+nrad(2):nrad(3))   =4.d-6
     r(1+nrad(3):nrad(4))   =5.d-6
     r(1+nrad(4):nrad(5))   =6.d-6
     r(1+nrad(5):nrad(6))   =7.d-6
     r(1+nrad(6):nrad(7))   =8.d-6
     r(1+nrad(7):nrad(8))   =10.d-6
     r(1+nrad(8):nrad(9))   =12.d-6
     r(1+nrad(9):nrad(10))  =14.d-6
     r(1+nrad(10):nrad(11)) =16.d-6
     r(1+nrad(11):nrad(12)) =8.d-6
     nrad(1)=341
     r_ccn(1        :nrad(1))   = 1.5119d-8
     nrad(2)=635
     nrad(3)=884
     r_ccn(1+nrad(1):nrad(2))   = 1.9049d-8
     r_ccn(1+nrad(2):nrad(3))   = 2.4000d-8
     nrad(4)=1096
     nrad(5)=1289
     nrad(6)=1507
     nrad(7)=1821
     nrad(8)=2312
     nrad(9)=3010
     nrad(10)=3843
     nrad(11)=4651
     nrad(12)=5285
     r_ccn(1+nrad(3):nrad(4))   = 3.0238d-8
     r_ccn(1+nrad(4):nrad(5))   = 3.8098d-8
     r_ccn(1+nrad(5):nrad(6))   = 4.80d-8
     r_ccn(1+nrad(6):nrad(7))   = 6.05d-8
     r_ccn(1+nrad(7):nrad(8))   = 7.62d-8
     r_ccn(1+nrad(8):nrad(9))   = 9.6d-8
     r_ccn(1+nrad(9):nrad(10))  = 1.21d-7
     r_ccn(1+nrad(10):nrad(11)) = 1.52d-7
     r_ccn(1+nrad(11):nrad(12)) = 1.92d-7
     nrad(13)=5687
     nrad(14)=5897
     nrad(15)=5993
     nrad(16)=6035
     r_ccn(1+nrad(12):nrad(13)) = 2.419d-7
     r_ccn(1+nrad(13):nrad(14)) = 3.048d-7
     r_ccn(1+nrad(14):nrad(15)) = 3.84d-7
     r_ccn(1+nrad(15):nrad(16)) = 4.84d-7
     nrad(17)=6056
     nrad(18)=6070
     r_ccn(1+nrad(16):nrad(17)) = 6.096d-7
     r_ccn(1+nrad(17):nrad(18)) = 7.68d-7
     nrad(19)=6078
     r_ccn(1+nrad(18):nrad(19)) = 9.676d-7
     nrad(20)=6084
     nrad(21)=6088
     r_ccn(1+nrad(19):nrad(20)) = 1.2191d-6
     r_ccn(1+nrad(20):nrad(21)) = 1.536d-6
     nrad(22)=6091
     r_ccn(1+nrad(21):nrad(22)) = 1.9352d-6
     nrad(23)=6092
     r_ccn(1+nrad(22):nrad(23)) = 2.4382d-6
     nrad(24)=6093
     nrad(25)=6094
     r_ccn(1+nrad(23):nrad(24)) = 3.072d-6
     r_ccn(1+nrad(24):nrad(25)) = 4.8765d-6
     nrad(26)=6796
     r_ccn(1+nrad(25):nrad(26)) = 1.d-6


  elseif (disp .eq. 7498) then
          !only work for N=128 eps=500 npe=64
          !double seeding
     allocate(nrad(26))
     nrad(1) = 341
     nrad(2) = 884
     nrad(3) = 5285
     nrad(4) = 6035
     nrad(5) = 6070!69	  !6
     nrad(6) = 6079!18   !7
     nrad(7) = 6088!20	  !8
     nrad(8) = 6091!5	  !10
     nrad(9) = 6092!3    !11
     nrad(10)= 6093
     nrad(11) =6094!1   !12
     nrad(12) =7498!seeding 
     r(1     : nrad(1))     = 2.d-6
     r(1+nrad(1):nrad(2))   = 3.d-6
     r(1+nrad(2):nrad(3))   =4.d-6
     r(1+nrad(3):nrad(4))   =5.d-6
     r(1+nrad(4):nrad(5))   =6.d-6
     r(1+nrad(5):nrad(6))   =7.d-6
     r(1+nrad(6):nrad(7))   =8.d-6
     r(1+nrad(7):nrad(8))   =10.d-6
     r(1+nrad(8):nrad(9))   =12.d-6
     r(1+nrad(9):nrad(10))  =14.d-6
     r(1+nrad(10):nrad(11)) =16.d-6
     r(1+nrad(11):nrad(12)) =4.d-6
     nrad(1)=341
     r_ccn(1        :nrad(1))   = 1.5119d-8
     nrad(2)=635
     nrad(3)=884
     r_ccn(1+nrad(1):nrad(2))   = 1.9049d-8
     r_ccn(1+nrad(2):nrad(3))   = 2.4000d-8
     nrad(4)=1096
     nrad(5)=1289
     nrad(6)=1507
     nrad(7)=1821
     nrad(8)=2312
     nrad(9)=3010
     nrad(10)=3843
     nrad(11)=4651
     nrad(12)=5285
     r_ccn(1+nrad(3):nrad(4))   = 3.0238d-8
     r_ccn(1+nrad(4):nrad(5))   = 3.8098d-8
     r_ccn(1+nrad(5):nrad(6))   = 4.80d-8
     r_ccn(1+nrad(6):nrad(7))   = 6.05d-8
     r_ccn(1+nrad(7):nrad(8))   = 7.62d-8
     r_ccn(1+nrad(8):nrad(9))   = 9.6d-8
     r_ccn(1+nrad(9):nrad(10))  = 1.21d-7
     r_ccn(1+nrad(10):nrad(11)) = 1.52d-7
     r_ccn(1+nrad(11):nrad(12)) = 1.92d-7
     nrad(13)=5687
     nrad(14)=5897
     nrad(15)=5993
     nrad(16)=6035
     r_ccn(1+nrad(12):nrad(13)) = 2.419d-7
     r_ccn(1+nrad(13):nrad(14)) = 3.048d-7
     r_ccn(1+nrad(14):nrad(15)) = 3.84d-7
     r_ccn(1+nrad(15):nrad(16)) = 4.84d-7
     nrad(17)=6056
     nrad(18)=6070
     r_ccn(1+nrad(16):nrad(17)) = 6.096d-7
     r_ccn(1+nrad(17):nrad(18)) = 7.68d-7
     nrad(19)=6078
     r_ccn(1+nrad(18):nrad(19)) = 9.676d-7
     nrad(20)=6084
     nrad(21)=6088
     r_ccn(1+nrad(19):nrad(20)) = 1.2191d-6
     r_ccn(1+nrad(20):nrad(21)) = 1.536d-6
     nrad(22)=6091
     r_ccn(1+nrad(21):nrad(22)) = 1.9352d-6
     nrad(23)=6092
     r_ccn(1+nrad(22):nrad(23)) = 2.4382d-6
     nrad(24)=6093
     nrad(25)=6094
     r_ccn(1+nrad(23):nrad(24)) = 3.072d-6
     r_ccn(1+nrad(24):nrad(25)) = 4.8765d-6
     nrad(26)=7498
     r_ccn(1+nrad(25):nrad(26)) = 9.6d-8
     


  elseif (disp .eq. 12188) then 
          !lulin2010 maritime case log-normal distribution after activation in up=2m/s, Spmax=1.59% doubled
          !number concentration
     allocate(nrad(26))
     nrad(1) =  682 !2 micron
     nrad(2) = 1769!1087 !3
     nrad(3) = 10569!8800 !4
     nrad(4) = 12070!1501 !5
     nrad(5) = 12139!69	  !6
     nrad(6) = 12157!18   !7
     nrad(7) = 12177!20	  !8
     nrad(8) = 12182!5	  !10
     nrad(9) = 12185!3    !11
     nrad(10) =12186!1   !12
     nrad(11) =12187!1   !14
     nrad(12) =12188!1   !16
     r(1     : nrad(1))     = 2.d-6
     r(1+nrad(1):nrad(2))   = 3.d-6
     r(1+nrad(2):nrad(3))   =4.d-6
     r(1+nrad(3):nrad(4))   =5.d-6
     r(1+nrad(4):nrad(5))   =6.d-6
     r(1+nrad(5):nrad(6))   =7.d-6
     r(1+nrad(6):nrad(7))   =8.d-6
     r(1+nrad(7):nrad(8))   =10.d-6
     r(1+nrad(8):nrad(9))   =11.d-6
     r(1+nrad(9):nrad(10))  =12.d-6
     r(1+nrad(10):nrad(11)) =14.d-6
     r(1+nrad(11):nrad(12)) =16.d-6
     nrad(1)=682
     r_ccn(1        :nrad(1))   = 1.5119d-8
     nrad(2)=1270
     nrad(3)=1769
     r_ccn(1+nrad(1):nrad(2))   = 1.9049d-8
     r_ccn(1+nrad(2):nrad(3))   = 2.4000d-8
     nrad(4)=2192
     nrad(5)=2579
     nrad(6)=3014
     nrad(7)=3642
     nrad(8)=4624
     nrad(9)=6021
     nrad(10)=7686
     nrad(11)=9302
     nrad(12)=10569
     r_ccn(1+nrad(3):nrad(4))   = 3.0238d-8
     r_ccn(1+nrad(4):nrad(5))   = 3.8098d-8
     r_ccn(1+nrad(5):nrad(6))   = 4.80d-8
     r_ccn(1+nrad(6):nrad(7))   = 6.05d-8
     r_ccn(1+nrad(7):nrad(8))   = 7.62d-8
     r_ccn(1+nrad(8):nrad(9))   = 9.6d-8
     r_ccn(1+nrad(9):nrad(10))  = 1.21d-7
     r_ccn(1+nrad(10):nrad(11)) = 1.52d-7
     r_ccn(1+nrad(11):nrad(12)) = 1.92d-7

     nrad(13)=11373
     nrad(14)=11794
     nrad(15)=11985
     nrad(16)=12070
     r_ccn(1+nrad(12):nrad(13)) = 2.419d-7
     r_ccn(1+nrad(13):nrad(14)) = 3.048d-7
     r_ccn(1+nrad(14):nrad(15)) = 3.84d-7
     r_ccn(1+nrad(15):nrad(16)) = 4.84d-7
     nrad(17)=12113
     nrad(18)=12139
     r_ccn(1+nrad(16):nrad(17)) = 6.096d-7
     r_ccn(1+nrad(17):nrad(18)) = 7.68d-7
     nrad(19)=12157
     r_ccn(1+nrad(18):nrad(19)) = 9.676d-7
     nrad(20)=12169
     nrad(21)=12177
     r_ccn(1+nrad(19):nrad(20)) = 1.2191d-6
     r_ccn(1+nrad(20):nrad(21)) = 1.536d-6
     nrad(22)=12182
     r_ccn(1+nrad(21):nrad(22)) = 1.9352d-6
     nrad(23)=12185
     r_ccn(1+nrad(22):nrad(23)) = 2.4382d-6
     nrad(24)=12186
     nrad(25)=12187
     nrad(26)=12188
     r_ccn(1+nrad(23):nrad(24)) = 3.072d-6
     r_ccn(1+nrad(24):nrad(25)) = 3.8705d-6
     r_ccn(1+nrad(25):nrad(26)) = 4.8765d-6

  elseif (disp .eq. 80) then !Raga1990continuous
     radmin = 1.d-6
     allocate(nrad(16))
     nsets = ndropp/disp
     nrad(1) = 1
     nrad(2) = 2
     nrad(3) = 3
     nrad(4) = 5
     nrad(5) = 8
     nrad(6) = 11
     nrad(7) = 15
     nrad(8) = 21
     nrad(9) = 28
     nrad(10) = 41
     nrad(11) = 63
     nrad(12) = 69
     nrad(13) = 73
     nrad(14) = 76
     nrad(15) = 78
     nrad(16) = 80
     do id = 1,nsets
        r( (id-1)*disp+1:(id-1)*disp+nrad(1)) = 5.*radmin
        r( (id-1)*disp+1+nrad(1):(id-1)*disp+nrad(2)) = 6.*radmin
        r( (id-1)*disp+1+nrad(2):(id-1)*disp+nrad(3)) = 7.*radmin
        r( (id-1)*disp+1+nrad(3):(id-1)*disp+nrad(4)) = 8.*radmin
        r( (id-1)*disp+1+nrad(4):(id-1)*disp+nrad(5)) = 9.*radmin
        r( (id-1)*disp+1+nrad(5):(id-1)*disp+nrad(6)) = 10.*radmin
        r( (id-1)*disp+1+nrad(6):(id-1)*disp+nrad(7)) = 11.*radmin
        r( (id-1)*disp+1+nrad(7):(id-1)*disp+nrad(8)) = 12.*radmin
        r( (id-1)*disp+1+nrad(8):(id-1)*disp+nrad(9)) = 13.*radmin
        r( (id-1)*disp+1+nrad(9):(id-1)*disp+nrad(10)) = 14.*radmin
        r( (id-1)*disp+1+nrad(10):(id-1)*disp+nrad(11)) = 15.*radmin
        r( (id-1)*disp+1+nrad(11):(id-1)*disp+nrad(12)) = 16.*radmin
        r( (id-1)*disp+1+nrad(12):(id-1)*disp+nrad(13)) = 17.*radmin
        r( (id-1)*disp+1+nrad(13):(id-1)*disp+nrad(14)) = 18.*radmin
        r( (id-1)*disp+1+nrad(14):(id-1)*disp+nrad(15)) = 19.*radmin
        r( (id-1)*disp+1+nrad(15):(id-1)*disp+nrad(16)) = 20.*radmin
     enddo

  elseif (disp .eq. 76) then
     radmin = 5.0d-6
     nsets = ndropp/disp                     ! number of droplets sets                    
     allocate(nrad(5))
     nrad(1) = 10                                 ! there are nrad1 drops with radius rad1 per set        
     nrad(2) = 30                                 ! there are nrad2-nrad1 drops radius 2*rad1 per set, etc       
     nrad(3) = 65
     nrad(4) = 75
     nrad(5) = 76 
     do id = 1,nsets
        r( (id-1)*disp+1:(id-1)*disp+nrad(1)) = radmin
        r( (id-1)*disp+1+nrad(1):(id-1)*disp+nrad(2)) = 2.*radmin
        r( (id-1)*disp+1+nrad(2):(id-1)*disp+nrad(3)) = 3.*radmin
        r( (id-1)*disp+1+nrad(3):(id-1)*disp+nrad(4)) = 4.*radmin
        r( (id-1)*disp+1+nrad(4):(id-1)*disp+nrad(5)) = 5.*radmin
     enddo
  endif ! disp
  deallocate(nrad)
!get dsd
  dsd = 0
!  dsd_log2=0
  ccn_log2=0
  radmin=minval(r(1:ndropp))
  do id = 1,ndropp
     entr = Nint(r(id)/1.d-6)
     if (entr .le. 1) then
	  dsd(1) = dsd(1)+1
     else
	  dsd(entr) = dsd(entr)+1
     endif
     if (isolu*thermo .eq. 1) then
      entr = Nint(log2(r_ccn(id)/r_ccnmin))
      if (entr .le. 1) then
         ccn_log2(1)=ccn_log2(1)+1
      else
         ccn_log2(entr)=ccn_log2(entr)+1
      endif
     endif !isolu
     
  enddo
  if(mype .eq. 0 )then
  dsdtot=dsd*npe
 ! dsdtot_log2=dsd_log2*npe 
  ccntot_log2=ccn_log2*npe
  endif
 ! Compute average radius and initialize inteS (INTESUP)
  rm=0.0d0
  do id = 1,ndropp
     rm = rm + r(id)**3
!         integsup(id)=0.0d0
!         integsup2(id)=0.0d0
  enddo
     rm=(rm/ndropp)**(1.d0/3.d0)
!deallocate(integsup,integsup2)
     if(mype .eq. 0) then
         write(45,145) 0.d0, (dsdtot(i), i=1,nbins) 
         write(46,145) 0.d0, (ccntot_log2(i),i=1,nbins)
     endif

 299   format(1x,(e12.6),3(f8.6))
 199   format(1x,3(f8.6))
 145   format(1x,101(i10,1x))

  end SUBROUTINE IDROPS  
 
     real FUNCTION rangen(seed)
        implicit none
        integer, parameter :: a = 16807
        integer, parameter :: m = 2147483647
        integer, parameter :: q = 127773
        integer, parameter :: r = 2836
        integer ::  seed, lo, hi, test

        hi = seed/q
        lo = mod(seed,q)
        test = a*lo - r*hi
        if (test .gt. 0) then
           seed = test
        else
           seed = test + m
        endif
        rangen = real(seed)/real(m) 
     end FUNCTION rangen 
 
   SUBROUTINE wetradius()
   ! this subroutine calculate the equilibrium radius at current thermodynamic environment
   use mic_mod
   implicit none
   include 'param.inc'
   include 'mpif.h'

   integer :: id
   real(8) :: esat,seq1,seq2
   real(8) :: curv, solu, sp_c
   real, parameter :: sp_init=-1.d-2
   integer, allocatable, dimension(:) :: nrad
   integer :: mype,ierror
   common/mpi/ mype
!/////////calculate necessary parameters& coefficients
   diffvnd1=1.d-5*(0.015*T0-1.9)
   ka1=1.5d-11*T0**3-4.8d-8*T0**2+1.d-4*T0-3.9d-4
   esat = 2.53d11*exp(-5.42d3/T0)
   !--------------first guess the radius of wet aerosol--------!
   !--------------start with r_wet=1.5*r_d-----------------------!
   r=r_ccn*1.5d0
   !--------------iterate until get equilibrium saturation status around each droplet surface----!
   if(mype .eq. 0) print*, 'start iteration to get wet radius'
   sp_c=min(sp0,sp_init)
   if(mype .eq. 0) print*, 'initial supersaturation for computing equilibrium wet radius=',sp_c
   do id=1,ndropp
      seq=sp_c+1.d-2
      seq1=seq+1.d-2
      seq2=seq1+1.d-2
      do while(abs(seq-sp_c) .gt. 1.d-7 .and. seq2 .ne. seq)
         seq2=seq1
         seq1=seq
         diffvnd2=diffvnd1/(r(id)/(r(id)+0.104d-6)+diffvnd1/(r(id)*0.036)*sqrt(2.d0*pi/(Rv*T0)))
         ka2=ka1/(r(id)/(r(id)+.216d-6)+ka1/(r(id)*.7*rhoa*Cp)*sqrt(2.d0*pi/(Ra*T0)))
         ks = 1.d0/(rhow*Rv*T0/(esat*diffvnd2)+rhow*Lat/(Ka2*T0)*(Lat/(Rv*T0)-1))
         curv=2.d0*sigma_sa/(Rv*rhow*T0*r(id))
         solu=(r(id)**3-r_ccn(id)**3)/(r(id)**3-(1-kappa(id))*r_ccn(id)**3)
         seq=solu*exp(curv)-1.d0!equillibrium supersat.
         if(seq .gt. sp_c .and. r(id) .gt. r_ccn(id)) then
	         r(id)=r(id)-r_ccn(id)*1.d-2
         elseif(seq .lt. sp_c) then 
   	      r(id)=r(id)+r_ccn(id)*1.d-2
         endif

      enddo !dowhile
   enddo !drop id
   if (mype .eq. 0) print*, 'get wet radius'
end SUBROUTINE wetradius
