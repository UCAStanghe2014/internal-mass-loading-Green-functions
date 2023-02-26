subroutine int_step(dep,step)
  implicit none
  real(kind=8),intent(in) ::  dep
  real(kind=8),intent(out) :: step
  if (dep<2.0) then
    step=0.1d0
  elseif (dep<25) then
    step=0.1d0
  else
    step=0.1d0
  endif
  IF (abs(dep-30.0)<1.0 ) step=0.01d0
  
  !step=1d0

end subroutine


subroutine yfile2green(yfile,depth)
  use  :: some_functions
  implicit none
  character(len=*) :: yfile
  real(kind=8),intent(in) :: depth

  real(kind=8),parameter :: pi=3.1415926535897932384626433832795E0_16
  integer,     parameter :: nangle_set=115
  integer,     parameter :: nmax_set=127000
  integer :: nmax,ii,jj,kk,npterm

  real(kind=8) :: angle(nangle_set)
  real(kind=8) ::    ns(0:nmax_set)
  real(kind=8) ::    y1(0:nmax_set)
  real(kind=8) ::    y3(0:nmax_set)
  real(kind=8) ::    y5(0:nmax_set)
  real(kind=16)::    pn(0:nmax_set)
  real(kind=16) ::  dpn(0:nmax_set)
  real(kind=8) :: green(0:nmax_set)
  real(kind=8) :: cumsumgreen(0:nmax_set)
  real(kind=8) :: green_sign(1:3+nmax_set)
  real(kind=8) :: green_sum(3,3),sn(3);
  integer :: flag_terms(3,4)=-1
  integer :: n_terms

  real(kind=8) :: yline(7)
  real(kind=8) :: ra,g0,G,tt

  real(kind=16) :: c
  logical :: flag
  integer :: flagid
  character(len=2) :: jj_str


! 156
!  angle=(/&
!      1d-5, 2d-5, 3d-5, 4d-5, 5d-5, 6d-5, 8d-5, 9d-5,             &
!      0.0001d0,  0.0002d0,  0.0003d0,  0.0004d0, 0.0005d0,        &
!      0.0006d0,  0.0007d0,  0.0008d0,  0.0009d0, 0.00095d0,       &
!      0.001d0,   0.0015d0,  0.0018d0,  0.002d0,  0.003d0, 0.005d0,&
!      0.006d0,   0.007d0,   0.008d0,   0.009d0,  0.010d0,&
!      0.012d0,   0.014d0,   0.016d0,   0.017d0,  0.018d0, 0.019d0,&
!      0.020d0,   0.025d0,   0.030d0,   0.040d0,  0.050d0, 0.060d0,&
!      0.070d0,   0.080d0,   0.090d0,   0.100d0,  0.110d0,         &
!      0.120d0,   0.130d0,   0.140d0,   0.150d0,  0.160d0,         &
!      0.170d0,   0.180d0,   0.190d0,   0.200d0,  0.210d0,         &
!      0.220d0,   0.230d0,   0.240d0,   0.250d0,  0.260d0,         &
!      0.270d0,   0.280d0,   0.290d0,   0.300d0,  0.310d0,         &
!      0.320d0,   0.330d0,   0.340d0,   0.350d0,  0.360d0,         &
!      0.370d0,   0.380d0,   0.390d0,   0.400d0,  0.410d0,         &
!      0.420d0,   0.430d0,   0.440d0,   0.450d0,  0.460d0,         &
!      0.470d0,   0.480d0,   0.490d0,   0.500d0,  0.550d0,         &
!      0.600d0,   0.650d0,   0.700d0,   0.750d0,  0.800d0,         &
!      0.850d0,   0.900d0,   0.950d0,   1.000d0,  1.100d0,         &
!      1.200d0,   1.300d0,   1.400d0,   1.500d0,  1.600d0,         &
!      2.000d0,   2.500d0,   3.000d0,   4.000d0,  5.000d0,         &
!      6.000d0,   7.000d0,   8.000d0,   9.000d0, 10.000d0,         &
!     12.000d0,  16.000d0,  20.000d0,  25.000d0, 30.000d0,         &
!     40.000d0,  50.000d0,  60.000d0,  70.000d0, 80.000d0,         &
!     90.000d0, 100.000d0, 110.000d0, 120.000d0,  130.000d0,       &
!    140.000d0, 150.000d0, 160.000d0, 170.000d0,  171.000d0,       &
!    172.000d0, 173.000d0, 174.000d0, 175.000d0,  176.000d0,       &
!    177.000d0, 178.000d0, 179.000d0, 179.100d0,  179.200d0,       &
!    179.300d0, 179.400d0, 179.500d0, 179.600d0,  179.700d0,       &
!    179.800d0, 179.900d0, 179.950d0, 179.980d0,  179.990d0,       &
!    179.992d0, 179.994d0, 179.996d0, 179.998d0,  179.999d0/)

! 115
      angle=(/&
         0.00010d0,&
         0.00020d0,&
         0.00030d0,&
         0.00040d0,&
         0.00050d0,&
         0.00060d0,&
         0.00070d0,&
         0.00080d0,&
         0.00090d0,&
         0.00095d0,&
         0.00100d0,&
         0.00200d0,&
         0.00300d0,&
         0.00400d0,&
         0.00500d0,&
         0.00600d0,&
         0.00700d0,&
         0.00800d0,&
         0.00900d0,&
         0.01000d0,&
         0.02000d0,&
         0.03000d0,&
         0.04000d0,&
         0.05000d0,&
         0.06000d0,&
         0.07000d0,&
         0.08000d0,&
         0.09000d0,&
         0.10000d0,&
         0.11000d0,&
         0.12000d0,&
         0.13000d0,&
         0.14000d0,&
         0.15000d0,&
         0.16000d0,&
         0.17000d0,&
         0.18000d0,&
         0.19000d0,&
         0.20000d0,&
         0.21000d0,&
         0.22000d0,&
         0.23000d0,&
         0.24000d0,&
         0.25000d0,&
         0.26000d0,&
         0.27000d0,&
         0.28000d0,&
         0.29000d0,&
         0.30000d0,&
         0.31000d0,&
         0.32000d0,&
         0.33000d0,&
         0.34000d0,&
         0.35000d0,&
         0.36000d0,&
         0.37000d0,&
         0.38000d0,&
         0.39000d0,&
         0.40000d0,&
         0.41000d0,&
         0.42000d0,&
         0.43000d0,&
         0.44000d0,&
         0.45000d0,&
         0.46000d0,&
         0.47000d0,&
         0.48000d0,&
         0.49000d0,&
         0.50000d0,&
         0.55000d0,&
         0.60000d0,&
         0.65000d0,&
         0.70000d0,&
         0.75000d0,&
         0.80000d0,&
         0.85000d0,&
         0.90000d0,&
         0.95000d0,&
         1.00000d0,&
         1.10000d0,&
         1.20000d0,&
         1.30000d0,&
         1.40000d0,&
         1.50000d0,&
         1.60000d0,&
         2.00000d0,&
         2.50000d0,&
         3.00000d0,&
         4.00000d0,&
         5.00000d0,&
         6.00000d0,&
         7.00000d0,&
         8.00000d0,&
         9.00000d0,&
        10.00000d0,&
        12.00000d0,&
        16.00000d0,&
        20.00000d0,&
        25.00000d0,&
        30.00000d0,&
        40.00000d0,&
        50.00000d0,&
        60.00000d0,&
        70.00000d0,&
        80.00000d0,&
        90.00000d0,&
       100.00000d0,&
       110.00000d0,&
       120.00000d0,&
       130.00000d0,&
       140.00000d0,&
       150.00000d0,&
       160.00000d0,&
       170.00000d0,&
       179.90000d0 &
      /)



  call interp1_love_file(yfile)



  flag_terms=-1
  g0 = 9.8156d0
  ra = 6371d3
  G  = 6.67d-11
  tt=1d0-(depth*1d3)/ra


  nmax=count_lines('DSLOV33.DAT')-1
  write(*,*) nmax

  if (nmax>nmax_set) then
    write (*,*) 'nmax_set is too small to read all the ylove file lines.'
    stop
  endif

  open(unit=34,file='DSLOV33.DAT',action='read')
  do ii=0,nmax,1
    read (34,*) yline
    ns(ii)=yline(1)
    y1(ii)=yline(2)
    y3(ii)=yline(4)
    y5(ii)=yline(6)
  enddo
  close(unit=34)

  !// insure the degree-0 degree-1 value of y5(a)
  y5(0)=G/ra
  y5(1)=G/ra*tt


  open(unit=60,file='GRNFN33.DAT',action='write')
  do kk=1,nangle_set

    c=cos(real(angle(kk),kind=16)/180E0_16*pi)
    call legender(c,pn(0:nmax),dpn(0:nmax))

    do jj=1,3

      !// for vertical displacement
      if (jj==1) then
        do ii=0,nmax
          green(ii) = y1(ii)*pn(ii)
        enddo
      endif

      if(jj==2) then
        !// for horizontal displacment
        do ii=0,nmax
          green(ii) = y3(ii)*dpn(ii)
        enddo
      endif

      if (jj==3) then
        !// for graivty change
        do ii=0,nmax
          green(ii)=(2*y1(ii)*g0/ra-(ii+1d0)/ra*(y5(ii)-(tt**ii)*G/ra))*pn(ii)
        enddo
      endif


      !// sum the series
      call is_alternating_series(green(0:nmax),flag,npterm)

      ! write(*,*) flag,npterm
      ! read(*,*)
      ! flag=.false.

      if (flag) then
        flagid=1
      else
        flagid=0
      endif

      if (flag) then !// approximate alternating series
        call toAlternatingSeries(green(0:nmax),green_sign(1:nmax+3))
        n_terms=int(green_sign(1))
        call vanWijngaarden_sum(green_sign(1:n_terms),sn)
      else           !// not an alternating series
        n_terms=0
        call kahan_sum(green(0:nmax),cumsumgreen(0:nmax))
        sn(1)=cumsumgreen(nmax)
        sn(2:3)=sn(1)
      end if

      green_sum(jj,:)=sn

      flag_terms(jj,1:4)=(/flagid,nmax,npterm,n_terms/)

    enddo



    write(60,'(F10.6,3(2X,1ES16.6),12(2X,I10))') &
      angle(kk),&
      green_sum(1:2,1)*(1.0d12*ra*angle(kk)/180d0*pi),&
      green_sum(3:3,1)*(1.0d18*ra*angle(kk)/180d0*pi),&
      flag_terms(1,1:4),flag_terms(2,1:4),flag_terms(3,1:4)

    write( *,'(F10.6,3(2X,1ES16.6),12(2X,I10))') &
      angle(kk),&
      green_sum(1:2,1)*(1.0d12*ra*angle(kk)/180d0*pi),&
      green_sum(3:3,1)*(1.0d18*ra*angle(kk)/180d0*pi),&
      flag_terms(1,1:4),flag_terms(2,1:4),flag_terms(3,1:4)

  enddo
  close(unit=60)

end subroutine yfile2green



