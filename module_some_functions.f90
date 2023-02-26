module some_functions
  implicit none
  integer,parameter :: qp=16

contains



!    ***************************************************************************
!    *   the cubic spline interpolation based on the three moment function     *
!    *                     free boundary condition by Dr. Zhou Xin             *
!    ***************************************************************************

Subroutine spline3(n, x, y, m, t, u)
  implicit none
  integer :: n,m,l,i,j,j1,j2
  real(kind=8) ::  x(n), y(n), t(m), u(m), a(n), b(n), c(n), d(n)
  real(kind=8) :: ss,tt,rr,e,f
  ! x(n),y(n) given values
  ! t(m),u(m) return values

  a(1) = 0.0D0
  a(n) = 1.0D0
  d(1) = 0.0D0
  d(n) = 0.0D0
  c(1) = -1.0D0
  c(n) = 0.0D0
  b(1) = 1.0D0
  b(n) = -1.0D0
  l = n - 1

  Do i = 2, l
    a(i) = (x(i)-x(i-1))/6.0D0
    c(i) = (x(i+1)-x(i))/6.0D0
    b(i) = 2.0D0*(a(i)+c(i))
    d(i) = (y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))
  End Do
  c(1) = c(1)/b(1)
  d(1) = d(1)/b(1)

  Do i = 2, n
    c(i) = c(i)/(b(i)-a(i)*c(i-1))
    d(i) = (d(i)-a(i)*d(i-1))/(b(i)-a(i)*c(i-1))
  End Do
  a(n) = d(n)

  Do i = 1, l
    j = n - i
    a(j) = d(j) - c(j)*a(j+1)
  End Do

  Do j1 = 1, m
    f = t(j1)
    Do j2 = 1, n - 1
      If (x(j2)<=f .And. f<=x(j2+1)) Goto 25
    End Do
    Goto 30
    25  e = x(j2+1) - x(j2)
    rr = (a(j2)*(x(j2+1)-f)**3+a(j2+1)*(f-x(j2))**3)/6.0D0/e
    ss = (x(j2+1)-f)*(y(j2)/e-a(j2)*e/6.0D0)
    tt = (f-x(j2))*(y(j2+1)/e-a(j2+1)*e/6.0D0)
    u(j1) = rr + ss + tt

    30 End Do

    Return
End Subroutine spline3

  logical function sammesign(num1,num2)
    implicit none
    real(kind=8),intent(in) :: num1,num2
    if ((num1>0 .and. num2>0) .or. (num1<0 .and. num2<0)) then
      sammesign=.true.
    else
      sammesign=.false.
    end if
  end function sammesign

  integer function count_lines(ifile)
    implicit none
    character(len=*),intent(in) :: ifile
    integer :: io

    open(unit=12,file=ifile,action='read',status='old')
    count_lines=0
    do
      read (12,*, iostat=io)
      if (io/=0) exit
      count_lines=count_lines+1
    end do
    close(unit=12)

  end function count_lines

  subroutine kahan_sum(x,insums)
    ! KAHAN_SUM Computes the sum of a given set of numbers using the Kahan summation algorithm.
    !
    ! Inputs:
    ! - x: An array of numbers to be summed.
    !
    ! Outputs:
    ! - s: The final sum of the numbers in x.
    ! - c: The final correction value used in the Kahan summation algorithm.
    ! - intermediate_sums: An array of the intermediate sums computed during the Kahan summation process.

    implicit none
    real(kind=8),intent(in)  :: x(:)
    real(kind=8),intent(out) :: insums(:)
    integer :: i
    real(kind=8) :: s,c,y,t


    s = 0d0
    c = 0d0
    insums=0d0
    do i = 1,size(x)
      y = x(i) - c
      t = s + y
      c = (t - s) - y
      s = t
      insums(i) = s
    enddo

  end subroutine

  subroutine toAlternatingSeries(an,bn)
    implicit none
    real(kind=8),intent(in) :: an(:)
    real(kind=8),intent(out):: bn(:) ! size of bn = size of an + 2
    integer :: n,i,j
    n=size(an)
    if (n<1) then
      return
    else
      bn(1)=+1.0d0
      if (an(1)>-tiny(0.0d0)) then
        bn(2)=+1.0d0
      else
        bn(2)=-1.0d0
      end if

      bn(3)=an(1)
      j=3
      do i=2,n
        if (sammesign(bn(j),an(i))) then
          bn(j)=bn(j)+an(i)
        else
          j=j+1
          bn(j)=an(i)
        end if
      end do
      bn(1)=real(j,kind=qp)
    end if

  end subroutine toAlternatingSeries

  subroutine vanWijngaarden_sum(bn,sn)
    !// https://en.wikipedia.org/wiki/Van_Wijngaarden_transformation
    implicit none
    !// with two element before the alternative sereies b(1)=n-terms+2, b(2)=sign(1st-term)
    real(kind=8),intent(in) :: bn(:)

    !// form rows of averages between neighbors
    real(kind=8),intent(out) :: sn(3)

    real(kind=8) :: an(size(bn)-2)
    integer :: n,nmax0,nmax

    nmax=int(bn(1))
    nmax=nmax-2

    !// partial sum in an without truncation
    ! an(1)=bn(3)
    ! do n=1,nmax-1
    !   an(n+1)=an(n)+bn(n+3)
    ! end do
    call kahan_sum(bn(3:nmax+2),an(1:nmax))
    sn(3)=an(nmax)

    !// less than 15 terms do not use Euler's transform,
    !// Otherwise, when some angles are small, the value jumps.
    if(nmax<=18) then
      sn(1:2)=sn(3)
      return
    end if


    !// keeping the total terms in nmax0
    nmax0=nmax

    !// average of neighbors
    !// Repeat this operation until you get the last number.
    if (nmax==1) then
      sn(2)=an(1)
    else
      do while(nmax>0)
        do n=1,nmax-1
          an(n)=(an(n)+an(n+1))/2.0E0_qp
          if (nmax-1==ceiling((1E0_qp/3E0_qp)*nmax0)) then
            sn(2)=an(nmax-1)
            !// return at 2/3 of nmax0
          end if
        end do
        nmax=nmax-1
      end do
    end if
    sn(1)=an(1)

  end subroutine vanWijngaarden_sum

  subroutine legender(ctheta,pn,dpn)
    implicit none
    real(kind=qp),intent(in) :: ctheta
    real(kind=qp),intent(inout),dimension(0:) :: pn,dpn
    real(kind=qp) :: p0,p1,dp0,dp1,dp2,xni
    integer :: nmax,ni

    nmax=size(pn)-1
    p0=1.0_qp
    p1=ctheta
    dp0=0.0_qp
    dp1=-sqrt(1.0_qp-p1**2)
    dp2=3*ctheta*dp1

    if (nmax==0) then
      pn(0)=p0
      dpn(0)=dp0
    elseif (nmax==1) then
      pn(0)=p0
      pn(1)=p1
      dpn(0)=dp0
      dpn(1)=dp1
    elseif (nmax>=2) then
      pn(0)=p0
      pn(1)=p1
      dpn(0)=dp0
      dpn(1)=dp1
      dpn(2)=dp2

      do ni=2,nmax
        xni=real(ni,kind=qp)
        pn(ni)=-(1.0_qp-1.0_qp/xni)*pn(ni-2)+&
          (2.0_qp-1.0_qp/xni)*pn(ni-1)*ctheta
      end do

      do ni=3,nmax
        xni=real(ni,kind=qp)
        dpn(ni)=((2d0*xni-1)/(xni-1))*ctheta*dpn(ni-1)-(xni/(xni-1))*dpn(ni-2)
      end do

    end if

  end subroutine

  subroutine write_list_data(data,file_name,out_formate)
    implicit none
    real(kind=8),intent(in) :: data(:)
    character(len=*),intent(in) :: file_name
    character(len=*),intent(in) :: out_formate
    integer :: i

    open(unit=211,file=file_name,action='write')
    do i=1,size(data,1)
      write(211,out_formate) i,data(i)
    end do
    close(unit=211)
  end subroutine write_list_data

  subroutine write_table_data(data,file_name,out_formate)
    implicit none
    real(kind=8),intent(in) :: data(:,:)
    character(len=*),intent(in) :: file_name
    character(len=*),intent(in) :: out_formate
    integer :: i

    open(unit=211,file=file_name,action='write')
    do i=1,size(data,1)
      write(211,out_formate) i,data(i,:)
    end do
    close(unit=211)
  end subroutine write_table_data

  subroutine is_alternating_series(xn,flag,npterm)
    implicit none
    real(kind=8),intent(in) :: xn(:)
    logical,intent(out)     :: flag
    integer,intent(out)     :: npterm
    integer :: ii,jj,n

    flag=.false.
    jj=0
    n=size(xn)
    do ii=1,n
      if(xn(ii)>0.0D0) then
        jj=jj+1
      endif
    end do
    npterm=jj
    if ( abs(real(jj)/real(n)-0.5d0)<0.30d0 ) flag=.true.


  end subroutine is_alternating_series


end module some_functions
