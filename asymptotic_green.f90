subroutine asymptotic_green(HL,delta,grs,theta,u)
  Use module_psum,only: pnm_sum,dpnm_sum
  implicit none
  real(kind=8),intent(in) :: HL(4,2) ! coefficient of asymptotic y1 and y3
  real(kind=8),intent(in) :: delta ! rs/ra
  real(kind=8),intent(in) :: grs   ! gravity at source [m/s^2]
  real(kind=8),intent(in) :: theta
  real(kind=8),intent(out):: u(2) ! ur, utheta

  real(kind=8) :: H(4),L(4),pi,c,s,ur,utheta
  real(kind=8) :: p1,p2,p3,p4,p5


  ! write(*,*) 'in asymptotic green ...'
  ! write(*,*) HL,delta,grs,theta
  ! stop
  ! write(*,*) 'in Green asymptotic: '
  ! write(*,*) HL,delta,grs,theta


  !! sum of {n, 1 ,1/n, 1/(n^2-1)}*pn for n = 2, ..., inf

  ! ur
  ! sum of legendre
  pi = 4d0*atan(1d0)
  c  = cos(theta/180D0*pi)
  H  = HL(:,1)
  L  = HL(:,2)

  ! write(*,*) 'pi = ',pi
  ! write(*,*) 'c = ',c

  ! 2阶以上求和
  call pnm_sum(2, delta, c, p1)
  p1 = p1-c
  call pnm_sum(1, delta, c, p2)
  p2 = p2-c
  call pnm_sum(5, delta, c, p3)
  p3 = p3-c
  call pnm_sum(6, delta, c, p4)
  call pnm_sum(7, delta, c, p5)
  p4 = (p4-p5)/2d0 + c/4d0

  !write(*,*) p1,p2,p3,p4

  ur = H(1)*p1 + H(2)*p2 + H(3)*p3 + H(4)*p4
  ur = (ur*delta)*grs

  ! utheta
  ! sum of legendre

  ! 2阶以上求和
  s = sqrt(1d0-c**2)
  call dpnm_sum(2, delta, c, p1)
  p1 = p1+s
  call dpnm_sum(1, delta, c, p2)
  p2 = p2+s
  call dpnm_sum(5, delta, c, p3)
  p3 = p3+s
  call dpnm_sum(6, delta, c, p4)
  call dpnm_sum(7, delta, c, p5)
  p4 = (p4-p5)/2d0 - s/4d0

  ! write(*,*) p1,p2,p3,p4

  utheta = L(1)*p1 + L(2)*p2 + L(3)*p3 + L(4)*p4
  utheta = (utheta*delta)*grs

  u(1) = ur
  u(2) = utheta



end subroutine asymptotic_green
