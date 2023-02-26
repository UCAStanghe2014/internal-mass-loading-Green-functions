subroutine asymptotic_love(la,mu,rs,a,HL)
implicit none
real(kind=8),intent(in) :: la,mu,rs,a
real(kind=8),intent(out) :: HL(4,2)

real(kind=8) :: H(4),L(4),pi,b,d

pi = 4d0*atan(1d0)
d = rs/a
b = mu/la

! coefficient of [n,1,1/n,1/n**2]
H(1) = (-1+d**2)/(8d0*a*pi*b*d*la)
H(2) = (-3-d**2+b*(-5-3*d**2))/(a*b*(16*pi*d*la+16*pi*b*d*la))
H(3) = (1+(-1-2*b)*d**2)/(a*(8*pi+b*(16*pi+8*pi*b))*d*la)
H(4) = (-9+3*d**2+b*(-26+18*d**2+b*(-15+29*d**2+b*(-2+10*d**2))))/&
       (a*b*(32*pi*d*la+b*(96*pi*d*la+b*(96*pi*d*la+32*pi*b*d*la))))

! coefficient of [n,1,1/n,1/n**2]
L(1) = 0d0
L(2) = (1-d**2)/(8*a*pi*b*d*la)
L(3) = (-3+3*d**2+b*(-1+5*d**2))/(a*b*(16*pi*d*la+16*pi*b*d*la))
L(4) = (-3-3*d**2+b*(-10-6*d**2+b*(-5-d**2)))/&
       (a*b*(16*pi*d*la+b*(32*pi*d*la+16*pi*b*d*la)))
       
HL(:,1) = H
HL(:,2) = L

end subroutine