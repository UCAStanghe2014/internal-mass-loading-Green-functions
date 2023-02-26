subroutine z2y_unit_mass_loading(YLOVE,YLOVET1,YLOVET2)
  !// convert z-numbes in SUN's code to y-numbers
  !// See Eq.(18) 孙文科, & 李瑞浩. (1986). 弹性地球的静态负荷响应. 大地测量与地球动力学, 6(4), 11.
  !// All constants and variables use the International System of Units
  implicit none
  real(kind=8),intent(inout) :: YLOVE(6)        !// z-numbers and y-numbers of Spheroid components
  real(kind=8),intent(inout) :: YLOVET1,YLOVET2 !// z-numbers and y-numbers of Toroidal components
  
  real(kind=8) :: factor
  
  factor = 6.67D-11/((6371d3**2)*9.8156d0)

  !// z-numbers to y-numbers of unit mass (1kg) internal mass loading by a factor 

  YLOVE(1) = (YLOVE(1)*factor) * 6371d3
  YLOVE(2) = (YLOVE(2)*factor) * 1308d9
  YLOVE(3) = (YLOVE(3)*factor) * 6371d3
  YLOVE(4) = (YLOVE(4)*factor) * 1308d9
  YLOVE(5) = (YLOVE(5)*factor) * 6371d3*9.8156d0 ! 9.8156 means surface gravity
  YLOVE(6) = (YLOVE(6)*factor) * 9.8156d0
  YLOVET1  = (YLOVET1 *factor) * 6371d3
  YLOVET2  = (YLOVET2 *factor) * 1308d9
 
end subroutine z2y_unit_mass_loading
