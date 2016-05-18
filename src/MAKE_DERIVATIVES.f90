!****************************************************************************
! MAKE DERIVATIVES
!****************************************************************************

subroutine MAKE_DERIVATIVES(alfa1,alfa2,beta1,beta2,gamma1,gamma2,nn,ct,dxdx,dxdy,dydx,dydy)
    !           
    implicit none
    ! intent IN
    integer*4, intent(in)               :: nn
    real*8, intent(in)                  :: alfa1,alfa2
    real*8, intent(in)                  :: beta1,beta2
    real*8, intent(in)                  :: gamma1,gamma2
    real*8, intent(in)   , dimension(nn):: ct
    ! intent INOUT
    real*8, intent(inout), dimension(nn):: dxdx,dxdy
    real*8, intent(inout), dimension(nn):: dydx,dydy
    ! 
    integer*4 :: i
    ! COMPUTE SHAPE-FUNCTION DERIVATIVES
    dxdx = 0.d0
    dxdy = 0.d0
    dydx = 0.d0
    dydy = 0.d0
    do i = 1,nn
        dxdy(i) = beta1 + gamma1 * ct(i)
        dydy(i) = beta2 + gamma2 * ct(i)
        dxdx(i) = alfa1 + gamma1 * ct(i)
        dydx(i) = alfa2 + gamma2 * ct(i)
    enddo
    ! 
    return
end subroutine MAKE_DERIVATIVES
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

