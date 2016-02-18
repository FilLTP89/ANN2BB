subroutine ALLOCATE_NL(nn,ct,ww,dd,dxdx_el,dxdy_el,dydx_el,dydy_el,     &
        dUxdx_el,dUxdy_el,dUydx_el,dUydy_el,Sxx_el,Syy_el,Sxy_el,Szz_el &
        lambda_el,mu_el,Syld_el,Ckin_el,kkin_el,Riso_el,Rinf_el,biso_el &
        Xkin_el,dEpl_el)
    implicit none

    integer*4, intent(in)                               :: nn
    real*8, intent(inout), dimension(:),    allocatable :: ct,ww
    real*8, intent(inout), dimension(:),    allocatable :: dxdx_el,dydy_el
    real*8, intent(inout), dimension(:),    allocatable :: dxdy_el,dydx_el
    real*8, intent(inout), dimension(:,:),  allocatable :: dd,det_j,fx_el,fy_el
    real*8, intent(inout), dimension(:,:),  allocatable :: dUxdx_el,dUydy_el
    real*8, intent(inout), dimension(:,:),  allocatable :: dUxdy_el,dUydx_el
    real*8, intent(inout), dimension(:,:),  allocatable :: sxx_el,syy_el,sxy_el,szz_el
    real*8, intent(inout), dimension(:,:),  allocatable :: lambda_el,mu_el,syld_el
    real*8, intent(inout), dimension(:,:),  allocatable :: Riso_el,biso_el,Rinf_el
    real*8, intent(inout), dimension(:,:),  allocatable :: Ckin_el,kkin_el
    real*8, intent(inout), dimension(:,:,:),allocatable :: Xkin_el,dEpl_el

    allocate(ct(nn))
    allocate(ww(nn))
    allocate(dd(nn,nn))
    allocate(dxdx_el(nn))
    allocate(dxdy_el(nn))
    allocate(dydx_el(nn))
    allocate(dydy_el(nn))
    allocate(duxdx_el(nn,nn))
    allocate(duxdy_el(nn,nn))
    allocate(duydx_el(nn,nn))
    allocate(duydy_el(nn,nn))
    allocate(sxx_el(nn,nn))
    allocate(syy_el(nn,nn))
    allocate(sxy_el(nn,nn))
    allocate(szz_el(nn,nn))
    allocate(Riso_el(nn,nn))
    allocate(fx_el(nn,nn))
    allocate(fy_el(nn,nn))
    allocate(det_j(nn,nn))
    allocate(mu_el(nn,nn))
    allocate(lambda_el(nn,nn))
    allocate(syld_el(nn,nn))
    allocate(Ckin_el(nn,nn))
    allocate(kkin_el(nn,nn))
    allocate(Rinf_el(nn,nn))
    allocate(biso_el(nn,nn))
    allocate(dEpl_el(4,nn,nn))
    allocate(Xkin_el(4,nn,nn))
    
    return
end subroutine ALLOCATE_NL
