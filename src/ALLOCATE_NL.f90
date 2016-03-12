!    Copyright (C) 2014 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with SPEED.  If not, see <http://wqw.gnu.org/licenses/>.

!> @brief Allocate all variables necessary for nonlinear calculations 
!! @author Filippo Gatti
!> @date February,2016
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] ct LGL nodes
!> @param[in] ww LGL weights
!> @param[in] dd matrix of spectral derivatives
!> @param[in] dxdx_el,dxdy_el,dydx_el,dydy_el coordinates spectral derivatives on element LGL
!> @param[in] Ux_el,Uy_el
!> @param[in] dUxdx,dUxdy,dUydx,dUydy 
!> @param[inout] sxx,syy,sxy,szz nodal values for the stress tensor (on element)
!> @param[inout] Xkin_el,Riso_el hardening variables on element LGL
!> @param[in] lambda_el,mu_el Lamè parameters on element LGL
!> @param[in] syld_el yield limit on element LGL
!> @param[in] Ckin_el,kkin_el kinematic hardening parameters on element LGL
!> @param[in] Rinf_el,biso_el isotropic hardening parameters on element LGL
!> @param[inout] dEpl_el plastic strain increment on element LGL
!> @param[inout] det_j determinant of jacobian
!> @param[inout] fx_el internal forces along x-direction on element LGL
!> @param[inout] fy_el internal forces along y-direction on element LGL
!> @param[in] nl_sism flag for seismic moment
!> @param[inout] fxs_el seismic moment equivalent forces in x-direction on element LGL
!> @param[inout] fys_el seismic moment equivalent forces in x-direction on element LGL

subroutine ALLOCATE_NL(nn,ct,ww,dd,ux_el,uy_el,      &
         dxdx_el, dxdy_el, dydx_el, dydy_el,det_j,   &
        duxdx_el,duxdy_el,duydx_el,duydy_el,         &
          sxx_el,  sxy_el,  syy_el,szz_el,           &
       lambda_el,   mu_el, syld_el, Ckin_el, kkin_el,&
         Rinf_el, biso_el, Riso_el, Xkin_el, dEpl_el,&
           fx_el,   fy_el,  fxs_el,  fys_el, nl_sism,&
         sxxs_el, sxys_el, syys_el,szzs_el)
    
    integer*4,                              intent(in ) :: nn, nl_sism
    real*8,     dimension(:),  allocatable, intent(out) :: ct,ww
    real*8,     dimension(:),  allocatable, intent(out) :: dxdx_el,dydy_el
    real*8,     dimension(:),  allocatable, intent(out) :: dxdy_el,dydx_el
    real*8,     dimension(:,:),allocatable, intent(out) :: dd,det_j,fx_el,fy_el
    real*8,     dimension(:,:),allocatable, intent(out) :: fxs_el,fys_el
    real*8,     dimension(:,:),allocatable, intent(out) :: ux_el,uy_el
    real*8,     dimension(:,:),allocatable, intent(out) :: duxdx_el,duydy_el
    real*8,     dimension(:,:),allocatable, intent(out) :: duxdy_el,duydx_el
    real*8,     dimension(:,:),allocatable, intent(out) :: sxx_el,syy_el,sxy_el,szz_el
    real*8,     dimension(:,:),allocatable, intent(out) :: sxxs_el,syys_el,sxys_el,szzs_el
    real*8,     dimension(:,:),allocatable, intent(out) :: lambda_el,mu_el,syld_el
    real*8,     dimension(:,:),allocatable, intent(out) :: Riso_el,biso_el,Rinf_el
    real*8,     dimension(:,:),allocatable, intent(out) :: Ckin_el,kkin_el
    real*8,     dimension(:,:,:),allocatable, intent(out) :: Xkin_el,dEpl_el
    
    allocate(ct(nn))
    allocate(ww(nn))
    allocate(dd(nn,nn))
    allocate(dxdx_el(nn))
    allocate(dxdy_el(nn))
    allocate(dydx_el(nn))
    allocate(dydy_el(nn))
    allocate(ux_el(nn,nn))
    allocate(uy_el(nn,nn))
    allocate(duxdx_el(nn,nn))
    allocate(duxdy_el(nn,nn))
    allocate(duydx_el(nn,nn))
    allocate(duydy_el(nn,nn))
    allocate(sxx_el(nn,nn))
    allocate(syy_el(nn,nn))
    allocate(szz_el(nn,nn))
    allocate(sxy_el(nn,nn))
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
    
    if (nl_sism.gt.0) then
        allocate(fxs_el(nn,nn))
        allocate(fys_el(nn,nn))
        allocate(sxxs_el(nn,nn))
        allocate(syys_el(nn,nn))
        allocate(szzs_el(nn,nn))
        allocate(sxys_el(nn,nn))
    endif
    return
end subroutine ALLOCATE_NL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

