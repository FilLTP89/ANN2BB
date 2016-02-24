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

!> @brief Deallocate all variables necessary for nonlinear calculations 
!! @author Filippo Gatti
!> @date February,2016
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] ct LGL nodes
!> @param[in] ww LGL weights
!> @param[in] dd matrix of spectral derivatives
!> @param[in] dxdx_el,dxdy_el,dydx_el,dydy_el coordinates spectral derivatives on element LGL
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

subroutine DEALLOCATE_NL(nn,ct,ww,dd,dxdx_el,dxdy_el,dydx_el,dydy_el,det_j,   &
    dUxdx_el,dUxdy_el,dUydx_el,dUydy_el,Sxx_el,Syy_el,Sxy_el,Szz_el,        &
    lambda_el,mu_el,Syld_el,Ckin_el,kkin_el,Riso_el,Rinf_el,biso_el,        &
    Xkin_el,dEpl_el,fx_el,fy_el,nl_sism,fxs_el,fys_el,Sxxs_el,Syys_el,      &
    Sxys_el,Szzs_el)
    
    implicit none

    integer*4,  intent(in)                             :: nn,nl_sism
    real*8,     intent(inout), dimension(:),    allocatable :: ct,ww
    real*8,     intent(inout), dimension(:),    allocatable :: dxdx_el,dydy_el
    real*8,     intent(inout), dimension(:),    allocatable :: dxdy_el,dydx_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: dd,det_j,fx_el,fy_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: fxs_el,fys_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: dUxdx_el,dUydy_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: dUxdy_el,dUydx_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: sxx_el,syy_el,sxy_el,szz_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: sxxs_el,syys_el,sxys_el,szzs_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: lambda_el,mu_el,syld_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: Riso_el,biso_el,Rinf_el
    real*8,     intent(inout), dimension(:,:),  allocatable :: Ckin_el,kkin_el
    real*8,     intent(inout), dimension(:,:,:),allocatable :: Xkin_el,dEpl_el
    
    deallocate(ct)
    deallocate(ww)
    deallocate(dd)
    deallocate(dxdx_el)
    deallocate(dxdy_el)
    deallocate(dydx_el)
    deallocate(dydy_el)
    deallocate(duxdx_el)
    deallocate(duxdy_el)
    deallocate(duydx_el)
    deallocate(duydy_el)
    deallocate(sxx_el)
    deallocate(syy_el)
    deallocate(szz_el)
    deallocate(sxy_el)
    deallocate(Riso_el)
    deallocate(fx_el)
    deallocate(fy_el)
    deallocate(det_j)
    deallocate(mu_el)
    deallocate(lambda_el)
    deallocate(syld_el)
    deallocate(Ckin_el)
    deallocate(kkin_el)
    deallocate(Rinf_el)
    deallocate(biso_el)
    deallocate(dEpl_el)
    deallocate(Xkin_el)
    
    if (nl_sism.gt.0) then
        deallocate(fxs_el)
        deallocate(fys_el)
        deallocate(sxxs_el)
        deallocate(syys_el)
        deallocate(szzs_el)
        deallocate(sxys_el)
    endif
    return
end subroutine DEALLOCATE_NL
