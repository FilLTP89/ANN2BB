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

subroutine DEALLOCATE_NL_EL(ct,ww,dd,dxdx_el,dxdy_el,dydx_el,dydy_el,det_j,   &
    dUxdx_el,dUxdy_el,dUydx_el,dUydy_el,Sxx_el,Syy_el,Sxy_el,Szz_el,        &
    lambda_el,mu_el,syld_el,Ckin_el,kkin_el,Riso_el,Rinf_el,biso_el,        &
    Xkin_el,dEpl_el,fx_el,fy_el,nl_sism,fxs_el,fys_el,Sxxs_el,Syys_el,      &
    Sxys_el,Szzs_el,ux_el,uy_el)
    
    implicit none    
    
    integer*4,  intent(in)                                    :: nl_sism
    real*8,     intent(inout), dimension(:),    allocatable   :: ct,ww
    real*8,     intent(inout), dimension(:),    allocatable   :: dxdx_el,dydy_el
    real*8,     intent(inout), dimension(:),    allocatable   :: dxdy_el,dydx_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: dd,det_j,fx_el,fy_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: fxs_el,fys_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: duxdx_el,duydy_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: duxdy_el,duydx_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: sxx_el,syy_el,sxy_el,szz_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: sxxs_el,syys_el,sxys_el,szzs_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: lambda_el,mu_el,syld_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: Riso_el,biso_el,Rinf_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: Ckin_el,kkin_el
    real*8,     intent(inout), dimension(:,:),  allocatable   :: ux_el,uy_el 
    real*8,     intent(inout), dimension(:,:,:),allocatable   :: Xkin_el,dEpl_el

    if (allocated(ct      )) deallocate(ct      )
    if (allocated(ww      )) deallocate(ww      )
    if (allocated(dd      )) deallocate(dd      )
    if (allocated(dxdx_el )) deallocate(dxdx_el )
    if (allocated(dxdy_el )) deallocate(dxdy_el )
    if (allocated(dydx_el )) deallocate(dydx_el )
    if (allocated(dydy_el )) deallocate(dydy_el )
    if (allocated(duxdx_el)) deallocate(duxdx_el)
    if (allocated(duxdy_el)) deallocate(duxdy_el)
    if (allocated(duydx_el)) deallocate(duydx_el)
    if (allocated(duydy_el)) deallocate(duydy_el)
    if (allocated(sxx_el  )) deallocate(sxx_el  )
    if (allocated(syy_el  )) deallocate(syy_el  )
    if (allocated(szz_el  )) deallocate(szz_el  )
    if (allocated(sxy_el  )) deallocate(sxy_el  )
    if (allocated(Riso_el )) deallocate(Riso_el )
    if (allocated(fx_el   )) deallocate(fx_el   )
    if (allocated(fy_el   )) deallocate(fy_el   )
    if (allocated(det_j   )) deallocate(det_j   )
    if (allocated(mu_el   )) deallocate(mu_el   )
    if (allocated(lambda_el)) deallocate(lambda_el)
    if (allocated(syld_el)) deallocate(syld_el)
    if (allocated(Ckin_el)) deallocate(Ckin_el)
    if (allocated(kkin_el)) deallocate(kkin_el)
    if (allocated(Rinf_el)) deallocate(Rinf_el)
    if (allocated(biso_el)) deallocate(biso_el)
    if (allocated(dEpl_el)) deallocate(dEpl_el)
    if (allocated(Xkin_el)) deallocate(Xkin_el)
    
    if (nl_sism.gt.0.0) then
        deallocate(fxs_el)
        deallocate(fys_el)
        deallocate(sxxs_el)
        deallocate(syys_el)
        deallocate(szzs_el)
        deallocate(sxys_el)
    endif
    return
end subroutine DEALLOCATE_NL_EL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

