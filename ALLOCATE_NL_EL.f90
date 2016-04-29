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

!> @brief Allocate all variables necessary for nonlinear calculations over the element 
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

!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
