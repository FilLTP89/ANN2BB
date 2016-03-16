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
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.

!> @brief Perform coordinates' jacobian derivatives.
!! @author Ilario Mazzieri and Filippo Gatti
!> @date February, 2016
!> @version 1.0

!> @param[in] nn number of nodes on element
!> @param[in] alfa,beta,gamma coefficients of bilinear map
!> @param[inout] dxdx_el,dydy_el,dxdy_el,dydx_el derivatives of coordinate functions

subroutine MAKE_DERIVATIVES(nn,alfa1,alfa2,beta1,beta2,gamma1,gamma2,ct,&
    dxdy_el,dydy_el,dxdx_el,dydx_el)
                
    implicit none
    integer*4, intent(in)               :: nn
    real*8, intent(in)                  :: alfa1,alfa2
    real*8, intent(in)                  :: beta1,beta2
    real*8, intent(in)                  :: gamma1,gamma2
    real*8, intent(in)   , dimension(nn):: ct
    real*8, intent(inout), dimension(nn):: dxdx_el,dxdy_el
    real*8, intent(inout), dimension(nn):: dydx_el,dydy_el
    
    integer*4 :: i

    do i = 1,nn
        dxdy_el(i) = beta1 + gamma1 * ct(i)
        dydy_el(i) = beta2 + gamma2 * ct(i)
        dxdx_el(i) = alfa1 + gamma1 * ct(i)
        dydx_el(i) = alfa2 + gamma2 * ct(i)
    enddo
    
    return
end subroutine MAKE_DERIVATIVES
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
