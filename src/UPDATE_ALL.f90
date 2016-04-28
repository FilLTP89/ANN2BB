!    Copyright (C) 2014 The SPEED FOuNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNu Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOuT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICuLAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNu Affero General Public License
!    along with SPEED.  If not, see <http://wqw.gnu.org/licenses/>.

!> @brief Update stress and hardening variables from nonlinear calculations.
!! @author Filippo Gatti
!> @date February,2016
!> @version 1.0
subroutine UPDATE_ALL(nnt,sxx,syy,szz,sxy,duxdx,duxdy,duydx,duydy,xkin,riso,epl,nodal_counter)
    
    implicit none
    integer*4, intent(in)                               :: nnt
    integer*4, dimension(nnt), intent(in)               :: nodal_counter
    real*8, dimension(nnt),   intent(inout)             :: sxx,syy,szz,sxy
    real*8, dimension(nnt),   intent(inout)             :: duxdx,duydy,duxdy,duydx
    real*8, dimension(nnt),   intent(inout)             :: riso
    real*8, dimension(4*nnt), intent(inout)             :: xkin
    real*8, dimension(4*nnt), intent(inout)             :: epl
    sxx   = sxx
    syy   = syy
    szz   = szz
    sxy   = sxy
    duxdx = duxdx
    duxdy = duxdy
    duydx = duydx
    duydy = duydy
    xkin  = xkin
    riso  = riso
    epl   = epl

    return
end subroutine UPDATE_ALL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

