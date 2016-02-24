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

!> @brief Computes spatial derivatives of displacement.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nn number of 1-D Legendre nodes
!> @param[in] dd matrix of spectral derivatives
!> @param[in] dxdx,dxdy,dydx,dydy component of the jacobian for coor. transf.
!> @param[in] ux x-displacement
!> @param[in] uy y-displacement
!> @param[out] duxdx nodal values for spatial derivatives of the displacement 
!> @param[out] duydx nodal values for spatial derivatives of the displacement
!> @param[out] duxdy nodal values for spatial derivatives of the displacement
!> @param[out] duydy nodal values for spatial derivatives of the displacement


subroutine MAKE_STRAIN(nn,dd,dxdx,dxdy,dydx,dydy,&
    ux,uy,duxdx,duxdy,duydx,duydy)
    real*8                                  :: t1ux,t1uy,t2ux,t2uy
    real*8                                  :: t1fx,t1fy,t2fx,t2fy,det_j
    integer*4                               :: ip,iq,il,im
    integer*4,               intent(in)     :: nn
    real*8, dimension(nn),   intent(in)     :: dxdx,dxdy,dydx,dydy
    real*8, dimension(nn,nn),intent(in)     :: dd,ux,uy
    real*8, dimension(nn,nn),intent(inout)  :: duxdx,duxdy,duydx,duydy

    !   DERIVATIVE CALCULATION

    do iq = 1,nn
        do ip = 1,nn
            t1ux = 0.d0; t1uy = 0.d0
            t2ux = 0.d0; t2uy = 0.d0
            det_j = dxdx(iq)*dydy(ip) - dxdy(ip)*dydx(iq)
            do il = 1,nn
               t1ux = t1ux + ux(il,iq) * dd(ip,il)
               t1uy = t1uy + uy(il,iq) * dd(ip,il)
            enddo

            do im = 1,nn
               t2ux = t2ux + ux(ip,im) * dd(iq,im)
               t2uy = t2uy + uy(ip,im) * dd(iq,im)
            enddo
            duxdx(ip,iq) = (1.0d0 / det_j)*((dydy(ip) * t1ux) - (dydx(iq) * t2ux))
            duydx(ip,iq) = (1.0d0 / det_j)*((dydy(ip) * t1uy) - (dydx(iq) * t2uy))
            duxdy(ip,iq) = (-1.0d0 / det_j)*((dxdy(ip) * t1ux) - (dxdx(iq) * t2ux))
            duydy(ip,iq) = (-1.0d0 / det_j)*((dxdy(ip) * t1uy) - (dxdx(iq) * t2uy))
        enddo
    enddo
    return
end subroutine MAKE_STRAIN
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

