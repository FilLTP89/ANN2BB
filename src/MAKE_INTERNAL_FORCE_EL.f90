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

!> @brief Makes internal forces.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights
!> @param[in] dd spectral derivative matrix
!> @param[in] dxdx,dxdy,dydx,dydy jacobian of the bilinear map
!> @param[in] sxx nodal values for the stress tensor
!> @param[in] syy nodal values for the stress tensor
!> @param[in] szz nodal values for the stress tensor
!> @param[in] sxy nodal values for the stress tensor
!> @param[out] fx x-componnent for internal forces
!> @param[out] fy y-componnent for internal forces


subroutine MAKE_INTERNAL_FORCE_EL(nn,ww,dd,dxdx,dxdy,dydx,dydy,&
    sxx,syy,sxy,fx,fy)                               
    !
    implicit none
    ! intent IN
    integer*4, intent(in)                   :: nn
    real*8, intent(in), dimension(nn)       :: ww,dxdx,dxdy,dydx,dydy
    real*8, intent(in), dimension(nn,nn)    :: dd,sxx,syy,sxy
    !intent INOUT
    real*8, intent(inout), dimension(nn,nn) :: fx,fy
    !
    integer*4                               :: ip,iq,il,im,i
    real*8                                  :: t1ux,t1uy,t2ux,t2uy
    real*8                                  :: t1fx,t1fy,t2fx,t2fy
    
    ! FORCE CALCULATION
    do iq = 1,nn
        do ip = 1,nn
            t1fx = 0.0d0
            t1fy = 0.0d0
            t2fx = 0.0d0
            t2fy = 0.0d0

            do il = 1,nn

               t1fx = t1fx + dd(il,ip) * ww(il)*ww(iq) &
                    * (sxx(il,iq)*dydy(il) - sxy(il,iq)*dxdy(il))
               t1fy = t1fy + dd(il,ip) * ww(il)*ww(iq) &
                    * (sxy(il,iq)*dydy(il) - syy(il,iq)*dxdy(il))

            enddo

            do im = 1,nn

               t2fx = t2fx + dd(im,iq) * ww(ip)*ww(im) &
                    * (sxx(ip,im)*dydx(im) - sxy(ip,im)*dxdx(im))
               t2fy = t2fy + dd(im,iq) * ww(ip)*ww(im) &
                    * (sxy(ip,im)*dydx(im) - syy(ip,im)*dxdx(im))
            
            enddo

            fx(ip,iq) = t1fx - t2fx
            fy(ip,iq) = t1fy - t2fy
        
        enddo
    enddo
    !
    return
end subroutine MAKE_INTERNAL_FORCE_EL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

