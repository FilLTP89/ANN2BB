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

!> @brief Computes the stress tensor.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0

!> @param[in] nn number of 1-D Legendre nodes
!> @param[in] lambda_el nodal values of Lame coefficient lambda 
!> @param[in] mu_el nodal values of Lame coefficient mu
!> @param[in] duxdx nodal values for spatial derivatives of the displacement 
!> @param[in] duydx nodal values for spatial derivatives of the displacement
!> @param[in] duxdy nodal values for spatial derivatives of the displacement
!> @param[in] duydy nodal values for spatial derivatives of the displacement
!> @param[out] sxx nodal values for the stress tensor
!> @param[out] syy nodal values for the stress tensor
!> @param[out] szz nodal values for the stress tensor
!> @param[out] sxy nodal values for the stress tensor

      subroutine MAKE_STRESS(lambda_el,mu_el,nn,&
			     duxdx,duxdy,duydx,duydy,&
			     sxx,syy,szz,sxy)
      
      implicit none
      
      real*8 :: lambda,mu

      integer*4 :: nn
      real*8, dimension(nn,nn) :: duxdx,duxdy,duydx,duydy
      real*8, dimension(nn,nn) :: sxx,syy,szz,sxy

      real*8, dimension(nn,nn) :: mu_el        
      real*8, dimension(nn,nn) :: lambda_el    
      
      integer*4 :: ip,iq
      
                 
      do iq = 1,nn   
         do ip = 1,nn

	    mu =mu_el(ip,iq)
            lambda = lambda_el(ip,iq)
			
            sxx(ip,iq) = (lambda +2.0d0*mu) * duxdx(ip,iq) &
                         + lambda * duydy(ip,iq)
            syy(ip,iq) = (lambda +2.0d0*mu) * duydy(ip,iq) &
                         + lambda * duxdx(ip,iq)
            sxy(ip,iq) = mu * (duxdy(ip,iq) + duydx(ip,iq))
			szz(ip,iq) = lambda * duxdx(ip,iq) + lambda * duydy(ip,iq)
 
         enddo
      enddo
      
      return
      
      end subroutine MAKE_STRESS
