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

!> @brief Computes ABC's. 
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] rho,lambda,mu material properties given node by node
!> @param[in] nn number of 1D GLL nodes
!> @param[in] ct GLL nodes
!> @param[in] wq GLL weights
!> @param[in] dd spectral derivative matrix
!> @param[in] dxdx,dxdy,dydx,dydy Jacobian of the coordinate transformation
!> @param[in] nnx x-component of the normal vector
!> @param[in] nny y-component of the normal vector
!> @param[in] ll lenght of the boundary edge
!> @param[in] ia index for selecting the element absorbing face
!> @param[in] ib index for selecting the element absorbing face
!> @param[in] ja index for selecting the element absorbing face
!> @param[in] jb index for selecting the element absorbing face
!> @param[in] ux x-displacement
!> @param[in] uy y-displacement
!> @param[in] vx x-velocity
!> @param[in] vy y-velocity
!> @param[out] fx x-absorbing forces
!> @param[out] fy y-absorbing forces


      subroutine MAKE_ABSORBING_FORCE(lambda,mu,rho,nn,ct,ww,dd,&
                                      dxdx,dxdy,dydx,dydy,&
                                      nnx,nny,ll,ia,ja,ib,jb,&
                                      ux,uy,vx,vy,fx,fy)
      
      implicit none
      
      real*8 :: lambda,mu,rho
      integer*4 :: nn
      real*8, dimension(nn) :: ct,ww
      real*8, dimension(nn,nn) :: dd
      real*8, dimension(nn) :: dxdx,dxdy,dydx,dydy
      real*8 :: nnx,nny,ll
      integer*4 :: ia,ja,ib,jb
      real*8, dimension(nn,nn) :: ux,uy,vx,vy
      real*8, dimension(nn,nn) :: fx,fy
      
      integer*4 :: ip,iq,il,im,i,j
      real*8 :: det_j,nnx2,nny2,term,term11,term12,term21,term22,term_k
      
      do j = 1,nn
         do i = 1,nn
            fx(i,j) = 0.0d0
            fy(i,j) = 0.0d0
         enddo
      enddo
      
      nnx2 = nnx * nnx
      nny2 = nny * nny
      
      term = 0.5d0*ll * (2.0d0*mu - dsqrt(mu * (lambda+2.0d0*mu)))
      
      term11 = nnx2 * dsqrt(lambda +2.0d0*mu) + nny2 * dsqrt(mu)
      term11 = -0.5d0 * ll * dsqrt(rho) * term11
      
      term12 = (dsqrt(mu) - dsqrt(lambda +2.0d0*mu)) * nnx * nny
      term12 = 0.5d0 * ll * dsqrt(rho) * term12
      
      term21 = nny2 * dsqrt(lambda +2.0d0*mu) + nnx2 * dsqrt(mu)
      term21 = -0.5d0 * ll * dsqrt(rho) * term21
      
      term22 = (dsqrt(mu) - dsqrt(lambda +2.0d0*mu)) * nnx * nny
      term22 = 0.5d0 * ll * dsqrt(rho) * term22
      
      j = 0
      
      do il = ia,ib
         do im = ja,jb
            j = j +1
            
            det_j = dxdx(im)*dydy(il) - dxdy(il)*dydx(im)
            
            term_k = term * (dydy(il)*nny + dxdy(il)*nnx) &
                   * ww(j) / det_j
            
            iq = im
            do ip = 1,nn
               fx(il,im) = fx(il,im) -(term_k * dd(il,ip) * uy(ip,iq))
               fy(il,im) = fy(il,im) +(term_k * dd(il,ip) * ux(ip,iq))
            enddo
            
            term_k = -1.0d0 * term * (dydx(im)*nny + dxdx(im)*nnx) &
                    * ww(j) / det_j
            
            ip = il
            do iq = 1,nn
               fx(il,im) = fx(il,im) -(term_k * dd(im,iq) * uy(ip,iq))
               fy(il,im) = fy(il,im) +(term_k * dd(im,iq) * ux(ip,iq))
            enddo
            
            fx(il,im) = fx(il,im) -term11*ww(j) *vx(il,im)
            fx(il,im) = fx(il,im) -term12*ww(j) *vy(il,im)
            
            fy(il,im) = fy(il,im) -term21*ww(j) *vy(il,im)
            fy(il,im) = fy(il,im) -term22*ww(j) *vx(il,im)
            
         enddo
      enddo
      
      return
      
      end subroutine MAKE_ABSORBING_FORCE
