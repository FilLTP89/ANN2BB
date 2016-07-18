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

!> @brief Make damping matrices.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nnode number of nodes
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vecotr
!> @param[in] nm number of materials
!> @param[in] tag_mat label of materials
!> @param[in] sdeg_mat polynomial degree
!> @param[in] prop_mat material properties
!> @param[in] ne number of elements
!> @param[in] alfa1 costant values for the bilinear map
!> @param[in] alfa2 costant values for the bilinear map
!> @param[in] beta1 costant values for the bilinear map 
!> @param[in] beta1 costant values for the bilinear map 
!> @param[in] gamma1 costant values for the bilinear map
!> @param[in] gamma2 costant values for the bilinear map
!> @param[in] delta1 costant values for the bilinear map
!> @param[in] delta2 costant values for the bilinear map
!> @param[out] cvec viscous forces proportional to the velocity field
!> @param[out] kcvec viscous forces proportional to the displacement field

      subroutine MAKE_CEL_KCEL(nnode,cs_nnz,cs,&
                               nm,tag_mat,sdeg_mat,prop_mat,&
                               ne,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                               cvec,kcvec)
      
      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,sdeg_mat
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(ne) :: alfa1,beta1,gamma1
      real*8, dimension(ne) :: alfa2,beta2,gamma2
      real*8, dimension(2*nnode) :: cvec
      real*8, dimension(2*nnode) :: kcvec
      
      real*8, dimension(:), allocatable :: dxdy,dydy,dxdx,dydx
      real*8, dimension(:,:), allocatable :: det_j
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8 :: rho,term_Cel,term_KCel,gamma
      integer*4 :: im,ie,i,j,is,in,id1,id2,nn,nn2
      
      
      nn = 2
      
      cvec = 0.d0;   kcvec = 0.d0
      ne = cs(0) -1
      
      do ie = 1,ne

         im = cs(cs(ie -1) +0);
         rho = prop_mat(im,1)
         gamma = prop_mat(im,4)
         nn = sdeg_mat(im) +1
  
         allocate(ct(nn),ww(nn),dd(nn,nn))
         call lgl(nn,ct,ww,dd)
         allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))

         do i = 1,nn
              dxdy(i) = beta1(ie) + gamma1(ie) * ct(i)
              dydy(i) = beta2(ie) + gamma2(ie) * ct(i)
         enddo
               
         do j = 1,nn
              dxdx(j) = alfa1(ie) + gamma1(ie) * ct(j)
              dydx(j) = alfa2(ie) + gamma2(ie) * ct(j)
         enddo
               
         do j = 1,nn
              do i = 1,nn
                 det_j(i,j) = dxdx(j)*dydy(i) - dxdy(i)*dydx(j)
              enddo
         enddo
               
         do j = 1,nn
            do i = 1,nn
               is = nn*(j -1) +i
               in = cs(cs(ie -1) +is)
                     
               id1 = in
               id2 = in + nnode
                        
               term_Cel = 2 * gamma * rho * det_j(i,j) * ww(i) * ww(j)
               cvec(id1) = cvec(id1) + term_Cel
               cvec(id2) = cvec(id2) + term_Cel
                                                
               term_KCel = (gamma**2) * rho * det_j(i,j) * ww(i) * ww(j)
               kcvec(id1) = kcvec(id1) + term_KCel
               kcvec(id2) = kcvec(id2) + term_KCel
                        
            enddo
         enddo
 
         deallocate(ct,ww,dd)
         deallocate(dxdx,dydx,dxdy,dydy,det_j)
        
      enddo
            
      return
      end subroutine MAKE_CEL_KCEL
