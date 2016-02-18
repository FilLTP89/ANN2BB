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

!> @brief Makes spectral connectivity vector for boundary elements.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] cs_nnz    length of cs
!> @param[in] cs        spectral connectivity vector
!> @param[in] ne_bc     number of boundary elements
!> @param[in] cm_bc     connectivity matrix for boundary elements
!> @param[in] nm        number of materials
!> @param[in] tm        material label
!> @param[in] sd        polynomial degree vector 
!> @param[in] Ennz      nnzero elements for Ebin vector
!> @param[in] Ebin      pointer for connectivity (see MAKE_GRID_NODES.f90) 
!> @param[in] cs_nnz_bc legth of cs_bc 
!> @param[out] cs_bc    spectral connectivity vector for boundary elements


      subroutine MAKE_SPECTRAL_BOUNDARY(cs_nnz,cs,ne_bc,cm_bc,nm,tm,sd,&
                                        Ennz,Ebin,cs_nnz_bc,cs_bc)
            
      implicit none
      
      integer*4 :: cs_nnz,ne_bc,nm,Ennz,cs_nnz_bc
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(ne_bc,*) :: cm_bc
      integer*4, dimension(nm) :: tm 
      integer*4, dimension(nm) :: sd 
      integer*4, dimension(0:Ennz) :: Ebin
      integer*4, dimension(0:cs_nnz_bc) :: cs_bc
      
      integer*4 :: im,ie,ielem,i,j,k,nn
      integer*4 :: an,bn,n1,n2,n3,n4
      
      
      cs_bc(0) = ne_bc +1
      
      do ie = 1,ne_bc
         an = cm_bc(ie,1 +1)
         bn = cm_bc(ie,2 +1)
         
         call get_edge_element(Ennz,Ebin,&
                               cm_bc(ie,2),cm_bc(ie,3),ielem)
         do im = 1,nm
            if (tm(im).eq.cs(cs(ielem -1) +0)) nn = sd(im) +1
         enddo
         cs_bc(ie) = cs_bc(ie -1) + nn +1
      enddo
      
      
      do ie = 1,ne_bc
         an = cm_bc(ie,1 +1)
         bn = cm_bc(ie,2 +1)
         
         cs_bc(cs_bc(ie -1) +0) = cm_bc(ie,1)
         
         call get_edge_element(Ennz,Ebin,&
                               cm_bc(ie,2),cm_bc(ie,3),ielem)
         do im = 1,nm
            if (tm(im).eq.cs(cs(ielem -1) +0)) nn = sd(im) +1
         enddo
         
         n1 = cs(cs(ielem -1) +1)
         n2 = cs(cs(ielem -1) +nn)
         n3 = cs(cs(ielem -1) +nn*nn)
         n4 = cs(cs(ielem -1) +nn*(nn -1) +1)
         
         if (((an.eq.n1).and.(bn.eq.n2)) &
              .or.((bn.eq.n1).and.(an.eq.n2))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +k)
            enddo
            
         else if (((an.eq.n2).and.(bn.eq.n3)) &
              .or.((bn.eq.n2).and.(an.eq.n3))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +nn*k)
            enddo
            
         else if (((an.eq.n3).and.(bn.eq.n4)) &
              .or.((bn.eq.n3).and.(an.eq.n4))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +nn*(nn -1) +k)
            enddo
            
         else if (((an.eq.n4).and.(bn.eq.n1)) &
              .or.((bn.eq.n4).and.(an.eq.n1))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +nn*(k -1) +1)
            enddo
            
         endif
         
      enddo
      
      return
      
      end subroutine MAKE_SPECTRAL_BOUNDARY
