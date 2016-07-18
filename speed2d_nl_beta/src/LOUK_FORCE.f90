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

!> @brief Calculates the Loukadidis forces px and py at DRM nodes
!! @author
!> @date
!> @version

!> @param[in] ns
!> @param[in] K_el stiffness matrix for the current DRM element
!> @param[in] nnode_BD number of DRM internal boundary nodes
!> @param[in] nnode_TOT number of DRM total nodes without duplicate
!> @param[in] node_PDRM
!> @param[in] node_TOT nodes of total DRM domain without duplicate
!> @param[in] nf_drm number of DRM time functions
!> @param[in] disp_t
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] ie
!> @param[out] px 
!> @param[out] py

      subroutine LOUK_FORCE(ns,K_el,nnode_BD,nnode_TOT,node_PDRM, &
						    node_TOT,nf_drm,disp_t,cs_nnz,cs,ie,px,py)
      
!     © POLIMI, 2006, All Rights Reserved
!     Author: Laura Scandella

     
      implicit none

	  integer*4 :: nnode_BD,nf,ie, nf_drm, ns

	  integer*4 :: cs_nnz,nnode_TOT
      integer*4, dimension(0:cs_nnz) :: cs
	  integer*4, dimension(nnode_TOT) :: node_TOT
      integer*4, dimension(nnode_TOT) :: node_PDRM
      real*8, dimension(nnode_TOT) ::px,py

	  real*8, dimension(nf_drm,2) :: disp_t
	  
	  real*8, dimension(2*(ns+1)*(ns+1),2*(ns+1)*(ns+1)) :: K_el

	  integer*4 ::nn,is,ip,n1,n2,n3,n4,g

	  integer*4 :: i,j,k,f,p,n,ii,check,imDRM,s,or

	  nn = ns+1


      do i=1,2
	     if(i.eq.1) then
		   n1 = 1
		   n2 = nnode_BD
		   n3 = nnode_BD + 1
		   n4 = nnode_TOT
		   g = +1
		 else
		   n1 = nnode_BD + 1
		   n2 = nnode_TOT
		   n3 = 1
		   n4 = nnode_BD
		   g = -1
		 endif

         do j = 1,((ns+1)**2)
		    check = 0
		    do n = n1,n2

			   if (cs(cs(ie-1)+j).eq.node_TOT(n)) then
                  or =0
			      do s = 1,nnode_TOT 
				    if (cs(cs(ie-1)+j).eq.node_PDRM(s)) then
					   or = s
					   exit
					endif
				  enddo	
			      check = 1		         
			      do k = 1,((ns+1)**2)
				     do ii = n3,n4
					    if(cs(cs(ie-1)+k).eq.node_TOT(ii)) then
						        px(ii) = px(ii) + (K_el(2*k-1,2*j-1)*disp_t(or,1)+ &
								         K_el(2*k-1,2*j)*disp_t(or,2))*g
						        py(ii) = py(ii) + (K_el(2*k,2*j-1)*disp_t(or,1)+ &
								         K_el(2*k,2*j)*disp_t(or,2))*g

						     exit  
						 endif
					 enddo
			      enddo
			   endif
			   if (check.eq.1) exit
		   	enddo
		 enddo
	  enddo  

	  return
	  									       
      end subroutine LOUK_FORCE    