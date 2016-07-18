!    Copyright (C) 2014 The SPEED FOUNDATION
!    Author: 
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

!> @brief Reads dimensions of DRM elements 
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] ns
!> @param[in] nquad number of quad elements
!> @param[in] con_quad connectivity for quad elements
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] nMDRM
!> @param[in] tag_MDRM labels for DRM blocks
!> @param[out] n_el_DRM number of DRM elements
!> @param[out] ii_node_DRM number of nodes of DRM elements


      subroutine DIME_DRM_EL(ns,nquad,con_quad,cs_nnz,cs, &
	                      nMDRM,tag_MDRM,n_el_DRM,ii_node_DRM)
      
      implicit none

	  integer*4 :: ns

      integer*4 :: cs_nnz,nquad
      integer*4, dimension(0:cs_nnz) :: cs
	  integer*4, dimension(nquad,*) :: con_quad

      integer*4 :: nMDRM
      integer*4, dimension(nMDRM) :: tag_MDRM

      integer*4 :: imMDRM,ie

	  integer*4 :: n_el_DRM,cont
	  integer*4 ::ii_node_DRM


!     Number of elements of DRM domain
      n_el_DRM=0
	  
	  do imMDRM = 1,nMDRM
	    do ie =1,nquad
		   if (con_quad(ie,1).eq.tag_MDRM(imMDRM)) then
           n_el_DRM = n_el_DRM+1
		   endif
		enddo
	  enddo

!     Number of nodes of DRM domain (with duplications)
	  ii_node_DRM = n_el_DRM*((ns+1)*(ns+1))

      return

      end subroutine DIME_DRM_EL 
