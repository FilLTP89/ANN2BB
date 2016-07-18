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

!> @brief Extracts number of DRM boundary lines and nodes
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] ns
!> @param[in] con_line connectivity for line elements
!> @param[in] line
!> @param[in] cs_nnz_bc legth of cs_bc 
!> @param[in] cs_bc spectral connectivity vector for boundary elements
!> @param[in] nBDRM
!> @param[in] tag_BDRM labels for BDRM blocks
!> @param[out] n_boun_DRM number of elements of DRM boundary domain
!> @param[out] ii_node_BDRM number of nodes of DRM boundary domain


      subroutine DIME_DRM_BOUN(ns,con_line,line,cs_nnz_bc,cs_bc,&
	                      nBDRM,tag_BDRM,n_boun_DRM,ii_node_BDRM)
      
      implicit none

      integer*4 :: ns
      
      integer*4 :: line,cs_nnz_bc
      integer*4, dimension(line,*) :: con_line
      integer*4, dimension(0:cs_nnz_bc) :: cs_bc

      integer*4 :: nBDRM
      integer*4, dimension(nBDRM) :: tag_BDRM

	  integer*4 :: n_boun_DRM,ii_node_BDRM

	  integer*4 :: imBDRM,il

!     Number of elements of DRM domain

      n_boun_DRM=0
	  
	  do imBDRM = 1,nBDRM
	    do il =1,line
		   if (con_line(il,1).eq.tag_BDRM(imBDRM)) then
           n_boun_DRM = n_boun_DRM+1
		   endif
		enddo
	  enddo

!     Number of nodes of DRM domain (with duplications)

	  ii_node_BDRM = n_boun_DRM*(ns+1)

      return

      end subroutine DIME_DRM_BOUN      