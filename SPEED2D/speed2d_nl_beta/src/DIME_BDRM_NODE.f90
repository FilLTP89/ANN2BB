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

!> @brief Extracts number of DRM internal boundary nodes without duplicates
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] ns
!> @param[in] line
!> @param[in] con_line connectivity for line elements
!> @param[in] cs_nnz_bc legth of cs_bc 
!> @param[in] cs_bc spectral connectivity vector for boundary elements
!> @param[in] nBDRM
!> @param[in] tag_BDRM labels for BDRM blocks
!> @param[in] n_boun_DRM number of elements of DRM boundary domain
!> @param[in] ii_node_BDRM number of nodes of DRM boundary domain
!> @param[out] node_BDRM node's numbers of DRM internal boundary
!> @param[out] nnode_BD size of boundary DRM nodes without duplicates


      subroutine DIME_BDRM_NODE(ns,line,con_line,cs_nnz_bc,cs_bc, &
	                           nBDRM,tag_BDRM,n_boun_DRM, &
						       node_BDRM,ii_node_BDRM,nnode_BD)
      
!     © POLIMI, 2006, All Rights Reserved
!     Author: Laura Scandella

      implicit none

      integer*4 :: ns
     
      integer*4 :: cs_nnz_bc,line,n_boun_DRM
      integer*4, dimension(0:cs_nnz_bc) :: cs_bc
      integer*4, dimension(line,*) :: con_line

      integer*4 :: nBDRM
      integer*4, dimension(nBDRM) :: tag_BDRM

      integer*4 :: nn
	  
	  integer*4 :: ii_node_BDRM
	  integer*4, dimension(ii_node_BDRM) :: node_BDRM
 
      integer*4 :: nnode_BD

	  integer*4 :: imBDRM,il,i,j,cont,cont1,min,nd,k   


!     Node's numbers of DRM internal boundary (node_BDRM) 
	  nn=ns+1 
      cont=1
	  cont1=nn

	   do imBDRM = 1,nBDRM
	    do il =1,line
		   if (cs_bc(cs_bc(il-1)+0).eq.tag_BDRM(imBDRM)) then 
               node_BDRM(cont:cont1) = cs_bc((cs_bc(il-1)+1):(cs_bc(il-1)+nn))
			   cont=cont1+1
			   cont1=cont1+nn
		   endif
		enddo
	  enddo

!    Check for duplicates 
     do j = 1,ii_node_BDRM-1
	    min = node_BDRM(j)
	    do i = ii_node_BDRM,j+1,-1
	        if  (node_BDRM(i).lt.min) then
		        nd = node_BDRM(i)
			    node_BDRM(i) = min
			    min = nd
		    endif
	    enddo
        node_BDRM(j) = min
	 enddo
!    Size of boundary DRM nodes without duplicates
	 k = 0
	 do i = 1,ii_node_BDRM-1
		if (node_BDRM(i+1).ne.node_BDRM(i)) then
		k = k+1
		endif
	 enddo
	 k = k+1 !In order to count the last node of the list

	 nnode_BD = k

	 return

     end subroutine DIME_BDRM_NODE