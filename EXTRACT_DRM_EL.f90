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

!> @brief Extracts DRM elements and their nodes 
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
!> @param[in] n_el_DRM number of DRM elements
!> @param[in] ii_node_DRM number of nodes of DRM elements
!> @param[out] el_DRM "Macro" nodes matrix of DRM domain
!> @param[out] node_DRM  Node's numbers of DRM domain


      subroutine EXTRACT_DRM_EL(ns,nquad,con_quad,cs_nnz,cs,&
	                      nMDRM,tag_MDRM,n_el_DRM,ii_node_DRM,&
						  el_DRM,node_DRM)
      
      implicit none
      
	  integer*4 :: ns

      integer*4 :: cs_nnz,nquad,nnode
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nquad,5) :: con_quad

      integer*4 :: nMDRM
      integer*4, dimension(nMDRM) :: tag_MDRM

	  integer*4 :: imMDRM,nnode_DRM

	  integer*4 :: n_el_DRM,ii_node_DRM
	  integer*4, dimension(n_el_DRM,5) :: el_DRM
	  integer*4, dimension(ii_node_DRM) :: node_DRM

	  integer*4 :: ie,nn,i,j,cont,cont1
      real*8 :: AA
	  
!    "Macro" nodes matrix of DRM domain (el_DRM: [ n. el|node1|node2|node3|node4 ]) 
	  nn=ns+1
      cont=0

      do imMDRM = 1,nMDRM
	    do ie =1,nquad
		   if (con_quad(ie,1).eq.tag_MDRM(imMDRM)) then
               cont = cont+1
               el_DRM(cont,1) = ie
               el_DRM(cont,2:5) = con_quad(ie,2:5)
		   endif
		enddo
	  enddo
     

!     Node's numbers of DRM domain (node_DRM)  
      cont=1
	  cont1=nn*nn
      
	  
	  do imMDRM = 1,nMDRM
	    do ie =1,nquad
		   if (cs(cs(ie-1)+0).eq.tag_MDRM(imMDRM)) then 
		       !AA=nn;
			   !write(*,*) AA,cont, cont1
			   !read(*,*)
               node_DRM(cont:cont1) = cs((cs(ie-1)+1):(cs(ie-1)+nn*nn))
               write(21,'(I6,100(I6,:))') ie, node_DRM(cont:cont1)  !Laura_Manu 15_03_2006
			   cont=cont1+1
			   cont1=cont1+nn*nn
		   endif
		enddo
	  enddo
	  
	   
      write(*,*) size(con_quad,2)

	  return

      end subroutine EXTRACT_DRM_EL
