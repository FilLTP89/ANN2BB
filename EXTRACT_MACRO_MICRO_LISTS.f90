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

!> @brief Extracts the list of macro and micro node coordinates
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] ns
!> @param[in] nquad number of quad elements
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] nMDRM
!> @param[in] tag_MDRM labels for DRM blocks
!> @param[in] n_el_DRM number of DRM elements
!> @param[in] nnode number of nodes
!> @param[out] xs x coordinate
!> @param[out] ys  y coordinate


      subroutine EXTRACT_MACRO_MICRO_LISTS(ns,nquad,cs_nnz,cs,&
	                      nMDRM,tag_MDRM,n_el_DRM,&
						  nnode,xs,ys)
      
      implicit none
      
	  integer*4 :: ns,nnode

      integer*4 :: cs_nnz,nquad,nnode_DRM_ma, nnode_DRM_ma_eff, nnode_DRM_mi, nnode_DRM_mi_eff    !DRM Scandella 09.05.2007
      integer*4, dimension(0:cs_nnz) :: cs

      integer*4 :: nMDRM
      integer*4, dimension(nMDRM) :: tag_MDRM
	  integer*4 :: imMDRM,nnode_DRM
	  integer*4 :: n_el_DRM

      integer*4, dimension(:), allocatable ::  node_DRM_ma, node_DRM_mi                           !DRM Scandella 09.05.2007
	  character*70 :: lista_macro, lista_macro_coor, lista_micro, lista_micro_coor                !DRM Scandella 09.05.2007
	  real*8, dimension(nnode) :: xs,ys                                                           !DRM Scandella 09.05.2007
      integer*4, dimension(:), allocatable :: node_DRM_ma_eff, node_DRM_mi_eff                    !DRM Scandella 09.05.2007 
	  real*8, dimension(:), allocatable :: x_DRM_ma_eff, y_DRM_ma_eff, x_DRM_mi_eff, y_DRM_mi_eff !DRM Scandella 09.05.2007 

	  integer*4 :: ie,nn,i,j,cont_ma, check, k, node_ma, nd, cont_mi, node_mi

	  nn=ns+1
	  cont_ma=1  !DRM Scandella 09.05.2007
	  cont_mi=1  !DRM Scandella 09.05.2007

      lista_macro='lista_macro.txt'                            !DRM Scandella 09.05.2007    
	  open(1,file=lista_macro)                                 !DRM Scandella 09.05.2007   
	  lista_macro_coor='lista_macro_coor.txt'                  !DRM Scandella 09.05.2007    
	  open(2,file=lista_macro_coor)                            !DRM Scandella 09.05.2007
      lista_micro='lista_micro.txt'                            !DRM Scandella 09.05.2007    
	  open(3,file=lista_micro)                                 !DRM Scandella 09.05.2007   
	  lista_micro_coor='lista_micro_coor.txt'                  !DRM Scandella 09.05.2007    
	  open(4,file=lista_micro_coor)                            !DRM Scandella 09.05.2007

	  nnode_DRM_ma=n_el_DRM*4                                  !DRM Scandella 09.05.2007 
	  nnode_DRM_mi=n_el_DRM*(nn*nn-4)                          !DRM Scandella 09.05.2007 
	  allocate(node_DRM_ma(nnode_DRM_ma))
	  allocate(node_DRM_mi(nnode_DRM_mi))

      

	  do imMDRM = 1,nMDRM
	    do ie =1,nquad
		   if (cs(cs(ie-1)+0).eq.tag_MDRM(imMDRM)) then 
               node_DRM_ma(cont_ma)=cs(cs(ie-1)+1)                 !DRM Scandella 09.05.2007  
			   node_DRM_ma(cont_ma+1)=cs(cs(ie -1) +nn)            !DRM Scandella 09.05.2007  
			   node_DRM_ma(cont_ma+2)=cs(cs(ie -1) +nn*nn)         !DRM Scandella 09.05.2007
			   node_DRM_ma(cont_ma+3)=cs(cs(ie -1) +nn*(nn -1) +1) !DRM Scandella 09.05.2007

			   node_DRM_mi(cont_mi:cont_mi+nn-3)=cs((cs(ie-1)+2):(cs(ie-1)+nn-1))                                  !DRM Scandella 09.05.2007  
			   node_DRM_mi(cont_mi+nn-2:cont_mi+nn*(nn -1) -3)=cs((cs(ie-1)+nn+1):(cs(ie-1)+nn*(nn -1)))           !DRM Scandella 09.05.2007  
			   node_DRM_mi(cont_mi+nn*(nn -1)-2:cont_mi+nn*nn-5)=cs((cs(ie-1)+nn*(nn -1)+2):(cs(ie-1)+nn*nn-1))    !DRM Scandella 09.05.2007

               cont_ma=cont_ma+4                                   !DRM Scandella 09.05.2007
               cont_mi=cont_mi+nn*nn-4   			                                   
		   endif
		enddo
	  enddo


    !  Size of DRM macro nodes without duplicates

    nnode_DRM_ma_eff = 0
	do i = 1,nnode_DRM_ma
	   check = 0
	   do j = i,nnode_DRM_ma
	      if (node_DRM_ma(i).eq.node_DRM_ma(j)) then
		     check = check+1 
		  endif
       enddo
	   if (check.eq.1) then
	   nnode_DRM_ma_eff = nnode_DRM_ma_eff+1
	   endif
	enddo

	!  Size of DRM micro nodes without duplicates

    nnode_DRM_mi_eff = 0
	do i = 1,nnode_DRM_mi
	   check = 0
	   do j = i,nnode_DRM_mi
	      if (node_DRM_mi(i).eq.node_DRM_mi(j)) then
		     check = check+1 
		  endif
       enddo
	   if (check.eq.1) then
	   nnode_DRM_mi_eff = nnode_DRM_mi_eff+1
	   endif
	enddo

	
	allocate(node_DRM_ma_eff(nnode_DRM_ma_eff))   !DRM Scandella 09.05.2007
	allocate(x_DRM_ma_eff(nnode_DRM_ma_eff))      !DRM Scandella 09.05.2007
	allocate(y_DRM_ma_eff(nnode_DRM_ma_eff))      !DRM Scandella 09.05.2007
	allocate(node_DRM_mi_eff(nnode_DRM_mi_eff))   !DRM Scandella 09.05.2007
	allocate(x_DRM_mi_eff(nnode_DRM_mi_eff))      !DRM Scandella 09.05.2007
	allocate(y_DRM_mi_eff(nnode_DRM_mi_eff))      !DRM Scandella 09.05.2007


    !  DRM macro nodes without duplicates
    k = 0
	do i = 1,nnode_DRM_ma
	   node_ma = node_DRM_ma(i)
	   if (k.ge.1) then
	      check = 0
	         do j = 1,k
	            if (node_ma.eq.node_DRM_ma_eff(j)) then
		           check = 1
				 exit
		        endif
	         enddo
	         if (check.ne.1) then
                k = k+1
	            node_DRM_ma_eff(k)=node_ma	   	              
             endif
		else
           k = k+1
	       node_DRM_ma_eff(k)=node_ma
		endif
	enddo

	do i = 1,nnode_DRM_ma_eff
	    write(*,*) node_DRM_ma_eff(i)
	    x_DRM_ma_eff(i)=xs(node_DRM_ma_eff(i))
		y_DRM_ma_eff(i)=ys(node_DRM_ma_eff(i))
    enddo

	  write(1,'(I6)')(node_DRM_ma_eff(i),i=1,nnode_DRM_ma_eff)                  !DRM Scandella 09.05.2007 
	  close(1)                                                                  !DRM Scandella 09.05.2007
      write(2,'(2f20.4)')(x_DRM_ma_eff(i),y_DRM_ma_eff(i),i=1,nnode_DRM_ma_eff) !DRM Scandella 09.05.2007
	  close(2)                                                                  !DRM Scandella 09.05.2007
	  deallocate(node_DRM_ma,node_DRM_ma_eff,x_DRM_ma_eff,y_DRM_ma_eff)

    !  DRM micro nodes without duplicates
    k = 0
	do i = 1,nnode_DRM_mi
	   node_mi = node_DRM_mi(i)
	   if (k.ge.1) then
	      check = 0
	         do j = 1,k
	            if (node_mi.eq.node_DRM_mi_eff(j)) then
		           check = 1
				 exit
		        endif
	         enddo
	         if (check.ne.1) then
                k = k+1
	            node_DRM_mi_eff(k)=node_mi	   	              
             endif
		else
           k = k+1
	       node_DRM_mi_eff(k)=node_mi
		endif
	enddo

	do i = 1,nnode_DRM_mi_eff
	    x_DRM_mi_eff(i)=xs(node_DRM_mi_eff(i))
		y_DRM_mi_eff(i)=ys(node_DRM_mi_eff(i))
    enddo

	  write(3,'(I6)')(node_DRM_mi_eff(i),i=1,nnode_DRM_mi_eff)                  !DRM Scandella 09.05.2007 
	  close(3)                                                                  !DRM Scandella 09.05.2007
      write(4,'(2f20.4)')(x_DRM_mi_eff(i),y_DRM_mi_eff(i),i=1,nnode_DRM_mi_eff) !DRM Scandella 09.05.2007
	  close(4)                                                                  !DRM Scandella 09.05.2007
	  deallocate(node_DRM_mi,node_DRM_mi_eff,x_DRM_mi_eff,y_DRM_mi_eff)

      return 
	                                                                          
      end subroutine EXTRACT_MACRO_MICRO_LISTS