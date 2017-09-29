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

!> @brief Reads dimensions in gridfile (*.mesh)
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] file_grid  file name (*.mesh)
!> @param[in] nmat  number of blocks (materials)
!> @param[in] tag_mat  tag for materials
!> @param[in] ndirX  number of Dirichlet b.c. (x-dir)
!> @param[in] ndirY  number of Dirichlet b.c. (y-dir)
!> @param[in] tag_dirX  label for Dirichlet b.c. (x-dir)
!> @param[in] tag_dirY  label for Dirichlet b.c. (y-dir)
!> @param[in] nneuX  number of Neumann boundary loads (x-dir)
!> @param[in] nneuY  number of Neumann boundary loads (y-dir)
!> @param[in] tag_neuX   label for Neumann boundary loads (x-dir)
!> @param[in] tag_neuY   label for Neumann boundary loads (y-dir)
!> @param[in] nabc  number of absorbing boundary conditions
!> @param[in] tag_abc   label for absorbing boundary condition
!> @param[in] nBDRM  number of DRM boundary conditions
!> @param[in] tag_BDRM   label for DRM boundary condition
!> @param[in] nb_dg  number of discontinuos interfaces (where DG scheme is applied)
!> @param[in] lab_dg   label for DG interface conditions
!> @param[out] nnode number of mesh node
!> @param[out] nquad number of quad elements
!> @param[out] nline number of edge elements

      subroutine READ_DIME_GRID_EL(file_grid,nmat,tag_mat,&
                                   ndirX,tag_dirX,ndirY,tag_dirY,&
                                   nneuX,tag_neuX,nneuY,tag_neuY,&
                                   nabc,tag_abc,&
								   nBDRM, tag_BDRM,& !Laura 27.09.2005
                                   nb_dg,lab_dg,&
                                   nnode,nquad,nline)
      
      implicit none
      
      character*70 :: file_grid
      integer*4 :: nmat,ndirX,ndirY,nneuX,nneuY,nabc,nb_dg
	  integer*4 :: nBDRM                        !Laura 27.09.2005

      integer*4, dimension(nmat) :: tag_mat
      integer*4, dimension(ndirX) :: tag_dirX
      integer*4, dimension(ndirY) :: tag_dirY
      integer*4, dimension(nneuX) :: tag_neuX
      integer*4, dimension(nneuY) :: tag_neuY
      integer*4, dimension(nabc) :: tag_abc
	  integer*4, dimension(nBDRM) :: tag_BDRM   !Laura 27.09.2005
      integer*4, dimension(nb_dg) :: lab_dg
      integer*4 :: nnode,nquad,nline
      
      integer*4 :: nelem,ie,i,mcode,ileft,iright,sl,trash,status,check
      character*80 :: input_line
      character*20 :: ecode
      
      
      nnode = 0
      nquad = 0
      nline = 0
      
      open(23,file=file_grid)
      
      do 
         read(23,'(A)') input_line
         if (input_line(1:1) .ne. '#') exit
      enddo
      
      read(input_line,*)nnode,nelem
      
      do i = 1,nnode
         read(23,'(A)')input_line
      enddo
      
      do ie = 1,nelem
         read(23,'(A)')input_line
         
         sl = len(input_line)
         ileft = 0
         iright = 0 
         do i = 1,sl
            if (input_line(i:i).ge.'A') exit
         enddo
         ileft = i
         do i = ileft,sl
            if (input_line(i:i).lt.'A') exit
         enddo
         iright = i
         
         ecode = input_line(ileft:iright)
         
         read(input_line(1:ileft),*)trash,mcode
         
         if ((ecode.eq.'quad').or.(ecode.eq.'QUAD')) then
            check = 0
            do i = 1,nmat
               if (tag_mat(i).eq.mcode) check = 1
            enddo
            !
            if (check.ne.0) nquad = nquad +1
         !
         elseif ((ecode.eq.'line').or.(ecode.eq.'LINE')) then
            check = 0
            do i = 1,ndirX
               if (tag_dirX(i).eq.mcode) check = 1
            enddo
            do i = 1,ndirY
               if (tag_dirY(i).eq.mcode) check = 1
            enddo
            do i = 1,nneuX
               if (tag_neuX(i).eq.mcode) check = 1
            enddo
            do i = 1,nneuY
               if (tag_neuY(i).eq.mcode) check = 1
            enddo
            do i = 1,nabc
               if (tag_abc(i).eq.mcode) check = 1
            enddo   
			do i = 1,nBDRM                             !DRM Scandella 27.09.2005
               if (tag_BDRM(i).eq.mcode) check = 1     !DRM Scandella 27.09.2005
            enddo                                      !DRM Scandella 27.09.2005
            do i = 1,nb_dg
               if (lab_dg(i).eq.mcode) check = 1
            enddo
                        
            if (check.ne.0) nline = nline +1
         endif
      enddo
      
      close(23)
      
      return
      end subroutine READ_DIME_GRID_EL
