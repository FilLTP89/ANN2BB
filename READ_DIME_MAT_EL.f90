!    Copyright (C) 2012 The SPEED FOUNDATION
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

!> @brief Reads dimensions in filemate (*.mate)
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] file_mat  file name (*.mate)
!> @param[out] nmat  number of blocks (materials)
!> @param[out] nload_dirX   number of Dirichlet b.c. (x-dir)
!> @param[out] nload_dirY  number of Dirichlet b.c. (y-dir)
!> @param[out] nload_neuX  number of Neumann boundary loads (x-dir)
!> @param[out] nload_neuY  number of Neumann boundary loads (y-dir)
!> @param[out] nload_poiX   number of point loads (x-dir)
!> @param[out] nload_poiY   number of point loads (y-dir)
!> @param[out] nload_plaX   number of plane wave loads (x-dir)
!> @param[out] nload_plaY   number of plane wave loads (y-dir)
!> @param[out] nload_abc   number of abc boundary conditions
!> @param[out] nload_MDRM
!> @param[out] nload_BDRM
!> @param[out] nload_PDRM
!> @param[out] nb_load_dg   number of DG interfaces
!> @param[out] nfunc  number of time functions
!> @param[out] nfunc_data  number of data for time functions
!> @param[out] nfunc_drm  number of DRM time functions
!> @param[out] nfunc_data_drm
!> @param[out] nload_sism  number of seismic loads
!> @param[out] n_test  1 for test case mode

      subroutine READ_DIME_MAT_EL(file_mat,nmat, &
                                  nload_dirX,nload_dirY, &
                                  nload_neuX,nload_neuY, & 
                                  nload_poiX,nload_poiY, &
                                  nload_plaX,nload_plaY, &
                                  nload_sism, &
                                  nload_abc, &
			                      nload_MDRM,nload_BDRM,nload_PDRM, &  !DRM Scandella 17.10.2005
                                  nfunc,nfunc_data,n_test, &
			                      nfunc_drm,nfunc_data_drm, &                   !DRM Scandella 11.04.2006
                                  nb_load_dg)
            
      implicit none
      
      character*70   :: file_mat

      integer*4 :: nmat
                         
      integer*4 :: nload_dirX,nload_dirY,nload_neuX,nload_neuY, nb_load_dg
      integer*4 :: nload_poiX,nload_poiY,nload_plaX,nload_plaY,nload_sism,nload_abc
	  integer*4 :: nload_MDRM,nload_BDRM         !DRM Scandella 27.09.2005
      integer*4 :: nload_PDRM                    !DRM Scandella 17.10.2005 

      integer*4 :: nfunc,nfunc_data, n_test
	  integer*4 :: nfunc_drm,nfunc_data_drm      !DRM Scandella 11.04.2006
      
      character*80 :: input_line
      character*4 :: keyword
      integer*4 :: status
      integer*4 :: tag_func,func_type,func_nd
	  integer*4 :: tag_func_drm,func_type_drm,func_nd_drm   !DRM Scandella 11.04.2006
      
      nmat = 0
       
      nload_dirX = 0
      nload_dirY = 0
      nload_neuX = 0 
      nload_neuY = 0 
      nload_poiX = 0 
      nload_poiY = 0
      nload_plaX = 0
      nload_plaY = 0 
      nload_sism = 0
      nload_abc = 0
	  nload_MDRM = 0       !DRM Scandella 27.09.2005
	  nload_BDRM = 0       !DRM Scandella 27.09.2005
	  nload_PDRM = 0       !DRM Scandella 17.10.2005
      nb_load_dg = 0
      nfunc = 0
      nfunc_data = 0
	  nfunc_drm = 0        !DRM Scandella 11.04.2006
      nfunc_data_drm = 0   !DRM Scandella 11.04.2006
      n_test = 0

      
      open(23,file=file_mat)
      
      do
         read(23,'(A)',IOSTAT = status) input_line
         
         if (status.ne.0) exit
         
         keyword = input_line(1:4)
         
         if (keyword.eq.'MATE') then
            nmat = nmat + 1				!
         elseif (keyword.eq.'DIRX') then
            nload_dirX = nload_dirX + 1
         elseif (keyword.eq.'DIRY') then
            nload_dirY = nload_dirY + 1
         elseif (keyword.eq.'NEUX') then
            nload_neuX = nload_neuX + 1
         elseif (keyword.eq.'NEUY') then
            nload_neuY = nload_neuY + 1
         elseif (keyword.eq.'PLOX') then
            nload_poiX = nload_poiX + 1
         elseif (keyword.eq.'PLOY') then
            nload_poiY = nload_poiY + 1
         elseif (keyword.eq.'PLAX') then
            nload_plaX = nload_plaX + 1
         elseif (keyword.eq.'PLAY') then
            nload_plaY = nload_plaY + 1
         elseif (keyword.eq.'SISM') then
            nload_sism = nload_sism + 1
         elseif (keyword.eq.'ABSO') then
            nload_abc = nload_abc + 1
		 elseif (keyword.eq.'MDRM') then     !DRM Scandella 27.09.2005 
            nload_MDRM = nload_MDRM + 1      !DRM Scandella 27.09.2005 
		 elseif (keyword.eq.'BDRM') then     !DRM Scandella 27.09.2005
            nload_BDRM = nload_BDRM + 1      !DRM Scandella 27.09.2005
		 elseif (keyword.eq.'PDRM') then     !DRM Scandella 27.09.2005
            nload_PDRM = nload_PDRM + 1      !DRM Scandella 17.10.2005
         elseif (keyword.eq.'DGIC') then
            nb_load_dg = nb_load_dg + 1            
         elseif (keyword.eq.'TEST') then
            n_test = n_test + 1

         elseif (keyword.eq.'FUNC') then
            nfunc = nfunc + 1
            read(input_line(5:),*)tag_func,func_type
            
            if (func_type.eq.0) then
               nfunc_data = nfunc_data + 0
            elseif (func_type.eq.1) then
               nfunc_data = nfunc_data + 2
            elseif (func_type.eq.2) then
               nfunc_data = nfunc_data + 2
            elseif (func_type.eq.3) then
               read(input_line(5:),*)tag_func,func_type,func_nd
               nfunc_data = nfunc_data + 2*func_nd
            elseif (func_type.eq.4) then
               nfunc_data = nfunc_data + 2
            elseif (func_type.eq.5) then
               nfunc_data = nfunc_data + 2
   	    elseif (func_type.eq.31) then
               nfunc_data = nfunc_data + 3
            
            ! Non-Linear Elasticity   
            elseif (func_type.eq.60) then                                    
               read(input_line(5:),*)tag_func,func_type,func_nd              
               nfunc_data = nfunc_data + 2*func_nd                           
	    elseif (func_type.eq.61) then                        
               read(input_line(5:),*)tag_func,func_type,func_nd              
               nfunc_data = nfunc_data + 2*func_nd                           
	    endif
		 elseif (keyword.eq.'FDRM') then                                     !DRM Scandella 11.04.2006 
            nfunc_drm = nfunc_drm + 1                                        !DRM Scandella 11.04.2006           
            read(input_line(5:),*)tag_func_drm,func_type_drm                 !DRM Scandella 11.04.2006
            if (func_type_drm.eq.50) then                                    !DRM Scandella 11.04.2006
               read(input_line(5:),*)tag_func_drm,func_type_drm,func_nd_drm  !DRM Scandella 11.04.2006
               nfunc_data_drm = nfunc_data_drm + 3*func_nd_drm               !DRM Scandella 11.04.2006
            endif
         endif
      enddo
      
      close(23)
      
      return
      end subroutine READ_DIME_MAT_EL

