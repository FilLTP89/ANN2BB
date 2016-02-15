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

!> @brief Calculates for time t the displacements for all the functions, so for all the DRM points.
!! @author 
!> @date August 2015
!> @version 1.0

!> @param[in] nf_drm number of DRM functions
!> @param[in] func_type_drm DRM function type
!> @param[in] func_indx_drm indices for the DRM data 
!> @param[in] func_data_drm DRM data for the calculation
!> @param[in] nfunc_data_drm
!> @param[in] t  instant time
!> @param[out] disp_t

      subroutine  GET_DISP_VALUEDRM(nf_drm,func_type_drm,func_indx_drm,nfunc_data_drm,func_data_drm,t,&
                                    disp_t)
      
!     © POLIMI, 2006, All Rights Reserved
!     Author: Laura Scandella

     
      implicit none
      
      integer*4 :: nf_drm,nfunc_data_drm
      integer*4, dimension(nf_drm) :: func_type_drm
      integer*4, dimension(nf_drm +1) :: func_indx_drm
      real*8, dimension(nfunc_data_drm) :: func_data_drm
      real*8 :: t
      
      integer*4 :: n,i
	  real*8, dimension(nf_drm,2) :: disp_t
      real*8 :: t0,t1,u0x,u1x,u0y,u1y,d_x,d_y
	  
      integer*4 :: number_of_threads

      number_of_threads = 1;

      call OMP_set_num_threads(number_of_threads)

      !call OMP_get_num_threads()

	  d_x = 0.0d0
	  d_y = 0.0d0 
 
!$OMP PARALLEL &
!$OMP PRIVATE(n,i,t0,t1,u0x,u1x,u0y,u1y,d_x,d_y )
 
!$OMP DO  
 
      do n = 1,nf_drm	   
         if (func_type_drm(n).eq.50) then
      
         ! - Displacement history DRM II step
             do i = func_indx_drm(n),func_indx_drm(n+1) -4,3
                t0 = func_data_drm(i)
                t1 = func_data_drm(i +3)
                u0x = func_data_drm(i +1)
                u1x = func_data_drm(i +4)
                u0y = func_data_drm(i +2)
                u1y = func_data_drm(i +5)
                if ((t.ge.t0).and.(t.le.t1)) then
                   d_x = (u1x - u0x) / (t1 - t0) * (t - t0)  + u0x
                   d_y = (u1y - u0y) / (t1 - t0) * (t - t0)  + u0y
                    exit 
                 endif
              enddo
			  disp_t(n,1) = d_x
			  disp_t(n,2) = d_y
          endif
      enddo

!$OMP END DO
!$OMP END PARALLEL		  
      
      return

      end subroutine GET_DISP_VALUEDRM
