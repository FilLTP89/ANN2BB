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

!> @brief Writes infos about monitored points.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in,out] filec file where writing the output
!> @param[in] nmon  number of monitors
!> @param[in] n_mon indices of the nodes for monitors
!> @param[in] xr_mon monitor x-coordinate
!> @param[in] yr_mon monitor y-coordinate

      subroutine WRITE_FILE_MPGM(filec,nmon,n_mon,xr_mon,yr_mon)

         
      implicit none

      character*70 :: filec
      character*100000 :: input_line

      integer*4 :: nmon, i, trash
      integer*4 :: ileft, iright, status

      integer*4,dimension(nmon) :: n_mon


      real*8,dimension(nmon) :: xr_mon
      real*8,dimension(nmon) :: yr_mon



      open(20,file=filec)
      write(20,'(I20)')nmon

      do i = 1,nmon
         write(20,'(1I20,1X,1I20,1X,1E20.12,1X,1E20.12)') &
                                        i,n_mon(i),xr_mon(i),yr_mon(i)

      enddo
                    
      close(20)

      return

      end subroutine WRITE_FILE_MPGM

