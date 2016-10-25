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

!> @brief Reads dimension of monitor files.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] filec file name
!> @param[out] num_nodes number of nodes 
 
      subroutine READ_DIME_FILEPG(filec,num_nodes)
         
      implicit none

      integer*4 :: num_nodes
      integer*4 :: ileft,iright
      integer*4 :: status

      character*70 :: filec
      character*100000 :: input_line

      open(20,file=filec)
      read(20,'(A)',IOSTAT = status) input_line
      ileft = 1
      iright = len(input_line)
      read(input_line(ileft:iright),*)num_nodes
      close(20)
       

      return

      end subroutine READ_DIME_FILEPG

