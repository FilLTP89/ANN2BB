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

!> @brief Assign Dirichlet boundary conditions in a matrix
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in] nnode_dirX number of Dirchlet-x nodes
!> @param[in] inode_dirX id of a Dirichlet-x node
!> @param[in] nnz_full_morse nnzero elements 
!> @param[in] IFULL_MORSE, JFULL_MORSE vector for morse format storage
!> @param[in] td pol. degree for time + 1
!> @param[inout] MFULL_MORSE matrix where Dirichlet-x boundary conditions are applied

       subroutine  MAKE_DIRIX_BC(nnode_dirX,inode_dirX, &
                          NNZ_FULL_MORSE, IFULL_MORSE, JFULL_MORSE, MFULL_MORSE, td )

       implicit none
       
       integer*4 :: NNZ_FULL_MORSE, nnode_dirX, td
       integer*4, dimension(nnode_dirX) :: inode_dirX
       integer*4, dimension(NNZ_FULL_MORSE) :: JFULL_MORSE, IFULL_MORSE
       real*8, dimension(NNZ_FULL_MORSE) :: MFULL_MORSE
       integer*4 :: i,j,in, ic
       

       do i = 1,nnode_dirX
          do j = 1, td
             in = (inode_dirX(i)-1)*td +j
             ic = 1
             do while (ic .le. NNZ_FULL_MORSE)
                if(IFULL_MORSE(ic) .eq. in) MFULL_MORSE(ic) = 0.d0
                if(IFULL_MORSE(ic) .eq. in .and. JFULL_MORSE(ic) .eq. in ) MFULL_MORSE(ic) = 1.d0
                ic = ic + 1
             enddo
          enddo
       enddo               


       end subroutine  MAKE_DIRIX_BC
