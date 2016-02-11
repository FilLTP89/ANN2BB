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

!> @brief Makes spectral connectivity vector.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nelem     number of elements
!> @param[in] con_mac   connectivity matrix for spectral nodes  
!> @param[in] nmat      number of materials
!> @param[in] tag_mat   material label
!> @param[in] sdeg      polynomial degree vector 
!> @param[in] Ennz      nnzero elements for Ebin vector
!> @param[in] Ebin      pointer for connectivity (see MAKE_GRID_NODES.f90) 
!> @param[out] con_spx  spectral connectivity vector 
!> @param[out] nnode    number of spectral nodes

      subroutine MAKE_SPECTRAL_CONNECTIVITY(nelem,con_mac,nmat,tag_mat,sdeg,&
                                            Ennz,Ebin,con_nnz,con_spx,nnode)
      
      implicit none
      
      integer*4 :: nelem,nmat,Ennz,con_nnz,nnode
      integer*4, dimension(nelem,5) :: con_mac
      integer*4, dimension(nmat) :: tag_mat
      integer*4, dimension(nmat) :: sdeg
      integer*4, dimension(0:Ennz) :: Ebin
      integer*4, dimension(0:con_nnz) :: con_spx
      
      integer*4 :: nnode_mac,nn,imat,ie,i,j,k,an,bn,check
      
!     Makes spectral connectivity for counter clockwise oriented elements
      
      
      con_spx(0) = nelem +1
      do ie = 1,nelem
         do imat = 1,nmat
            if (con_mac(ie,1).eq.tag_mat(imat)) then
               nn = sdeg(imat) +1
               con_spx(ie) = con_spx(ie -1) + nn*nn +1
            endif
         enddo
      enddo
      
      nnode_mac = Ebin(0) -1
      
      nnode = nnode_mac
      
      do ie = 1,nelem
         do imat = 1,nmat
            if (con_mac(ie,1).eq.tag_mat(imat)) then
               nn = sdeg(imat) +1
               
! First put material index
               con_spx(con_spx(ie -1) +0) = con_mac(ie,1)
               
! Then put vertices
               con_spx(con_spx(ie -1) +1) = con_mac(ie,2)
               con_spx(con_spx(ie -1) + nn) = con_mac(ie,3)
               con_spx(con_spx(ie -1) + nn*(nn-1) +1) = con_mac(ie,5)
               con_spx(con_spx(ie -1) + nn*nn) = con_mac(ie,4)
               
! Then construct edge connectivity
               
! First edge (down)
               
               an = con_spx(con_spx(ie -1) +1)
               bn = con_spx(con_spx(ie -1) +nn)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(k -1) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn -k +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -1) +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k +1))
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) + k) = nnode
                  enddo
               endif
               
               
! Fourth edge (left)
               
               an = con_spx(con_spx(ie -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*(nn -1) +1)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*nn -k +1)
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*(k -1) +1) = nnode
                  enddo
               endif
               
               
! Third edge (up)
               
               an = con_spx(con_spx(ie -1) +nn*(nn -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*nn -k +1)
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*(nn -1) +k) = nnode
                  enddo
               endif
               
               
! Second edge (rigth)
               
               an = con_spx(con_spx(ie -1) +nn)
               bn = con_spx(con_spx(ie -1) +nn*nn)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(k -1) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn -k +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -1) +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k +1))
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*k) = nnode
                  enddo
               endif
               
               
! End of edge connectivity
               
! Finally fill inside the elements
               
               do j = 2,(nn -1)
                  do i = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*(j -1) +i) = nnode
                  enddo
               enddo
            endif
         enddo
      enddo
      
      return
      
      end subroutine MAKE_SPECTRAL_CONNECTIVITY
