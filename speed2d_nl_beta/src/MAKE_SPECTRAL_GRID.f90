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

!> @brief Makes local spectral grid.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nn_m number of grid nodes
!> @param[in,out] xx_m x-coordinate of grid nodes 
!> @param[in,out] yy_m y-coordinate of grid nodes
!> @param[in] cs_nnz length of cs
!> @param[in] cs  local spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tm material labels
!> @param[in] sd polynomial degree vector
!> @param[in] ne number of local elements
!> @param[out] alfa1 costant values for the bilinear map
!> @param[out] alfa2 costant values for the bilinear map
!> @param[out] beta1 costant values for the bilinear map 
!> @param[out] beta2 costant values for the bilinear map
!> @param[out] gamma1 costant values for the bilinear map
!> @param[out] gamma2 costant values for the bilinear map
!> @param[out] delta1 costant values for the bilinear map
!> @param[out] delta2 costant values for the bilinear map
!> @param[in] nn_s number of GLL nodes
!> @param[out] xx_s x-coordinate of GLL nodes 
!> @param[out] yy_s y-coordinate of GLL nodes

      subroutine MAKE_SPECTRAL_GRID(nn_m,xx_m,yy_m,cs_nnz,cs,nm,tm,sd,ne,&
                                    alfa1,beta1,gamma1,delta1,&
                                    alfa2,beta2,gamma2,delta2,&
                                    nn_s,xx_s,yy_s)
      
      implicit none
      
      integer*4 :: nn_m,cs_nnz,nm,ne,nn_s
      real*8, dimension(nn_m) :: xx_m,yy_m
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tm
      integer*4, dimension(nm) :: sd
      real*8, dimension(ne) :: alfa1,beta1,gamma1,delta1
      real*8, dimension(ne) :: alfa2,beta2,gamma2,delta2
      real*8, dimension(nn_s) :: xx_s,yy_s
      
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8 :: xp,yp,x1,x2,x3,x4,y1,y2,y3,y4,jac
      integer*4 :: im,ie,i,j,nn,ip
      
      
      nn = 2
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      
      do im = 1,nm
         if ((sd(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            nn = sd(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
         endif
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tm(im)) then
               x1 = xx_m(cs(cs(ie -1) +1))
               x2 = xx_m(cs(cs(ie -1) +nn))
               x3 = xx_m(cs(cs(ie -1) +nn*nn))
               x4 = xx_m(cs(cs(ie -1) +nn*(nn -1) +1))
               
               y1 = yy_m(cs(cs(ie -1) +1))
               y2 = yy_m(cs(cs(ie -1) +nn))
               y3 = yy_m(cs(cs(ie -1) +nn*nn))
               y4 = yy_m(cs(cs(ie -1) +nn*(nn -1) +1))
               
               alfa1(ie) = 0.25d0*(-x1 +x2 +x3 -x4)
               beta1(ie) = 0.25d0*(-x1 -x2 +x3 +x4)
               gamma1(ie) = 0.25d0*(+x1 -x2 +x3 -x4)
               delta1(ie) = 0.25d0*(+x1 +x2 +x3 +x4)
               
               alfa2(ie) = 0.25d0*(-y1 +y2 +y3 -y4)
               beta2(ie) = 0.25d0*(-y1 -y2 +y3 +y4)
               gamma2(ie) = 0.25d0*(+y1 -y2 +y3 -y4)
               delta2(ie) = 0.25d0*(+y1 +y2 +y3 +y4)
               
               jac = alfa1(ie)*beta2(ie) - alfa2(ie)*beta1(ie)
               
               if (jac.le.0.d0) then
                  write(*,*)'Error ! Orientation non-conforming !'
                  write(*,*)'(element ',ie,' is clockwise oriented)'
                  stop
               endif
               
               do j = 1,nn
                  do i = 1,nn
                     xp = alfa1(ie) * ct(i) + beta1(ie) * ct(j) &
                        + gamma1(ie) * ct(i) * ct(j) + delta1(ie)
                     yp = alfa2(ie) * ct(i) + beta2(ie) * ct(j) &
                        + gamma2(ie) * ct(i) * ct(j) + delta2(ie)
                     
                     ip = nn*(j -1) +i
                     xx_s(cs(cs(ie -1) + ip)) = xp
                     yy_s(cs(cs(ie -1) + ip)) = yp
                  enddo
               enddo
            endif
         enddo
      enddo
      
      return
      
      end subroutine MAKE_SPECTRAL_GRID
