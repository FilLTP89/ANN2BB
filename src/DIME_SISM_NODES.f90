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

!> @brief Find the position for storing the column index for the i-th row
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] Xipo,Yipo coordinate of the Hypocenter 
!> @param[in] X1,Y1 first end point of the fault
!> @param[in] X1,Y1 second end point of the fault
!> @param[in] slip1,slip2 direction of the slip
!> @param[in] xs,ys coordinates of the spectral nodes
!> @param[out] node_sism number of nodes where seismic force is applied 

      subroutine DIME_SISM_NODES(Xipo,Yipo,X1,Y1,&
				  X2,Y2,slip1,slip2,nnod,xs,ys,&
                                  node_sism)

      implicit none
      
      integer*4 :: nnod,node_sism
      real*8, dimension(*) :: xs,ys
      integer*4 :: isn
      real*8 :: Xipo,Yipo,X1,Y1,X2,Y2,slip1,slip2,tol
      real*8 :: Xmax,Xmin,Ymax,Ymin,cost

      node_sism = 0
      tol = 1.0d-2
      cost = xipo*slip2-yipo*slip1

	 !Check on the relative position of fault end points
	 !In this preliminary part we check if P1 point is the lower-left point or not
	 !If this check would fail the two points are swapped
  
	 if (X1.gt.X2) then
	    Xmax=X1
	    Xmin=X2
	 else
	    Xmax=X2
            Xmin=X1
	 endif
         if (Y1.gt.Y2) then
            Ymax=Y1
            Ymin=Y2
         else
            Ymax=Y2 
            Ymin=Y1
         endif


  
	 do isn = 1,nnod

	    if (dabs(ys(isn)*slip1-xs(isn)*slip2+cost).le.tol) then

		if ((xs(isn).ge.Xmin).and.(xs(isn).le.Xmax)) then
		   if ((ys(isn).ge.Ymin).and.(ys(isn).le.Ymax)) then
 			node_sism = node_sism+1
	          endif
	        endif
            endif
             
	 enddo

      return
      end subroutine DIME_SISM_NODES
      
