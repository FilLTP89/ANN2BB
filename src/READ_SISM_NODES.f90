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

!> @brief Generates seismic line faults.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] Xipo Hipocenter x-coordinate
!> @param[in] Yipo Hipocenter y-coordinate
!> @param[in] X1 x-coordinate of the line fault 1-node
!> @param[in] Y1 y-coordinate of the line fault 1-node
!> @param[in] X2 x-coordinate of the line fault 2-node
!> @param[in] Y2 y-coordinate of the line fault 2-node
!> @param[in] slip1,slip2 component of the slip vector
!> @param[in] xs x-coord local spectral nodes
!> @param[in] ys y-coord local spectral nodes
!> @param[in] num_node_sism  number of seismic nodes
!> @param[in] i  column index for filling sour_ns
!> @param[in] nl_sism  number of seismic loads
!> @param[in] nnod number of local nodes
!> @param[in] max_num_node_sism max number of nodes for the specific line fault 
!> @param[out] sour_node_sism info about seismic nodes 
!> @param[out] dist_sour_node_sism distance of seismic nodes in sour_ns

      subroutine READ_SISM_NODES(Xipo,Yipo,X1,Y1,&
				  X2,Y2,slip1,slip2,nnod,xs,ys,&
                                  num_node_sism,sour_node_sism,i,&
				  dist_sour_node_sism,nl_sism,&
				  max_num_node_sism)
         
      implicit none
      
      integer*4 :: nnod,node_sism,num_node_sism,nl_sism,max_num_node_sism
      real*8, dimension(*) :: xs,ys
      integer*4, dimension(max_num_node_sism,nl_sism) :: sour_node_sism
      real*8, dimension(max_num_node_sism,nl_sism) :: dist_sour_node_sism
      integer*4 :: i,isn
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
                      node_sism = node_sism + 1	
                      sour_node_sism(node_sism,i) = isn
		      dist_sour_node_sism(node_sism,i) = sqrt((Xipo - xs(isn))**2 +&
			 (Yipo - ys(isn))**2)
                      

	           endif
	        endif
            endif
             
	 enddo

	!write(47,*)'node_sism =',node_sism
     
      return
      end subroutine READ_SISM_NODES
