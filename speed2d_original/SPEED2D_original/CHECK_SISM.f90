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

!> @brief Fills array check_ns for seismic force.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] cs_nnz length of cs  
!> @param[in] cs spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tag_mat label for materias 
!> @param[in] sdeg_mat polynomial degree vector
!> @param[in] ne number of elements
!> @param[in] nl_sism  number of seismic loads
!> @param[in] num_node_sism  number of seismic nodes
!> @param[in] max_num_node_sism  max number of seismic nodes
!> @param[in] sour_node_sism seismic source
!> @param[in] dist_sour_node_sism  distance of the node from the seismic source
!> @param[in] length_cns  length check seismic nodes (useless)
!> @param[in] fun_sism  function associeted with the seismic load
!> @param[in] nf  number of functions
!> @param[in] tag_func tag for seismic functions
!> @param[in] valsism values for the seismic loads
!> @param[out] check_node_sism  check_node_sism(i,1)  seismic source node, 
!!      check_node_sism(i,2)  seismic function for the node i,
!!      check_node_sism(i,3)  i for seismic load,
!!      check_node_sism(i,4)  number of the element
!!      check_node_sism(i,4)  func number       
!> @param[out]  check_dist_node_sism check_dist_node_sism(i,1)  dist_sour_node_sism / valsism

      subroutine CHECK_SISM(cs_nnz,cs,&
                          nm,tag_mat,sdeg_mat,&
                          ne,&
                          nl_sism,&
                          num_node_sism,max_num_node_sism,&
                          sour_node_sism,dist_sour_node_sism,&
                          check_node_sism,check_dist_node_sism,&
                          length_cns,&
                          fun_sism,nf,tag_func,valsism)


      implicit none

      integer*4 :: cs_nnz,nm,ne,nl_sism
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,sdeg_mat
      integer*4, dimension(nl_sism) :: num_node_sism
      integer*4, dimension(nl_sism) :: fun_sism
      real*8, dimension(nl_sism,12) :: valsism
      real*8 :: vel_rup
      integer*4 :: nf 
      integer*4, dimension(nf) :: tag_func
      integer*4 :: max_num_node_sism
      integer*4, dimension(max_num_node_sism,nl_sism) :: sour_node_sism
      real*8, dimension(max_num_node_sism,nl_sism) :: dist_sour_node_sism
      integer*4 :: length_cns
      integer*4, dimension(length_cns,5) :: check_node_sism
      real*8, dimension(length_cns,1) :: check_dist_node_sism
      integer*4 :: im,ie,isism,nn
      integer*4 :: is,in
      integer*4 :: i,j,k,h
      integer*4 :: fn,fun_sism_k
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:), allocatable :: dxdx,dydx,dxdy,dydy
      real*8, dimension(:,:), allocatable :: det_j

      !write(50,*),length_cns
      nn = 2
      h = 1

      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
   
      ne = cs(0) -1  
         
      do im = 1,nm
         if ((sdeg_mat(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            deallocate(dxdx,dydx,dxdy,dydy,det_j)
         
            nn = sdeg_mat(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
            allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))   
         endif

         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
            
            !h = 1
                                                
            do isism = 1,nl_sism
               do fn = 1,nf
                  if (tag_func(fn).eq.fun_sism(isism)) then
                     vel_rup = valsism(isism,12)
                     fun_sism_k = fn
                  endif
               enddo
               do k = 1,num_node_sism(isism)
                  do j = 1,nn   
                     do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) +is)
                        if (in.eq.sour_node_sism(k,isism)) then
                            check_node_sism(h,1) = ie
                            check_node_sism(h,2) = i
                            check_node_sism(h,3) = j
                            check_node_sism(h,4) = isism
                            check_node_sism(h,5) = fun_sism_k
                            !check_node_sism(h,5) = 5
                            check_dist_node_sism(h,1) = dist_sour_node_sism(k,isism) &
                                                        / vel_rup

                            !write(50,*),check_node_sism(h,1),' | ',&
                            !            check_node_sism(h,2),' | ',&
                            !            check_node_sism(h,3),' | ',&
                            !            check_node_sism(h,4),' | ',&
                            !            check_node_sism(h,5),' | ',& 
                            !            check_dist_node_sism(h,1),' | ',&
                            !            vel_rup
                            h = h + 1
                        endif
                     enddo
                  enddo

               enddo
            enddo                           

! Seismic moment scale factor - end
                                    
      endif
    enddo
  enddo

      !write(50,*),h
      !do isism =  1,h-1
      !	                             write(52,*),check_node_sism(isism,1),' | ',&
      !                                  check_node_sism(isism,2),' | ',&
      !                                  check_node_sism(isism,3),' | ',&
      !                                  check_node_sism(isism,4),' | ',&
      !                                  check_dist_node_sism(isism,1),' | '
      !enddo
       
      return
      end subroutine CHECK_SISM
