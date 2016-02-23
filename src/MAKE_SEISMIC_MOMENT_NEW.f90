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

!> @brief Creates the seismic moment tensor.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] nl_sism  number of sismic loads
!> @param[in] length_cns length of check_ns 
!> @param[in] check_node_sism   info about sismic load 
!> @param[in] check_dist_node_sism distance from epicenter
!> @param[in] ielem  element index
!> @param[in] facsmom seismic moment factor 
!> @param[in] nf number of fuctions
!> @param[in] func_type function types
!> @param[in] func_indx function indeces for function data
!> @param[in] func_data  function data
!> @param[in] tt2 time 
!> @param[in] nfunc_data number of data for functions
!> @param[in] tag_func functin label
!> @param[in,out] sxx nodal values for the stress tensor
!> @param[in,out] syy nodal values for the stress tensor
!> @param[in,out] szz nodal values for the stress tensor
!> @param[in,out] sxy nodal values for the stress tensor
!> @param[in] NLFLAG flag to run nonlinear calculations

subroutine MAKE_SEISMIC_MOMENT_NEW(nn,sxx,syy,szz,sxy,&
    check_node_sism,check_dist_node_sism,&
    length_cns,ielem,facsmom,nl_sism,& 
    func_type,func_indx,func_data,nf,tt2, &
    nfunc_data,tag_func) 

    implicit none

    integer*4 :: nn
    real*8, dimension(nn,nn) :: sxx,syy,szz,sxy

    integer*4 :: ip,iq,ielem,i
    integer*4 :: length_cns
    integer*4, dimension(length_cns,5) :: check_node_sism
    real*8, dimension(length_cns,1) :: check_dist_node_sism
    integer*4 :: nl_sism,nf
    real*8, dimension(nl_sism,3) :: facsmom
    real*8 :: tt2
    integer*4, dimension(nf) :: func_type
    integer*4, dimension(nf +1) :: func_indx
    real*8, dimension(*) :: func_data

    integer*4 :: nfunc_data                    
    integer*4, dimension(nf) :: tag_func       

    real*8 :: get_func_value
    
    sxx=0d0
    syy=0d0
    szz=0d0
    sxy=0d0

    if (nl_sism.gt.0) then
        if ((ielem.ge.check_node_sism(1,1)).and. &
          (ielem.le.check_node_sism(length_cns,1))) then
              do i = 1,length_cns
                  if (ielem.eq.check_node_sism(i,1)) then
                      do iq = 1,nn
                          do ip = 1,nn
                              if ((check_node_sism(i,3).eq.iq).and. &
                                  (check_node_sism(i,2).eq.ip)) then

                                  sxx(ip,iq) = sxx(ip,iq) &
                                      + get_func_value(nf,func_type,func_indx,func_data, &  
                                      check_node_sism(i,5),tt2,check_dist_node_sism(i,1)) &
                                      * facsmom(check_node_sism(i,4),1)
                                  syy(ip,iq) = syy(ip,iq) &
                                      + get_func_value(nf,func_type,func_indx,func_data, &  
                                      check_node_sism(i,5),tt2,check_dist_node_sism(i,1)) &
                                      * facsmom(check_node_sism(i,4),2)
                                  sxy(ip,iq) = sxy(ip,iq) &
                                      + get_func_value(nf,func_type,func_indx,func_data, &  
                                      check_node_sism(i,5),tt2,check_dist_node_sism(i,1)) &
                                      * facsmom(check_node_sism(i,4),3)
                              endif
                          enddo
                      enddo

                  endif
              enddo
        endif  
    endif 

    return

end subroutine MAKE_SEISMIC_MOMENT_NEW
