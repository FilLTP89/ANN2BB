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
module seismic
    use fields
    !    
    implicit none
    !
    contains
        
        subroutine MAKE_SEISMIC_FORCES(nnt,nm,ne,nf,cs_nnz,cs,sdeg_mat,nfunc_data,nl_sism,length_cns,&
            check_node_sism,check_dist_node_sism,func_data,func_type,tag_func,func_indx,facsmom,tt1,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,displ,mvec,sism,fe)
            !
            implicit none
            ! INTENT IN
            integer*4, intent(in)                               :: nnt,nm,ne,nf,cs_nnz
            integer*4, intent(in)                               :: nfunc_data,nl_sism,length_cns
            integer*4, intent(in), dimension(nm)                :: sdeg_mat
            integer*4, intent(in), dimension(nf)                :: func_type,tag_func
            integer*4, intent(in), dimension(nf+1)              :: func_indx
            integer*4, intent(in), dimension(0:cs_nnz)          :: cs
            integer*4, intent(in), dimension(length_cns,5)      :: check_node_sism
            real*8,    intent(in)                               :: tt1
            real*8,    intent(in), dimension(ne)                :: alfa1,beta1,gamma1
            real*8,    intent(in), dimension(ne)                :: alfa2,beta2,gamma2
            real*8,    intent(in), dimension(2*nnt)             :: displ,mvec
            real*8,    intent(in), dimension(nl_sism,3)         :: facsmom       
            real*8,    intent(in), dimension(nfunc_data)        :: func_data         
            real*8,    intent(in), dimension(length_cns,1)      :: check_dist_node_sism
            ! INTENT INOUT
            real*8,           intent(inout), dimension(2*nnt)   :: sism,fe
            ! 
            real*8,     dimension(:),  allocatable      :: ct,ww
            real*8,     dimension(:),  allocatable      :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable      :: dd,det_j,fxs,fys
            real*8,     dimension(:,:),allocatable      :: sxxs,syys,szzs,sxys
            real*8,     dimension(:,:,:), allocatable   :: dstrain,dstrial
            real*8                                      :: t1ux,t1uy,t2ux
            real*8                                      :: t2uy,t1fx,t1fy,t2fx,t2fy
            integer*4                                   :: ie,ip,iq,il,im,nn,is,in
            !
            fe = 0.d0
            do ie = 1,ne 
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                
                ! ALLOCATION/INITIALIZATION LOCAL VARIABLES 
                call ALLOINIT_LOC_EL(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain, &
                    fxs,fys,displ,alfa1(ie),alfa2(ie),beta1(ie),beta2(ie),gamma1(ie),gamma2(ie))
                allocate(sxxs(nn,nn))
                allocate(syys(nn,nn))
                allocate(szzs(nn,nn))
                allocate(sxys(nn,nn))
                sxxs = 0.d0
                syys = 0.d0
                sxys = 0.d0
                szzs = 0.d0
                ! LOOP OVER GLL
                do iq = 1,nn
                    do ip = 1,nn
                        call MAKE_SEISMIC_MOMENT_NEW(nn,sxxs,syys,szzs,sxys,&
                            check_node_sism,check_dist_node_sism,length_cns,ie,facsmom, &
                            nl_sism,func_type,func_indx,func_data,nf,tt1, &
                            nfunc_data,tag_func)
                        call MAKE_INTERNAL_FORCE(nn,ww,dd,dxdx,dxdy,dydx,dydy,&
                            sxxs,syys,sxys,fxs,fys)
                    enddo
                enddo
                ! ASSIGN TO GLOBAL NODAL EXTERNAL FORCES
                do iq = 1,nn
                    do ip = 1,nn
                        is = nn*(iq-1)+ip
                        in = cs(cs(ie-1)+is)
                        sism(in)    = sism(in)      + fxs(ip,iq)
                        sism(in+nnt)= sism(in+nnt)  + fys(ip,iq)
                        fe(in)      = fe(in)        + sism(in)/mvec(in)
                        fe(in+nnt)  = fe(in+nnt)    + sism(in+nnt)/mvec(in+nnt)  
                    enddo
                enddo
            enddo
    end subroutine MAKE_SEISMIC_FORCES 
    
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

end module seismic
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
