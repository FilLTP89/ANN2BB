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

!> @brief Perform time integration with leap-frog method.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nnt number of nodes
!> @param[in] xs,ys
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tag_mat label for materials
!> @param[in] sdeg_mat pol. degrees
!> @param[in] prop_mat material properties
!> @param[in] alfa_el,beta_el,gamma_el constatn for bilinear map
!> @param[in] cs_nnz_bc length of cs_bc
!> @param[in] cs_bc spectral connectivity vector for boundary
!> @param[in] nl_dirX   number of Dirichlet b.c. (x-dir)
!> @param[in] nl_dirY   number of Dirichlet b.c. (y-dir)
!> @param[in] tag_dirX  label for Dirichlet b.c. (x-dir)
!> @param[in] tag_dirY  label for Dirichlet b.c. (y-dir)
!> @param[in] nelem_abc number of ab elements
!> @param[in] nedge_abc number of ab edges
!> @param[in] ielem_abc index of ab elements 
!> @param[in] iedge_abc index of ab edges
!> @param[in] nf  number of time functions
!> @param[in] func_type type for functions
!> @param[in] func_indx indeces for functions
!> @param[in] nfunc_data number of data for functions
!> @param[in] func_data data for functions
!> @param[in] tag_func label for functions
!> @param[in] nf_drm  number of DRM time functions
!> @param[in] func_type_drm type for DRM functions
!> @param[in] func_indx_drm indeces for DRM functions
!> @param[in] nfunc_data_drm number of data for DRM functions
!> @param[in] func_data_drm data for DRM functions
!> @param[in] ndt_monitor time step for monitors
!> @param[in] N_TOT,IN_TOT,JN_TOT damping matrix in morse format
!> @param[in] NNZ_N nnzero el for damping matrix
!> @param[in] mvec mass matrix 
!> @param[in] Fmat external loads
!> @param[in] u0,v1 initial conditions
!> @param[in] nts number of time steps
!> @param[in] dt time step
!> @param[in] nmonit number of monitors
!> @param[in] node_m id for monitors
!> @param[in] nsnap number of snapshots (dummy)
!> @param[in] itersnap dummy
!> @param[in] check_node_sism array for seismic sources
!> @param[in] check_dist_node_sism distance of the node from the seismic source
!> @param[in] length_cns ength check seismic nodes (useless)
!> @param[in] facsmom factor for seismic momentum
!> @param[in] nl_sism number of seismic loads
!> @param[in] make_damping_yes_or_not damping variable
!> @param[in] option_out_var options for output
!> @param[in] nnode_TOT_eff number of DRM total nodes without duplicate
!> @param[in] node_TOT_eff nodes of total DRM domain without duplicate
!> @param[in] tagstep label of the step of DRM Analysis
!> @param[in] ns
!> @param[in] n_el_DRM number of DRM elements
!> @param[in] el_DRM_eff "Macro" nodes matrix of DRM domain
!> @param[in] K_DRM stiffnes matrix of DRM elements
!> @param[in] nnode_BD_eff number of DRM internal boundary nodes
!> @param[in] nload_MDRM_el
!> @param[in] tag_MDRM_el labels for DRM blocks
!> @param[in] val_PDRM_el
!> @param[in] fn_ord order of PDRM function
!> @param[in] node_PDRM_el
!> @param[in] glob_drm_x
!> @param[in] glob_drm_y
subroutine TIME_LOOP_NL(nnt,xs,ys,cs_nnz,cs,nm,tag_mat,sdeg_mat,prop_mat,ne,    &
    alfa1,beta1,gamma1,alfa2,beta2,gamma2,delta1,delta2,cs_nnz_bc,cs_bc,        &
    nl_dirX,tag_dirX,nl_dirY,tag_dirY,nl_abc,tag_abc,nelem_abc,nedge_abc,       &
    ielem_abc,iedge_abc,nf,func_type,func_indx,nfunc_data,func_data,tag_func,   &
    nf_drm,func_type_drm,func_indx_drm,nfunc_data_drm,func_data_drm,ndt_monitor,&
    N_TOT,IN_TOT,JN_TOT,NNZ_N,mvec,Fmat,u0,v1,nts,dt,nmonit,node_m,nsnap,       &
    itersnap,check_node_sism,check_dist_node_sism,length_cns,facsmom,nl_sism,   &
    make_damping_yes_or_not,nnode_TOT,node_TOT,tagstep,ns,n_el_DRM,el_DRM,K_DRM,&
    nnode_BD,nMDRM,tag_MDRM,val_PDRM,fun_ord,node_PDRM,glob_x,glob_y,           &
    option_out_var,test,nelem_dg,IDG_only_uv,JDG_only_uv,MDG_only_uv,nnz_dg_only_uv)    
    
    use fields 
    use write_output
    use seismic
    use nonlinear2d
    !
    implicit none
    !
    integer*4                               :: NNZ_K,NNZ_N,error,test,nelem_dg,nnz_dg_only_uv
    integer*4                               :: i,is,in,id1,id2,idt,j,jn,k,kn,iaz,jaz,kaz
    integer*4                               :: clock_start,clock_finish,start,finish
    integer*4                               :: nnt,cs_nnz,nm,ne,cs_nnz_bc,nf
    integer*4                               :: isnap,its,fn,ie,nn,ip,im
    integer*4                               :: nl_dirX,nl_dirY,nl_abc
    integer*4                               :: nnd,nts,nmonit,nsnap
    integer*4                               :: number_of_threads 
    integer*4                               :: nfunc_data
    integer*4                               :: i4t
    
    integer*4, dimension(3)                 :: clock
    integer*4, dimension(nm)                :: tag_mat,sdeg_mat
    integer*4, dimension(nf)                :: func_type,tag_func
    integer*4, dimension(nf+1)              :: func_indx
    integer*4, dimension(NNZ_N)             :: JN_TOT
    integer*4, dimension(nsnap)             :: itersnap
    integer*4, dimension(nmonit)            :: node_m
    integer*4, dimension(nl_abc)            :: tag_abc
    integer*4, dimension(nl_dirX)           :: tag_dirX
    integer*4, dimension(nl_dirY)           :: tag_dirY
    integer*4, dimension(0:2*nnt)           :: IN_TOT,IDG_only_uv
    integer*4, dimension(0:cs_nnz)          :: cs
    integer*4, dimension(0:cs_nnz_bc)       :: cs_bc
    integer*4, dimension(nnz_dg_only_uv)    :: JDG_only_uv
    
    real*8                                  :: dt,r8t,eps,pi
    real*8                                  :: time_disp_DRM,total_disp_DRM
    real*8                                  :: get_func_value, tt0,tt1,tt2,dt2
    real*8                                  :: time_in_seconds,time_total,time_u,total_u
    real*8                                  :: time_fk,total_fk,time_fd,total_fd,time_fe
    real*8                                  :: total_fe,time_force_DRM,total_force_DRM
    
    real*8, dimension(ne)                   :: alfa1,beta1,gamma1,delta1
    real*8, dimension(ne)                   :: alfa2,beta2,gamma2,delta2 
    real*8, dimension(nnt)                  :: xs,ys
    real*8, dimension(nm,9)                 :: prop_mat
    real*8, dimension(NNZ_N)                :: N_TOT
    real*8, dimension(nf,2*nnt)             :: Fmat
    real*8, dimension(nfunc_data)           :: func_data         
    real*8, dimension(nnz_dg_only_uv)       :: MDG_only_uv

    real*8, dimension(:), allocatable       :: func_value

    !*************************************************
    !
    !            DRM VARIABLES
    !
    !*************************************************

    integer*4                               :: p,ns,tagstep
    integer*4                               :: nnode_TOT,nnode_BD
    integer*4                               :: n_el_DRM,nf_drm 
    integer*4                               :: nMDRM,imDRM
    integer*4                               :: nfunc_data_drm
    integer*4, dimension(n_el_DRM,5 )       :: el_DRM        
    integer*4, dimension(nMDRM      )       :: tag_MDRM
    integer*4, dimension(nnode_TOT  )       :: node_TOT,node_PDRM,fun_ord
    integer*4, dimension(nf_drm     )       :: func_type_drm 
    integer*4, dimension(nf_drm +1  )       :: func_indx_drm 
    integer*4, dimension(nf_drm,2   )       :: glob_x,glob_y        
    real*8, dimension(nnode_TOT,3)          :: val_PDRM           
    real*8, dimension(nnode_TOT)            :: px_t,py_t         
    real*8, dimension(nfunc_data_drm)       :: func_data_drm     
    real*8, dimension(:,:), allocatable     :: disp_PDRM_t
    real*8, dimension(2*(ns+1)*(ns+1)*n_el_DRM,2*(ns+1)*(ns+1)) :: K_DRM 
    real*8, dimension(2*(ns+1)**2,2*(ns+1)**2)                  :: K_el 


    !************************************************
    !
    !             OUTPUT VARIABLES
    !
    !************************************************

    integer*4, dimension(nmonit) :: unit_disp
    integer*4, dimension(nmonit) :: unit_vel
    integer*4, dimension(nmonit) :: unit_acc
    integer*4, dimension(nmonit) :: unit_stress
    integer*4, dimension(nmonit) :: unit_strain
    integer*4, dimension(nmonit) :: unit_omega
    integer*4, dimension(nmonit) :: unit_uDRM 
    integer*4, dimension (6)     :: option_out_var           
    type(nodepatched)            :: disout
    logical                      :: condition

    !************************************************
    !
    !            BOUNDARY CONDITIONS
    !
    !************************************************

    integer*4                               :: nnode_dirX,nnode_dirY,nnode_neuX,nnode_neuY
    integer*4                               :: nelem_abc,nedge_abc
    integer*4                               :: iedge,nedge,ied1,ied2,iel1,iel2,iel3,iel4
    integer*4                               :: edge_ia,edge_ja,edge_ib,edge_jb
    integer*4, dimension(nelem_abc)         :: ielem_abc
    integer*4, dimension(nedge_abc)         :: iedge_abc
    integer*4, dimension(:), allocatable    :: i4count
    integer*4, dimension(:), allocatable    :: inode_dirX
    integer*4, dimension(:), allocatable    :: inode_dirY
    integer*4, dimension(:), allocatable    :: inode_neuX
    integer*4, dimension(:), allocatable    :: inode_neuY
    real*8                                  :: edge_lx,edge_ly,edge_ll,edge_nx,edge_ny

    !************************************************
    !
    !           NONLINEAR ELEMENT-WISE VECTORS 
    !
    !************************************************
    type(nl_element), dimension(:), allocatable :: snl

    !************************************************
    !
    !           NODE-WISE VECTORS 
    !
    !************************************************
    
    integer*4,  dimension(:), allocatable   :: update_index_el_az
    real*8,     dimension(:), allocatable   :: u1,u2,vel,acc,u_predictor
    real*8,     dimension(:), allocatable   :: fk,fe,fd
    real*8,     dimension(:), allocatable   :: sism
    real*8, dimension(2*nnt), intent(in)    :: mvec
    real*8, dimension(2*nnt), intent(inout) :: u0,v1
    
    !************************************************
    !
    !           SEISMIC MOMENT 
    !
    !************************************************

    integer*4                               :: nl_sism,length_cns
    integer*4, dimension(length_cns,5)      :: check_node_sism
    real*8, dimension(length_cns,1)         :: check_dist_node_sism
    real*8, dimension(nl_sism,3)            :: facsmom       

    !************************************************
    !
    !           DAMPING 
    !
    !************************************************

    integer*4                               :: make_damping_yes_or_not
    real*8                                  :: ndt_monitor                                   
    integer*4, dimension(:), allocatable    :: nodal_counter

    !************************************************
    !
    !           NONLINEAR FUNCTIONS - INTERFACE 
    !
    !************************************************
    
    interface
        ! ALLOCATE/INITIALIZE ALL NL VARIABLES
        subroutine ALLOINIT_NL_ALL(ne,sdeg_mat,nm,nnt,cs_nnz,cs,prop_mat,u1,u2,vel,acc,v1,fk,fe,fd,&
            snl,option_out_var,disout,update_index_el_az,nodal_counter)  
            !
            use fields
            !
            implicit none
            ! intent IN
            integer*4,                              intent(in)      :: ne,nm,nnt,cs_nnz
            integer*4,  dimension(6),               intent(in)      :: option_out_var
            integer*4,  dimension(nm),              intent(in)      :: sdeg_mat
            integer*4,  dimension(0:cs_nnz),        intent(in)      :: cs
            real*8,     dimension(nm,9),            intent(in)      :: prop_mat
            real*8,     dimension(2*nnt),           intent(in)      :: v1
            ! intent INOUT
            real*8,     dimension(:), allocatable,  intent(inout)   :: u1,u2,vel,acc,fk,fe,fd
            integer*4,  dimension(:), allocatable,  intent(inout)   :: update_index_el_az,nodal_counter
            type(nl_element), dimension(:), allocatable, intent(inout)   :: snl
            type(nodepatched), intent(inout)                        :: disout
            ! counters
            integer*4                                               :: nn,im,iaz,ie,in,is,i,j
        end subroutine ALLOINIT_NL_ALL

        subroutine DEALLOCATE_ALL(ne,u1,u2,vel,acc,fk,fe,fd,snl,disout,update_index_el_az,nodal_counter)  
            ! 
            use fields 
            !
            implicit none
            ! intent IN
            integer*4, intent(in)                                      :: ne
            ! intent INOUT
            real*8,     dimension(:), allocatable,  intent(inout)      :: u1,u2,vel,acc,fk,fe,fd
            integer*4,  dimension(:), allocatable,  intent(inout)      :: update_index_el_az,nodal_counter
            type(nl_element), dimension(:), allocatable, intent(inout) :: snl
            type(nodepatched), intent(inout)                           :: disout
            ! counters
            integer*4                                                  :: ie 
        end subroutine DEALLOCATE_ALL 
    end interface
   
    !tagstep = 0
    pi = 4.d0*datan(1.d0)
    IN_TOT = IN_TOT -1

    ne = cs(0) -1
    eps = 1.0d3 * dabs(epsilon(mvec(1)))

    if (nf.gt.0) allocate(func_value(nf)) 

    if (tagstep.eq.2) then               !DRM Scandella 21.10.2005
        allocate(disp_PDRM_t(nf_drm,2))  !DRM Scandella 11.04.2006 
    endif                                !DRM Scandella 21.10.2005 
    !********************************************************************************************
    ! Set the number of points where Neumann load is applied (nnode_neuX,nnode_neuY)
    ! Set the point-id where Neumann load is applied (inode_neuX,inode_neuY)

    !if (nMDRM.eq.0 .or. tagstep.eq.3 ) then !! Kiana and Ali 
    allocate(i4count(nnt));       i4count = 0

    nnode_neuX = 0;       
    do in = 1,nnt
    id1 = in
        if (nf.gt.0) then
            r8t = 0.0d0
            do fn = 1,nf
                r8t = r8t + dabs(Fmat(fn,id1))
            enddo
            if (r8t.gt.eps) then
                nnode_neuX = nnode_neuX + 1
                i4count(in) = nnode_neuX
            endif
        endif
    enddo


    if (nnode_neuX.gt.0) then
        allocate(inode_neuX(nnode_neuX))
        do in = 1,nnt
            if (i4count(in).ne.0) then
                inode_neuX(i4count(in)) = in
            endif
        enddo
    endif

    nnode_neuY = 0;         i4count = 0

    do in = 1,nnt
        id2 = in + nnt

        if (nf.gt.0) then
            r8t = 0.0d0
            do fn = 1,nf
                r8t = r8t + dabs(Fmat(fn,id2))
            enddo
            if (r8t.gt.eps) then
                nnode_neuY = nnode_neuY +1
                i4count(in) = nnode_neuY
            endif
        endif
    enddo

    if (nnode_neuY.gt.0) then
        allocate(inode_neuY(nnode_neuY))
        do in = 1,nnt
            if (i4count(in).ne.0) then
                inode_neuY(i4count(in)) = in
            endif
        enddo
    endif

    !endif

    !********************************************************************************************
    ! Set the number of points where Dirichlet load is applied (nnode_dirX,nnode_dirY)
    ! Set the point-id where Dirichlet load is applied (inode_dirX,inode_dirY)

    !if (nMDRM.eq.0 .or. tagstep.eq.3 ) then !! Kiana and Ali

    nnode_dirX = 0
    i4count = 0

    call GET_EDGE_NODES(nnt,cs_nnz_bc,cs_bc,nl_dirX,tag_dirX,nnode_dirX,i4count)

    if (nnode_dirX.gt.0) then
        allocate(inode_dirX(nnode_dirX))
        do in = 1,nnt
            if (i4count(in).ne.0) then
                inode_dirX(i4count(in)) = in
            endif
        enddo
    endif

    nnode_dirY = 0;     i4count = 0

    call GET_EDGE_NODES(nnt,cs_nnz_bc,cs_bc,nl_dirY,tag_dirY,nnode_dirY,i4count)

    if (nnode_dirY .gt. 0) then
        allocate(inode_dirY(nnode_dirY))
        do in = 1,nnt
            if (i4count(in).ne.0) then
                inode_dirY(i4count(in)) = in
            endif
        enddo
    endif
    deallocate(i4count)

    ! endif


    !********************************************************************************************
    ! SET MONITOR FILES FOR OUTPUT
    !********************************************************************************************

    call MAKE_MONITOR_FILES(nnt,nmonit,option_out_var,nnode_TOT,tagstep,& 
        unit_disp,unit_vel,unit_acc,unit_stress,unit_strain,unit_omega,unit_uDRM)

    !********************************************************************************************
    ! INITIALIZATION
    !********************************************************************************************
    call ALLOINIT_NL_ALL(ne,sdeg_mat,nm,nnt,cs_nnz,cs,prop_mat,u1,u2,vel,acc,v1,fk,fe,fd,&
        snl,option_out_var,disout,update_index_el_az,nodal_counter)  
    dt2 = dt*dt
    number_of_threads = 1                                                !PARALLEL Kiana 06.10.2015
    call OMP_set_num_threads(number_of_threads)                          !PARALLEL Kiana 06.10.2015
    
    !********************************************************************************************
    !     FIRST STEP
    !********************************************************************************************	 

    if (tagstep.eq.2) then !DRM Scandella 28-10-05   

        !-----DRM----------------------------------------------------------------------------------
        ! Calculation of displacements at all DRM points at t=0

        call get_disp_valueDRM(nf_drm,func_type_drm,func_indx_drm,nfunc_data_drm,func_data_drm,&     !DRM Scandella 11.04.2006
            0.0d0,disp_PDRM_t)                                              !DRM Scandella 16.11.2005

        if (maxval(disp_PDRM_t).gt.0.d0) then                                                        !DRM Scandella 16.11.2005 
            write(*,'(A,E12.4)')'Error!!Initial condition not zero, but',maxval(disp_PDRM_t)         !DRM Scandella 16.11.2005
        endif                                                                                        !DRM Scandella 16.11.2005 

        !Initialization of Loukakis Forces at t=0

        do i = 1,nf_drm                     !DRM Scandella 12.04.2006
            px_t(i)=0.0d0                   !DRM Scandella 16.11.2005
            py_t(i)=0.0d0                   !DRM Scandella 16.11.2005
        enddo                               !DRM Scandella 16.11.2005 


        do fn = 1,nf_drm  
            id1= glob_x(fn,1)               !DRM Scandella 12.04.2006
            iaz = update_index_el_az(id1)   !DRM Scandella 12.04.2006 !change here "+1" 24.11.2015
            fe(iaz) = 0.0d0                 !DRM Scandella 12.04.2006 
            fe(iaz) = fe(iaz) + glob_x(fn,2) * px_t(fn) !DRM Scandella 12.04.2006
            id2 = glob_y(fn,1)              !DRM Scandella 12.04.2006
            iaz = update_index_el_az(id2)   !DRM Scandella 12.04.2006
            fe(iaz) = 0.0d0                 !DRM Scandella 12.04.2006 !change here "+1" 24.11.2015
            fe(iaz) = fe(iaz) + glob_y(fn,2) * py_t(fn)                   !DRM Scandella 12.04.2006 
        enddo                               !DRM Scandella 12.04.2006  		 

    else                                    !DRM Scandella 16.11.2005 
        !-----------------------------------------------------------------------------------------------

        do fn = 1,nf
            func_value(fn) = get_func_value(nf,func_type,func_indx,func_data,fn,0.0d0,0.0d0)
        enddo

        if (nnode_neuX.gt.0) then
            do i = 1,nnode_neuX
                in = inode_neuX(i)                 
                fe(in) = 0.0d0
                do fn = 1,nf
                    fe(in) = fe(in) + Fmat(fn,in) * func_value(fn)
                enddo
            enddo
        endif

        if (nnode_neuY.gt.0) then
            do i = 1,nnode_neuY
                in = inode_neuY(i) + nnt
                fe(in) = 0.0d0
                do fn = 1,nf
                    fe(in) = fe(in) + Fmat(fn,in) * func_value(fn)
                enddo
            enddo
        endif
    endif

    if (tagstep.eq.2) then 
        call get_disp_valueDRM(nf_drm,func_type_drm,func_indx_drm,nfunc_data_drm,func_data_drm, &     !DRM Scandella 12.04.2006  
            dt,disp_PDRM_t)                                                   !DRM Scandella 02.11.2005 
    else
        do fn = 1,nf
            func_value(fn) = get_func_value(nf,func_type,func_indx,func_data, &
                fn,dt,0.d0)
        enddo                                                                                          !DRM Scandella 02.11.2005 
    endif                                                                                            !DRM Scandella 02.11.2005  
    allocate(u_predictor(2*nnt))
    
    !********************************************************************************************
    !     ALL STEPS
    !********************************************************************************************
    v1 = 0.0d0
    do its = 0,nts
        fd = 0.d0 
        if (nl_sism.gt.0) sism=0.0d0

        call system_clock(COUNT=start,COUNT_RATE=clock(2))

        tt0 = dfloat(its -1) * dt
        tt1 = dfloat(its) * dt
        tt2 = dfloat(its +1) * dt

        write(*,'(A,F7.4)') 'TIME: ', tt1

        !-----DRM-----------------------------------------------------------------------------------------		       
        if (tagstep.eq.2) then     !DRM Scandella 28-10-05   

            ! Loukakis forces initialization
            do i = 1,nf_drm          !DRM Scandella 12.04.2006 
                px_t(i) = 0.0d0        !DRM Scandella 16.11.2005
                py_t(i) = 0.0d0        !DRM Scandella 16.11.2005 
            enddo                   !DRM Scandella 16.11.2005 

            call system_clock(COUNT=clock_start,COUNT_RATE=clock(2))  

            ! Selection of stiffness matrix for the current DRM element

            !$OMP PARALLEL &                                                                  
            !$OMP PRIVATE(i,ie,K_el)                                                          

            !$OMP DO 
            do i = 1,n_el_DRM                                                      !DRM Scandella 16.11.2005
                ie = el_DRM(i,1)                                                    !DRM Scandella 16.11.2005
                K_el = K_DRM(1+(i-1)*2*(ns+1)**2:i*2*(ns+1)**2,1:2*(ns+1)**2)       !DRM Scandella 16.11.2005

                ! Calculation of Loukakis forces
                call Louk_force(ns,K_el,nnode_BD,nnode_TOT,node_PDRM,&              !DRM Scandella 16.11.2005
                    node_TOT,nf_drm,disp_PDRM_t,cs_nnz,cs,ie,px_t,py_t) !DRM Scandella 16.11.2005 
            enddo                                                                   !DRM Scandella 16.11.2005

            !$OMP END DO
            !$OMP END PARALLEL		

            call system_clock(COUNT=clock_finish)
            time_force_DRM = float(clock_finish - clock_start) / float(clock(2))

            ! open(21,file=out_file)
            ! write(21,'(E12.4,E12.4,E12.4,E12.4)')(disp_PDRM_t(i,1),disp_PDRM_t(i,2),px_t(i),py_t(i),i =1,nnode_TOT)

            call system_clock(COUNT=clock_start,COUNT_RATE=clock(2))  

            do fn = 1,nf_drm  
                id1 = glob_x(fn,1)                                            !DRM Scandella 12.04.2006
                iaz = update_index_el_az(id1)                                 !DRM Scandella 12.04.2006 !change here "+1" 24.11.2015
                fe(iaz) = 0.0d0                                               !DRM Scandella 12.04.2006
                fe(iaz) = fe(iaz) + glob_x(fn,2) * px_t(fn)                   !DRM Scandella 12.04.2006

                id2 = glob_y(fn,1)                                            !DRM Scandella 12.04.2006 
                iaz = update_index_el_az(id2)                                 !DRM Scandella 12.04.2006 !change here "+1" 24.11.2015
                fe(iaz) = 0.0d0                                               !DRM Scandella 12.04.2006
                fe(iaz) = fe(iaz) + glob_y(fn,2) * py_t(fn)                   !DRM Scandella 12.04.2006
                !DRM Scandella 16.11.2005
            enddo                                                            !DRM Scandella 12.04.2006

            call system_clock(COUNT=clock_finish)
            time_fe = float(clock_finish - clock_start) / float(clock(2))

            fe = fe/mvec

                !-----------------------------------------------------------------------------------------------------
        else

            if (nnode_neuX.gt.0) then
                do i = 1,nnode_neuX
                    in = inode_neuX(i)                 
                    fe(in) = 0.0d0
                    do fn = 1,nf
                        fe(in) = fe(in) + Fmat(fn,in) * func_value(fn)
                    enddo
                enddo
            endif
!
            if (nnode_neuY.gt.0) then
                do i = 1,nnode_neuY
                    in = inode_neuY(i) + nnt
                    fe(in) = 0.0d0
                    do fn = 1,nf
                        fe(in) = fe(in) + Fmat(fn,in) * func_value(fn)
                    enddo
                enddo
            endif
        endif

        !********************************************************************************************
        ! LEAP-FROG        
        !********************************************************************************************
        !Damping
        !(u2 - 2*u1 + u0)/dt^2 + N_TOT*(u2-u0)/2dt= - Fint(u1) + Fel

        !NO Damping
        !(u2 - 2*u1 + u0)/dt^2  = - Fint(u1) + Fel
        ! COMPUTE NL INTERNAL FORCES fk(u)
        !write(*,'(A)')
        !write(*,'(A)') '*************************************************************'
        !write(*,'(A)') '----------COMPUTING NON LINEAR INTERNAL FORCES---------------'
        !
        !
        !write(*,'(A)') 'Non Linear Internal Forces: OK'
        ! COMPUTE EXTERNAL SEISMIC FORCES sism(fe) 
        if (nl_sism.gt.0) then 
            !write(*,'(A)')
            !write(*,'(A)') '*************************************************************'
            !write(*,'(A)') '----------COMPUTING SEISMIC EXTERNAL FORCES---------------'
            !
            call system_clock(COUNT=clock_start,COUNT_RATE=clock(2))   
            !
            call MAKE_SEISMIC_FORCES(nnt,nm,ne,nf,cs_nnz,cs,sdeg_mat,nfunc_data,nl_sism,length_cns,&
            check_node_sism,check_dist_node_sism,func_data,func_type,tag_func,func_indx,facsmom,tt1,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,u1,mvec,sism,fe)
            ! 
            call system_clock(COUNT=clock_finish)
            time_fe = time_fe+float(clock_finish - clock_start) / float(clock(2))
            !
            !write(*,'(A)') 'Seismic External Forces: OK'
        else
            !write(*,*) "NO EXTERNAL SEISMIC FORCES"
        endif
        
        !*******************************************************************************************
        !COMPUTE VISCID FORCES fd = N_TOT*v1
        !*******************************************************************************************
        
        call system_clock(COUNT=clock_start,COUNT_RATE=clock(2)) 
        call MATMUL_SPARSE(N_TOT, NNZ_N, JN_TOT, IN_TOT, fd, 2*nnt, v1, 2*nnt, error)
        call system_clock(COUNT=clock_finish)
        time_fd = float(clock_finish - clock_start) / float(clock(2))
        
        !*******************************************************************************************
        ! COMPUTE NONLINEAR INTERNAL FORCES
        !*******************************************************************************************
        call system_clock(COUNT=clock_start,COUNT_RATE=clock(2)) 
        
        
        u_predictor=v1*dt+acc*0.5d0*dt2
        
        fk(:) = 0.0d0
        call MAKE_INTERNAL_FORCES_NL(nnt,ne,nm,cs_nnz,cs,sdeg_mat,snl,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,u_predictor,fk,mvec,dt)
        call system_clock(COUNT=clock_finish)
        time_fk = float(clock_finish - clock_start) / float(clock(2))

        !*******************************************************************************************
        ! SOLVE SYSTEM
        !*******************************************************************************************

        call system_clock(COUNT=clock_start,COUNT_RATE=clock(2)) 
        u2 = 2.0d0 * u1 - u0 + dt2*(fe - fk - fd)
        call system_clock(COUNT=clock_finish)
        time_u = float(clock_finish - clock_start) / float(clock(2))
        ! 
        
        !-----DRM---------------------------------------------------------------------------------------------------        
        if (tagstep.eq.2) then                                                                      !DRM Scandella 28.10.2005 
            ! Calculation of displacements at all DRM points at a fixed time step                               !DRM Scandella 12.04.2006  

            call system_clock(COUNT=clock_start,COUNT_RATE=clock(2)) 

            call get_disp_valueDRM(nf_drm,func_type_drm,func_indx_drm,nfunc_data_drm,func_data_drm,&  !DRM Scandella 12.04.2006 
                tt2,disp_PDRM_t)                                                    !DRM Scandella 12.04.2006

            call system_clock(COUNT=clock_finish)
            time_disp_DRM = float(clock_finish - clock_start) / float(clock(2))

            !-----------------------------------------------------------------------------------------------------------                                 
        else
            do fn = 1,nf
            func_value(fn) = GET_FUNC_VALUE(nf,func_type,func_indx,func_data, &
                fn,tt2,0.0d0)
            enddo
        endif                                                                                       !DRM Scandella 28.10.2005 

        !********************************************************************************************
        ! DIRICHLET BOUNDARY CONDITIONS
        !********************************************************************************************

        if (nnode_dirX.gt.0) then
            do i = 1,nnode_dirX
                in = inode_dirX(i)
                u2(in) = 0.0d0;
                u1(in) = 0.d0; 
                v1(in) = 0.d0;
            enddo
        endif

        if (nnode_dirY.gt.0) then
            do i = 1,nnode_dirY
                in = inode_dirY(i) + nnt
                u1(in) = 0.0d0; 
                u2(in) = 0.d0; 
                v1(in) = 0.d0;
            enddo
        endif

        call system_clock(COUNT=finish)

        time_in_seconds = float(finish - start) / float(clock(2))
        time_total = time_total + time_in_seconds

        total_force_DRM = total_force_DRM + time_force_DRM
        total_disp_DRM = total_disp_DRM + time_disp_DRM
        total_fe = total_fe + time_fe
        total_fk = total_fk + time_fk
        total_fd = total_fd + time_fd
        total_u = total_u + time_u

        !********************************************************************************************
        !     WRITE OUTPUT FILE
        !********************************************************************************************
        condition = (nmonit.ge.1).and.(int(real(its)/ndt_monitor).eq.(real(its)/ndt_monitor)) 
        if (condition) then
            
            call WRITE_MONITOR_NL(unit_disp,unit_vel,unit_acc,unit_strain,unit_stress,unit_omega,&
                unit_uDRM,option_out_var,nmonit,ndt_monitor,node_m,nm,ne,nnt,cs,cs_nnz,sdeg_mat,snl,&
                its,tt1,node_TOT,nnode_TOT,tagstep,u1,vel,acc,nodal_counter,disout)
        endif


        !********************************************************************************************
        ! UPDATE
        !********************************************************************************************

        if(its .gt. 0) then
            vel = (u2 - u0) / (2*dt) 
            acc = (u2 - 2*u1 + u0) / (dt*dt)
            v1 = vel
            u0 = u1
            u1 = u2
        endif    

    enddo
    !write(*,'(A)') '*************************************************************'


    if (nmonit.ge.1) then
        call CLOSE_OUTPUT_FILES(option_out_var,nmonit,&
            unit_disp,unit_vel,unit_acc,unit_stress,unit_strain,unit_omega)
    endif
    !-----DRM---------------------------------------------------------------------------------------------------
    !DRM out files of I step closed

    if ((nnode_TOT.ne.0).and.(tagstep.eq.1)) then  !DRM Scandella 16.11.2005 
        if (nnode_TOT.ne.0) then                    !DRM Scandella 16.11.2005
            do i = 1,nnode_TOT                       !DRM Scandella 16.11.2005
                close(unit_uDRM(i))                   !DRM Scandella 16.11.2005
            enddo                                    !DRM Scandella 16.11.2005
        endif                                      !DRM Scandella 16.11.2005
    endif                                          !DRM Scandella 16.11.2005
    !-----------------------------------------------------------------------------------------------------------

    write(*,'(A,F16.8,A)')'total force_DRM time = ',total_force_DRM,' s' 
    write(*,'(A,F16.8,A)')'total disp_DRM time = ',total_disp_DRM,' s' 
    write(*,'(A,F16.8,A)')'total fe time = ',total_fe,' s' 
    write(*,'(A,F16.8,A)')'total fk time = ',total_fk,' s' 
    write(*,'(A,F16.8,A)')'total fd time = ',total_fd,' s' 
    write(*,'(A,F16.8,A)')'total u time = ',total_u,' s' 

    write(*,'(A,F16.8,A)')'time loop total time= ',time_total,' s'

    write(*,'(A,F8.4,A)')'Mean time-step CPU time= ', &
        time_total / dfloat(nts - 1),' s'

    call DEALLOCATE_ALL(ne,u1,u2,vel,acc,fk,fe,fd,snl,disout,update_index_el_az,nodal_counter)   
    if (nf.gt.0) deallocate(func_value) 
    if (nnode_dirX.gt.0) deallocate(inode_dirX)
    if (nnode_dirY.gt.0) deallocate(inode_dirY)
    if (nnode_neuX.gt.0) deallocate(inode_neuX)
    if (nnode_neuY.gt.0) deallocate(inode_neuY)

    return      
end subroutine TIME_LOOP_NL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
