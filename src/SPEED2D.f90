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


!> @brief SPEED (SPectral Elements in Elastodynamics with discontinuous Galerkin) 
!! is an open-source code for the simulation of seismic wave propagation in 
!! three-dimensional complex media. SPEED is jointly developed by MOX (The Laboratory for Modeling and Scientific 
!! Computing, Department of Mathematics) and DICA (Department of Civil and Environmental Engineering)
!! at Politecnico di Milano.
!> @see Website http://mox.polimi.it/speed/SPEED/Home.html
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0


! Here starts the code SPEED2D

program SPEED2D

    use speed_par_dg

    implicit none

    logical :: filefound

    !     INPUT FILES
    character*70 :: head_file,grid_file,mat_file,out_file

    !   SYSTEM CLOCK
    integer*4 :: COUNTCLOCK(1:1), COUNTRATE, COUNTMAX
    integer*4 :: START_TIME(1:8),END_TIME(1:8), OMP_GET_NUM_PROCS    
    integer*4 :: start, finish, clock_start,clock_finish
    integer*4, dimension(3)  :: clock

    !     OPTIONS OUTPUT      
    integer*4, dimension (6) :: opt_out_var

    !     TIME VARIABLES
    real*8 :: time_in_seconds, deltat, xtime, &
        deltat_cfl, fmax, fpeak, ndt_mon
    real*8 :: time_MAKE_EBW_MACRO, time_MAKE_EBIN_MACRO, &
        time_MAKE_SPECTRAL_CONNECTIVITY, time_MAKE_SPECTRAL_GRID, &
        time_MAKE_SPECTRAL_BOUNDARY, time_MAKE_MEL, &
        time_MAKE_CEL_KCEL, time_MAKE_FEL, &
        time_MAKE_PATTERN_STIFF_MATRIX, time_MAKE_STIFF_MATRIX, &
        time_MAKE_PATTERN_ABC_MATRIX, time_MAKE_ABCS_MATRIX, &
        time_SETUP_DG, time_SETUP_DG_ELEM, &
        time_MAKE_DG_INTERFACE, time_DELTAT_MAX


    !     MODEL VARIABLES
    real*8, dimension (:), allocatable :: xx_macro,yy_macro
    real*8, dimension (:), allocatable :: xx_spx,yy_spx
    integer*4 :: nnod_macro,nnod, nedge_abc, nelem_abc
    integer*4, dimension (:,:), allocatable :: con
    integer*4, dimension (:), allocatable :: con_spx
    integer*4 :: nelem,con_nnz
    integer*4, dimension (:,:), allocatable :: con_bc
    integer*4, dimension (:), allocatable :: con_spx_bc, i4count
    integer*4, dimension (:), allocatable :: iedge_abc, ielem_abc
    integer*4 ::nedge,con_nnz_bc
    integer*4, dimension (:), allocatable :: Ebw,Ebin,Nbw,Nbin
    integer*4 :: Ennz,Nnnz, iedge,ied1,ied2,iel1,iel2,iel3,iel4
    real*8, dimension (:), allocatable :: alfa1,beta1,gamma1,delta1
    real*8, dimension (:), allocatable :: alfa2,beta2,gamma2,delta2

    !     MECHANICAL VARIABLES
    integer*4 :: nmat
    integer*4, dimension (:), allocatable :: sdeg_mat
    integer*4, dimension(:), allocatable :: tag_mat
    real*8, dimension (:,:), allocatable :: prop_mat
    real*8, dimension (:), allocatable :: u0,v0,Mel,Cel,KCel
    real*8, dimension (:,:), allocatable :: Fel
    real*8, dimension (:,:), allocatable :: val_dirX_el,val_dirY_el
    real*8, dimension (:,:), allocatable :: val_neuX_el,val_neuY_el
    real*8, dimension (:,:), allocatable :: val_poiX_el,val_poiY_el
    real*8, dimension (:,:), allocatable :: val_plaX_el,val_plaY_el
    real*8, dimension (:,:), allocatable :: val_sism_el
    integer*4, dimension (:), allocatable :: fun_dirX_el,fun_neuX_el
    integer*4, dimension (:), allocatable :: fun_dirY_el,fun_neuY_el
    integer*4, dimension (:), allocatable :: fun_poiX_el,fun_poiY_el
    integer*4, dimension (:), allocatable :: fun_plaX_el,fun_plaY_el
    integer*4, dimension (:), allocatable :: fun_sism_el
    integer*4, dimension(:), allocatable :: tag_dirX_el,tag_neuX_el,tag_plaX_el
    integer*4, dimension(:), allocatable :: tag_dirY_el,tag_neuY_el,tag_plaY_el
    integer*4, dimension(:), allocatable :: tag_sism_el
    integer*4, dimension(:), allocatable :: tag_abc_el
    integer*4 :: nload_dirX_el,nload_dirY_el
    integer*4 :: nload_neuX_el,nload_neuY_el
    integer*4 :: nload_poiX_el,nload_poiY_el
    integer*4 :: nload_plaX_el,nload_plaY_el
    integer*4 :: nload_sism_el
    integer*4 :: nload_abc_el
    integer*4, dimension (:), allocatable :: fun_test
    integer*4, dimension (:), allocatable :: tag_func
    integer*4, dimension (:), allocatable :: func_type
    integer*4, dimension (:), allocatable :: func_indx
    real*8, dimension (:), allocatable :: func_data
    integer*4 :: nfunc,nfunc_data
    real*8, dimension (:), allocatable :: tsnap
    real*8, dimension (:), allocatable :: x_monitor,y_monitor
    integer*4, dimension (:),   allocatable :: itersnap
    integer*4, dimension (:),   allocatable :: n_monitor   
    integer*4 :: ns,nn,nn2,im,ie,i,j,in,ic,id,nts,nsnaps,nmonitors,trash
    real*8 :: eps
    integer*4 :: nnode_dom

    !     DRM VARIABLES

    integer*4 :: tagstep                                               !DRM Scandella 2005
    integer*4, dimension(:), allocatable :: tag_MDRM_el,tag_BDRM_el    !DRM Scandella 2005
    integer*4 :: nload_MDRM_el,n_el_DRM ,i_node_DRM                    !DRM Scandella 2005
    integer*4 :: nload_BDRM_el,n_boun_DRM,i_node_BDRM                  !DRM Scandella 2005 
    integer*4 :: nload_PDRM_el                                         !DRM Scandella 2005

    integer*4, dimension (:), allocatable :: node_PDRM_el              !DRM Scandella 2005
    real*8, dimension (:,:), allocatable :: val_PDRM_el                !DRM Scandella 2005 
    integer*4, dimension (:), allocatable :: fun_PDRM_el               !DRM Scandella 2005  

    integer*4, dimension(:,:), allocatable :: el_DRM_eff               !DRM Scandella 2005
    integer*4, dimension(:), allocatable :: node_BDRM_eff              !DRM Scandella 2005	  
    integer*4, dimension(:), allocatable :: node_DRM_eff               !DRM Scandella 2005

    integer*4 :: nnode_MDRM                                            !DRM Scandella 2005
    integer*4, dimension(:), allocatable :: node_MDRM                  !DRM Scandella 2005 

    integer*4 :: nnode_BD_eff,nnode_EL_eff,nnode_TOT_eff               !DRM Scandella 2005 
    integer*4, dimension(:), allocatable :: node_BD_eff                !DRM Scandella 2005
    integer*4, dimension(:), allocatable :: node_TOT_eff               !DRM Scandella 2005
    real*8, dimension(:), allocatable :: val_TOT_eff                   !DRM Scandella 12.04.2006
    real*8, dimension(:), allocatable :: xx_TOT_eff,yy_TOT_eff         !DRM Scandella 2005

    integer*4 :: nfunc_drm,nfunc_data_drm                              !DRM Scandella 11.04.2006
    integer*4, dimension (:), allocatable :: tag_func_drm              !DRM Scandella 11.04.2006
    integer*4, dimension (:), allocatable :: func_type_drm             !DRM Scandella 11.04.2006 
    integer*4, dimension (:), allocatable :: func_indx_drm             !DRM Scandella 11.04.2006
    real*8, dimension (:), allocatable :: func_data_drm                !DRM Scandella 11.04.2006 
    integer*4, dimension(:,:), allocatable :: glob_drm_x,glob_drm_y    !DRM Scandella 11.04.2006

    integer*4, dimension (:), allocatable :: fn_ord                    !DRM Scandella 2005

    real*8, dimension(:,:), allocatable ::  K_DRM                      !DRM Scandella 2005


    !     SEISMIC MOMENT VARIABLE
    real*8, dimension(:,:), allocatable :: facsmom
    integer*4, dimension(:), allocatable :: num_node_sism
    integer*4, dimension(:,:), allocatable :: sour_node_sism  
    real*8, dimension(:,:), allocatable :: dist_sour_node_sism
    integer*4, dimension(:,:), allocatable :: check_node_sism
    real*8, dimension(:,:), allocatable :: check_dist_node_sism
    integer*4 :: max_num_node_sism,length_check_node_sism
    integer*4 :: conta

    !   DAMPING MATRIX VARIABLES
    integer*4 :: make_damping_yes_or_not 

    !   MLST VARIABLES
    real*8 :: depth_search_mon_lst
    integer*4 ::num_lst, file_mon_lst      
    real*8, dimension(:), allocatable :: dist_monitor_lst, x_monitor_lst, y_monitor_lst
    character*70 :: file_LS, file_MLST

    !   TIME APPROXIMATION
    integer*4 :: test, n_test, time_degree
    real*8 :: pi  

    !   SPARSE MATRIX VARIABLES
    integer*4, parameter :: max_nodes    = 10000
    integer*4, dimension(:,:), allocatable :: STIFF_PATTERN, ABC_PATTERN     
    integer*4, dimension(:), allocatable :: I_STIFF, J_STIFF, I_ABC, J_ABC, J_MASS, I_MASS
    integer*4, dimension(:), allocatable :: NDEGR,IW, I_SUM, J_SUM, IN_SUM, JN_SUM,&
        IC_SUM, JC_SUM, ID_SUM, JD_SUM, IE_SUM, JE_SUM, &
        IN_TOT, JN_TOT, IK_TOT, JK_TOT, IK_MORSE, JK_MORSE, &
        IN_MORSE, JN_MORSE, I_FULL, J_FULL, IFULL_MORSE, JFULL_MORSE

    real*8, dimension(:), allocatable :: M_STIFF, M_ABC_U, M_ABC_V, C_MASS, D_MASS, M_MASS, C_SUM, &
        D_SUM, E_SUM, N_TOT, K_TOT, K_MORSE, N_MORSE, M_FULL, MFULL_MORSE

    integer*4 :: length, length_abc, NNZ_AB, NNZ_K, NNZ_N, ierr, NNZ_N_MORSE, NNZ_K_MORSE, NNZ_FULL, NNZ_FULL_MORSE

    !   DG VARIABLES
    character*70 :: file_face
    integer*4 :: nload_dg_el, nelem_dg, nnode_dg, nnz_dg_total, nnz_dg, &
        nnz_dg_total_only_uv, nnz_dg_only_uv

    integer*4, dimension(:), allocatable ::  tag_dg_el, tag_dg_yn
    integer*4, dimension (:,:), allocatable :: faces
    real*8 :: dg_const, dg_pen
    real*8, dimension (:,:), allocatable :: area_nodes

    !   SPARSE MATRIX DG VARIABLES    
    integer*4, dimension(:), allocatable :: IDG_TOTAL, JDG_TOTAL, IDG, JDG, IDG_MORSE, JDG_MORSE, &
        IDG_SUM, JDG_SUM, IDG_only_uv, JDG_only_uv
    real*8, dimension(:), allocatable :: MDG_TOTAL, MDG, MDG_MORSE, MDG_SUM, MDG_only_uv
    logical :: NLFLAG

    !*****************************************************************************************      
    !  START
    !*****************************************************************************************      
    call system_clock(COUNT=clock_start,COUNT_RATE=clock(2))

    write(*,'(A)')''
    write(*,'(A)')'*******************************************************'
    write(*,'(A)')'*                                                     *'
    write(*,'(A)')'*                        SPEED                        *'
    write(*,'(A)')'*          SPectral Elements in Elastodynamics        *'
    write(*,'(A)')'*              with Discontinuous Galerkin            *'
    write(*,'(A)')'*                                                     *'
    write(*,'(A)')'*            PoliMi, 2012, All Rights Reserved        *'
    write(*,'(A)')'*                                                     *'
    write(*,'(A)')'*******************************************************'
    write(*,'(A)')'*                                                     *'
    write(*,'(A)')'*                  2D - DG Space-Time SEM             *'  
    write(*,'(A)')'*                      Serial version                 *'
    write(*,'(A)')'*                                                     *'
    write(*,'(A)')'*******************************************************'


    !*****************************************************************************************      
    !  READ HEADER FILE
    !*****************************************************************************************      

    head_file = 'SPEED.input'

    write(*,'(A)') 
    write(*,'(A)')'*******************************************************'
    write(*,'(A)')'------------------READING HEADER FILE------------------'
    write(*,'(A,A36)') 'Header File: ',head_file

    inquire(file=head_file,exist=filefound); 
    if(filefound .eqv. .FALSE.) then
        read(*,*)
        stop
    endif 

    opt_out_var = 0   !do not write output files 


    call READ_DIME_HEADER(head_file,nsnaps)

    if (nsnaps.gt.0) allocate(tsnap(nsnaps),itersnap(nsnaps));

    time_degree = 0;
    call READ_HEADER(head_file, grid_file, mat_file, out_file,&
        deltat, xtime, opt_out_var, &
        nsnaps, tsnap, ndt_mon, depth_search_mon_lst, num_lst,&        
        file_mon_lst,time_degree, test, dg_const, dg_pen, NLFLAG)     

    !*****************************************************************************************      
    !  READ MATE FILE
    !*****************************************************************************************      

    write(*,'(A)')'HEADER FILE: OK'
    write(*,'(A)')'*******************************************************'

    mat_file = mat_file(1:len_trim(mat_file)) // '.mate'

    write(*,'(A)')    
    write(*,'(A)')'*******************************************************'
    write(*,'(A,A20)')'-----------------READING MATERIAL FILE-----------------'
    write(*,'(A,A20)')'Material File : ',mat_file

    inquire(file=mat_file,exist=filefound); 
    if(filefound .eqv. .FALSE.) then
        read(*,*);
        stop
    endif


    call READ_DIME_MAT_EL(mat_file,nmat, &
        nload_dirX_el,nload_dirY_el, &
        nload_neuX_el,nload_neuY_el, &
        nload_poiX_el,nload_poiY_el, &
        nload_plaX_el,nload_plaY_el, &
        nload_sism_el,nload_abc_el,nload_MDRM_el,  &
        nload_BDRM_el,nload_PDRM_el, &  !DRM Scandella 17.10.2005
        nfunc,nfunc_data,n_test,&
        nfunc_drm,nfunc_data_drm,&                     !DRM Scandella 11.04.2006 
        nload_dg_el)

    write(*,'(A,I8)')'Materials      : ',nmat
    write(*,'(A,I8)')'Dichlet X B.C. : ',nload_dirX_el
    write(*,'(A,I8)')'Dichlet Y B.C. : ',nload_dirY_el
    write(*,'(A,I8)')'Neumann X B.C. : ',nload_neuX_el
    write(*,'(A,I8)')'Neumann Y B.C. : ',nload_neuY_el
    write(*,'(A,I8)')'DG Interfaces  : ',nload_dg_el
    write(*,'(A,I8)')'Absorbing B.C. : ',nload_abc_el
    if (nload_MDRM_el.ne.0) then                                        ! DRM Scandella 10.05.2007  
        write(*,'(A,I8)')'Domain DRM     : ',nload_MDRM_el                  ! DRM Scandella 27.09.2005
        write(*,'(A,I8)')'Internal B. DRM: ',nload_BDRM_el                  ! DRM Scandella 27.09.2005
        write(*,'(A,I8)')'DRM loaded points: ',nload_PDRM_el                ! DRM Scandella 17.10.2005
    endif                                                               ! DRM Scandella 10.05.2007
    write(*,'(A,I8)')'Point Loads X  : ',nload_poiX_el
    write(*,'(A,I8)')'Point Loads Y  : ',nload_poiY_el
    write(*,'(A,I8)')'Plane Loads X  : ',nload_plaX_el
    write(*,'(A,I8)')'Plane Loads Y  : ',nload_plaY_el
    write(*,'(A,I8)')'Moment Loads  :  ',nload_sism_el
    if (nload_PDRM_el.ne.0) then                                        ! DRM Scandella 10.05.2007  
        write(*,'(A,I8)')'Functions DRM:      ',nfunc_drm                   ! DRM Scandella 11.04.2006 
    endif                                                               ! DRM Scandella 10.05.2007
    if (nmat.le.0) then
        write(*,*)'Error ! nmat = 0';   
        read(*,*)
        stop
    endif

    allocate (sdeg_mat(nmat), tag_mat(nmat))
    allocate(prop_mat(nmat,9))
    prop_mat=-1d0

    if (nload_dirX_el.gt.0) &
        allocate (val_dirX_el(nload_dirX_el,2),fun_dirX_el(nload_dirX_el),tag_dirX_el(nload_dirX_el))
    if (nload_dirY_el.gt.0) & 
        allocate (val_dirY_el(nload_dirY_el,2),fun_dirY_el(nload_dirY_el),tag_dirY_el(nload_dirY_el))

    if (nload_neuX_el.gt.0) &
        allocate (val_neuX_el(nload_neuX_el,2),fun_neuX_el(nload_neuX_el),tag_neuX_el(nload_neuX_el))     
    if (nload_neuY_el.gt.0) &
        allocate (val_neuY_el(nload_neuY_el,2),fun_neuY_el(nload_neuY_el),tag_neuY_el(nload_neuY_el))

    if (nload_poiX_el.gt.0) allocate (val_poiX_el(nload_poiX_el,3), fun_poiX_el(nload_poiX_el))
    if (nload_poiY_el.gt.0) allocate (val_poiY_el(nload_poiY_el,3), fun_poiY_el(nload_poiY_el))

    if (nload_plaX_el.gt.0) &
        allocate (val_plaX_el(nload_plaX_el,1),fun_plaX_el(nload_plaX_el),tag_plaX_el(nload_plaX_el))
    if (nload_plaY_el.gt.0) &
        allocate (val_plaY_el(nload_plaY_el,1),fun_plaY_el(nload_plaY_el),tag_plaY_el(nload_plaY_el))


    if (nload_sism_el.gt.0) &
        allocate (val_sism_el(nload_sism_el,12),fun_sism_el(nload_sism_el),tag_sism_el(nload_sism_el))  

    if (nload_abc_el.gt.0) allocate (tag_abc_el(nload_abc_el))
    if (nload_dg_el.gt.0) allocate (tag_dg_el(nload_dg_el), tag_dg_yn(nload_dg_el))      

    if (nfunc.gt.0) &
        allocate (tag_func(nfunc),func_type(nfunc),func_indx(nfunc +1),func_data(nfunc_data))

    !----DRM---------------------------------------------------------------------------
    if (nload_MDRM_el.gt.0) then                      ! DRM Scandella 27.09.2005
        allocate (tag_MDRM_el(nload_MDRM_el))          ! DRM Scandella 27.09.2005
    endif                                             ! DRM Scandella 27.09.2005

    if (nload_BDRM_el.gt.0) then                      ! DRM Scandella 27.09.2005 
        allocate (tag_BDRM_el(nload_BDRM_el))          ! DRM Scandella 27.09.2005
    endif                                             ! DRM Scandella 27.09.2005

    if (nload_PDRM_el.gt.0) then                      ! DRM Scandella 27.09.2005
        allocate (val_PDRM_el(nload_PDRM_el,3))        ! DRM Scandella 20.10.2005 
        allocate (node_PDRM_el(nload_PDRM_el))         ! DRM Scandella 16.11.2005 
        allocate (fun_PDRM_el(nload_PDRM_el))          ! DRM Scandella 20.10.2005 
        allocate (fn_ord(nload_PDRM_el))               ! DRM Scandella 16.11.2005
        allocate (glob_drm_x(nload_PDRM_el,2))         ! DRM Scandella 11.04.2006
        allocate (glob_drm_y(nload_PDRM_el,2))         ! DRM Scandella 11.04.2006
    endif                                             ! DRM Scandella 16.11.2005

    if (nfunc_drm.gt.0) then                          ! DRM Scandella 11.04.2006 
        allocate (tag_func_drm(nfunc_drm))             ! DRM Scandella 11.04.2006   
        allocate (func_type_drm(nfunc_drm))            ! DRM Scandella 11.04.2006 
        allocate (func_indx_drm(nfunc_drm +1))         ! DRM Scandella 11.04.2006 
        allocate (func_data_drm(nfunc_data_drm))       ! DRM Scandella 11.04.2006 
    endif                                             ! DRM Scandella 11.04.2006 

    if (n_test.gt.0) allocate (fun_test(n_test))


    call READ_MATERIAL_EL(mat_file,nmat,prop_mat,sdeg_mat,tag_mat,&
        nload_dirX_el,val_dirX_el,fun_dirX_el,tag_dirX_el, &
        nload_dirY_el,val_dirY_el,fun_dirY_el,tag_dirY_el, &
        nload_neuX_el,val_neuX_el,fun_neuX_el,tag_neuX_el, &
        nload_neuY_el,val_neuY_el,fun_neuY_el,tag_neuY_el, &
        nload_poiX_el,val_poiX_el,fun_poiX_el, &
        nload_poiY_el,val_poiY_el,fun_poiY_el, &
        nload_plaX_el,val_plaX_el,fun_plaX_el,tag_plaX_el, &
        nload_plaY_el,val_plaY_el,fun_plaY_el,tag_plaY_el, &
        nload_sism_el,val_sism_el,fun_sism_el,tag_sism_el, &
        nload_abc_el,tag_abc_el, &
        nload_MDRM_el,nload_BDRM_el,tag_MDRM_el,tag_BDRM_el,tagstep, &     !DRM Scandella 17.10.2005
        nload_PDRM_el,val_PDRM_el,fun_PDRM_el,&                            !DRM Scandella 20.10.2005
        nfunc,func_type,func_indx,func_data,tag_func, &
        nfunc_drm,func_type_drm,func_indx_drm,func_data_drm,tag_func_drm,& !DRM Scandella 11.04.2006
        fmax, n_test, fun_test, &
        nload_dg_el,tag_dg_el,tag_dg_yn,NLFLAG)

    do im = 1,nmat
        write(*,'(A,I8)')    'MATERIAL : ',tag_mat(im)
        write(*,'(A,I8)')    'DEGREE   : ',sdeg_mat(im)
        write(*,'(A,E12.4)') 'rho      : ',prop_mat(im,1)
        write(*,'(A,E12.4)') 'Vp       : ',((prop_mat(im,2) + 2*prop_mat(im,3))/prop_mat(im,1))**0.5
        write(*,'(A,E12.4)') 'Vs       : ',(prop_mat(im,3)/prop_mat(im,1))**0.5
        write(*,'(A,E12.4)') 'gamma    : ',prop_mat(im,4)
        if (NLFLAG) then
            write(*,'(A,E12.4)') 'sigma_yld    : ',prop_mat(im,5)
            write(*,'(A,E12.4)') 'Ckin         : ',prop_mat(im,6)
            write(*,'(A,E12.4)') 'kkin         : ',prop_mat(im,7)  
            write(*,'(A,E12.4)') 'Rinf         : ',prop_mat(im,8)
            write(*,'(A,E12.4)') 'biso         : ',prop_mat(im,9)  
        endif
    enddo
    write(*,'(A)') 'MATERIAL FILE:OK'
    write(*,'(A)') '****************************************************'

    !----DRM-----------------------------------------------------------------------------------------
    if (nload_MDRM_el.ne.0) then                                        !DRM Scandella 10.05.2007
        do im = 1,nmat                                                  !DRM Kiana 19.10.2015
            if (im.eq.tag_MDRM_el(1)) then                              !DRM Kiana 19.10.2015
            ns = sdeg_mat(im)                                           !DRM Kiana 19.10.2015
            endif                                                       !DRM Kiana 19.10.2015
        enddo  
        write(*,'(A,I8)') 'DRMSTEP : ',tagstep                          !DRM Scandella 10.05.2007
        if ((tagstep.eq.2).and.(nfunc_drm.eq.0)) then                   !DRM Scandella 11.05.2007
            write(*,'(A)')'ERROR: Not equivalent DRM forces applied!'   !DRM Scandella 11.05.2007
            read(*,*)
            stop                                                        !DRM Scandella 11.05.2007
        endif                                                           !DRM Scandella 11.05.2007
        if ((tagstep.eq.2).and.(nload_PDRM_el.eq.0)) then               !DRM Scandella 11.05.2007
            write(*,'(A)')'ERROR: Not equivalent DRM points defined!'   !DRM Scandella 11.05.2007
            read(*,*)
            stop                                                        !DRM Scandella 11.05.2007
        endif                                                           !DRM Scandella 11.05.2007
        write(*,*)                                                      !DRM Scandella 10.05.2007 
        do j = 1,nload_MDRM_el                                          !DRM Scandella 27.09.2005
            write(*,*)                                                  !DRM Scandella 27.09.2005 
            write(*,'(A,I8)') 'DRM BLOCKS : ',tag_MDRM_el(j)            !DRM Scandella 27.09.2005
        enddo                                                           !DRM Scandella 27.09.2005
        write(*,*)                                                      !DRM Scandella 27.09.2005  
        do j = 1,nload_BDRM_el                                          !DRM Scandella 27.09.2005
            write(*,*)                                                  !DRM Scandella 27.09.2005 
            write(*,'(A,I8)') 'DRM BOUNDARIES : ',tag_BDRM_el(j)        !DRM Scandella 27.09.2005 
        enddo                                                           !DRM Scandella 27.09.2005
        write(*,*)                                                      !DRM Scandella 27.09.2005  
        write(*,'(A)') '****************************************************'
    endif                                                               !DRM Scandella 10.05.2007 

    !----------------------------------------------------------------------------------

    !*****************************************************************************************      
    !  READ GRID FILE
    !*****************************************************************************************      

    grid_file = grid_file(1:len_trim(grid_file)) // '.mesh'    
    write(*,'(A)') ''
    write(*,'(A)') '*******************************************************'
    write(*,'(A)') '-------------------READING MESH FILE-------------------'
    write(*,'(A,A20)') 'Grid File : ',grid_file

    inquire(file=grid_file,exist=filefound); 

    if (filefound .eqv. .FALSE.) then
        read(*,*)
        stop
    endif

    call READ_DIME_GRID_EL(grid_file,nmat,tag_mat,&
        nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&
        nload_neuX_el,tag_neuX_el,nload_neuY_el,tag_neuY_el,&
        nload_abc_el,tag_abc_el,&
        nload_BDRM_el,tag_BDRM_el,&   !DRM Scandella 27.09.2005
        nload_dg_el,tag_dg_el, &
        nnod_macro,nelem,nedge)

    write(*,'(A,I8)')'Nodes : ',nnod_macro
    write(*,'(A,I8)')'Elements : ',nelem
    write(*,'(A,I8)')'Edges : ',nedge

    if (nnod_macro.gt.0) then
        allocate (xx_macro(nnod_macro),yy_macro(nnod_macro))
    else
        write(*,*)'Error ! Vertex number = 0'; 
        read(*,*);
        stop
    endif

    if (nelem.gt.0) then
        allocate (con(nelem,5))
    else
        write(*,*)'Error ! Element number = 0'; 
        read(*,*);
        stop
    endif

    if (nedge.gt.0) allocate (con_bc(nedge,3))

    call READ_GRID_EL(grid_file,nmat,tag_mat,prop_mat,       &
        nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el, &
        nload_neuX_el,tag_neuX_el,nload_neuY_el,tag_neuY_el, &
        nload_abc_el,tag_abc_el,nload_BDRM_el,tag_BDRM_el,   &
        nload_dg_el,tag_dg_el,nnod_macro,xx_macro,yy_macro,  &
        nelem,con,nedge,con_bc)

    write(*,'(A)')'MESH FILE: OK'    
    write(*,'(A)') '****************************************************'
    write(*,'(A)')

    !*****************************************************************************************      
    !  MAKING SPECTRAL CONNECTIVITY
    !*****************************************************************************************      
    write(*,'(A)') '*******************************************************'
    write(*,'(A)') '-------------MAKING SPECTRAL CONNECTIVITIES------------'

    call system_clock(COUNT=start,COUNT_RATE=clock(2))

    allocate(Ebw(nnod_macro))    
    call MAKE_EBW_MACRO(nnod_macro,nelem,con,Ebw,Ennz)

    call system_clock(COUNT=finish)
    time_MAKE_EBW_MACRO = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'MAKE_EBW_MACRO time = ',time_MAKE_EBW_MACRO,' s'

    call system_clock(COUNT=start,COUNT_RATE=clock(2))

    allocate(Ebin(0:Ennz))
    call MAKE_EBIN_MACRO(nnod_macro,nelem,con,Ebw,Ennz,Ebin)

    call system_clock(COUNT=finish)
    time_MAKE_EBIN_MACRO = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'MAKE_EBIN_MACRO time = ',time_MAKE_EBIN_MACRO,' s'

    deallocate(Ebw)

    con_nnz = nelem +1
    do ie = 1,nelem
        do j = 1,nmat
            if (tag_mat(j).eq.con(ie,1)) nn = sdeg_mat(j) +1
        enddo
        con_nnz = con_nnz + nn*nn +1
    enddo

    call system_clock(COUNT=start,COUNT_RATE=clock(2))
    allocate(con_spx(0:con_nnz))
    call MAKE_SPECTRAL_CONNECTIVITY(nelem,con,nmat,tag_mat,sdeg_mat,&
        Ennz,Ebin,con_nnz,con_spx,nnod)
    call system_clock(COUNT=finish)
    time_MAKE_SPECTRAL_CONNECTIVITY = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'MAKE_SPECTRAL_CONNECTIVITY time = ',time_MAKE_SPECTRAL_CONNECTIVITY,' s'


    deallocate(Ebin)

    ! Dual connectivity for spectral nodes

    allocate(Ebw(nnod))
    call MAKE_EBW(nnod,con_nnz,con_spx,Ebw,Ennz)

    allocate(Ebin(0:Ennz))
    call MAKE_EBIN(nnod,con_nnz,con_spx,Ebw,Ennz,Ebin)

    deallocate(Ebw)

    allocate(xx_spx(nnod),yy_spx(nnod))
    allocate(alfa1(nelem),beta1(nelem),gamma1(nelem),delta1(nelem))
    allocate(alfa2(nelem),beta2(nelem),gamma2(nelem),delta2(nelem))

    write(*,'(A,I8)')'Spectral Nodes : ',nnod

    call system_clock(COUNT=start,COUNT_RATE=clock(2))

    call MAKE_SPECTRAL_GRID(nnod_macro,xx_macro,yy_macro,con_nnz,con_spx,&
        nmat,tag_mat,sdeg_mat,nelem,&
        alfa1,beta1,gamma1,delta1,&
        alfa2,beta2,gamma2,delta2,&
        nnod,xx_spx,yy_spx)

    call system_clock(COUNT=finish)
    time_MAKE_SPECTRAL_GRID = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'MAKE_SPECTRAL_GRID time = ',time_MAKE_SPECTRAL_GRID,' s'


    !     Make the spectral connectivities for the boundary

    con_nnz_bc = 0

    call system_clock(COUNT=start,COUNT_RATE=clock(2))

    if (nedge.gt.0) then
        con_nnz_bc = nedge +1
        do i = 1,nedge
            call GET_EDGE_ELEMENT(Ennz,Ebin,con_bc(i,2),con_bc(i,3),ie)
            do j = 1,nmat
                if (tag_mat(j).eq.con(ie,1)) nn = sdeg_mat(j) +1
            enddo
            con_nnz_bc = con_nnz_bc +nn +1
        enddo

        allocate(con_spx_bc(0:con_nnz_bc))

        call MAKE_SPECTRAL_BOUNDARY(con_nnz,con_spx,nedge,con_bc,&
            nmat,tag_mat,sdeg_mat,Ennz,Ebin,&
            con_nnz_bc,con_spx_bc)
    endif

    write(*,'(A)')'SPECTRAL CONNECTIVITY: OK'

    call system_clock(COUNT=finish)
    time_MAKE_SPECTRAL_BOUNDARY = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'MAKE_SPECTRAL_BOUNDARY time = ',time_MAKE_SPECTRAL_BOUNDARY,' s'
    write(*,'(A)') '****************************************************'

    nnode_dom = nnod

    !*****************************************************************************************      
    !  MAKE MASS MATRIX
    !*****************************************************************************************      

    write(*,'(A)')
    write(*,'(A)') '****************************************************'
    write(*,'(A)') '--------------COMPUTING MASS MATRIX-----------------'
    write(*,'(A)')
    
    call system_clock(COUNT=start,COUNT_RATE=clock(2))

    allocate (Mel(2*nnod))
    call MAKE_MEL(nnod, con_nnz,con_spx,&
        nmat,tag_mat,sdeg_mat,prop_mat,&
        nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,Mel)

    write(*,'(A)') 'Mass matrix: OK'

    call system_clock(COUNT=finish)
    time_MAKE_MEL = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'MAKE_MEL time = ',time_MAKE_MEL,' s'
    write(*,'(A)') '****************************************************'


    !*****************************************************************************************      
    !  MAKE DAMPING MATRIX
    !*****************************************************************************************      

    !Check if there is any damping factor between materials characteristics    
    make_damping_yes_or_not = 0

    do im = 1,nmat
    if (abs(prop_mat(im,4)) .gt. 10e-10) make_damping_yes_or_not = 1
    enddo

    if (make_damping_yes_or_not.eq.1) then

        write(*,'(A)') '****************************************************'
        write(*,'(A)') '-------------COMPUTING DAMPING MATRIX---------------'

        call system_clock(COUNT=start,COUNT_RATE=clock(2))

        allocate (Cel(2*nnod),KCel(2*nnod))

        call MAKE_CEL_KCEL(nnod,con_nnz,con_spx,&
            nmat,tag_mat,sdeg_mat,prop_mat,&
            nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
            Cel,KCel)

        write(*,'(A)') 'Damping matrix: OK'

        call system_clock(COUNT=finish)
        time_MAKE_CEL_KCEL = float(finish-start)/float(clock(2))
        write(*,'(A,F8.4,A)')'MAKE_CEL_KCEL time = ',time_MAKE_CEL_KCEL,' s'
        write(*,'(A)') '****************************************************'

    else

        write(*,'(A)') '****************************************************'
        write(*,'(A)') 'NO Damping Matrix!!!'
        write(*,'(A)') 'There are no materials with damping defined on'
        write(*,'(A)') '****************************************************'

    endif

    !*****************************************************************************************      
    !  MAKE THE LOAD MATRIX --- EXTERNAL LOADS --- SEISmIC MOMENT
    !*****************************************************************************************      
    write(*,'(A)') 
    write(*,'(A)') '****************************************************'

    ! Dimensioning vector 'num_node_sism'(nodes number generating each single fault)
    if (nload_sism_el.gt.0) then

        write(*,'(A)') '---------------COMPUTING SEISMIC MOMENT-----------------'

        allocate (num_node_sism(nload_sism_el))

        do i = 1,nload_sism_el
        if (((val_sism_el(i,1).eq.val_sism_el(i,3)).and.(val_sism_el(i,3).eq.val_sism_el(i,5))) &
            .and.((val_sism_el(i,1).eq.val_sism_el(i,3)).and.(val_sism_el(i,3).eq.val_sism_el(i,5))))  then

            num_node_sism(i)=1
        else  
            call DIME_SISM_NODES(val_sism_el(i,1),val_sism_el(i,2),val_sism_el(i,3),val_sism_el(i,4),&
                val_sism_el(i,5),val_sism_el(i,6),val_sism_el(i,7),val_sism_el(i,8),&
                nnod,xx_spx,yy_spx,num_node_sism(i))
        endif
        enddo

        !Checking the maximum number of fault nodes
        max_num_node_sism = num_node_sism(1)
        do i = 1,nload_sism_el
        if (num_node_sism(i).gt.max_num_node_sism) then
            max_num_node_sism = num_node_sism(i)
        endif
        enddo

        allocate (sour_node_sism(max_num_node_sism,nload_sism_el))
        allocate (dist_sour_node_sism(max_num_node_sism,nload_sism_el))

        !Searching the node 'id' in the global numeration for each fault.
        !sour_node_sism = node id (global numeration) generating the fault 'i'
        do i = 1,nload_sism_el
        if (num_node_sism(i).eq.1) then

            call FIND_NEAREST_NODE(nnod,xx_spx,yy_spx,val_sism_el(i,1),val_sism_el(i,2),sour_node_sism(1,i))
            dist_sour_node_sism(1,i) = 0

        else 
            call READ_SISM_NODES(val_sism_el(i,1),val_sism_el(i,2),val_sism_el(i,3),val_sism_el(i,4),&
                val_sism_el(i,5),val_sism_el(i,6),val_sism_el(i,7),val_sism_el(i,8),&
                nnod,xx_spx,yy_spx,num_node_sism(i),sour_node_sism,i,&
                dist_sour_node_sism,nload_sism_el,&
                max_num_node_sism)
        endif

        if (i.eq.1) then
            write(*,'(A)')'Sesmic moment & fault'
            write(*,'(A)')
        endif
        if (num_node_sism(i).eq.1) then
            write(*,'(A,I6,A)')'Sism ',i,' is located on:'
            write(*,'(I6,2E14.5,I6)')(j,xx_spx(sour_node_sism(j,i)),&
                yy_spx(sour_node_sism(j,i)),sour_node_sism(j,i),&
                j=1,num_node_sism(i))
        else
            write(*,'(A,I6,A,I6,A)')'Fault ',i,' is generated by ',num_node_sism(i),' nodes'
            write(*,'(I6,2E14.5,I6)')(j,xx_spx(sour_node_sism(j,i)),&
                yy_spx(sour_node_sism(j,i)),sour_node_sism(j,i),&
                j=1,num_node_sism(i)) 
        endif
        write(*,'(A)')
        enddo

        write(*,'(A)')
        write(*,'(A)') 'SEISMIC MOMENT:OK'
        write(*,'(A)') '****************************************************'
    endif

    !---DRM---------------------------------------------------------------------------------------------

    !     DRM effective elements and nodes 
    if (nload_MDRM_el.ne.0) then

        !     Number of DRM elements (n_el_DRM)                                       !DRM Scandella 16.11.2005 
        call DIME_DRM_EL(ns,nelem,con,con_nnz,con_spx,&                         !DRM Scandella 16.11.2005 
            nload_MDRM_el,tag_MDRM_el,n_el_DRM,i_node_DRM)         !DRM Scandella 16.11.2005 

        write(*,'(A)')                                                          !DRM Scandella 16.11.2005 
        write(*,'(A,I8)')'Number of DRM elements =',n_el_DRM                    !DRM Scandella 16.11.2005 
        write(*,'(A)')                                                          !DRM Scandella 16.11.2005 

        allocate(el_DRM_eff(n_el_DRM,5))                                        !DRM Scandella 16.11.2005 
        allocate(node_DRM_eff(i_node_DRM))                                      !DRM Scandella 16.11.2005

        !     DRM elements and their vertex (el_DRM_eff)                              !DRM Scandella 15-03-2006 
        open(21,file=out_file)                                                  !DRM Scandella 15-03-2006 
        write(21,'(A,I6)')'# Number of DRM domain elements =', n_el_DRM         !DRM Scandella 16.11.2005 
        write(21,'(A)')'# DRM domain elements connectivity'  
        call EXTRACT_DRM_EL(ns,nelem,con,con_nnz,con_spx,&                      !DRM Scandella 16.11.2005 
            nload_MDRM_el,tag_MDRM_el,n_el_DRM,i_node_DRM, &    !DRM Scandella 16.11.2005 
            el_DRM_eff,node_DRM_eff)                            !DRM Scandella 16.11.2005 
        if (tagstep.eq.3) then 
            call EXTRACT_MACRO_MICRO_LISTS(ns,nelem,con_nnz,con_spx,&            !DRM Scandella 09.05.2007 
                nload_MDRM_el,tag_MDRM_el,n_el_DRM, &               !DRM Scandella 09.05.2007 
                nnode_dom,xx_spx,yy_spx)                            !DRM Scandella 09.05.2007

            write(*,'(A)')
            write(*,'(A)')'Macro and micro node lists created.'      
            write(*,'(A)')'**************************************'
            write(*,'(A)')
            write(*,'(A)')'Bye.'
            read(*,*)
            stop 

        endif


        !     Number of DRM boundary lines and nodes (n_boun_DRM)                     !DRM Scandella 16.11.2005 
        call DIME_DRM_BOUN(ns,con_bc,nedge,con_nnz_bc,con_spx_bc,&              !DRM Scandella 16.11.2005 
            nload_BDRM_el,tag_BDRM_el,n_boun_DRM,i_node_BDRM)    !DRM Scandella 16.11.2005 

        write(*,'(A)')                                                          !DRM Scandella 16.11.2005 
        write(*,'(A,I8)')'Number of DRM boundaries =',n_boun_DRM                !DRM Scandella 16.11.2005 
        write(*,'(A,I8)')'tag DRM boundaries =',tag_BDRM_el                     !DRM Scandella 16.11.2005 
        write(*,'(A)')

        allocate(node_BDRM_eff(i_node_BDRM))                                    !DRM Scandella 16.11.2005 

        !     Number of DRM boundary nodes without duplicates (nnode_BD_eff)          !DRM Scandella 16.11.2005 
        call DIME_BDRM_NODE(ns,nedge,con_bc,con_nnz_bc,con_spx_bc, &            !DRM Scandella 16.11.2005 
            nload_BDRM_el,tag_BDRM_el,n_boun_DRM, &             !DRM Scandella 16.11.2005 
            node_BDRM_eff,i_node_BDRM,nnode_BD_eff)             !DRM Scandella 16.11.2005 

        allocate (node_BD_eff(nnode_BD_eff))                                    !DRM Scandella 16.11.2005 

        !     DRM boundary nodes without duplicates (node_BD_eff)                     !DRM Scandella 16.11.2005 
        call EXTRACT_BDRM_NODE(node_BDRM_eff,i_node_BDRM,&                      !DRM Scandella 16.11.2005 
            nnode_BD_eff,node_BD_eff)                        !DRM Scandella 16.11.2005 

        !     Number of DRM elements nodes without DRM boundary nodes (nnode_MDRM)    !DRM Scandella 16.11.2005 
        call DIME_MDRM_NODE(i_node_DRM,node_DRM_eff,nnode_BD_eff,node_BD_eff,&  !DRM Scandella 16.11.2005 
            nnode_MDRM)                                         !DRM Scandella 16.11.2005 

        allocate (node_MDRM(nnode_MDRM))                                        !DRM Scandella 16.11.2005 

        !     DRM elements nodes without DRM boundary nodes (node_MDRM)               !DRM Scandella 16.11.2005 
        call MDRM_NODE(i_node_DRM,node_DRM_eff,nnode_BD_eff,&                   !DRM Scandella 16.11.2005 
            node_BD_eff,nnode_MDRM,node_MDRM)              !DRM Scandella 16.11.2005 


        !     Number of DRM elements nodes without DRM boundary nodes and duplicates (nnode_EL_eff) !DRM Scandella 16.11.2005 
        call DIME_EL_DRM_NODE(nnode_MDRM,node_MDRM,nnode_EL_eff)                !DRM Scandella 16.11.2005 

        nnode_TOT_eff = nnode_BD_eff+nnode_EL_eff                               !DRM Scandella 16.11.2005 

        allocate (node_TOT_eff(nnode_TOT_eff))                                  !DRM Scandella 16.11.2005 
        allocate (val_TOT_eff(nnode_TOT_eff))                                   !DRM Scandella 12.04.2006 
        allocate (xx_TOT_eff(nnode_TOT_eff))                                    !DRM Scandella 16.11.2005 
        allocate (yy_TOT_eff(nnode_TOT_eff))                                    !DRM Scandella 16.11.2005 



        !     DRM elements nodes without DRM boundary nodes and duplicates (node_TOT_eff)         !DRM Scandella 16.11.2005 
        call EXTRACT_TOT_DRM_NODE(nnode_dom,nnode_MDRM,node_MDRM,nnode_BD_eff,node_BD_eff,& !DRM Scandella 16.11.2005 
            xx_spx,yy_spx,&                                           !DRM Scandella 16.11.2005 
            nnode_EL_eff,nnode_TOT_eff,node_TOT_eff,&                 !DRM Scandella 16.11.2005 
            xx_TOT_eff,yy_TOT_eff)                                    !DRM Scandella 16.11.2005 


        write(21,'(A,I6)')'# Number of DRM boundary lines =', n_boun_DRM                    !DRM Scandella 16.11.2005 
        write(21,'(A,I6)')'# Number of DRM element nodes =', nnode_EL_eff                   !DRM Scandella 16.11.2005 
        write(21,'(A,I6)')'# Number of DRM internal boundary nodes =', nnode_BD_eff         !DRM Scandella 16.11.2005 
        write(21,'(A,I6)')'# Number of DRM total nodes =', nnode_TOT_eff                    !DRM Scandella 16.11.2005 
        write(21,'(A)')'# DRM domain nodes'                                                 !DRM Scandella 16.11.2005 
        !     write(21,'(I6,2f12.4)')(node_TOT_eff(i),xx_TOT_eff(i),yy_TOT_eff(i),i=1,nnode_TOT_eff)      !DRM Scandella 20.02.2006 
        write(21,'(A,I6,3f12.4)')('PDRM  ',i, xx_TOT_eff(i), yy_TOT_eff(i),1.0,i=1,nnode_TOT_eff)   !DRM Scandella 20.02.2006 
        write(21,'(A)')'  '                                                                         !DRM Scandella 16.11.2005 



        write(*,'(A)')                                                                      !DRM Scandella 16.11.2005 
        write(*,'(A,I8)')'Number of DRM element nodes =',nnode_EL_eff                       !DRM Scandella 16.11.2005 
        write(*,'(A)')                                                                      !DRM Scandella 16.11.2005 
        write(*,'(A,I8)')'Number of DRM internal boundary nodes =',nnode_BD_eff             !DRM Scandella 16.11.2005 
        write(*,'(A)')                                                                      !DRM Scandella 16.11.2005 
        write(*,'(A,I8)')'Number of DRM total nodes =',nnode_TOT_eff                        !DRM Scandella 16.11.2005 


        deallocate (node_DRM_eff)                                        !DRM Scandella 19-04-2006 
        deallocate (node_BDRM_eff)                                       !DRM Scandella 19-04-2006
        deallocate (node_BD_eff)                                         !DRM Scandella 19-04-2006
        deallocate (node_MDRM)                                           !DRM Scandella 19-04-2006
        deallocate (xx_TOT_eff)                                          !DRM Scandella 19-04-2006
        deallocate (yy_TOT_eff)                                          !DRM Scandella 19-04-2006

        !Node number association to load DRM nodes 
        if (tagstep.eq.2) then                                                 !DRM Scandella 16.11.2005 
            do i = 1,nload_PDRM_el                                              !DRM Scandella 16.11.2005  
            call FIND_NEAREST_NODE(nnod,xx_spx,yy_spx,val_PDRM_el(i,1),&     !DRM Scandella 16.11.2005 
                val_PDRM_el(i,2),node_PDRM_el(i))         !DRM Scandella 16.11.2005                                    
            enddo	                                                             !DRM Scandella 16.11.2005  

            !List of DRM loaded functions associated with DRM nodes ordered for Loukakis forces calculation

            do i = 1,nnode_TOT_eff                                              !DRM Scandella 16.11.2005 
            do j = 1,nnode_TOT_eff                                           !DRM Scandella 16.11.2005
            if (node_TOT_eff(i).eq.node_PDRM_el(j)) then                   !DRM Scandella 16.11.2005
                fn_ord(i)=fun_PDRM_el(j)                                    !DRM Scandella 16.11.2005
                val_TOT_eff(i)=val_PDRM_el(j,3)                             !DRM Scandella 12.04.2006
                exit                                                        !DRM Scandella 16.11.2005
            endif                                                          !DRM Scandella 16.11.2005
            enddo                                                            !DRM Scandella 16.11.2005                                                                        
            enddo                                                               !DRM Scandella 16.11.2005


            ! Stiffnes matrix generation, only for DRM elements and without assembling   
            allocate (K_DRM(2*nn*nn*n_el_DRM,2*nn*nn))                         !DRM Scandella 16.11.2005
            call MAKE_K(nn,con_nnz,con_spx,nelem,nload_MDRM_el, &              !DRM Scandella 16.11.2005
                tag_MDRM_el,n_el_DRM,el_DRM_eff,&                      !DRM Scandella 16.11.2005
                alfa1,alfa2,beta1,beta2,gamma1,gamma2,&                !DRM Scandella 16.11.2005
                nmat,prop_mat,tag_mat,sdeg_mat,K_DRM)                  !DRM Scandella 16.11.2005

            !   do i =1,2*nn*nn*n_el_DRM                                                               !DRM Scandella 16.11.2005                                                       
            !   write(21,'(E12.4,E12.4,E12.4,E12.4,E12.4,E12.4,E12.4,E12.4)')(K_DRM(i,j),j =1,2*nn*nn) !DRM Scandella 16.11.2005
            !   write(21,'(a)')''                                                                      !DRM Scandella 16.11.2005
            !   enddo                                                                                  !DRM Scandella 16.11.2005
        endif                                                                   !DRM Scandella 16.11.2005
    endif                                                                    !DRM Scandella 16.11.2005

    !---DRM---------------------------------------------------------------------------------------------      

    if (nfunc .le. 0) nfunc = 1
    if (nload_sism_el.gt.0) allocate (facsmom(nload_sism_el,3))

    write(*,'(A)')
    write(*,'(A)') '-----------------Make the load vector------------------'

    call system_clock(COUNT=start,COUNT_RATE=clock(2))

    allocate (Fel(nfunc,2*nnod))

    call MAKE_FEL(nnod,xx_spx,yy_spx,con_nnz,con_spx,&
        nmat,tag_mat,sdeg_mat,prop_mat,&
        nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
        con_nnz_bc,con_spx_bc,&
        nload_dirX_el,val_dirX_el,fun_dirX_el,tag_dirX_el,&
        nload_dirY_el,val_dirY_el,fun_dirY_el,tag_dirY_el,&
        nload_neuX_el,val_neuX_el,fun_neuX_el,tag_neuX_el,&
        nload_neuY_el,val_neuY_el,fun_neuY_el,tag_neuY_el,&
        nload_poiX_el,val_poiX_el,fun_poiX_el,&
        nload_poiY_el,val_poiY_el,fun_poiY_el,&
        nload_plaX_el,val_plaX_el,fun_plaX_el,tag_plaX_el,&
        nload_plaY_el,val_plaY_el,fun_plaY_el,tag_plaY_el,&
        nload_sism_el,val_sism_el,fun_sism_el,tag_sism_el,&
        nload_PDRM_el,val_TOT_eff,node_TOT_eff,&              !DRM Scandella 12.04.2006
        nload_MDRM_el,tag_MDRM_el,&                           !DRM Scandella 26.10.2005
        nfunc,tag_func,Fel,&
        xx_macro,yy_macro,con,nelem,con_bc,nedge,&
        num_node_sism,max_num_node_sism,&  
        sour_node_sism,dist_sour_node_sism,&
        length_check_node_sism,facsmom,&
        tagstep,&                              !DRM Scandella 21.10.2005
        nfunc_drm,glob_drm_x,glob_drm_y, &     !DRM Scandella 11.04.2006
        test, n_test, fun_test)  

    call system_clock(COUNT=finish)
    time_MAKE_FEL = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'MAKE_FEL time = ',time_MAKE_FEL,' s'              

    if (nload_sism_el.gt.0) then
        allocate (check_node_sism(length_check_node_sism,5)) 
        allocate (check_dist_node_sism(length_check_node_sism,1))


        call CHECK_SISM(con_nnz,con_spx,&
            nmat,tag_mat,sdeg_mat,&   
            nelem,&
            nload_sism_el,&
            num_node_sism,max_num_node_sism,&
            sour_node_sism,dist_sour_node_sism,&
            check_node_sism,check_dist_node_sism,&
            length_check_node_sism,&
            fun_sism_el,nfunc,tag_func,val_sism_el)

        deallocate (sour_node_sism)
        deallocate (dist_sour_node_sism)
    endif


    write(*,'(A)')'LOAD MATRIX:OK'
     

    deallocate (Ebin)

    write(*,'(A)')
    write(*,'(A)') '*******************************************************'
    write(*,'(A)') '--------------COMPUTING TIME-INDEPENDENT MATRICES-----------------'

    !********************************************************************************************************
    !     REWRITE MASS MATRIX & DUMPING MATRIX IN CRS SPARSE FORMAT
    !********************************************************************************************************

    write(*,'(A)')
    write(*,'(A)') '------------COMPUTING SPARSE MASS MATRIX------------'

    allocate(I_MASS(0:2*nnod));     I_MASS = 0;
    do i = 1, 2*nnod
        I_MASS(i) = I_MASS(i-1) + 1;
    enddo

    allocate(J_MASS(2*nnod), M_MASS(2*nnod));   J_MASS = 0; M_MASS = 0.d0
    allocate(C_MASS(2*nnod), D_MASS(2*nnod));   C_MASS = 0.d0; D_MASS = 0.d0

    do i = 1, 2*nnod
        J_MASS(i) = i;
        M_MASS(i) = Mel(i);
        if (make_damping_yes_or_not .eq. 1) C_MASS(i) = Cel(i);
        if (make_damping_yes_or_not .eq. 1) D_MASS(i) = KCel(i);                
    enddo   

    write(*,'(A)') 'SPARSE MATRIX: OK'

    !********************************************************************************************************
    !     BUILD THE STIFFNESS MATRIX  IN CRS SPARSE FORMAT
    !********************************************************************************************************
    if (.not.NLFLAG) then

        write(*,'(A)')     
        write(*,'(A)') '----------COMPUTING STIFFNESS MATRIX PATTERN----------'

        call system_clock(COUNT=start,COUNT_RATE=clock(2))

        allocate(STIFF_PATTERN(2*nnod, max_nodes));    STIFF_PATTERN = 0;
        allocate(I_STIFF(0:2*nnod));                   I_STIFF = 0;

        call MAKE_PATTERN_STIFF_MATRIX(STIFF_PATTERN, nnod, max_nodes, nelem, nmat, sdeg_mat, &
            con_spx, con_nnz, I_STIFF, length)

        allocate(J_STIFF(1:length),M_STIFF(1:length))
        J_STIFF = 0;     M_STIFF = 0.d0;

        j = 0
        do i = 1 , 2*nnod
        call COUNT_NNZ_EL(STIFF_PATTERN, 2*nnod, max_nodes, i,ic)
        J_STIFF(j+1 : j + ic) = STIFF_PATTERN(i,1:ic)
        j = j + ic
        enddo

        deallocate(STIFF_PATTERN)       
        write(*,'(A)') 'STIFFNESS MATRIX PATTERN:OK'

        call system_clock(COUNT=finish)
        time_MAKE_PATTERN_STIFF_MATRIX = float(finish-start)/float(clock(2))
        write(*,'(A,F8.4,A)')'MAKE_PATTERN_STIFF_MATRIX time = ',time_MAKE_PATTERN_STIFF_MATRIX,' s'        

        write(*,'(A)')
        write(*,'(A)') '-------------COMPUTING STIFFNESS MATRIX-------------'

        call system_clock(COUNT=start,COUNT_RATE=clock(2))
        call MAKE_STIFF_MATRIX(nnod,length,I_STIFF,J_STIFF,M_STIFF, &
                nelem, con_nnz, con_spx, nmat, tag_mat, prop_mat, sdeg_mat, &
                alfa1, alfa2, beta1, beta2, gamma1, gamma2)
        write(*,'(A)') 'Done.'

        call system_clock(COUNT=finish)
        time_MAKE_STIFF_MATRIX = float(finish-start)/float(clock(2))
        write(*,'(A,F8.4,A)')'MAKE_STIFF_MATRIX time = ',time_MAKE_STIFF_MATRIX,' s'  
    endif
    write(*,'(A)') '*******************************************************'
    write(*,'(A)')
    write(*,'(A)') '*******************************************************'

    !********************************************************************************************************
    !     BUILD THE ABCs MATRIX IN CRS SPARSE FORMAT
    !********************************************************************************************************

    nelem_abc = 0; nedge_abc= 0;

    if(nload_abc_el .ge. 1) then

        nedge_abc = 0; 
        if (con_nnz_bc.gt.0) then
            nedge = con_spx_bc(0) -1
            allocate(i4count(nedge)); i4count = 0

            call MAKE_ABC(nedge_abc, nedge, i4count,con_nnz_bc,con_spx_bc,&
                con_nnz,con_spx, nload_abc_el,tag_abc_el)

            if (nedge_abc.gt.0) then
                allocate(iedge_abc(nedge_abc))
                do iedge = 1,nedge
                if (i4count(iedge).ne.0) then
                    iedge_abc(i4count(iedge)) = iedge
                endif
                enddo

                deallocate(i4count)
            endif    
        endif

        nelem_abc = nedge_abc

        if (nelem_abc.gt.0) then
            allocate(ielem_abc(nelem_abc))

            do i = 1,nedge_abc
                iedge = iedge_abc(i)
                ied1 = con_spx_bc(con_spx_bc(iedge -1) +1)
                ied2 = con_spx_bc(con_spx_bc(iedge) -1)

                do ie = 1, nelem
                    nn = con_spx_bc(iedge) - con_spx_bc(iedge -1) -1
                    iel1 = con_spx(con_spx(ie -1) +1)
                    iel2 = con_spx(con_spx(ie -1) +nn)
                    iel3 = con_spx(con_spx(ie -1) +nn*nn)
                    iel4 = con_spx(con_spx(ie -1) +nn*(nn -1) +1)
                    if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. (ied1.eq.iel3).or.(ied1.eq.iel4)).and. &
                        ((ied2.eq.iel1).or.(ied2.eq.iel2).or. (ied2.eq.iel3).or.(ied2.eq.iel4))) then
                        ielem_abc(i) = ie
                    endif
                enddo
            enddo

            write(*,'(A)') '----------COMPUTING ABC PATTERN------------'

            call system_clock(COUNT=start,COUNT_RATE=clock(2))

            allocate(ABC_PATTERN(2*nnod, max_nodes));    ABC_PATTERN = 0;
            allocate(I_ABC(0:2*nnod));                   I_ABC = 0;

            call MAKE_PATTERN_ABC_MATRIX(ABC_PATTERN, nnod, max_nodes, nelem_abc, ielem_abc, nmat, sdeg_mat, &
                con_spx, con_nnz, I_ABC, length_abc)

            allocate(J_ABC(1:length_abc),M_ABC_U(1:length_abc),M_ABC_V(1:length_abc))
            J_ABC = 0;     M_ABC_U = 0.d0;  M_ABC_V = 0.d0; 

            j = 0
            do i = 1 , 2*nnod
                call COUNT_NNZ_EL(ABC_PATTERN, 2*nnod, max_nodes, i,ic)
                J_ABC(j+1 : j + ic) = ABC_PATTERN(i,1:ic)
                j = j + ic
            enddo

            deallocate(ABC_PATTERN)       
        endif

        write(*,'(A)') 'ABC PATTERN:OK'

        call system_clock(COUNT=finish)
        time_MAKE_PATTERN_ABC_MATRIX = float(finish-start)/float(clock(2))
        write(*,'(A,F8.4,A)')'MAKE_PATTERN_ABC_MATRIX time = ',time_MAKE_PATTERN_ABC_MATRIX,' s' 

        write(*,'(A)')
        write(*,'(A)') '----------------COMPUTING ABC MATRIX--------------'

        call system_clock(COUNT=start,COUNT_RATE=clock(2))

        call MAKE_ABCS_MATRIX(nnod,length_abc,I_ABC,J_ABC,M_ABC_U,M_ABC_V, &
            nelem, nelem_abc, ielem_abc, nedge_abc, iedge_abc, & 
            con_nnz_bc, con_spx_bc, con_nnz, con_spx, nmat, tag_mat, prop_mat, sdeg_mat, &
            alfa1, alfa2, beta1, beta2, gamma1, gamma2, xx_spx, yy_spx)

        write(*,'(A)') 'ABC MATRIX:OK'

        call system_clock(COUNT=finish)
        time_MAKE_ABCS_MATRIX = float(finish-start)/float(clock(2))
        write(*,'(A,F8.4,A)')'MAKE_ABCS_MATRIX time = ',time_MAKE_ABCS_MATRIX,' s' 

    endif
    write(*,*)'*********************************************'
    !********************************************************************************************************
    !     MAKE LOCAL DG MATRIX
    !********************************************************************************************************
    nelem_dg = 0;     nnode_dg = 0;
    allocate(i4count(nnod));       i4count = 0;


    call GET_NODE_FROM_FACE(nnod, con_nnz_bc, con_spx_bc, nload_dg_el, tag_dg_el,&
        nnode_dg, i4count)

    ! output: - number of dg local elements 
    !         - number of dg global elements 

    call GET_DIME_DG(nmat, sdeg_mat, tag_mat, con_nnz, con_spx, &
        nnod, xx_spx, yy_spx, nelem_dg, i4count)

    allocate(faces(3,nelem_dg), area_nodes(9,nelem_dg))  
    faces = 0; area_nodes = 0.d0

    if(nelem_dg .gt. 0) then

        file_face = 'FACES.input'

        call system_clock(COUNT=start,COUNT_RATE=clock(2))

        call SETUP_DG(nmat, sdeg_mat, tag_mat, con_nnz, con_spx, &
            nnod, nelem, xx_spx, yy_spx, nelem_dg, i4count,&
            alfa1, alfa2, beta1, beta2, gamma1, gamma2, delta1, delta2,&
            faces, area_nodes)

        call system_clock(COUNT=finish)
        time_SETUP_DG = float(finish-start)/float(clock(2))
        write(*,'(A,F8.4,A)')'SETUP_DG time = ',time_SETUP_DG,' s' 

        inquire(file=file_face, exist=filefound)
        if(filefound .eqv. .FALSE.) then

            write(*,'(A)') 'Writing FACES.input'

            open(400,file=file_face)
            do j = 1, nelem_dg
                write(400,"(1I2,1X,1I12,1X,1I2,9(2X,ES16.9))") &
                    faces(1,j), faces(2,j), faces(3,j), &
                    area_nodes(1,j), area_nodes(2,j), area_nodes(3,j), &
                    area_nodes(4,j), area_nodes(5,j), area_nodes(6,j), &
                    area_nodes(7,j), area_nodes(8,j), area_nodes(9,j)
            enddo
            close(400)

        endif  


        head_file = 'DGFS.input'
        inquire(file=head_file, exist=filefound)
        if(filefound .eqv. .FALSE.) then     

            write(*,'(A)') 'Writing DGFS.input'

            ! output: - write DGFS.input file

            call system_clock(COUNT=start,COUNT_RATE=clock(2))

            allocate(dg_els(nelem_dg), scratch_dg_els(nelem_dg)) 

            call SETUP_DG_ELEM(nmat, sdeg_mat, tag_mat, con_nnz, con_spx, &
                nnod, nelem, xx_spx,yy_spx,&
                nelem_dg, i4count, &
                alfa1, alfa2,beta1, beta2, gamma1, gamma2, delta1, delta2,&
                dg_els, scratch_dg_els, &
                tag_dg_el, tag_dg_yn, nload_dg_el, &
                con_bc, nedge)

            call system_clock(COUNT=finish)
            time_SETUP_DG_ELEM = float(finish-start)/float(clock(2))
            write(*,'(A,F8.4,A)')'SETUP_DG_ELEM time = ',time_SETUP_DG_ELEM,' s'  

            call WRITE_FILE_DGFS(nmat, sdeg_mat, tag_mat, con_nnz, con_spx, &
                nnod, nelem, xx_spx, yy_spx,&
                nelem_dg, &
                alfa1, alfa2, beta1, beta2, gamma1, gamma2, delta1, delta2, &
                faces, area_nodes, dg_els, scratch_dg_els, &
                head_file)

            deallocate(dg_els, scratch_dg_els) 


        endif


        allocate(el_new(nelem_dg))  

        write(*,'(A)') 
        write(*,'(A)') '-----------------Making DG interfaces------------------' 

        call system_clock(COUNT=start,COUNT_RATE=clock(2))

        call MAKE_DG_INTERFACE(nmat, sdeg_mat, tag_mat, prop_mat, con_nnz, con_spx, &
            nnod, nelem, xx_spx, yy_spx, &
            nelem_dg, i4count, &
            alfa1, alfa2, beta1, beta2, gamma1, gamma2, &
            delta1, delta2, dg_const, dg_pen, & 
            faces, area_nodes, el_new, head_file, test)


        write(*,'(A)') 'Done.'      

        call system_clock(COUNT=finish)
        time_MAKE_DG_INTERFACE = float(finish-start)/float(clock(2))
        write(*,'(A,F8.4,A)')'MAKE_DG_INTERFACE time = ',time_MAKE_DG_INTERFACE,' s' 	

    endif


    deallocate(faces, area_nodes,i4count)

    !********************************************************************************************************
    !     MAKE GLOBAL DG MATRIX IN CRS SPARSE FORMAT
    !********************************************************************************************************
    if(nelem_dg .gt. 0) then

        nnz_dg_total = 0; 

        do i = 1, nelem_dg
        nnz_dg_total = nnz_dg_total + el_new(i)%nnz_plus + el_new(i)%nnz_minus
        enddo

        allocate(IDG_TOTAL(nnz_dg_total),JDG_TOTAL(nnz_dg_total),MDG_TOTAL(nnz_dg_total))


        call MAKE_DG_SPARSE_TOTAL(nelem_dg, el_new, nnz_dg_total, con_nnz, con_spx, &
            nnod, nmat, sdeg_mat, IDG_TOTAL, &
            JDG_TOTAL, MDG_TOTAL)

        call COUNT_EFFECTIVE_LENGTH(nnz_dg_total,IDG_TOTAL,JDG_TOTAL,MDG_TOTAL,nnz_dg)

        allocate(IDG_MORSE(nnz_dg), JDG_MORSE(nnz_dg), MDG_MORSE(nnz_dg))                          

        call MAKE_DG_SPARSE(nnz_dg_total,IDG_TOTAL,JDG_TOTAL,MDG_TOTAL,nnz_dg,IDG_MORSE,JDG_MORSE,MDG_MORSE)

        !Convert sparse matrix into RCS matrix
        allocate(IDG(0:2*nnod), JDG(nnz_dg), MDG(nnz_dg))                  


        call MORSE_2_RCS(nnod, nnz_dg, IDG_MORSE, JDG_MORSE, MDG_MORSE, &
            IDG, JDG, MDG)

        deallocate(IDG_TOTAL,JDG_TOTAL,MDG_TOTAL,IDG_MORSE,JDG_MORSE,MDG_MORSE)


        if(test .eq. 1) then
            nnz_dg_total_only_uv = 0; 

            do i = 1, nelem_dg
            nnz_dg_total_only_uv = nnz_dg_total_only_uv + &
                el_new(i)%nnz_plus_only_uv + el_new(i)%nnz_minus_only_uv
            enddo

            allocate(IDG_TOTAL(nnz_dg_total_only_uv),JDG_TOTAL(nnz_dg_total_only_uv),&
                MDG_TOTAL(nnz_dg_total_only_uv))


            call MAKE_DG_SPARSE_TOTAL_ONLY_UV(nelem_dg, el_new, nnz_dg_total_only_uv, con_nnz, con_spx, &
                nnod, nmat, sdeg_mat, IDG_TOTAL, &
                JDG_TOTAL, MDG_TOTAL)

            call COUNT_EFFECTIVE_LENGTH(nnz_dg_total_only_uv,IDG_TOTAL,JDG_TOTAL,&
                MDG_TOTAL, nnz_dg_only_uv)

            allocate(IDG_MORSE(nnz_dg_only_uv), JDG_MORSE(nnz_dg_only_uv), MDG_MORSE(nnz_dg_only_uv))                          

            call MAKE_DG_SPARSE(nnz_dg_total_only_uv,IDG_TOTAL,JDG_TOTAL,MDG_TOTAL,&
                nnz_dg_only_uv,IDG_MORSE,JDG_MORSE,MDG_MORSE)

            !Convert sparse matrix into RCS matrix
            allocate(IDG_only_uv(0:2*nnod), JDG_only_uv(nnz_dg_only_uv), MDG_only_uv(nnz_dg_only_uv))                  


            call MORSE_2_RCS(nnod, nnz_dg_only_uv, IDG_MORSE, JDG_MORSE, MDG_MORSE, &
                IDG_only_uv, JDG_only_uv, MDG_only_uv)

            deallocate(IDG_TOTAL,JDG_TOTAL,MDG_TOTAL,IDG_MORSE,JDG_MORSE,MDG_MORSE)

        endif


    endif   



    !********************************************************************************************************
    !     BUILDING GLOBAL MATRICES -- LEGEND:
    !     M = mass matrix;  C,D = damping matrix; R,S = abc matrix; A = stiffness matrix; B = dg matrix 
    !                                     MU''  + (C-S)U' + (A+B+D-R)U = F 
    !     M = M_MASS
    !     C = C_MASS
    !     D = D_MASS
    !     S = M_ABC_V
    !     R = M_ABC_U
    !     A = M_STIFF
    !     B = MDG
    !********************************************************************************************************

    if(nload_abc_el .eq. 0) then 
        allocate(I_ABC(0:2*nnod), J_ABC(2*nnod),M_ABC_U(2*nnod),M_ABC_V(2*nnod))
        I_ABC = 0;   J_ABC = 0;     M_ABC_U = 0.d0;  M_ABC_V = 0.d0; 
    endif

    write(*,'(A)') '*******************************************************'

    write(*,'(A)')
    write(*,'(A)') '*******************************************************'
    write(*,'(A)') '-------------------COMPUTE MATRIX SUM-------------------'

    allocate(NDEGR(2*nnod),IW(2*nnod))

    M_ABC_V = - M_ABC_V;      M_ABC_U = - M_ABC_U;
    I_ABC = I_ABC + 1;        I_MASS = I_MASS + 1;

    call aplbdg ( 2*nnod, 2*nnod,  J_MASS, I_MASS, J_ABC, I_ABC, NDEGR, NNZ_AB, IW )

    allocate(IC_SUM(0:2*nnod), JC_SUM(NNZ_AB), C_SUM(NNZ_AB))
    allocate(ID_SUM(0:2*nnod), JD_SUM(NNZ_AB), D_SUM(NNZ_AB))

    !Computing the nonzero elements for the sum
    call aplb ( 2*nnod, 2*nnod, 0, C_MASS, J_MASS, I_MASS, M_ABC_V, J_ABC, I_ABC, &
        C_SUM, JC_SUM, IC_SUM, NNZ_AB, IW, ierr)

    call aplb ( 2*nnod, 2*nnod, 0, D_MASS, J_MASS, I_MASS, M_ABC_U, J_ABC, I_ABC, &
        D_SUM, JD_SUM, ID_SUM, NNZ_AB, IW, ierr)

    !Summing C_SUM = C_MASS - M_ABC_V 
    call aplb ( 2*nnod, 2*nnod, 1, C_MASS, J_MASS, I_MASS, M_ABC_V, J_ABC, I_ABC, &
        C_SUM, JC_SUM, IC_SUM, NNZ_AB, IW, ierr)

    !Summing D_SUM = D_MASS - M_ABC_U 
    call aplb ( 2*nnod, 2*nnod, 1, D_MASS, J_MASS, I_MASS, M_ABC_U, J_ABC, I_ABC, &
        D_SUM, JD_SUM, ID_SUM, NNZ_AB, IW, ierr)
    

    if (NLFLAG) then
        !Summing E_SUM = D_SUM = D_MASS - M_ABC_U
        allocate(IE_SUM(0:2*nnod), JE_SUM(NNZ_AB), E_SUM(NNZ_AB))
        E_SUM  = D_SUM
        IE_SUM = ID_SUM
        JE_SUM = JD_SUM
    else
        I_STIFF = I_STIFF + 1
        !Computing the nonzero elements for the sum
        call aplbdg ( 2*nnod, 2*nnod,  J_STIFF, I_STIFF, JD_SUM, ID_SUM, NDEGR, NNZ_AB, IW )
        allocate(IE_SUM(0:2*nnod), JE_SUM(NNZ_AB), E_SUM(NNZ_AB))
        !Summing E_SUM = M_STIFF + D_SUM = M_STIFF + D_MASS - M_ABC_U
        call aplb ( 2*nnod, 2*nnod, 0, M_STIFF, J_STIFF, I_STIFF, D_SUM, JD_SUM, ID_SUM, &
            E_SUM, JE_SUM, IE_SUM, NNZ_AB, IW, ierr)
        call aplb ( 2*nnod, 2*nnod, 1, M_STIFF, J_STIFF, I_STIFF, D_SUM, JD_SUM, ID_SUM, &
            E_SUM, JE_SUM, IE_SUM, NNZ_AB, IW, ierr)
        deallocate(I_STIFF,J_STIFF,M_STIFF)
    endif

    if(nelem_dg .gt. 0) then
        IDG = IDG + 1
        !Computing the nonzero elements for the sum
        call aplbdg ( 2*nnod, 2*nnod,  JE_SUM, IE_SUM, JDG, IDG, NDEGR, NNZ_AB, IW )
        allocate(IDG_SUM(0:2*nnod), JDG_SUM(NNZ_AB), MDG_SUM(NNZ_AB))
        !Summing MDG_SUM = E_SUM + MDG = M_STIFF + MDG + D_MASS - M_ABC_U 
        call aplb ( 2*nnod, 2*nnod, 0, E_SUM, JE_SUM, IE_SUM, MDG, JDG, IDG, &
            MDG_SUM, JDG_SUM, IDG_SUM, NNZ_AB, IW, ierr)
        call aplb ( 2*nnod, 2*nnod, 1, E_SUM, JE_SUM, IE_SUM, MDG, JDG, IDG, &
            MDG_SUM, JDG_SUM, IDG_SUM, NNZ_AB, IW, ierr)
        deallocate(IE_SUM,JE_SUM,E_SUM)
        allocate(IE_SUM(0:2*nnod), JE_SUM(NNZ_AB), E_SUM(NNZ_AB))
        IE_SUM = IDG_SUM; JE_SUM = JDG_SUM; E_SUM = MDG_SUM;
        deallocate(IDG_SUM,JDG_SUM,MDG_SUM)
    endif

    deallocate(C_MASS, D_MASS, I_ABC, J_ABC, M_ABC_U, M_ABC_V)
    write(*,'(A)') 'MATRIX SUM:OK'
    write(*,'(A)') '*******************************************************'

    !********************************************************************************************************
    !     COMPUTE M^-1*(C-S) AND M^-1*(A+B+D-R)
    !********************************************************************************************************

    write(*,'(A)')
    write(*,'(A)') '*******************************************************'
    write(*,'(A)') '----------------MATRIX MULTIPLICATION------------------'

    M_MASS = 1./M_MASS;

    !Computing the nonzero elements for the mul      
    call amubdg ( 2*nnod, 2*nnod, 2*nnod, J_MASS, I_MASS, JC_SUM, IC_SUM, NDEGR, NNZ_AB, IW )
    allocate(IN_TOT(0:2*nnod), JN_TOT(NNZ_AB), N_TOT(NNZ_AB))

    !Multiplying N = M^-1*(C-S)      
    call amub ( 2*nnod, 2*nnod, 0, M_MASS, J_MASS, I_MASS, C_SUM, JC_SUM, IC_SUM, &
        N_TOT, JN_TOT, IN_TOT, NNZ_AB, IW, ierr )

    call amub ( 2*nnod, 2*nnod, 1, M_MASS, J_MASS, I_MASS, C_SUM, JC_SUM, IC_SUM, &
        N_TOT, JN_TOT, IN_TOT, NNZ_AB, IW, ierr )

    NNZ_N = NNZ_AB

    !Computing the nonzero elements for the mul      
    call amubdg ( 2*nnod, 2*nnod, 2*nnod, J_MASS, I_MASS, JE_SUM, IE_SUM, NDEGR, NNZ_AB, IW )
    
    allocate(IK_TOT(0:2*nnod), JK_TOT(NNZ_AB), K_TOT(NNZ_AB))

    !Multiplying K = M^-1*(A+B+D-R)
    call amub ( 2*nnod, 2*nnod, 0, M_MASS, J_MASS, I_MASS, E_SUM, JE_SUM, IE_SUM, &
        K_TOT, JK_TOT, IK_TOT, NNZ_AB, IW, ierr )

    call amub ( 2*nnod, 2*nnod, 1, M_MASS, J_MASS, I_MASS, E_SUM, JE_SUM, IE_SUM, &
        K_TOT, JK_TOT, IK_TOT, NNZ_AB, IW, ierr )

    NNZ_K = NNZ_AB
    deallocate(M_MASS, I_MASS, J_MASS, C_SUM, JC_SUM, IC_SUM, E_SUM, JE_SUM, IE_SUM, NDEGR, IW)
    write(*,'(A)') 'MATRIX MULTIPLICATION:OK'
    write(*,'(A)') '*******************************************************'

    !********************************************************************************************************
    !     BUILDING THE RHS M^-1*F
    !********************************************************************************************************

    write(*,'(A)')
    write(*,'(A)') '*******************************************************'
    write(*,'(A)') '-------------------COMPUTING RHS--------------------'

    do i = 1, nfunc
        Fel(i,:) = Fel(i,:)/Mel
    enddo

    write(*,'(A)') 'RHS:OK'
    write(*,'(A)') '*******************************************************'

    !********************************************************************************************************
    !     SETTING CALCULATION PARAMETERS
    !********************************************************************************************************

    write(*,'(A)')
    write(*,'(A)') '*******************************************************'
    write(*,'(A)') '-----------Setting calculation parameters--------------'



    !     INITIAL CONDITIONS

    allocate (u0(2*nnod),v0(2*nnod));   u0 = 0.0d0; v0 = 0.0d0

    if (test .ge. 1) then
        pi = 4.d0*datan(1.d0)
        do i = 1, nnod 
        !u0(i) = (xx_spx(i)**2-1.d0)*(yy_spx(i)**2-1.d0)
        !u0(i+nnod) = 0.d0 
        u0(i) = 0.d0
        u0(i+nnod) = 0.d0 
        v0(i) = - dsqrt(2.d0) * pi * dsin(pi*xx_spx(i))**2 * dsin(2.d0*pi*yy_spx(i)) 
        v0(i+nnod) = dsqrt(2.d0) * pi * dsin(2.d0*pi*xx_spx(i)) * dsin(pi*yy_spx(i))**2

        enddo  
        !v0 = 0.0d0
    endif

    !     CFL CONDITION

    write(*,'(A,E12.4)')'Time step = ',deltat

    call system_clock(COUNT=start,COUNT_RATE=clock(2))

    call DELTAT_MAX(deltat,nnod,nmat,tag_mat,prop_mat,sdeg_mat,&
        xx_macro,yy_macro,nelem,con,deltat_cfl,fmax)

    call system_clock(COUNT=finish)
    time_DELTAT_MAX = float(finish-start)/float(clock(2))
    write(*,'(A,F8.4,A)')'DELTAT_MAX time = ',time_DELTAT_MAX,' s'

    nts = int(xtime / deltat)
    write(*,'(A,I8)')'Number of time-steps = ',nts

    xtime = dfloat(nts) * deltat
    write(*,'(A,E12.4)')'Final time = ',xtime

    !     SNAPSHOTS

    if (nsnaps.gt.0) then
        do i = 1,nsnaps
        itersnap(i) = int(tsnap(i)/deltat)
        if (itersnap(i).gt.nts) itersnap(i) = nts
        enddo
    endif


    !     MONITORED NODES

    nmonitors = 0
    if (num_lst .eq. 1) then

        file_LS = 'LS.input'
        call READ_DIME_FILEPG(file_LS,nmonitors)     

        allocate(n_monitor(nmonitors),dist_monitor_lst(nmonitors))
        allocate(x_monitor_lst(nmonitors),y_monitor_lst(nmonitors))
        allocate(x_monitor(nmonitors),y_monitor(nmonitors))

        call READ_FILEPG(file_LS,nmonitors,x_monitor_lst,y_monitor_lst)

        if (file_mon_lst.eq.0) then ! NO input file with the position of LST monitors

            do i = 1,nmonitors
            call GET_NEAREST_NODE_PGM(nnod, xx_spx, yy_spx, &
                x_monitor_lst(i), y_monitor_lst(i),&
                n_monitor(i), dist_monitor_lst(i), depth_search_mon_lst)

            x_monitor(i) = xx_spx(n_monitor(i))
            y_monitor(i) = yy_spx(n_monitor(i))
            enddo

            file_MLST = 'MLST.input'
            call WRITE_FILE_MPGM(file_MLST, nmonitors, n_monitor, x_monitor, y_monitor)

        else ! YES, it exists an input file with the position of LST monitors

            file_MLST = 'MLST.input'
            call READ_FILE_MPGM(file_MLST, nmonitors, n_monitor, &
                x_monitor, y_monitor)

        endif

        deallocate(dist_monitor_lst,x_monitor_lst,y_monitor_lst)

    else
        write(*,'(A)') 'MLST key not found!'
    endif

    write(*,'(A)')'Monitored nodes'
    write(*,'(I6,2E14.5,I6)')(i,x_monitor(i),y_monitor(i),n_monitor(i), i=1,nmonitors)


    call system_clock(COUNT=clock_finish)
    time_in_seconds = float(clock_finish-clock_start)/float(clock(2))

    write(*,'(A)')
    write(*,'(A)') '-------------------------------------------------------'
    write(*,'(A,F8.4,A)')'Set-up time = ',time_in_seconds,' s'
    write(*,'(A)') '-------------------------------------------------------'
    write(*,'(A)')
    write(*,'(A)')
    write(*,'(A)') '-------------------------------------------------------'
    write(*,'(A)')'             Beginning of the time-loop                 '
    write(*,'(A)')

    if (.not.NLFLAG) then
        call TIME_LOOP_EL(nnod,xx_spx,yy_spx,con_nnz,con_spx,&                                                         ! 5
            nmat,tag_mat,sdeg_mat,prop_mat,&                                                                 ! 4
            nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,delta1,delta2,&                                      ! 9
            con_nnz_bc,con_spx_bc,&                                                                          ! 2
            nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&                                            ! 4
            nload_abc_el,tag_abc_el,&                                                                        ! 2
            nelem_abc,nedge_abc,ielem_abc,iedge_abc,&                                                        ! 4
            nfunc,func_type,func_indx,nfunc_data,func_data,tag_func,&                                        ! 6
            nfunc_drm,func_type_drm,func_indx_drm,nfunc_data_drm,func_data_drm,&  !DRM Scandella 11.04.2006  ! 5
            ndt_mon, &                                                            !DRM Scandella 24.01.2006  ! 1          
            K_TOT, IK_TOT, JK_TOT, NNZ_K, &                                                                  ! 4
            N_TOT, IN_TOT, JN_TOT, NNZ_N, Mel,&                                                              ! 5
            Fel,u0,v0,&                                                                                      ! 3
            nts,deltat,nmonitors,n_monitor,nsnaps,itersnap,&                                                 ! 6
            check_node_sism,check_dist_node_sism,&                                                           ! 2
            length_check_node_sism,facsmom,&                                                                 ! 2
            nload_sism_el,&                                                                                  ! 1
            make_damping_yes_or_not,&                                                                        ! 1
            nnode_TOT_eff,node_TOT_eff,tagstep,&                                  !DRM Scandella 17.10.2005  ! 3 
            ns,n_el_DRM,el_DRM_eff,K_DRM,nnode_BD_eff,&                           !DRM Scandella 21.10.2005  ! 5
            nload_MDRM_el,tag_MDRM_el,val_PDRM_el,&                               !DRM Scandella 25.10.2005  ! 3
            fn_ord,node_PDRM_el,&                                                 !DRM Scandella 16.11.2005  ! 2
            glob_drm_x,glob_drm_y, &                                              !DRM Scandella 11.04.2006  ! 2
            opt_out_var,test,nelem_dg,&                                                                      ! 3
            IDG_only_uv, JDG_only_uv, MDG_only_uv,nnz_dg_only_uv)
    else
        call TIME_LOOP_NL(nnod,xx_spx,yy_spx,con_nnz,con_spx,nmat,tag_mat,sdeg_mat,   &
            prop_mat,nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,delta1,delta2,       &
            con_nnz_bc,con_spx_bc,nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&
            nload_abc_el,tag_abc_el,nelem_abc,nedge_abc,ielem_abc,iedge_abc,          &
            nfunc,func_type,func_indx,nfunc_data,func_data,tag_func,                  &
            nfunc_drm,func_type_drm,func_indx_drm,nfunc_data_drm,func_data_drm,       & 
            ndt_mon, N_TOT, IN_TOT, JN_TOT, NNZ_N, Mel,Fel,u0,v0,nts,deltat,nmonitors,&
            n_monitor,nsnaps,itersnap,check_node_sism,check_dist_node_sism,           &
            length_check_node_sism,facsmom,nload_sism_el,make_damping_yes_or_not,     &
            nnode_TOT_eff,node_TOT_eff,tagstep,ns,n_el_DRM,el_DRM_eff,K_DRM,          &
            nnode_BD_eff,nload_MDRM_el,tag_MDRM_el,val_PDRM_el,fn_ord,node_PDRM_el,   &
            glob_drm_x,glob_drm_y,opt_out_var,test,nelem_dg,IDG_only_uv,JDG_only_uv,  &
            MDG_only_uv,nnz_dg_only_uv)
    endif
    write(*,'(A)')
    write(*,'(A)')'Bye.'


    deallocate (tag_mat,sdeg_mat,prop_mat)
    deallocate (xx_macro,yy_macro,con)
    deallocate (xx_spx,yy_spx,con_spx)
    deallocate (alfa1,beta1,gamma1,delta1,alfa2,beta2,gamma2,delta2)
    deallocate (Fel,u0,v0)
    if (make_damping_yes_or_not .eq. 1) deallocate(Cel,KCel)
    if (nedge_abc.gt.0) deallocate(iedge_abc)
    if (nelem_abc.gt.0) deallocate(ielem_abc)

    if (nmonitors.gt.0) deallocate(x_monitor,y_monitor,n_monitor)
    if (nsnaps.gt.0) deallocate(tsnap,itersnap)

    if (nload_dirX_el.gt.0) deallocate (val_dirX_el,fun_dirX_el,tag_dirX_el)
    if (nload_dirY_el.gt.0) deallocate (val_dirY_el,fun_dirY_el,tag_dirY_el)
    if (nload_neuX_el.gt.0) deallocate (val_neuX_el,fun_neuX_el,tag_neuX_el)
    if (nload_neuY_el.gt.0) deallocate (val_neuY_el,fun_neuY_el,tag_neuY_el)
    if (nload_poiX_el.gt.0) deallocate (val_poiX_el,fun_poiX_el)
    if (nload_poiY_el.gt.0) deallocate (val_poiY_el,fun_poiY_el)
    if (nedge.gt.0) deallocate (con_bc,con_spx_bc)
    if (n_test.gt.0) deallocate (fun_test)      

    !--------DRM-----------------------------------------------------------------------------
    if (nload_PDRM_el.gt.0) deallocate (val_PDRM_el,fun_PDRM_el,K_DRM,fn_ord,& !DRM Scandella 21.10.2005 
        node_PDRM_el,glob_DRM_x,glob_DRM_y)    !DRM Scandella 21.10.2005
    if (nload_MDRM_el.gt.0) deallocate (tag_MDRM_el)                           !DRM Scandella 09.08.2007 
    if (nload_BDRM_el.gt.0) deallocate (tag_BDRM_el)                           !DRM Scandella 09.08.2007  
    if (n_el_DRM.gt.0) deallocate (el_DRM_eff)                                 !DRM Scandella 19.04.2006
    if (nfunc.gt.1) deallocate (tag_func,func_type,func_indx,func_data)
    if (nfunc_drm.gt.0) deallocate (tag_func_drm,func_type_drm,func_indx_drm,func_data_drm) !DRM Scandella 19.04.2006     
    !----------------------------------------------------------------------------------------  

    deallocate(IN_TOT, JN_TOT, N_TOT)
    if(.not.NLFLAG) then
        deallocate(IK_TOT, JK_TOT, K_TOT)
    endif
end program SPEED2D
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

