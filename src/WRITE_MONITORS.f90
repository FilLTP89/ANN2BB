module write_output
    !
    use fields
    !
    implicit none
    !
    contains
        
        !****************************************************************************
        ! OPEN MONITOR FILE
        !****************************************************************************
        
        subroutine MAKE_MONITOR_FILES(nnt,nmonit,option_out_var,nnode_TOT,tagstep,&
                unit_disp,unit_vel,unit_acc,unit_stress,unit_strain,unit_omega,unit_uDRM)
            !
            implicit none
            ! intent IN
            integer*4 ,                   intent(in)    :: nnt,nmonit,nnode_TOT,tagstep
            integer*4, dimension(6)     , intent(in)    :: option_out_var    
            ! intent OUT
            integer*4, dimension(nmonit), intent(inout) :: unit_disp
            integer*4, dimension(nmonit), intent(inout) :: unit_vel
            integer*4, dimension(nmonit), intent(inout) :: unit_acc
            integer*4, dimension(nmonit), intent(inout) :: unit_stress
            integer*4, dimension(nmonit), intent(inout) :: unit_strain
            integer*4, dimension(nmonit), intent(inout) :: unit_omega
            integer*4, dimension(nnode_TOT), intent(inout) :: unit_uDRM
            !
            character*70                                :: file_disp
            character*70                                :: file_vel
            character*70                                :: file_acc
            character*70                                :: file_stress 
            character*70                                :: file_strain 
            character*70                                :: file_omega 
            character*70                                :: file_uDRM
            integer*4                                   :: i,last_number
            !
            if (nmonit.ge.1) then
                last_number = 0
                ! DISPLACEMENT
                if (option_out_var(1).eq.1) then  
                    file_disp = 'monitorXXXXX.d'  
                    call MAKE_MONITOR_NAME(nmonit,file_disp,unit_disp,last_number)
                endif
                ! VELOCITY
                if (option_out_var(2).eq.1) then  
                    file_vel = 'monitorXXXXX.v' 
                    call MAKE_MONITOR_NAME(nmonit,file_vel,unit_vel,last_number)
                endif
                ! ACCELERATION
                if (option_out_var(3).eq.1) then  
                    file_acc = 'monitorXXXXX.a' 
                    call MAKE_MONITOR_NAME(nmonit,file_acc,unit_acc,last_number)
                endif
                ! STRESS
                if (option_out_var(4).eq.1) then
                    file_stress = 'monitorXXXXX.s'
                    call MAKE_MONITOR_NAME(nmonit,file_stress,unit_stress,last_number)
                endif
                ! STRAIN
                if (option_out_var(5).gt.0) then  
                    file_strain = 'monitorXXXXX.e'
                    call MAKE_MONITOR_NAME(nmonit,file_strain,unit_strain,last_number)
                endif
                ! OMEGA
                if (option_out_var(6).eq.1) then  
                    file_omega = 'monitorXXXXX.w'
                    call MAKE_MONITOR_NAME(nmonit,file_omega,unit_omega,last_number)
                endif
            endif

            !---DRM---------------------------------------------------------------------
            !Open output files for DRM I step                     !DRM Scandella 25.11.2005
            if ((nnode_TOT.ne.0).and.(tagstep.eq.1)) then   !DRM Scandella 25.11.2005 
                file_uDRM = 'monDRMXXXXX.d'                  !DRM Scandella 25.11.2005 
                do i = 1,nnode_TOT                           !DRM Scandella 25.11.2005 
                    unit_uDRM(i)=700000+i                    !DRM Scandella 25.11.2005
                    if (i.lt.10) then                      !DRM Scandella 25.11.2005 
                        write(file_uDRM(7:10),'(a4)')'0000'     !DRM Scandella 25.11.2005
                        write(file_uDRM(11:11),'(i1)')i     !DRM Scandella 25.11.2005
                    else if (i.le.99) then                 !DRM Scandella 25.11.2005
                        write(file_uDRM(7:9),'(a3)')'000'     !DRM Scandella 25.11.2005
                        write(file_uDRM(10:11),'(i2)')i     !DRM Scandella 25.11.2005
                    else if (i.le.999) then                !DRM Scandella 25.11.2005
                        write(file_uDRM(7:8),'(a2)')'00'     !DRM Scandella 25.11.2005
                        write(file_uDRM(9:11),'(i3)')i      !DRM Scandella 25.11.2005
                    else if (i.le.9999) then               !DRM Scandella 25.11.2005
                        write(file_uDRM(7:7),'(a1)')'0'     !DRM Scandella 25.11.2005 
                        write(file_uDRM(8:11),'(i4)')i      !DRM Scandella 25.11.2005
                    else if (i.le.99999) then              !DRM Scandella 25.11.2005
                        write(file_uDRM(7:11),'(i5)')i      !DRM Scandella 25.11.2005  
                    endif                                  !DRM Scandella 25.11.2005              
                    open(unit_uDRM(i),file=file_uDRM)         !DRM Scandella 25.11.2005
                    !close(unit_uDRM)                      !DRM Scandella 13.12.2005 			       
                enddo                                        !DRM Scandella 25.11.2005
            endif                                           !DRM Scandella 25.11.2005
            return
        end  subroutine MAKE_MONITOR_FILES
       
        subroutine MAKE_MONITOR_NAME(nmonit,file_name,file_unit,last_number)
            !
            implicit none
            ! intent IN
            integer*4, intent(in) :: nmonit
            ! integer INOUT 
            integer*4, intent(inout) :: last_number
            character*70, intent(inout) :: file_name
            ! intent OUT
            integer*4, dimension(nmonit), intent(out) :: file_unit
            !
            integer*4 :: i
            !
            last_number = last_number+100000
            
            do i = 1,nmonit
                file_unit(i)=last_number+i
                if (i.lt.10) then
                    write(file_name(8:11),'(a4)')'0000'
                    write(file_name(12:12),'(i1)')i
                else if (i.le.99) then
                    write(file_name(8:10),'(a3)')'000'
                    write(file_name(11:12),'(i2)')i
                else if (i.le.999) then
                    write(file_name(8:9),'(a2)')'00'
                    write(file_name(10:12),'(i3)')i
                else if (i.le.9999) then
                    write(file_name(8:8),'(a1)')'0'
                    write(file_name(9:12),'(i4)')i    
                else if (i.le.99999) then
                    write(file_name(8:12),'(i5)')'i'     
                endif
                open(file_unit(i),file=file_name,status="REPLACE")
            enddo
            !
            return
        end subroutine MAKE_MONITOR_NAME
        !****************************************************************************
        ! WRITE MONITORS - ELASTIC CASE
        !****************************************************************************
        
        subroutine WRITE_MONITOR_EL(cs_nnz,cs,ne,nm,sdeg_mat,prop_mat,nmonit,option_out_var, &
            unit_disp,unit_vel,unit_acc,unit_strain,unit_stress,unit_omega,unit_uDRM,node_m,nnt,         &
            tt1,node_TOT,nnode_TOT,tagstep,dis,vel,acc,nodal_counter,alfa1,alfa2,beta1,beta2, &
            gamma1,gamma2)
            ! 
            implicit none
            ! intent IN
            real*8,                        intent(in)  :: tt1
            integer*4,                     intent(in)  :: nmonit,nnt,nm,ne,tagstep
            integer*4,                     intent(in)  :: nnode_TOT,cs_nnz
            integer*4, dimension(nm),      intent(in)  :: sdeg_mat
            integer*4, dimension(0:cs_nnz),intent(in)  :: cs
            real*8,    dimension(nm,8)  ,  intent(in)  :: prop_mat
            real*8,    dimension(ne)    ,  intent(in)  :: alfa1,beta1,gamma1
            real*8,    dimension(ne)    ,  intent(in)  :: alfa2,beta2,gamma2 
            integer*4, dimension(6),       intent(in)  :: option_out_var
            integer*4, dimension(nnt),     intent(in)  :: nodal_counter
            integer*4, dimension(nmonit),  intent(in)  :: node_m
            integer*4, dimension(nnode_TOT),intent(in) :: node_TOT
            ! intent INOUT
            real*8,    dimension(2*nnt) ,  intent(inout) :: dis,vel,acc
            integer*4, dimension(nmonit),  intent(inout) :: unit_disp
            integer*4, dimension(nmonit),  intent(inout) :: unit_vel
            integer*4, dimension(nmonit),  intent(inout) :: unit_acc
            integer*4, dimension(nmonit),  intent(inout) :: unit_stress
            integer*4, dimension(nmonit),  intent(inout) :: unit_strain
            integer*4, dimension(nmonit),  intent(inout) :: unit_omega
            integer*4, dimension(nmonit),  intent(inout) :: unit_uDRM
            integer*4                                    :: ii,jj,i,ie,im,is,in,nn
            !
            real*8, dimension(:),   allocatable :: ct,ww,dxdx_el,dxdy_el,dydx_el,dydy_el
            real*8, dimension(:,:), allocatable :: dd,ux_el,uy_el
            real*8, dimension(:,:), allocatable :: duxdx_el,duxdy_el,duydx_el,duydy_el
            real*8, dimension(:,:), allocatable :: sxx_el,syy_el,szz_el,sxy_el
            real*8, dimension(:,:), allocatable :: det_j,mu_el,lambda_el 
            real*8, dimension(nnt)              :: sxx,syy,sxy,szz 
            real*8, dimension(nnt)              :: duxdx,duydy,duxdy,duydx
            
            ! Displacements
            if (option_out_var(1).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(dis(in)).lt.(1.0d-99))     dis(in)= 0.d0
                    if (dabs(dis(in+nnt)).lt.(1.0d-99)) dis(in+nnt)=0.d0
                    write(unit_disp(i),'(3E16.6)') tt1,dis(in),dis(in+nnt)
                enddo
            endif
            ! Velocity
            if (option_out_var(2).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(vel(in)).lt.(1.0d-99))     vel(in)= 0.d0
                    if (dabs(vel(in+nnt)).lt.(1.0d-99)) vel(in+nnt)=0.d0
                    write(unit_vel(i),'(3E16.6)') tt1,vel(in),vel(in+nnt)
                enddo
            endif
            ! Acceleration
            if (option_out_var(3).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(acc(in)).lt.(1.0d-99))     acc(in)= 0.d0
                    if (dabs(acc(in+nnt)).lt.(1.0d-99)) acc(in+nnt)=0.d0
                    write(unit_acc(i),'(3E16.6)') tt1,acc(in),acc(in+nnt)
                enddo
            endif
            ! Compute strain
            if (sum(option_out_var(4:6)).gt.0) then
                duxdx = 0.d0
                duydy = 0.d0
                duxdy = 0.d0
                duydx = 0.d0
                
                if (option_out_var(4).eq.1) then
                    sxx = 0.d0
                    syy = 0.d0 
                    sxy = 0.d0
                    szz = 0.d0  
                endif

                do ie=1,ne
                    im = cs(cs(ie -1) +0)   
                    nn = sdeg_mat(im) +1

                    allocate(ct(nn))
                    allocate(ww(nn))
                    allocate(dd(nn,nn))
                    allocate(dxdx_el(nn))
                    allocate(dxdy_el(nn))
                    allocate(dydx_el(nn))
                    allocate(dydy_el(nn))
                    allocate(ux_el(nn,nn))
                    allocate(uy_el(nn,nn))
                    allocate(duxdx_el(nn,nn))
                    allocate(duxdy_el(nn,nn))
                    allocate(duydx_el(nn,nn))
                    allocate(duydy_el(nn,nn))
                    
                    
                    call LGL(nn,ct,ww,dd)
                    
                    call MAKE_DERIVATIVES(alfa1(ie),alfa2(ie),beta1(ie),beta2(ie),gamma1(ie),gamma2(ie),nn,ct,&
                        dxdx_el,dxdy_el,dydx_el,dydy_el)
                    
                    do jj=1,nn
                        do ii=1,nn
                            is = nn*(jj -1) + ii
                            in = cs(cs(ie -1) + is)
                            
                            ux_el(ii,jj)=dis(in)
                            uy_el(ii,jj)=dis(in+nnt)
                        enddo
                    enddo

                    call MAKE_STRAIN(nn,dd,dxdx_el,dxdy_el,dydx_el,dydy_el,ux_el,uy_el, &
                        duxdx_el,duxdy_el,duydx_el,duydy_el)
                    
                    do jj=1,nn
                        do ii=1,nn
                            is = nn*(jj -1) + ii
                            in = cs(cs(ie -1) + is)
                            
                            duxdx(in) = duxdx(in) + duxdx_el(ii,jj)
                            duydy(in) = duydy(in) + duydy_el(ii,jj)
                            duxdy(in) = duxdy(in) + duxdy_el(ii,jj)
                            duydx(in) = duydx(in) + duydx_el(ii,jj)
                        enddo
                    enddo
                    

                    if (option_out_var(4).eq.1) then
                        allocate(sxx_el(nn,nn))
                        allocate(syy_el(nn,nn))
                        allocate(szz_el(nn,nn))
                        allocate(sxy_el(nn,nn))
                        allocate(mu_el(nn,nn))
                        allocate(lambda_el(nn,nn))
                        lambda_el(:,:) = prop_mat(im,2)
                        mu_el(:,:)     = prop_mat(im,3)
                        
                        call MAKE_STRESS(nn,lambda_el,mu_el,duxdx_el,duxdy_el,duydx_el,duydy_el, &
                            sxx_el,syy_el,szz_el,sxy_el)
                        do jj=1,nn
                            do ii=1,nn
                                is = nn*(jj -1) + ii
                                in = cs(cs(ie -1) + is)
                                sxx(in) = sxx(in) + sxx_el(ii,jj)
                                syy(in) = syy(in) + syy_el(ii,jj)
                                szz(in) = szz(in) + szz_el(ii,jj)
                                sxy(in) = sxy(in) + sxy_el(ii,jj)
                            enddo
                        enddo
                        
                    endif
                    if(allocated(ct)) deallocate(ct)
                    if(allocated(ww)) deallocate(ww)
                    if(allocated(dd)) deallocate(dd)
                    if(allocated(dxdx_el)) deallocate(dxdx_el)
                    if(allocated(dxdy_el)) deallocate(dxdy_el)
                    if(allocated(dydx_el)) deallocate(dydx_el)
                    if(allocated(dydy_el)) deallocate(dydy_el)
                    if(allocated(ux_el   )) deallocate(ux_el   )   
                    if(allocated(uy_el   )) deallocate(uy_el   )
                    if(allocated(duxdx_el)) deallocate(duxdx_el)
                    if(allocated(duxdy_el)) deallocate(duxdy_el)
                    if(allocated(duydx_el)) deallocate(duydx_el)
                    if(allocated(duydy_el)) deallocate(duydy_el)
                    if(allocated(sxx_el   )) deallocate(sxx_el   )
                    if(allocated(syy_el   )) deallocate(syy_el   )  
                    if(allocated(szz_el   )) deallocate(szz_el   )   
                    if(allocated(sxy_el   )) deallocate(sxy_el   )
                    if(allocated(mu_el    )) deallocate(mu_el    )
                    if(allocated(lambda_el)) deallocate(lambda_el)

                enddo            
            endif

            if (option_out_var(4) .eq. 1) then

                do in = 1, nnt  
                    sxx(in) = sxx(in) / nodal_counter(in)
                    syy(in) = syy(in) / nodal_counter(in)
                    sxy(in) = sxy(in) / nodal_counter(in)
                    szz(in) = szz(in) / nodal_counter(in)
                enddo
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(sxx(in)).lt.(1.0d-99)) sxx(in)=0.0
                    if (dabs(syy(in)).lt.(1.0d-99)) syy(in)=0.0
                    if (dabs(sxy(in)).lt.(1.0d-99)) sxy(in)=0.0
                    if (dabs(szz(in)).lt.(1.0d-99)) szz(in)=0.0

                    write(unit_stress(i),'(5E16.8)') tt1,sxx(in),syy(in),sxy(in),szz(in)
                enddo
            endif


            if (option_out_var(5) .eq. 1) then  !Out Options Scandella 02.07.2007
                do in = 1, nnt  
                    duxdx(in) = duxdx(in) / nodal_counter(in)
                    duydy(in) = duydy(in) / nodal_counter(in)
                    duxdy(in) = duxdy(in) / nodal_counter(in)
                    duydx(in) = duydx(in) / nodal_counter(in) 
                enddo 

                do i = 1,nmonit 
                    in = node_m(i) 
                    if (dabs(duxdx(in)).lt.(1.0d-99)) duxdx(in)=0.0
                    if (dabs(duydy(in)).lt.(1.0d-99)) duydy(in)=0.0
                    if (dabs(duxdy(in)).lt.(1.0d-99)) duxdy(in)=0.0
                    if (dabs(duydx(in)).lt.(1.0d-99)) duydx(in)=0.0

                    write(unit_strain(i),'(4E16.8)') tt1,duxdx(in),duydy(in),0.5*(duxdy(in)+duydx(in)) 
                enddo
            endif

            if (option_out_var(6) .eq. 1) then  !Out Options Scandella 02.07.2007
                if (option_out_var(5) .ne. 1) then
                    do in = 1, nnt  
                        duxdx(in) = duxdx(in) / nodal_counter(in)
                        duydy(in) = duydy(in) / nodal_counter(in)
                        duxdy(in) = duxdy(in) / nodal_counter(in)
                        duydx(in) = duydx(in) / nodal_counter(in) 
                    enddo 
                endif

                do i = 1,nmonit 
                in = node_m(i) 
                if (dabs(duxdy(in)).lt.(1.0d-99)) duxdy(in)=0.0
                if (dabs(duydx(in)).lt.(1.0d-99)) duydx(in)=0.0

                write(unit_omega(i),'(2E16.8)') tt1, 0.5*(duxdy(in)-duydx(in)) 
                enddo
            endif

            !***********************************************************************************
            ! DRM
            !***********************************************************************************
            
            if ((nnode_TOT.ne.0).and.(tagstep.eq.1)) then  
                do i = 1,nnode_TOT                                  
                    in = node_TOT(i)                                
                    if (dabs(dis(in)).lt.(1.0d-99)) then
                        dis(in)=0.0
                    endif
                    if (dabs(dis(in+nnt)).lt.(1.0d-99)) then
                        dis(in+nnt)=0.0
                    endif
                    write(unit_uDRM(i),'(3E16.8)') &                
                        tt1,dis(in),dis(in+nnt)                       
                enddo 
            endif
            !
            return
            !
        end subroutine WRITE_MONITOR_EL
       
        !******************************************************************************************
        ! NONLINEAR OUTPUT 
        !******************************************************************************************

        subroutine WRITE_MONITOR_NL(unit_disp,unit_vel,unit_acc,unit_strain,unit_stress,unit_omega,&
            unit_uDRM,option_out_var,nmonit,ndt_monitor,node_m,nm,ne,nnt,cs,cs_nnz,sdeg_mat,snl,&
            its,tt1,node_TOT,nnode_TOT,tagstep,dis,vel,acc,nodal_counter,disout)
            !
            implicit none
            ! intent IN
            real*8,                         intent(in)  :: tt1,ndt_monitor
            integer*4,                      intent(in)  :: nm,ne,nnt,cs_nnz,tagstep
            integer*4,                      intent(in)  :: nmonit,its,nnode_TOT
            integer*4, dimension(6),        intent(in)  :: option_out_var
            integer*4, dimension(nm),       intent(in)  :: sdeg_mat
            integer*4, dimension(nnt),      intent(in)  :: nodal_counter
            integer*4, dimension(nmonit),   intent(in)  :: node_m
            integer*4, dimension(0:cs_nnz), intent(in)  :: cs
            integer*4, dimension(nnode_TOT),intent(in)  :: node_TOT
            type(nl_element), dimension(ne),intent(in)  :: snl
            ! intent INOUT
            integer*4, dimension(nmonit), intent(inout) :: unit_disp
            integer*4, dimension(nmonit), intent(inout) :: unit_vel
            integer*4, dimension(nmonit), intent(inout) :: unit_acc
            integer*4, dimension(nmonit), intent(inout) :: unit_stress
            integer*4, dimension(nmonit), intent(inout) :: unit_strain
            integer*4, dimension(nmonit), intent(inout) :: unit_omega
            integer*4, dimension(nmonit), intent(inout) :: unit_uDRM
            real*8,    dimension(2*nnt) , intent(inout) :: dis,vel,acc
            type(nodepatched),            intent(inout) :: disout
            !
            integer*4                                     ::  in,i
            real*8                                        :: sxx_out,syy_out,szz_out,sxy_out
            real*8                                        :: exx_out,eyy_out,gxy_out
            real*8                                        :: epxx_out,epyy_out,epzz_out,gpxy_out
           

            !***********************************************************************************
            ! DISPLACEMENT
            !***********************************************************************************
            
            if (option_out_var(1).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(dis(in)).lt.(1.0d-99))     dis(in)= 0.d0
                    if (dabs(dis(in+nnt)).lt.(1.0d-99)) dis(in+nnt)=0.d0
                    write(unit_disp(i),'(3E16.6)') tt1,dis(in),dis(in+nnt)
                enddo
            endif

            !***********************************************************************************
            ! VELOCITY
            !***********************************************************************************
            
            if (option_out_var(2).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(vel(in)).lt.(1.0d-99))     vel(in)= 0.d0
                    if (dabs(vel(in+nnt)).lt.(1.0d-99)) vel(in+nnt)=0.d0
                    write(unit_vel(i),'(3E16.6)') tt1,vel(in),vel(in+nnt)
                enddo
            endif
            
            !***********************************************************************************
            ! ACCELERATION
            !***********************************************************************************

            if (option_out_var(3).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(acc(in)).lt.(1.0d-99))     acc(in)= 0.d0
                    if (dabs(acc(in+nnt)).lt.(1.0d-99)) acc(in+nnt)=0.d0
                    write(unit_acc(i),'(3E16.8)') tt1,acc(in),acc(in+nnt)
                enddo
            endif

            !***********************************************************************************
            ! STRESS 
            !***********************************************************************************
            
            if (option_out_var(4).eq.1) then
                disout%stress(:,:) = 0.0d0

                call UPDATE_OUT_STRESS(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
                do i = 1,nmonit
                    in = node_m(i)
                    sxx_out = disout%stress(1,in)/nodal_counter(in)
                    syy_out = disout%stress(2,in)/nodal_counter(in)
                    szz_out = disout%stress(3,in)/nodal_counter(in)
                    sxy_out = disout%stress(4,in)/nodal_counter(in)
                    if (dabs(disout%stress(1,in)).lt.(1.0d-99)) sxx_out=0.d0
                    if (dabs(disout%stress(2,in)).lt.(1.0d-99)) syy_out=0.d0
                    if (dabs(disout%stress(3,in)).lt.(1.0d-99)) szz_out=0.d0
                    if (dabs(disout%stress(4,in)).lt.(1.0d-99)) sxy_out=0.d0
                    write(unit_stress(i),'(5E16.8)') tt1,sxx_out,syy_out,sxy_out,szz_out
                enddo
            endif
            
            
            !***********************************************************************************
            ! STRAIN 
            !***********************************************************************************
            
            if (option_out_var(5).eq.1) then
                disout%strain(:,:) = 0.0d0

                call UPDATE_OUT_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
                do i = 1,nmonit
                    in = node_m(i) 
                    exx_out = disout%strain(1,in) / nodal_counter(in)
                    eyy_out = disout%strain(2,in) / nodal_counter(in)
                    gxy_out = disout%strain(3,in) / nodal_counter(in)
                    if (dabs(disout%strain(1,in)).lt.(1.0d-99)) exx_out=0.d0
                    if (dabs(disout%strain(2,in)).lt.(1.0d-99)) eyy_out=0.d0
                    if (dabs(disout%strain(3,in)).lt.(1.0d-99)) gxy_out=0.d0
                    write(unit_strain(i),'(4E16.8)') tt1,exx_out,eyy_out,gxy_out 
                enddo
            endif

            !***********************************************************************************
            ! PLASTIC STRAIN 
            !***********************************************************************************
            
            if (option_out_var(5).eq.2) then
                disout%pstrain(:,:) = 0.0d0 

                call UPDATE_OUT_PLASTIC_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
                do i = 1,nmonit
                    in = node_m(i) 
                    epxx_out = disout%pstrain(1,in) / nodal_counter(in)
                    epyy_out = disout%pstrain(2,in) / nodal_counter(in)
                    gpxy_out = disout%pstrain(4,in) / nodal_counter(in)
                    epzz_out = disout%pstrain(3,in) / nodal_counter(in)
                    if (dabs(disout%pstrain(1,in)).lt.(1.0d-99)) epxx_out=0.d0
                    if (dabs(disout%pstrain(2,in)).lt.(1.0d-99)) epyy_out=0.d0
                    if (dabs(disout%pstrain(4,in)).lt.(1.0d-99)) gpxy_out=0.d0
                    if (dabs(disout%pstrain(3,in)).lt.(1.0d-99)) epzz_out=0.d0
                    write(unit_strain(i),'(5E16.8)') tt1,epxx_out,epyy_out,gpxy_out,epzz_out 
                enddo
            endif
            
            !***********************************************************************************
            ! STRAIN & PLASTIC STRAIN 
            !***********************************************************************************
            
            if (option_out_var(5).eq.3) then
                disout%strain(:,:) = 0.0d0
                disout%pstrain(:,:) = 0.0d0 

                call UPDATE_OUT_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
                call UPDATE_OUT_PLASTIC_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
                do i = 1,nmonit
                    in = node_m(i) 
                    exx_out = disout%strain(1,in) / nodal_counter(in)
                    eyy_out = disout%strain(2,in) / nodal_counter(in)
                    gxy_out = disout%strain(3,in) / nodal_counter(in)
                    epxx_out = disout%pstrain(1,in) / nodal_counter(in)
                    epyy_out = disout%pstrain(2,in) / nodal_counter(in)
                    gpxy_out = disout%pstrain(4,in) / nodal_counter(in)
                    epzz_out = disout%pstrain(3,in) / nodal_counter(in)
                    if (dabs(disout%strain(1,in)).lt.(1.0d-99)) exx_out=0.d0
                    if (dabs(disout%strain(2,in)).lt.(1.0d-99)) eyy_out=0.d0
                    if (dabs(disout%strain(3,in)).lt.(1.0d-99)) gxy_out=0.d0
                    if (dabs(disout%pstrain(1,in)).lt.(1.0d-99)) epxx_out=0.d0
                    if (dabs(disout%pstrain(2,in)).lt.(1.0d-99)) epyy_out=0.d0
                    if (dabs(disout%pstrain(4,in)).lt.(1.0d-99)) gpxy_out=0.d0
                    if (dabs(disout%pstrain(3,in)).lt.(1.0d-99)) epzz_out=0.d0
                    write(unit_strain(i),'(8E16.8)') tt1,exx_out,eyy_out,gxy_out,&
                        epxx_out,epyy_out,gpxy_out,epzz_out 
                enddo
            endif

            !***********************************************************************************
            ! DRM
            !***********************************************************************************
            
            if ((nnode_TOT.gt.0).and.(tagstep.eq.1)) then  
                do i = 1,nnode_TOT                                  
                    in = node_TOT(i)                                
                    if (dabs(dis(in)).lt.(1.0d-99)) then
                        dis(in)=0.0
                    endif
                    if (dabs(dis(in+nnt)).lt.(1.0d-99)) then
                        dis(in+nnt)=0.0
                    endif
                    write(unit_uDRM(i),'(3E16.8)') &                
                        tt1,dis(in),dis(in+nnt)                       
                enddo 
            endif                                           
            !
            return
            !
        end subroutine WRITE_MONITOR_NL
        !
        subroutine UPDATE_OUT_STRESS(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
            !
            implicit none
            ! intent IN 
            integer*4, intent(in)                       :: ne,nnt,cs_nnz,nm
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            type(nl_element), intent(in), dimension(ne) :: snl
            ! intent INOUT
            type(nodepatched), intent(inout)            :: disout
            !
            integer*4                                   :: ie,im,nn,in,is,i,j 
            
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                do j = 1,nn
                    do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) + is)
                        
                        disout%stress(:,in) = disout%stress(:,in) + &
                            snl(ie)%stress(:,i,j)

                    enddo
                enddo
            enddo
            return
        end subroutine UPDATE_OUT_STRESS
        !
        subroutine UPDATE_OUT_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
            !
            implicit none
            ! intent IN 
            integer*4, intent(in)                       :: ne,nnt,cs_nnz,nm
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            type(nl_element), intent(in), dimension(ne) :: snl
            ! intent INOUT
            type(nodepatched), intent(inout)            :: disout
            !
            integer*4                                   :: ie,im,nn,in,is,i,j 
            
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                do j = 1,nn
                    do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) + is)
                        
                        disout%strain(:,in) = disout%strain(:,in) + &
                            snl(ie)%strain(:,i,j)

                    enddo
                enddo
            enddo
            return
        end subroutine UPDATE_OUT_STRAIN
        !
        subroutine UPDATE_OUT_PLASTIC_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
            !
            implicit none
            ! intent IN 
            integer*4, intent(in)                       :: ne,nnt,cs_nnz,nm
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            type(nl_element), intent(in), dimension(ne) :: snl
            ! intent INOUT
            type(nodepatched), intent(inout)            :: disout
            !
            integer*4                                   :: ie,im,nn,in,is,i,j 
            
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                do j = 1,nn
                    do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) + is)
                        
                        disout%pstrain(:,in) = disout%pstrain(:,in) + &
                            snl(ie)%pstrain(:,i,j)
                    enddo
                enddo
            enddo
            return
    end subroutine UPDATE_OUT_PLASTIC_STRAIN
end module write_output
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
