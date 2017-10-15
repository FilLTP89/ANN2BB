module nonlinear2d
    !
    use fields
    !
    implicit none
    !
    real*8, parameter :: zero=0.d0,one=1.0d0
    real*8, parameter :: half=0.5d0,two=2.0d0,three=3.0d0
    !
    real*8, parameter :: FTOL = 0.0010D0
    real*8, parameter :: LTOL = 0.001D0
    real*8, parameter :: STOL = 0.00010D0
    real*8, parameter :: PSI  = one!5.0d0!one
    real*8, parameter :: OMEGA= zero!1.0d6!zero
    !
    real*8, parameter, dimension(4,4) :: MM = reshape((/ &
        one , zero, zero, zero, &
        zero, one , zero, zero, &
        zero, zero, one , zero, &
        zero, zero, zero, two   &
        /), (/4,4/))
    
    real*8, parameter, dimension(4,4) :: MM1 = reshape((/ &   
        one , zero, zero, zero, &
        zero, one , zero, zero, &
        zero, zero, one , zero, &
        zero, zero, zero, half  &
        /), (/4,4/))
    
    real*8, parameter, dimension(4) :: m = (/one,one,one,zero/)

    contains
        
        !****************************************************************************
        ! MAKE NL INTERNAL FORCES
        !****************************************************************************
        
        subroutine MAKE_INTERNAL_FORCES_NL(nnt,ne,nm,cs_nnz,cs,sdeg_mat,snl,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,displ,fk,mvec,dt)
            !
            implicit none
            ! intent IN
            integer*4, intent(in)                       :: nnt,ne,nm,cs_nnz
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            real*8,    intent(in), dimension(ne)        :: alfa1,beta1,gamma1
            real*8,    intent(in), dimension(ne)        :: alfa2,beta2,gamma2
            real*8,    intent(in), dimension(2*nnt)     :: displ,mvec
            real*8,    intent(in)                       :: dt
            ! intent INOUT
            real*8,    intent(inout), dimension(2*nnt)  :: fk
            type(nl_element), intent(inout), dimension(ne) :: snl 
            ! 
            real*8,     dimension(:),  allocatable      :: ct,ww
            real*8,     dimension(:),  allocatable      :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable      :: dd,fx,fy
            real*8,     dimension(:,:,:), allocatable   :: dstrain,dstrial
            !
            real*8, dimension(4)                        :: stress_,dstrial_,center_
            real*8, dimension(4)                        :: dstrain_,dpstrain_,pstrain_
            real*8                                      :: radius_,syld_ 
            real*8                                      :: lambda_,mu_     
            real*8                                      :: ckin_,kkin_  
            real*8                                      :: biso_,rinf_   
            !
            real*8                                      :: alpha_epl,t1ux,t1uy,t2ux
            real*8                                      :: t2uy,t1fx,t1fy,t2fx,t2fy
            integer*4                                   :: ie,ip,iq,il,im,nn,is,in
            logical                                     :: st_epl
            !
            integer*4                                   :: number_of_threads
            real*8:: stat
            !
!            number_of_threads =1;
!
!            call OMP_set_num_threads(number_of_threads)

            !call OMP_get_num_threads()


!!$OMP PARALLEL &
!!$OMP PRIVATE(ie, im, nn, iq, ip, is, in, ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial) &
!!$OMP PRIVATE(stress_,dstrial_,dstrain_,pstrain_,center_,radius_,lambda_,mu_,ckin_,kkin_,rinf_,biso_,syld_) &
!!$OMP PRIVATE(alpha_epl,st_epl,fx,fy)


!!$OMP DO 
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                !*********************************************************************************
                ! COMPUTE STRAIN 
                !*********************************************************************************
                
                call ALLOINIT_LOC_NL(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial, &
                    fx,fy,displ,alfa1(ie),alfa2(ie),beta1(ie),beta2(ie),gamma1(ie),gamma2(ie),&
                    snl(ie)%lambda,snl(ie)%mu)
                ! VELOCITY FORMULATION
!                dstrain = dstrain*dt
!                dstrial = dstrial*dt
                !snl(ie)%strain(:,:,:) = snl(ie)%strain(:,:,:)+dstrain(:,:,:)
                !snl(ie)%stress(:,:,:) = snl(ie)%stress(:,:,:)+dstrial(:,:,:) 
                
!                ! DISPLACEMENT FORMULATION (ELASTIC)
!                snl(ie)%strain(:,:,:) = dstrain
!                snl(ie)%stress(:,:,:) = dstrial

                ! DISPLACEMENT FORMULATION (PLASTIC)
                dstrain(:,:,:) = dstrain(:,:,:) - snl(ie)%strain(:,:,:)
                dstrial(:,:,:) = dstrial(:,:,:) - snl(ie)%stress(:,:,:)!
                snl(ie)%strain(:,:,:) = snl(ie)%strain(:,:,:)+dstrain(:,:,:)
!                snl(ie)%stress(:,:,:) = snl(ie)%stress(:,:,:)+dstrial(:,:,:) 
                
                !*********************************************************************************
                ! COMPUTE STRESS
                !*********************************************************************************
                
                do iq = 1,nn
                    do ip = 1,nn
                        is = nn*(iq -1) + ip
                        ! STARTING POINT
                        stress_  = snl(ie)%stress(:,ip,iq)
                        center_  = snl(ie)%center(:,ip,iq)
                        radius_  = snl(ie)%radius(ip,iq)
                        pstrain_ = snl(ie)%pstrain(:,ip,iq)
                        lambda_  = snl(ie)%lambda(ip,iq)
                        syld_    = snl(ie)%syld(ip,iq)
                        ckin_    = snl(ie)%ckin(ip,iq)
                        kkin_    = snl(ie)%kkin(ip,iq)
                        biso_    = snl(ie)%biso(ip,iq)
                        rinf_    = snl(ie)%rinf(ip,iq)
                        mu_      = snl(ie)%mu(ip,iq)
                        
                        ! STRAIN INCREMENT
                        dstrain_(:)   = zero
                        dstrain_(1:2) = dstrain(1:2,ip,iq)
                        dstrain_(4)   = dstrain(3,ip,iq)
                        ! TRIAL STRESS INCREMENT
                        dstrial_(:)   = zero
                        dstrial_(:)   = dstrial(:,ip,iq)
                        
                        ! CHECK PLASTICITY
                        call check_plasticity(dstrial_,stress_,center_,radius_,&
                            syld_,st_epl,alpha_epl)
                        ! PLASTIC CORRECTION
                        if (st_epl) then
                            dstrain_(:) = (one-alpha_epl)*dstrain_(:)
                            call plastic_corrector(dstrain_,dstrial_,center_,radius_,syld_,&
                                biso_,rinf_,ckin_,kkin_,mu_,lambda_,pstrain_)
                        endif
                        ! STRESS VECTOR
                        snl(ie)%stress(:,ip,iq) = dstrial_(:)
                        ! CENTER
                        snl(ie)%center(:,ip,iq) = center_(:)
                        ! RADIUS
                        snl(ie)%radius(ip,iq)   = radius_
                        ! PLASTIC STRAIN VECTOR
                        snl(ie)%pstrain(:,ip,iq) = pstrain_(:)
                        
                    enddo
                enddo
              
               !*********************************************************************************
               ! COMPUTE LOCAL INTERNAL FORCES
               !*********************************************************************************
                ! DISPLACEMENT FORMULATION (ELASTIC)
                call MAKE_INTERNAL_FORCE(nn,ww,dd,dxdx,dxdy,dydx,dydy,dstrial(1,:,:),&
                    dstrial(2,:,:),dstrial(4,:,:),fx,fy)
               ! DISPLACEMENT FORMULATION (PLASTIC)
                call MAKE_INTERNAL_FORCE(nn,ww,dd,dxdx,dxdy,dydx,dydy,snl(ie)%stress(1,:,:),&
                    snl(ie)%stress(2,:,:),snl(ie)%stress(4,:,:),fx,fy)
                
                !*********************************************************************************
                ! COMPUTE LOCAL INTERNAL FORCES
                !*********************************************************************************
                do iq = 1,nn
                    do ip = 1,nn
                        is = nn*(iq-1)+ip
                        in = cs(cs(ie-1)+is)
                        fk(in)     = fk(in)     + fx(ip,iq)/mvec(in)
                        fk(in+nnt) = fk(in+nnt) + fy(ip,iq)/mvec(in+nnt)
                    enddo
                enddo
                ! DEALLOCATE ELEMENT-WISE VARIABLES
                call DEALLOCATE_LOC(ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial,fx,fy)
            enddo
!!$OMP END DO
!!$OMP END PARALLEL 
            return
            !
        end subroutine MAKE_INTERNAL_FORCES_NL
        
        !*********************************************************************************
        ! STIFFNESS MATRIX 
        !*********************************************************************************
        
        subroutine STIFF_MATRIX(lambda,mu,DEL)
            implicit none
            ! intent IN
            real*8, intent(in) :: lambda,mu
            ! intent OUT
            real*8, intent(out), dimension(4,4) :: DEL
            !
            DEL(:,:) = zero
            DEL(1,1) = lambda+2*mu
            DEL(2,2) = lambda+2*mu
            DEL(3,3) = lambda+2*mu
            DEL(4,4) = mu
            DEL(1,2) = lambda
            DEL(1,3) = lambda
            DEL(2,1) = lambda
            DEL(2,3) = lambda
            DEL(3,1) = lambda
            DEL(3,2) = lambda
            !
            return
        end subroutine STIFF_MATRIX
        
        !*********************************************************************************
        ! CRITICAL STATE STIFFNESS MATRIX 
        !*********************************************************************************
        
        subroutine STIFF_MATRIX_CRITICAL(stress0,dstrain,lambda,mu,DEL)
            implicit none
            ! intent IN
            real*8, intent(in) :: lambda,mu
            real*8, dimension(4) :: stress0,dstrain
            ! intent OUT
            real*8, intent(out), dimension(4,4) :: DEL
            !
            real*8 :: p0,devol,B,k,nu,spec_vol,mu_crit,lambda_crit
           
            spec_vol = one+0.7d0
            ! pressure 
            p0 = dot_product(stress0,m)/three
            ! volumetric strain increment
            devol = dot_product(dstrain,m)
            ! Bulk's modulus
            B = p0/devol*(exp(spec_vol*devol/k)-one)
            ! Poisson's ratio
            nu = half*lambda/(lambda+mu)
            ! NEW shear modulus
            mu_crit = three*half*(one-two*nu)*B/(one+nu)
            ! NEW lambda 
            lambda_crit = two*mu_crit*nu/(one-two*nu)
            call STIFF_MATRIX(lambda_crit,mu_crit,DEL)
            !
            return
        end subroutine STIFF_MATRIX_CRITICAL
        
        !****************************************************************************
        ! MISES YIELDING LOCUS AND GRADIENT
        !****************************************************************************

        subroutine mises_yld_locus(stress, center, radius, syld, FM, gradFM)
            ! 
            implicit none
            ! intent IN
            real*8,               intent(in)    :: radius,syld
            real*8, dimension(4), intent(in)    :: stress,center
            ! intent OUT
            real*8, dimension(4), intent(out)   :: gradFM
            real*8,               intent(out)   :: FM
            !
            real*8, dimension(4)                :: dev
            real*8                              :: tau_eq
            ! COMPUTE TENSOR COMPONENTS 
            call tensor_components (stress, dev)
            ! COMPUTE MISES FUNCTION
            dev(:) = dev(:) - center(:)
            call tau_mises(dev,tau_eq)
            FM = tau_eq - syld - radius
            ! COMPUTE MISES FUNCTION GRADIENT
            gradFM(1) = three*half*dev(1)/tau_eq
            gradFM(2) = three*half*dev(2)/tau_eq
            gradFM(3) = three*half*dev(2)/tau_eq
            gradFM(4) = three*dev(4)/tau_eq
            !  
            return
            !
        end subroutine mises_yld_locus
        
        !****************************************************************************
        ! OCTAHEDRAL SHEAR STRESS
        !****************************************************************************
        
        subroutine tau_mises(dev,J2M)
            ! correspondence rule:
            !  2-d case     3-d case  
            ! (11) <=> sxx   (11) <=> sxx 
            ! (22) <=> syy   (22) <=> syy 
            ! (33) <=> szz   (33) <=> szz 
            ! (12) <=> sxy   (12) <=> sxy
            !                (23) <=> syz
            !                (31) <=> szx
            ! and continuum mechanics sign conventions
            !
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: dev
            ! intent OUT
            real*8,               intent(out)   :: J2M
            !
            real*8, dimension(4)                :: temp
            real*8                              :: J2M2
            
            ! COMPUTE OCTAHEDRAL SHEAR STRESS 
            temp = three*half*matmul(MM,dev)
            J2M2 = dot_product(dev,temp)
            J2M  = sqrt(J2M2)
            !
            return
            ! 
        end subroutine
        
        !****************************************************************************
        ! TENSOR COMPONENTS (SPHERICAL & DEVIATORIC)
        !****************************************************************************
        
        subroutine tensor_components(stress, dev)
            !     
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: stress
            ! intent OUT
            real*8, dimension(4), intent(out)   :: dev
            !
            real*8                              :: I1 
            !
            I1 = dot_product(stress,m)/three
            !
            dev = stress - I1*m
            !
            return
            !
        end subroutine tensor_components
            
        !****************************************************************************
        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
        !****************************************************************************
        
        subroutine check_plasticity(dtrial, stress0, center, radius, syld, &
            st_epl, alpha_epl) 
            ! ***** CRITICAL STATE EXTENSION *****
            ! st_epl,alpha_epl,dstrain,lambda,mu)    
            implicit none
            ! intent IN
            real*8,                 intent(in)      :: radius,syld
            real*8, dimension(4),   intent(in)      :: center,stress0
            ! intent INOUT
            real*8, dimension(4),   intent(inout)   :: dtrial
            ! intent OUT
            real*8,                 intent(out)     :: alpha_epl
            logical,                intent(out)     :: st_epl
            !
            logical                                 :: flagxit
            real*8                                  :: FS,FT,checkload
            real*8, dimension(4)                    :: gradFS, gradFT, stress1
            ! ***** CRITICAL STATE EXTENSION *****
            ! real*8, dimension(4), intent(in) :: dstrain
            ! real*8,               intent(in) :: lambda,mu
            !
            ! PREDICTION STRESS
            call update_stress(stress0,stress1,dtrial)
            ! ***** CRITICAL STATE EXTENSION *****
            ! call update_stress(stress0,stress1,dstrain,lambda,mu)
            ! 
            ! CHECK MISES FUNCTION
            call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
            call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)
           
            alpha_epl = zero
            if (FT.le.FTOL) then
                alpha_epl = one
                st_epl = .false.
                flagxit = .true.
            endif

            if ((FS.lt.-FTOL).and.(FT.gt.FTOL)) then
                st_epl = .true.
                call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_epl)
                ! ***** CRITICAL STATE EXTENSION *****
                ! call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_epl,dstrain,lambda,mu)
                flagxit = .true.
            endif

            if ((abs(FS).le.FTOL).and.(FT.gt.FTOL)) then
                ! CHECK LOAD DIRECTION 
                checkload = dot_product(gradFS,dtrial)/&
                    sqrt(dot_product(gradFS,gradFS)*dot_product(dtrial,dtrial))

                if (checkload.ge.-LTOL) then! PLASTIC LOADING  
                    alpha_epl = zero
                elseif (checkload.lt.-LTOL)then! PLASTIC UNLOADING  
                    call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_epl)
                    ! ***** CRITICAL STATE EXTENSION *****
                    ! call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_epl,dstrain,lambda,mu)
                endif
                st_epl = .true.
                flagxit = .true.
            endif
            
            if (.not.flagxit)then
                write(*,*) "ERROR IN FINDING INTERSECTION!!  F = ",FS,FT
                alpha_epl=zero
            endif

            ! ON-LOCUS STRESS STATE
            call update_stress(stress0,stress1,alpha_epl*dtrial)
            dtrial = stress1
            call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
            !
            return
        end subroutine check_plasticity
        
        !****************************************************************************
        ! CORRECT STRESS STATE
        !****************************************************************************
        
        subroutine plastic_corrector(dEps_alpha,stress,center,syld, &
            radius,biso,Rinf,Ckin,kkin,mu,lambda,pstrain)
            !
            implicit none
            !
            real*8,               intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
            !
            real*8, dimension(4), intent(inout) :: dEps_alpha,pstrain 
            real*8, dimension(4), intent(inout) :: stress,center
            real*8,               intent(inout) :: radius 
            real*8, dimension(4,4)              :: DEL
            real*8                              :: Ttot,deltaTk,qq,R1,R2,dR1,dR2
            real*8                              :: err0,err1,err2,err3
            real*8                              :: FM,hard0,hard1,hard2,deltaTmin
            real(8)                             :: Resk
            logical                             :: flag_fail
            real*8, dimension(4)                :: gradFM,S1,S2,X1,X2,Epl1
            real*8, dimension(4)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2
            integer                             :: counter
            call stiff_matrix(lambda,mu,DEL)
            deltaTk = one
            Ttot    = zero
            deltaTmin = 0.001D0
            flag_fail =.false.
            counter = 1
            do while ((Ttot.lt.one).and.counter.le.10)
                Resk  = zero
                dS1   = zero
                dX1   = zero
                dS2   = zero
                dX2   = zero
                dR1   = zero
                dR2   = zero
                dEpl1 = zero
                dEpl2 = zero
                ! FIRST ORDER COMPUTATION
                call ep_integration(dEps_alpha*deltaTk,stress,center,radius,syld,&
                    mu,lambda,biso,Rinf,Ckin,kkin,dS1,dX1,dR1,dEpl1,hard1,pstrain)

                S1 = stress + dS1
                X1 = center + dX1 
                R1 = radius + dR1
                Epl1 = pstrain + dEpl1 
                
                ! SECOND ORDER COMPUTATION
                call ep_integration(dEps_alpha*deltaTk,S1,X1,R1,syld,mu,lambda,&
                    biso,Rinf,Ckin,kkin,dS2,dX2,dR2,dEpl2,hard2,Epl1)
                
                ! TEMPORARY VARIABLES
                S1 = stress + half*(dS1+dS2)
                X1 = center + half*(dX1+dX2)
                R1 = radius + half*(dR1+dR2)
                Epl1 = pstrain + half*(dEpl1+dEpl2)
                
                ! ERROR
                call tau_mises(dS2-dS1,err0)
                call tau_mises(S1,err1)
                call tau_mises(dX1-dX2,err2)
                call tau_mises(X1,err3)
                Resk = max(epsilon(Resk),half*err0/err1,half*err2/err3)
               
                ! CHECK CONVERGENCE
                if (Resk.le.STOL) then
                    stress = S1
                    center = X1
                    radius = R1
                    pstrain = Epl1

                    call mises_yld_locus (stress, center,radius,syld,FM,gradFM)
                    if (abs(FM).gt.FTOL) then
                        call drift_corr(stress,center,radius,syld,&
                                biso,Rinf,Ckin,kkin,lambda,mu,pstrain)
                    endif
                    qq = min(0.9d0*sqrt(STOL/Resk),1.1d0)
                    if (flag_fail) then
                        qq = min(qq,one)
                    endif
                    flag_fail=.false.
                    counter = 1
                    Ttot=Ttot+deltaTk
                    deltaTk=qq*deltaTk
                    deltaTk=max(deltaTk,deltaTmin)
                    deltaTk=min(deltaTk,one-Ttot)
                else
                    counter = counter+1
                    qq=max(0.9d0*sqrt(STOL/Resk),0.1d0)
                    deltaTk=max(qq*deltaTk,deltaTmin)
                    flag_fail=.true.
                end if
            end do
            if (counter.eq.10)then
                write(*,*) "FAILED CORRECTION"
                stop
            endif
            !
            return
            !
        end subroutine plastic_corrector

        !****************************************************************************
        ! ELASTO-PLASTIC INTEGRATOR
        !****************************************************************************
        
        subroutine ep_integration(dstrain,stress,center,radius,syld,mu,lambda,biso,rinf,&
            ckin,kkin,dStress,dcenter,dradius,dpstrain,hard,pstrain)
            !
            implicit none
            ! intent IN
            real*8              , intent(in)    :: radius,syld,mu,lambda,biso,rinf,ckin,kkin
            real*8, dimension(4), intent(in)    :: dstrain,stress,center,pstrain
            ! intent INOUT
            real*8, dimension(4), intent(inout) :: dstress,dcenter,dpstrain
            real*8,               intent(inout) :: dradius
            ! intent OUT
            real*8,               intent(out)   :: hard
            !
            real*8, dimension(4)                :: gradF
            real*8, dimension(4,4)              :: DEL
            real*8                              :: FM,dPlast
            
            ! PREDICTION
            call mises_yld_locus (stress,center,radius,syld,FM,gradF)
            call stiff_matrix(lambda,mu,DEL)
            ! ***** CRITICAL STATE EXTENSION *****
            ! call stiff_matrix_critical(stress,dstrain,lambda,mu,DEL)

            ! PLASTIC MULTIPLIER
            call compute_plastic_modulus(dstrain,stress,center,radius,mu,lambda,syld, &
                biso,rinf,ckin,kkin,dPlast,hard,pstrain)
            
            ! INCREMENTS
            call hardening_increments(stress,radius,center,syld, &
                biso,Rinf,ckin,kkin,dradius,dcenter,pstrain)
            
            dradius  = dPlast*dradius
            dcenter  = dPlast*dcenter
            dpstrain = dPlast*gradF
            dstress = matmul(DEL,dstrain-dpstrain) 
            ! 
            return
            !
        end subroutine ep_integration
        
        !****************************************************************************
        ! PLASTIC MULTIPLIER
        !****************************************************************************

        subroutine compute_plastic_modulus(dstrain,stress,center,radius,mu,lambda, &
            syld,biso,rinf,ckin,kkin,dPlastM,hard,pstrain)
            !
            implicit none
            ! intent IN
            real*8,                 intent(in) :: mu,lambda,syld   
            real*8,                 intent(in) :: radius,biso,Rinf,Ckin,kkin
            real*8, dimension(4),   intent(in) :: dstrain,stress,center
            real*8, dimension(4),   intent(in) :: pstrain
            ! intent OUT
            real*8,                 intent(out):: dPlastM,hard
            !
            real*8                             :: Ah,FM
            real*8, dimension(4)               :: gradF,tempv
            real*8, dimension(4,4)             :: DEL
            real*8                             :: PHI,PlastM

            call stiff_matrix(lambda,mu,DEL)
            ! ***** CRITICAL STATE EXTENSION *****
            ! call stiff_matrix_critical(stress,dEps,lambda,mu,DEL)
            call mises_yld_locus(stress,center,radius,syld,FM,gradF)
            
            PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
            PHI  = one+(PSI-one)*exp(-OMEGA*PlastM)
            hard = PHI*Ckin+biso*(rinf-radius)
            hard = hard-kkin*dot_product(center,gradF)
            
            dPlastM = zero 
            
            tempv = matmul(DEL,tempv)
            Ah    = dot_product(tempv,gradF)
            tempv = matmul(DEL,dstrain)
            dPlastM = dot_product(gradF,tempv)
            dPlastM = max(zero,dPlastM/(hard+Ah))
            !
            return
        end subroutine compute_plastic_modulus

        !****************************************************************************
        ! HARDENING INCREMENTS
        !****************************************************************************

        subroutine hardening_increments(stress,radius,center,syld, &
            biso,rinf,ckin,kkin,dradius,dcenter,pstrain)

            ! INCREMENTS OF INTRINSIC STATIC VARIABLES
            ! intent IN
            real*8,               intent(in) :: syld,radius
            real*8,               intent(in) :: biso,rinf,ckin,kkin
            real*8, dimension(4), intent(in) :: stress,center
            real*8, dimension(4), intent(in) :: pstrain
            ! intent OUT
            real*8,               intent(out):: dradius
            real*8, dimension(4), intent(out):: dcenter
            ! 
            real*8                           :: FM,PlastM,PHI
            real*8, dimension(4)             :: gradFM

            ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
            dradius = biso*(rinf-radius)

            ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
            PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
            PHI = 1+(PSI-1)*exp(-OMEGA*PlastM)
            call mises_yld_locus (stress,center,radius,syld,FM,gradFM)
            dcenter = (PHI*ckin*two/three)*matmul(MM1,gradFM)-center*kkin
            !
            return
            !
        end subroutine hardening_increments

        !****************************************************************************
        ! DRIFT CORRECTION
        !****************************************************************************
        
        subroutine drift_corr(stress,center,radius,syld, &
            biso,Rinf,Ckin,kkin,lambda,mu,pstrain)

            ! DRIFT CORRECTION (RADIAL RETURN)
            real*8, dimension(4), intent(inout) :: stress,center,pstrain
            real*8,               intent(inout) :: radius
            real*8,               intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
            real*8                              :: F0,F1,beta,hard,radiust,PlastM,PHI
            real*8, dimension(4)                :: tempv,gradF0,gradF1,dstress,stresst,centert,pstraint
            real*8, dimension(4,4)              :: DEL
            integer*4                           :: counter,k,j
            real*8, parameter :: FTOL_DRIFT =   0.000001D0 
            ! INITIAL PLASTIC CONDITION
            call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
            call stiff_matrix(lambda,mu,DEL)
            ! ***** CRITICAL STATE EXTENSION *****
            ! call STIFF_MATRIX_CRITICAL(stress0,dincrement,lambda,mu,DEL)
            do counter=1,10 
                ! MISES FUNCTION
                call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
                ! COMPUTE HARDENING INCREMENTS
                PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
                PHI = one+(PSI-one)*exp(-OMEGA*PlastM)
                hard = biso*(Rinf-radius)
                hard = hard + PHI*Ckin
                hard = hard - kkin*dot_product(gradF0,center)
                ! COMPUTE BETA FOR DRIFT CORRECTION
                beta  = zero
                tempv = matmul(DEL,tempv)
                beta  = dot_product(tempv,gradF0)
                beta = F0/(hard+beta)
                ! STRESS-STRAIN-HARDENING CORRECTION
                dstress = zero
                dstress = -beta*tempv
                call update_stress(stress,stresst,dstress)
                centert = center+beta*((two*PHI*ckin/three)*matmul(MM1,gradF0)-center*kkin)
                radiust = radius+beta*(Rinf-radius)*biso
                ! CHECK DRIFT
                call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
                if (abs(F1).gt.abs(F0)) then
                    beta    = F0/dot_product(gradF0,gradF0)
                    dstress = -beta*gradF0
                    !stresst = stress+dstress
                    call update_stress(stress,stresst,dstress)
                    centert = center
                    radiust = radius
                    pstraint = pstrain
                    call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
                endif
                stress = stresst
                center = centert
                radius = radiust
                pstrain = pstrain+beta*gradF0
                
                if (abs(F1).le.FTOL_DRIFT) then
                    exit
                endif
            enddo
            !
            return
            !
        end subroutine drift_corr

        !****************************************************************************
        ! UPDATE STRESS
        !****************************************************************************
        
        subroutine update_stress(stress0,stress1,dincrement,lambda,mu)
            implicit none
            ! intent IN
            real*8, optional,     intent(in) :: lambda,mu
            real*8, dimension(4), intent(in) :: stress0,dincrement
            ! intent INOUT
            real*8, dimension(4), intent(inout) :: stress1
            !
            real*8, dimension(4,4)           :: DEL
            !
            stress1 = zero
            if (present(mu).and.present(lambda)) then 
                ! dincrement = strain increment
                call STIFF_MATRIX_CRITICAL(stress0,dincrement,lambda,mu,DEL)
                stress1 = stress0 + matmul(DEL,dincrement)
            else
                ! dincrement = stress increment
                stress1 = zero
                stress1 = stress0 + dincrement
            endif
            return
        end subroutine update_stress
        
        !****************************************************************************
        ! FIND INTERSECTION
        !****************************************************************************
        
        subroutine gotoFpegasus(start0,dtrial,center,radius,s0,nsub,alpha)
            ! ***** CRITICAL STATE EXTENSION *****
            ! dstrain,lambda,mu)
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: start0,dtrial,center
            real*8,               intent(in)    :: radius,s0
            integer*4,            intent(in)    :: nsub
            ! intent OUT
            real*8,               intent(out)   :: alpha
            real*8, dimension(4)                :: stress0,stress1,stress,gradF
            real*8                              :: dalpha,alpha0,alpha1,F0,F1,FM,Fsave
            integer*4                           :: counter0,counter1
            logical                             :: flagxit
            ! ***** CRITICAL STATE EXTENSION *****
            ! real*8, dimension(4), intent(in) :: dstrain
            ! real*8, intent(in) :: lambda,mu
            alpha0  = zero
            alpha1  = one
            call update_stress(start0,stress0,alpha0*dtrial)
            call update_stress(start0,stress1,alpha1*dtrial)
            ! ***** CRITICAL STATE EXTENSION *****
            !call update_stress(start0,stress0,alpha0*dstrain,lambda,mu)
            !call update_stress(start0,stress1,alpha1*dstrain,lambda,mu)
            
            call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
            call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
            stress = zero
            ! LOAD REVERSAL
            if (nsub.gt.1)then
                Fsave=F0
                do counter0=1,3
                    dalpha = (alpha1-alpha0)/nsub
                    flagxit=.false.
                    do counter1=1,nsub
                        alpha  = alpha0+dalpha
                        call update_stress(start0,stress,alpha*dtrial)
                        !stress = start0+alpha*dtrial
                        ! ***** CRITICAL STATE EXTENSION *****
                        ! call update_stress(start0,stress,alpha*dstrain,lambda,mu)
                        call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                        if (FM.gt.FTOL) then
                            alpha1=alpha
                            if (F0.lt.-FTOL) then
                                F1=FM
                                flagxit=.true.
                            else
                                alpha0=zero
                                F0=Fsave
                            endif
                            exit ! exit loop counter1=1,nsub
                        else
                            alpha0=alpha
                            F0=FM
                        endif
                    end do
                    if (flagxit) then
                        exit ! exit loop counter0=1,3
                    endif
                end do
                call update_stress(start0,stress0,alpha0*dtrial)
                call update_stress(start0,stress1,alpha1*dtrial)  
                ! ***** CRITICAL STATE EXTENSION *****
                ! call update_stress(start0,stress0,alpha0*dstrain,lambda,mu)
                ! call update_stress(start0,stress1,alpha1*dstrain,lambda,mu)
                call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
                call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
            end if
            
            ! ORIGINAL PEGASUS ALGORITHM
            do counter0=1,10
                alpha  = alpha1 - F1*(alpha1 - alpha0)/(F1-F0)
                call update_stress(start0,stress,alpha*dtrial)
                ! ***** CRITICAL STATE EXTENSION *****
                !call update_stress(start0,stress,alpha*dstrain,lambda,mu)
                
                call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                if (abs(FM).le.FTOL) then
                    exit
                else
                    if(FM*F0.lt.zero) then
                        alpha1=alpha0
                        F1=F0
                    else
                        F1=F1*F0/(F0+FM)
                    endif
                    F0=FM
                    alpha0=alpha
                endif
            end do
            ! 
            return
            !
        end subroutine gotoFpegasus
        !
!        subroutine gotoFtan(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
!            real*8, dimension(0:5), intent(in)    :: start0,dtrial0,center
!            real*8,                 intent(in)    :: radius,s0
!            real*8,                 intent(inout) :: F0,Ftrial
!            real*8,                 intent(out)   :: alpha
!            real*8, dimension(0:5)                :: start,gradF,dev,dev_temp,temp,dev0
!            real*8                                :: Fstart,err0,err1,beta
!            integer                             :: counter
!            call tensor_deviator(start0,dev0)
!            alpha=F0/(F0-Ftrial)
!            start(0:5)=start0(0:5)
!            temp=start+alpha*dtrial0
!            do counter=0,9
!                call tensor_deviator(temp, dev_temp)
!                call tensor_deviator(start,dev)
!                call tau_mises(-dev+dev_temp,err0)
!                call tau_mises(dev-dev0,err1)
!                err0=err0/err1
!                start(0:5)=temp(0:5)
!                call mises_yld_locus(start,center,radius,s0,Fstart,gradF)
!                if (abs(Fstart).le.FTOL .or. err0.lt.0.000001) then
!                    exit
!                end if
!                beta=sum(gradF(0:5)*dtrial0(0:5))
!                alpha=alpha-Fstart/beta
!                temp(0:5)=start(0:5)-(Fstart/beta)*dtrial0(0:5)
!            end do
!        end subroutine gotoFtan
!
!        subroutine gotoFsec(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
!            real*8, dimension(0:5), intent(in)    :: start0,dtrial0,center
!            real*8,                 intent(in)    :: radius,s0
!            real*8,                 intent(inout) :: F0,Ftrial
!            real*8,                 intent(out)   :: alpha
!            real*8, dimension(0:5)                :: start,gradF,dev,dev_temp,temp,dev0
!            real*8                                :: alphanew,err0,err1,beta
!            integer                             :: counter
!            call tensor_deviator(start0,dev0)
!            alpha=F0/(F0-Ftrial)
!            beta=0d0
!            alphanew=0d0
!            start(0:5)=start0(0:5)
!            temp(0:5)=start(0:5)+alpha*dtrial0(0:5)
!            do counter=0,9
!                call tensor_deviator(temp, dev_temp)
!                call tensor_deviator(start,dev)
!                call tau_mises(-dev+dev_temp,err0)
!                call tau_mises(dev-dev0,err1)
!                err0=err0/err1
!                call mises_yld_locus(start,center,radius,s0,F0,gradF)
!                call mises_yld_locus(temp,center,radius,s0,Ftrial,gradF)
!                start(0:5)=temp(0:5)
!                if (abs(Ftrial).le.FTOL .or. err0.lt.0.0000001) exit
!                alphanew=alpha-(Ftrial/(Ftrial-F0))*(alpha-beta)
!                beta=alpha
!                alpha=alphanew
!                temp(0:5)=start(0:5)-(alpha-beta)*dtrial0(0:5)
!
!            enddo
!
!        end subroutine gotoFsec
end module nonlinear2d
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
