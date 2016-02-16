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
!    along with SPEED.  If not, see <http://wqw.gnu.org/licenses/>.

!> @brief Module to perform nonlinear integration
!! @author Filippo Gatti
!> @date February,2016
!> @version 1.0

module nonlinear2d

contains
    
    ! MISES YIELDING LOCUS AND GRADIENT
    subroutine mises_yld_locus(Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
        implicit none
        real, dimension(0:2), intent(in)  :: Sigma_ij,X_ij
        real,                 intent(in)  :: R,sigma_yld
        real, dimension(0:2), intent(out) :: gradF_mises ! Mises yield locus gradient
        real,                 intent(out) :: F_mises     ! Mises yield locus
        real, dimension(0:2), parameter   :: A = (/1.0,1.0,2.0/)
        real, dimension(0:2)              :: Sigma_ij_dev
        real                              :: tau_eq_mises
        integer                           :: k
        
        call tensor_components (Sigma_ij, Sigma_ij_dev)
        call tau_mises(Sigma_ij_dev-X_ij, tau_eq_mises)
        F_mises             =   tau_eq_mises-sigma_yld-R
        gradF_mises(0:2)    =   0d0
        do k=0,2
            gradF_mises(k) = 0.5*A(k)*(Sigma_ij_dev(k)-X_ij(k))/tau_eq_mises
        end do
    
    end subroutine mises_yld_locus
    
    ! OCTAHEDRAL SHEAR STRESS
    subroutine tau_mises(tensor,tau_eq_mises)
        implicit none
        real, dimension(0:2), intent(in) :: tensor
        real,                 intent(out):: tau_eq_mises
        real, dimension(0:2), parameter  :: A = (/1.0,1.0,2.0/) 
        integer                          :: k
        
        tau_eq_mises = 0.0d0
        do k=0,2
            tau_eq_mises = tau_eq_mises+A(k)*(tensor(k)**2)
        end do
        tau_eq_mises = sqrt(0.5*tau_eq_mises)

    end subroutine
    
    ! TENSOR COMPONENTS (SPHERICAL & DEVIATORIC)
    subroutine tensor_components (Sigma_ij, Sigma_ij_dev)
        implicit none
        real, dimension(0:2), intent(in)  :: Sigma_ij
        real, dimension(0:2), intent(out) :: Sigma_ij_dev
        real                              :: Sigma_P
        integer :: k

        Sigma_ij_dev(0:2) = Sigma_ij(0:2)
        Sigma_P=sum(Sigma_ij(0:1))*0.5
        do k=0,1
            Sigma_ij_dev(k) = Sigma_ij_dev(k)-Sigma_P
        end do

    end subroutine tensor_components
        
    ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
    subroutine check_plasticity (dSigma_ij_trial, Sigma_ij_start, X_ij, R, sigma_yld, &
        st_elp, alpha_elp, nx,ny,nz,nelement)
        implicit none
        real,                 intent(in)    :: R,sigma_yld
        real, dimension(0:2), intent(in)    :: X_ij,Sigma_ij_start
        integer, optional,    intent(in)    :: nx,ny,nz,nelement
        real, dimension(0:2), intent(inout) :: dSigma_ij_trial
        real,                 intent(out)   :: alpha_elp
        integer,              intent(out)   :: st_elp
        real, dimension(0:2), parameter     :: A=(/1.0,1.0,2.0/)
        real, dimension(0:2)                :: gradFstart, gradFtrial, Sigma_ij_trial
        real                                :: Fstart, Ftrial,checkload
        integer                             :: GLLx,GLLy,GLLz,k

        ! Yield function at Sigma_ij_start
        call mises_yld_locus (Sigma_ij_start, X_ij, R, sigma_yld, &
            Fstart, gradFstart)
        ! Yield function at Sigma_ij_trial
        Sigma_ij_trial=Sigma_ij_start+dSigma_ij_trial
        call mises_yld_locus(Sigma_ij_trial, X_ij, R, sigma_yld, &
            Ftrial, gradFtrial)

        write(*,*) "*********************************"
        write(*,*) "Fstart:",Fstart,"Ftrial:",Ftrial
        ! LOADING CONDITION    
        checkload=0d0
        do k=0,2
            checkload=checkload+10*gradFstart(k)*dSigma_ij_trial(k)
        end do
        ! KKT CONDITION
        if (abs(Fstart) .le. tol_nl) then
            if (checkload .ge. 0) then
                alpha_elp = 0d0
                st_elp    = 1
                write(*,*) "PURE-PLASTIC STEP"
            else
                if (Ftrial .gt. tol_nl) then
                    alpha_elp = 2*sigma_yld/(2*sigma_yld+Ftrial)
                    st_elp    = 1
                    write(*,*) "ELASTO-PLASTIC WITH REVERSAL STEP"
                else
                    alpha_elp = 1d0
                    st_elp    = 2 
                    write(*,*) "PURE ELASTIC WITH REVERSAL STEP"
                endif
            end if
        elseif (Fstart .lt. -tol_nl) then
            if (Ftrial .lt. tol_nl) then
                alpha_elp = 1d0
                st_elp    = 2
                write(*,*) "PURE ELASTIC STEP"
            else
                call gotoFtan(Sigma_ij_start,dSigma_ij_trial,Fstart,Ftrial,X_ij,R,sigma_yld,alpha_elp)
                st_elp    = 1
                write(*,*) "ELASTO-PLASTIC STEP"
            end if
        elseif (Fstart .gt. tol_nl) then
            write(*,*) "ERROR!"
            write(*,*) "Fstart: ",Fstart,">",tol_nl,"!!!!"
            write(*,*) "ERROR!"
            stop
        end if
        dSigma_ij_trial(0:2)=Sigma_ij_start(0:2)+dSigma_ij_trial(0:2)*alpha_elp
        call mises_yld_locus(dSigma_ij_trial,X_ij,R,sigma_yld,Fstart,gradFstart)
        write(*,*) "check plasticity:",Fstart,"<=",tol_nl
        write(*,*) "ALPHA:",alpha_elp
        write(*,*) "*********************************"
        write(*,*) ""
    end subroutine check_plasticity
    
    ! NR ALGORITHM AND DRIFT CORRECTION
    subroutine plastic_corrector (dEpsilon_ij_alpha, Sigma_ij, X_ij, sigma_yld, &
        R, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, mu, lambda, dEpsilon_ij_pl)
        implicit none
        real,                 intent(in)    :: sigma_yld            ! first yield limit
        real,                 intent(in)    :: b_lmc, Rinf_lmc      ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in)    :: C_lmc, kapa_lmc      ! Lamaitre and Chaboche parameters (kinematic hardening)
        real,                 intent(in)    :: mu, lambda           ! elastic parameters
        real,                 intent(inout) :: R                    ! starting mises radius
        real, dimension(0:2), intent(inout) :: Sigma_ij             ! starting stress state
        real, dimension(0:2), intent(inout) :: X_ij                 ! starting back stress
        real, dimension(0:2), intent(inout) :: dEpsilon_ij_alpha    ! percentage of elastic-plastic strain
        real, dimension(0:2), intent(inout) :: dEpsilon_ij_pl       ! percentage of plastic strain

        real, dimension(0:2)                :: dX_ij, gradF_0, gradF_mises
        real, dimension(0:2)                :: temp_vec
        real, dimension(0:2),     parameter :: A = (/1.0,1.0,0.5/)
        real, dimension(0:2,0:2)            :: DEL_ijhk
        integer,                  parameter :: N_incr = 10
        integer                             :: i,j,k
        real                                :: dR, dPlastMult, F_mises_0, F_mises,Ffinal

        ! COMPUTE ELASTIC STIFFNESS MATRIX
        DEL_ijhk(:,:) = 0d0
        DEL_ijhk(0:1,0:1)   = DEL_ijhk(0:1,0:1) + lambda
        DEL_ijhk(0,0)       = DEL_ijhk(0,0) + 2*mu
        DEL_ijhk(0,0)       = DEL_ijhk(0,0) + 2*mu
        DEL_ijhk(2,2)       = DEL_ijhk(2,2) + mu

        ! ELASTO-PLASTIC SUB-STEPPING
        dEpsilon_ij_alpha(0:2) = dEpsilon_ij_alpha(0:2)/N_incr

        do i = 0,N_incr-1

            ! PREDICTION
            call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises_0, gradF_0)
            write(*,*) "prediction",F_mises_0
            write(*,*) "grad",gradF_0
            ! COMPUTE PLASTIC MULTIPLIER
            call compute_plastic_modulus(dEpsilon_ij_alpha, Sigma_ij, X_ij, R, mu, lambda, sigma_yld, &
                b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult)
            
            ! HARDENING INCREMENTS
            call hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
                b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)
            dR = dPlastMult*dR
            dX_ij(0:2) = dPlastMult*dX_ij(0:2)

            ! VARIABLE UPDATE
            do k=0,2 ! plastic strains
                dEpsilon_ij_pl(k)=dEpsilon_ij_pl(k)+dPlastMult*gradF_0(k)*A(k)
            end do
            R=R+dR          ! isotropic hardening update
            X_ij(0:2)=X_ij(0:2)+dX_ij(0:2) ! back-stress update
            ! stress update
            do j = 0,2
                do k = 0,2
                    Sigma_ij(j)=Sigma_ij(j)+DEL_ijhk(k,j)*(dEpsilon_ij_alpha(k)-A(k)*dPlastMult*gradF_0(k))
                end do
            end do
            ! DRIFT CORRECTION
            call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
            write(*,*) "Ffinal (BD): ",Ffinal
            if (Ffinal .gt. tol_nl) then
                call drift_corr(.true.,Sigma_ij, X_ij, R, sigma_yld,&
                    b_lmc, Rinf_lmc, C_lmc, kapa_lmc, lambda, mu, dEpsilon_ij_pl)
                call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
            end if
            write(*,*) "Ffinal (AD): ",Ffinal
            write(*,*) ""
        end do
    end subroutine plastic_corrector
    
    ! HARDENING MODULUS AND PLASTIC MULTIPLIER INCREMENT
    subroutine compute_plastic_modulus(dEpsilon_ij, Sigma_ij, X_ij, R, mu, lambda, &
        sigma_yld, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult)
        implicit none
        real, dimension(0:2), intent(in) :: dEpsilon_ij         ! percentage of elastic-plastic strain
        real, dimension(0:2), intent(in) :: Sigma_ij            ! starting stress state
        real, dimension(0:2), intent(in) :: X_ij                ! starting back stress
        real,                 intent(in) :: R                   ! starting mises radius
        real,                 intent(in) :: b_lmc, Rinf_lmc     ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in) :: C_lmc, kapa_lmc     ! Lamaitre and Chaboche parameters (kinematic hardening)
        real,                 intent(in) :: mu, lambda          ! elastic parameters
        real,                 intent(in) :: sigma_yld           ! first yielding limit
        real,                 intent(out):: dPlastMult          ! plastic multiplier increment
        real, dimension(0:2)                :: gradF_mises
        real, dimension(0:2),     parameter :: A =(/1.0,1.0,0.5/)
        real, dimension(0:2,0:2)            :: DEL_ijhk
        real                             :: h_iso, h_kin        ! isotropic and kinematic hardening modula
        real                             :: F_mises
        real                             :: h_lmc               ! total hardening modulus
        real                             :: temp_vec
        integer                          :: j,k
        
        ! COMPUTE ELASTIC STIFFNESS MATRIX
        DEL_ijhk(:,:) = 0d0
        DEL_ijhk(0:1,0:1)   = DEL_ijhk(0:1,0:1) + lambda
        DEL_ijhk(0,0)       = DEL_ijhk(0,0) + 2*mu
        DEL_ijhk(0,0)       = DEL_ijhk(0,0) + 2*mu
        DEL_ijhk(2,2)       = DEL_ijhk(2,2) + mu

        ! COMPUTE HARDENING MODULUS
        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)

        h_iso = b_lmc*(Rinf_lmc-R)
        h_kin = 0d0
        do k=0,2
            h_kin=h_kin+kapa_lmc*X_ij(k)*gradF_mises(k)
        end do
        h_kin=C_lmc-h_kin
        h_lmc = h_iso+h_kin
        ! COMPUTE PLASTIC MULTIPLIER
        temp_vec   = 0d0
        dPlastMult = 0d0

        do k = 0,2
            do j = 0,2
                temp_vec   = temp_vec+gradF_mises(j)*DEL_ijhk(j,k)*gradF_mises(k)*A(k)
                dPlastMult = dPlastMult+gradF_mises(j)*DEL_ijhk(j,k)*dEpsilon_ij(k)
            end do
        end do
        write(*,*) "nominatore",dPlastMult
        write(*,*) "deps",dEpsilon_ij
        write(*,*) "grad",gradF_mises
        dPlastMult = dPlastMult/(h_lmc+temp_vec)
        if (dPlastMult .lt. 0) then
            write(*,*) "DPLAST NEGATIVE"
            stop
        end if

    end subroutine compute_plastic_modulus

    ! INCREMENTS OF INTRINSIC STATIC VARIABLES
    subroutine hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)
        implicit none
        real, dimension(0:2), intent(in) :: Sigma_ij        ! actual stress state
        real, dimension(0:2), intent(in) :: X_ij            ! actual back stress state
        real,                 intent(in) :: sigma_yld       ! first yielding limit
        real,                 intent(in) :: R               ! actual mises radius
        real,                 intent(in) :: b_lmc, Rinf_lmc ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in) :: C_lmc, kapa_lmc ! Lamaitre and Chaboche parameters (kinematic hardening)
        real, dimension(0:2), intent(out):: dX_ij           ! back stress increment
        real,                 intent(out):: dR              ! mises radius increment
        real, dimension(0:2), parameter  :: A = (/1.0,1.0,0.5/)
        real, dimension(0:2)             :: gradF_mises
        real                             :: F_mises
        integer                          :: k

        ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
        dR = b_lmc*(Rinf_lmc-R)
        ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
        dX_ij(0:2)=2*A(0:2)*gradF_mises(0:2)*C_lmc/3-X_ij(0:2)*kapa_lmc

    end subroutine hardening_increments
    
    ! DRIFT CORRECTION (RADIAL RETURN)
    subroutine drift_corr(drift,Sigma_ij, X_ij, R, sigma_yld, &
        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, lambda, mu, dEpsilon_ij_pl)

        logical, intent(in)                 :: drift
        real, dimension(0:2), intent(inout) :: Sigma_ij, X_ij, dEpsilon_ij_pl
        real,                 intent(inout) :: R
        real,                 intent(in)    :: lambda, mu, sigma_yld, b_lmc, Rinf_lmc, C_lmc, kapa_lmc
        real, dimension(0:2)                :: gradF_mises,gradF0,Sigma_temp,Sigma_dev_temp,Sigma_dev_ij
        real, dimension(0:2), parameter     :: A = (/1.0,1.0,0.5/)
        real, dimension(0:2,0:2)            :: DEL_ijhk
        real, dimension(0:1,0:1), parameter :: M = 1d0
        real :: Fmises,Fmises0,beta,h_kin,h_iso,h_lmc,dbeta
        integer :: k,j
        
        ! INITIAL PLASTIC CONDITION
        call mises_yld_locus(Sigma_ij, X_ij, R, sigma_yld, Fmises0, gradF0)
        if (drift) then ! B METHOD
            beta=0d0
            do 
                call mises_yld_locus(Sigma_ij,X_ij,R,sigma_yld,Fmises,gradF_mises)
                ! COMPUTE BETA FOR DRIFT CORRECTION
                dbeta=0d0
                do k=0,2
                    dbeta = dbeta+gradF0(k)*gradF_mises(k)
                end do
                beta=beta-Fmises/dbeta
                ! STRESS CORRECTION (RADIAL RETURN)
                Sigma_temp(0:2)=Sigma_ij(0:2)+beta*gradF0(0:2)*A(0:2)
                call mises_yld_locus(Sigma_temp,X_ij,R,sigma_yld,Fmises,gradF_mises)
                call tensor_components(Sigma_temp, Sigma_dev_temp)
                call tensor_components(Sigma_ij,Sigma_dev_ij)
                call tau_mises(Sigma_dev_ij-Sigma_dev_temp,err0)
                call tau_mises(Sigma_dev_ij-X_ij,err1)
                err0=err0/err1
                Sigma_ij(0:2)=Sigma_temp(0:2)
                !if (abs(Fmises) .le. tol_nl) then
                if (err0 .le. 0.0000001 .or. abs(Fmises) .le. tol_nl) then
                    exit
                end if
            end do
        else ! E METHOD
            ! COMPUTE ELASTIC STIFFNESS MATRIX
            DEL_ijhk(:,:) = 0d0
            DEL_ijhk(0:1,0:1) = DEL_ijhk(0:1,0:1) + lambda * M + id_matrix *2*mu
            DEL_ijhk(2,2) = DEL_ijhk(2,2) + mu
            
            do 
                ! COMPUTE HARDENING INCREMENTS
                h_iso = b_lmc*(Rinf_lmc-R)
                h_kin = C_lmc
                do k=0,2
                    h_kin = h_kin-kapa_lmc*X_ij(k)*gradF0(k)
                end do
                h_lmc=h_iso+h_kin
                write(*,*) "H: ", h_lmc
               ! COMPUTE BETA FOR DRIFT CORRECTION
                beta = 0d0
                do j=0,3
                    do k=0,3
                        beta=beta+gradF0(k)*DEL_ijhk(k,j)*A(j)*gradF0(j)
                    end do
                end do
                write(*,*) "beta",beta
                 
                beta=Fmises/(-h_lmc+beta)
                write(*,*) "beta",beta
                ! STRESS-STRAIN-HARDENING CORRECTION
                do k=0,3
                    do j=0,3
                        Sigma_temp(k)=Sigma_ij(k)-beta*DEL_ijhk(j,k)*A(j)*gradF0(j)
                    end do
                end do
                
                call mises_yld_locus(Sigma_temp,X_ij,R,sigma_yld,Fmises,gradF_mises)
                call tensor_components(Sigma_temp, Sigma_dev_temp)
                call tensor_components(Sigma_ij,Sigma_dev_ij)
                call tau_mises(Sigma_dev_ij-Sigma_dev_temp,err0)
                call tau_mises(Sigma_dev_ij-X_ij,err1)
                err0=err0/err1
                Sigma_ij(0:2)=Sigma_temp(0:2)

                dEpsilon_ij_pl(0:2)=dEpsilon_ij_pl(0:3)+beta*A(0:3)*gradF0(0:2)
                X_ij(0:2)=X_ij(0:2)+beta*(2*A(0:2)*gradF0(0:2)*C_lmc/3-X_ij(0:2)*kapa_lmc)
                R=R+beta*(Rinf_lmc-R)*b_lmc
                call mises_yld_locus(Sigma_ij, X_ij, R, sigma_yld, Fmises0, gradF0)

                if (abs(Fmises0) .lt. tol_nl .or. err0/err1.lt.tol_nl/100) exit
            end do
        end if
    end subroutine drift_corr

    subroutine gotoFtan(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
        real, dimension(0:2), intent(in) :: start0,dtrial0,center
        real,                 intent(in) :: radius,s0
        real,                 intent(inout):: F0,Ftrial
        real,                 intent(out):: alpha
        real, dimension(0:2)             :: start,gradF,dev,dev_temp,temp,dev0
        real                             :: Fstart,err0,err1
        integer                          :: counter
        call tensor_components(start0,dev0)
        alpha=F0/(F0-Ftrial)
        start(0:2)=start0(0:2)
        temp=start+alpha*dtrial0
        do counter=0,9 
            call tensor_components(temp, dev_temp)
            call tensor_components(start,dev)
            call tau_mises(-dev+dev_temp,err0)
            call tau_mises(dev-dev0,err1)
            write(*,*) "err0",err0
            write(*,*) "err1",err1
            err0=err0/err1
            start(0:2)=temp(0:2)
            call mises_yld_locus(start,center,radius,s0,Fstart,gradF)
            if (abs(Fstart).le.tol_nl .or. err0.lt.0.000001) then
                exit
            end if
            beta=sum(gradF(0:2)*dtrial0(0:2))
            alpha=alpha-Fstart/beta
            temp(0:2)=start(0:2)-(Fstart/beta)*dtrial0(0:2)
        end do
        write(*,*) "F:",Fstart
        write(*,*) "err0",err0
        write(*,*) "gotoF: stress:",start
        write(*,*) ""
    end subroutine gotoFtan
     
    subroutine gotoFsec(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
        real, dimension(0:2), intent(in) :: start0,dtrial0,center
        real,                 intent(in) :: radius,s0
        real,                 intent(inout):: F0,Ftrial
        real,                 intent(out):: alpha
        real, dimension(0:2)             :: start,gradF,dev,dev_temp,temp,dev0
        real                             :: alphanew,err0,err1
        integer                          :: counter
        call tensor_components(start0,dev0)
        alpha=F0/(F0-Ftrial)
        beta=0d0
        alphanew=0d0
        start(0:2)=start0(0:2)
        temp(0:2)=start(0:2)+alpha*dtrial0(0:2)
        do counter=0,9 
            call tensor_components(temp, dev_temp)
            call tensor_components(start,dev)
            call tau_mises(-dev+dev_temp,err0)
            call tau_mises(dev-dev0,err1)
            err0=err0/err1
            call mises_yld_locus(start,center,radius,s0,F0,gradF)
            call mises_yld_locus(temp,center,radius,s0,Ftrial,gradF)
            start(0:2)=temp(0:2)
            if (abs(Ftrial).le.tol_nl .or. err0.lt.0.0000001) exit
            alphanew=alpha-(Ftrial/(Ftrial-F0))*(alpha-beta)
            beta=alpha
            alpha=alphanew
            temp(0:2)=start(0:2)-(alpha-beta)*dtrial0(0:2)
            
        enddo

    end subroutine gotoFsec
end module nonlinear

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
