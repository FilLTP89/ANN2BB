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

!> @brief Computes max time step allowed according to the CFL condition
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] time_step time step for time integration 
!> @param[in] nnode number of spectral nodes
!> @param[in] nm number of materials
!> @param[in] tm label for materials
!> @param[in] pm material properties
!> @param[in] sdeg polynomial degrees
!> @param[in] xx,yy coordinates of spectral nodes
!> @param[in] nquad number of elements
!> @param[in] con_quad connectivity for elements
!> @param[in] fmax max frequency of the signal (dummy)
!> @param[out] time_step_cfl max time step allowed for avoiding instability

        subroutine DELTAT_MAX(time_step,nnode,nm,tm,pm,sdeg,&
                            xx,yy,nquad,con_quad,time_step_cfl,fmax)   
      
        implicit none
         
        integer*4 :: nm,im
        integer*4, dimension(nm) :: tm,sdeg
        real*8, dimension(nm,4) :: pm
        integer*4 :: nnode,nquad
        real*8, dimension(nnode) :: xx,yy
        integer*4, dimension(nquad,5) :: con_quad
        integer*4, dimension(12):: vet_node
        real*8 :: length_min,length,vs_length_min,vs_length,length_vp_min,length_vp
        real*8 :: time_step,time_step_cfl,percent_deltat
        real*8 :: num_of_points,fmax,vs_npoints,vp_deltat_cfl
        real*8 :: rho,lambda,mu,vs,vp
        
        integer*4 :: iquad,conta
        integer*4 :: mat_deltat_cfl,sdeg_deltat_cfl
        integer*4 :: mat__npoints,sdeg_npoints
        
        integer*4 :: ie,i,j,mcode,smcode,nn

        real*8, dimension(:), allocatable :: ct,ww
        real*8, dimension(:,:), allocatable :: dd
	
	character*3 :: deltat_fixed
	      
        vet_node= (/2,2,2,3,3,4,3,4,5,4,5,5/)
        conta=1
        iquad=1

       do j = 1,nm
               if (tm(j).eq.con_quad(iquad,1)) im=j
        enddo
        
	smcode=sdeg(im)        
        rho = pm(im,1)   
        lambda = pm(im,2)
        mu = pm(im,3)
        vs = dsqrt(mu/rho)
        vp = dsqrt((lambda+2*mu)/rho)


       !The program assumes the first value of length_min the distance between
       !       nodes 1 2 belonging to the first quad

       length=dsqrt((xx(con_quad(iquad,vet_node(conta))) &
               -xx(con_quad(iquad,vet_node(conta+6))))**2 &
                +(yy(con_quad(iquad,vet_node(conta))) &
                -yy(con_quad(iquad,vet_node(conta+6))))**2)
                        
        length_vp_min = length/vp
	vp_deltat_cfl = vp
	sdeg_deltat_cfl = smcode

        vs_length_min = vs/length
	vs_npoints = vs
	sdeg_npoints = smcode


        do iquad=1,nquad
                mcode=con_quad(iquad,1)
 
                do j = 1,nm
                    if (tm(j).eq.con_quad(iquad,1)) im=j
                enddo
                
                smcode = sdeg(im)
                rho = pm(im,1)
                lambda = pm(im,2)
                mu = pm(im,3)
                vs = dsqrt(mu/rho)
                vp = dsqrt((lambda+2*mu)/rho)
 
                !Min Length check

                do conta=1,6

                        length=dsqrt((xx(con_quad(iquad,vet_node(conta))) &
                                -xx(con_quad(iquad,vet_node(conta+6))))**2 &
                                +(yy(con_quad(iquad,vet_node(conta))) &
                                -yy(con_quad(iquad,vet_node(conta+6))))**2)
                        
                        
                        length_vp = length/vp

                        vs_length = vs/length

                        !Check on the minimum length/vp of the element                  

                        if (length_vp_min.gt.length_vp) then
                                length_vp_min = length_vp
                                vp_deltat_cfl = vp
                                sdeg_deltat_cfl = smcode
                        endif

                        !Check on the number of points per element
                        !vs_length_min = Minimun ratio between the maximum element length (length) and
                        !                minimum shear velocity (vs)
                        !mat_length_npoints = label of the material associated with the vs_lenght_min
                        !sdeg_mat_length_npoints = spectral degree of the material associated with the vs_lenght_min
                        
                        if (vs_length_min.gt.vs_length) then
                                vs_length_min = vs_length
                                vs_npoints = vs
                                sdeg_npoints = smcode
                        endif
                                
               enddo

        enddo

        
        
        nn=sdeg_deltat_cfl+1

        allocate(ct(nn),ww(nn),dd(nn,nn))
        call lgl(nn,ct,ww,dd)   
        
        time_step_cfl=length_vp_min*0.5d0*(ct(2)-ct(1))
        write(*,'(A)')' '
        write(*,'(A)')'Stability for time advancing algorithm '
        write(*,'(A,E12.4)')'Min. element length :',length_vp_min*vp_deltat_cfl
        write(*,'(A,E12.4)')'Min. pressure wave velocity :',vp_deltat_cfl
        write(*,'(A)') '*******************************************************'
        if (time_step.le.time_step_cfl) then
                percent_deltat=time_step/time_step_cfl*100
                if (percent_deltat.le.1) then
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This time step is excedingly lower the deltat CFL'
			if (deltat_fixed.eq.'not') then	 
                        	write(*,'(A)')'deltat chosen will be substituted with 1% of deltat CFL'
                        	time_step=time_step_cfl*0.01
                        	write(*,'(A,E12.4)')'deltat chosen :',time_step
			endif
                elseif (percent_deltat.le.25) then
                        write(*,'(A,E12.4)')'OK!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                else
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This could be not enough for the correct working of ABC'
			if (deltat_fixed.eq.'not') then
                        	write(*,'(A)')'deltat chosen will be substituted with 25% of deltat CFL'
                        	time_step=time_step_cfl*0.25
                        	write(*,'(A,E12.4)')'deltat chosen :',time_step
			endif
                endif
        elseif (time_step.gt.time_step_cfl) then
                write(*,'(2X,A,E12.4)')'ERROR!!! deltat CFL = ',&
                time_step_cfl,' < deltat = ',time_step
                write(*,'(A)')'The time advancing scheme will be unstable!'
		if (deltat_fixed.eq.'not') then
                	write(*,'(A)')'deltat chosen will be substituted with 25% of deltat CFL'
                	time_step=time_step_cfl*0.25
                	write(*,'(A,E12.4)')'deltat chosen :',time_step
		endif
        endif
        write(*,'(A)')'*******************************************************'
        write(*,'(A)')' '

        !Writing of the results for the maximum number of points per wave length
        
        nn=sdeg_npoints+1
        deallocate(ct,ww,dd)
	allocate(ct(nn),ww(nn),dd(nn,nn))
        call lgl(nn,ct,ww,dd)

        num_of_points = vs_length_min*(ct(1)-ct(nn))/(ct(int(nn/2))-ct(int(nn/2)+1))/fmax

        write(*,'(A)')'Number of points per wave length '
        write(*,'(A,E12.4)')'Max. element length :',1/vs_length_min*vs_npoints
        write(*,'(A,E12.4)')'Max. shear wave velocity :',vs_npoints
        write(*,'(A,E12.4)')'Number of points per wave length :',num_of_points
        write(*,'(A)')' '


        return
               
        end subroutine DELTAT_MAX

