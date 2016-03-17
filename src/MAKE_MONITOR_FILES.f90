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

!> @brief Create monitor files. 
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nmonit number of monitored points
!> @param[in] option_out_var options for output
!> @param[in] nnode_TOT nodes of total DRM domain without duplicate
!> @param[in] tagstep label of the step of DRM Analysis

subroutine MAKE_MONITOR_FILES(nnt,nmonit,option_out_var,nnode_TOT,tagstep,&
        unit_disp,unit_vel,unit_acc,unit_stress,unit_strain,unit_omega,unit_uDRM)

    implicit none

    integer*4 ,                   intent(in)    :: nnt,nmonit,nnode_TOT,tagstep
    integer*4, dimension(nmonit), intent(inout) :: unit_disp
    integer*4, dimension(nmonit), intent(inout) :: unit_vel
    integer*4, dimension(nmonit), intent(inout) :: unit_acc
    integer*4, dimension(nmonit), intent(inout) :: unit_stress
    integer*4, dimension(nmonit), intent(inout) :: unit_strain
    integer*4, dimension(nmonit), intent(inout) :: unit_omega
    integer*4, dimension(nnode_TOT), intent(inout) :: unit_uDRM
    integer*4, dimension(6)     , intent(in)    :: option_out_var    
    character*70                                :: file_disp
    character*70                                :: file_vel
    character*70                                :: file_acc
    character*70                                :: file_stress 
    character*70                                :: file_strain 
    character*70                                :: file_omega 
    character*70                                :: file_uDRM
    integer*4                                   :: i 
    if (nmonit.ge.1) then
        if (option_out_var(1).eq.1) then  
            file_disp = 'monitorXXXXX.d'  
            do i = 1,nmonit
                unit_disp(i)=40+i
                if (i.lt.10) then
                    write(file_disp(8:11),'(a4)')'0000'
                    write(file_disp(12:12),'(i1)')i
                else if (i.le.99) then
                    write(file_disp(8:10),'(a3)')'000'
                    write(file_disp(11:12),'(i2)')i
                else if (i.le.999) then
                    write(file_disp(8:9),'(a2)')'00'
                    write(file_disp(10:12),'(i3)')i
                else if (i.le.9999) then
                    write(file_disp(8:8),'(a1)')'0'
                    write(file_disp(9:12),'(i4)')i    
                else if (i.le.99999) then
                    write(file_disp(8:12),'(i5)')'i'     
                endif
                open(unit_disp(i),file=file_disp)
            enddo
        endif

        if (option_out_var(2).eq.1) then  
            file_vel = 'monitorXXXXX.v' 
            do i = 1,nmonit
                unit_vel(i)=100000+i
                if (i.lt.10) then
                    write(file_vel(8:11),'(a4)')'0000'
                    write(file_vel(12:12),'(i1)')i
                else if (i.le.99) then
                    write(file_vel(8:10),'(a3)')'000'
                    write(file_vel(11:12),'(i2)')i
                else if (i.le.999) then
                    write(file_vel(8:9),'(a2)')'00'
                    write(file_vel(10:12),'(i3)')i
                else if (i.le.9999) then
                    write(file_vel(8:8),'(a1)')'0'
                    write(file_vel(9:12),'(i4)')i    
                else if (i.le.99999) then
                    write(file_vel(8:12),'(i5)')'i'     
                endif
                open(unit_vel(i),file=file_vel)
            enddo
        endif


        if (option_out_var(3).eq.1) then  

            file_acc = 'monitorXXXXX.a' 
            do i = 1,nmonit
                unit_acc(i)=200000+i
                if (i.lt.10) then
                    write(file_acc(8:11),'(a4)')'0000'
                    write(file_acc(12:12),'(i1)')i
                else if (i.le.99) then
                    write(file_acc(8:10),'(a3)')'000'
                    write(file_acc(11:12),'(i2)')i
                else if (i.le.999) then
                    write(file_acc(8:9),'(a2)')'00'
                    write(file_acc(10:12),'(i3)')i
                else if (i.le.9999) then
                    write(file_acc(8:8),'(a1)')'0'
                    write(file_acc(9:12),'(i4)')i    
                else if (i.le.99999) then
                    write(file_acc(8:12),'(i5)')'i'     
                endif
                open(unit_acc(i),file=file_acc)
            enddo 
        endif

        if (option_out_var(4).eq.1) then
            file_stress = 'monitorXXXXX.s'
            do i = 1,nmonit
                unit_stress(i)=300000+i
                if (i.lt.10) then
                    write(file_stress(8:11),'(a4)')'0000'
                    write(file_stress(12:12),'(i1)')i
                else if (i.le.99) then
                    write(file_stress(8:10),'(a3)')'000'
                    write(file_stress(11:12),'(i2)')i
                else if (i.le.999) then
                    write(file_stress(8:9),'(a2)')'00'
                    write(file_stress(10:12),'(i3)')i
                else if (i.le.9999) then
                    write(file_stress(8:8),'(a1)')'0'
                    write(file_stress(9:12),'(i4)')i    
                else if (i.le.99999) then
                    write(file_stress(8:12),'(i5)')'i'     
                endif         
                open(unit_stress(i),file=file_stress)
            enddo
        endif

        if (option_out_var(5).eq.1) then  
            file_strain = 'monitorXXXXX.e'
            do i = 1,nmonit
                unit_strain(i)=400000+i
                if (i.lt.10) then
                    write(file_strain(8:11),'(a4)')'0000'
                    write(file_strain(12:12),'(i1)')i
                else if (i.le.99) then
                    write(file_strain(8:10),'(a3)')'000'
                    write(file_strain(11:12),'(i2)')i
                else if (i.le.999) then
                    write(file_strain(8:9),'(a2)')'00'
                    write(file_strain(10:12),'(i3)')i
                else if (i.le.9999) then
                    write(file_strain(8:8),'(a1)')'0'
                    write(file_strain(9:12),'(i4)')i    
                else if (i.le.99999) then
                    write(file_strain(8:12),'(i5)')'i'     
                endif                   
                open(unit_strain(i),file=file_strain)
            enddo
        endif


        if (option_out_var(6).eq.1) then  
            file_omega = 'monitorXXXXX.w'
            do i = 1,nmonit
                unit_omega(i)=500000+i
                if (i.lt.10) then
                    write(file_omega(8:11),'(a4)')'0000'
                    write(file_omega(12:12),'(i1)')i
                else if (i.le.99) then
                    write(file_omega(8:10),'(a3)')'000'
                    write(file_omega(11:12),'(i2)')i
                else if (i.le.999) then
                    write(file_omega(8:9),'(a2)')'00'
                    write(file_omega(10:12),'(i3)')i
                else if (i.le.9999) then
                    write(file_omega(8:8),'(a1)')'0'
                    write(file_omega(9:12),'(i4)')i    
                else if (i.le.99999) then
                    write(file_omega(8:12),'(i5)')'i'     
                endif                   
                open(unit_omega(i),file=file_omega)
            enddo
        endif

    endif

    !---DRM---------------------------------------------------------------------
    !Open output files for DRM I step                     !DRM Scandella 25.11.2005

    if ((nnode_TOT.ne.0).and.(tagstep.eq.1)) then   !DRM Scandella 25.11.2005 
        file_uDRM = 'monDRMXXXXX.d'                  !DRM Scandella 25.11.2005 
        do i = 1,nnode_TOT                           !DRM Scandella 25.11.2005 
            unit_uDRM(i)=600000+i                    !DRM Scandella 25.11.2005
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
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
