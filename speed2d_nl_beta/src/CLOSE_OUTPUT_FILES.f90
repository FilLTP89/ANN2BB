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

!> @brief Close output files.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] option_out_var options for output
!> @param[in] nmonit number of monitors
!> @param[in] unit_dips/vel/acc/stress/strain/omega logical units for files

subroutine CLOSE_OUTPUT_FILES(option_out_var, nmonit, &
    unit_disp, unit_vel, unit_acc,unit_stress, unit_strain, unit_omega)

    implicit none

    integer*4 :: i,nmonit

    integer*4, dimension (6)    , intent(in)    :: option_out_var           
    integer*4, dimension(nmonit), intent(inout) :: unit_disp
    integer*4, dimension(nmonit), intent(inout) :: unit_vel
    integer*4, dimension(nmonit), intent(inout) :: unit_acc
    integer*4, dimension(nmonit), intent(inout) :: unit_stress
    integer*4, dimension(nmonit), intent(inout) :: unit_strain
    integer*4, dimension(nmonit), intent(inout) :: unit_omega

    do i = 1,nmonit  
        if (option_out_var(1).eq.1)  then                   
            close(unit_disp(i))
        endif
        if (option_out_var(2).eq.1)  then                     
            close(unit_vel(i))
        endif
        if (option_out_var(3).eq.1)  then                     
            close(unit_acc(i))
        endif
        if (option_out_var(4).eq.1)  then                     
            close(unit_stress(i))
        endif
        if (option_out_var(5).eq.1)  then                     
            close(unit_strain(i))
        endif
        if (option_out_var(6).eq.1)  then                     
            close(unit_omega(i))
        endif
    enddo
    return
end subroutine CLOSE_OUTPUT_FILES
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

