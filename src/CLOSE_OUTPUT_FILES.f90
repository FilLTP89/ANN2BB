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
                                    unit_disp, unit_vel, unit_acc,&
                                    unit_stress, unit_strain, unit_omega)

      implicit none
      
      integer*4 :: i,nmonit
      integer*4 :: unit_disp, unit_vel, unit_acc, unit_stress, unit_strain, unit_omega
      integer*4, dimension (6) :: option_out_var           
      
      do i = 1,nmonit  
         if (option_out_var(1) .eq. 1)  then                   
             unit_disp = 40 + i; close(unit_disp)
         endif
         if (option_out_var(2) .eq. 1)  then                     
             unit_vel = 100000 + i; close(unit_vel)
         endif
         if (option_out_var(3).eq.1)  then                     
             unit_acc = 200000 + i; close(unit_acc)
         endif
         if (option_out_var(4).eq.1)  then                     
             unit_stress = 300000 + i; close(unit_stress)
         endif
         if (option_out_var(5).eq.1)  then                     
             unit_strain = 400000 + i; close(unit_strain)
         endif
         if (option_out_var(6).eq.1)  then                     
             unit_omega = 500000 + i; close(unit_omega)
         endif


      enddo

      end subroutine CLOSE_OUTPUT_FILES
