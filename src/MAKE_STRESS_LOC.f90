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

!> @brief Compute Stress Locally
!! @author Ilario Mazzieri and Filippo Gatti
!> @date February,2016
!> @version 1.0
!> @param[in] lambda,mu elastic parameters
!> @param[in] dE strain tensor
!> @param[inout] dS stress tensor


subroutine MAKE_STRESS_LOC(lambda,mu,dE,dS)
    implicit none
    
    real*8, intent(in)                  :: lambda,mu
    real*8, intent(in), dimension(3)    :: dE
    real*8, intent(inout), dimension(4) :: dS

    dS(1)=(lambda+2.0d0*mu)*dE(1)+lambda*(dE(2)+dE(3))
    dS(2)=(lambda+2.0d0*mu)*dE(2)+lambda*(dE(1)+dE(3))
    dS(3)=lambda*(dE(1)+dE(2))
    dS(4)=mu*dE(4)
    
    return
end subroutine MAKE_STRESS_LOC
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

