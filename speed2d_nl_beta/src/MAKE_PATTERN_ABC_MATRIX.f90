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

!> @brief Computes sparsitiy pattern for abc matrix. 
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nnod number of nodes
!> @param[in] nnod_max max number of elements to be stored for each row of MM
!> @param[in] nelem_abc number of abc elements
!> @param[in] ielem_abc list of element having an abc edge
!> @param[in] nmat number of materials
!> @param[in] sd poly. degrees
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[out] MM matrix containing indeces for sparsity pattern
!> @param[out] I_abc vector for CRS format
!> @param[out] length number of nnzero elements for CRS format

       subroutine MAKE_PATTERN_ABC_MATRIX(MM,nnod,nnod_max,nelem_abc, ielem_abc,nmat, sd, &
                                            cs, cs_nnz, I_abc, length_abc)

       implicit none
       integer*4 :: nnod, nnod_max, nelem_abc,cs_nnz
       integer*4 :: im, nmat, nn, ie, j,i, n, m, it, jm 
       integer*4 :: is, in, irow1, irow2, icol1, icol2
       integer*4 :: ind_c, tt, jk, ji, length_abc,k
       integer*4, dimension(0:cs_nnz) :: cs
       integer*4, dimension(nmat) :: sd
       integer*4, dimension(nelem_abc) ::ielem_abc
  
       integer*4, dimension(2*nnod,nnod_max) :: MM
       integer*4, dimension(0:2*nnod) :: I_abc
	   
       integer*4 :: number_of_threads
       
	   number_of_threads = 1;

       call OMP_set_num_threads(number_of_threads)

       !call OMP_get_num_threads()


!$OMP PARALLEL &
!$OMP PRIVATE(k,ie, im, nn, j, i, is, in, irow1, irow2, n, m, it, jm) &
!$OMP PRIVATE(icol1, icol2, ind_c, tt)
 
!$OMP DO  	   
             
             
       do k = 1,nelem_abc
           ie = ielem_abc(k)      
           im = cs(cs(ie -1) +0);
           nn = sd(im) + 1;      
 
           do j = 1,nn
              do i = 1,nn
                 is = nn*(j -1) + i
                 in = cs(cs(ie -1) + is)

                 irow1 = in;    irow2 = in + nnod


                  do n = 1,nn
                     do m = 1,nn
                        it = nn*(n -1) + m
                        jm = cs(cs(ie -1) + it)
                        icol1 = jm;    icol2 = jm + nnod;
                               
                        call FIND_INDEX(MM, 2*nnod, nnod_max, irow1, ind_c)  
                        call FIND_INT(MM, 2*nnod, nnod_max, irow1, icol1, tt)
                               
                        if(tt .eq. 0)  MM(irow1,ind_c) = icol1

                        call FIND_INDEX(MM, 2*nnod, nnod_max, irow1, ind_c)  
                        call FIND_INT(MM, 2*nnod, nnod_max, irow1, icol2, tt)

                        if(tt .eq. 0) MM(irow1,ind_c) = icol2 

                        call FIND_INDEX(MM, 2*nnod, nnod_max, irow2, ind_c)    
                        call FIND_INT(MM, 2*nnod, nnod_max, irow2, icol1, tt)

                        if(tt .eq. 0) MM(irow2,ind_c) = icol1

                        call FIND_INDEX(MM, 2*nnod, nnod_max, irow2, ind_c)  
                        call FIND_INT(MM, 2*nnod, nnod_max, irow2, icol2, tt)

                        if(tt .eq. 0) MM(irow2,ind_c) = icol2 
                                      
                       enddo
                   enddo
                enddo
             enddo
         enddo 

!$OMP END DO
!$OMP END PARALLEL		 

       jk = 0
       do ji = 1 , 2*nnod
           call COUNT_NNZ_EL(MM, 2*nnod, nnod_max, ji, ind_c)

           if (ind_c .ne. 0) then
              I_abc(ji) = I_abc(ji-1) + ind_c
           else 
              I_abc(ji) = I_abc(ji-1)
           endif
           jk = jk + ind_c
       enddo



       length_abc = jk
       
       
       end subroutine MAKE_PATTERN_ABC_MATRIX
