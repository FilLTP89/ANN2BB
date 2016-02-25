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

!> @brief Read header file.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] file_head  header file
!> @param[out] file_grid  mesh file
!> @param[out] file_mat   mate file
!> @param[out] file_out   snapshot file
!> @param[out] time_step  time step
!> @param[out] stop_time  final time
!> @param[out] option_out_var options for output
!> @param[out] nsnapshots   number of snapshots
!> @param[out] t_snapshot   time for snapshots
!> @param[out] ndt_monitor  number of time steps to complete before writing the solution
!> @param[out] depth_search_mon_lst depth for serching monitors
!> @param[out] n_lst  monitor list
!> @param[out] monfile_lst  0/1=write/read monitor list
!> @param[out] time_deg polynomial degree in time
!> @param[out] test  flag for test case
!> @param[out] dg_c  constant (-1,0,1) for DG method
!> @param[out] pen_c penalty constant for DG method

    subroutine READ_HEADER(file_head,file_grid,file_mat,file_out,&
      time_step,stop_time,option_out_var,nsnapshots,t_snapshot,ndt_monitor, &
      depth_search_mon_lst,n_lst, monfile_lst,time_deg,test,dg_c,pen_c,NLFLAG)   

      implicit none

      integer*4 :: nsnapshots

      integer*4, dimension (6) :: option_out_var  
      real*8 :: time_step,stop_time
      real*8, dimension(*) :: t_snapshot
      character*70 :: file_head,file_grid,file_mat,file_out

      character*80 :: input_line
      character*8 :: keyword
      integer*4 :: status
      integer*4 :: ileft,iright
      integer*4 :: i,j,im,is
      real*8 :: val,dg_c, pen_c

      real*8 :: ndt_monitor     
      real*8 :: depth_search_mon_lst 
      integer*4 :: n_lst
      integer*4 :: monfile_lst, time_deg, test,temp                
      logical, intent(out) :: NLFLAG

      ndt_monitor = 1; n_lst = 0;
      im = 0
      is = 0

      open(20,file=file_head)

      do
          read(20,'(A)',IOSTAT = status) input_line
          if (status.ne.0) exit
              keyword = input_line(1:8)
              ileft = 0
              iright = len(input_line)
              do i = 1,iright
                  if (input_line(i:i).eq.' ') exit
              enddo
              ileft = i

              if (keyword(1:8).eq.'GRIDFILE') then
                  read(input_line(ileft:iright),*) file_grid
                  write(*,*) 'Grid File: ',file_grid
              elseif (keyword(1:7).eq.'MATFILE') then
                  read(input_line(ileft:iright),*) file_mat
                  write(*,*) 'Material File: ',file_mat
              elseif (keyword(1:7).eq.'OUTFILE') then
                  read(input_line(ileft:iright),*) file_out
                  write(*,*) 'Output File: ',file_out
              elseif (keyword(1:8).eq.'TIMESTEP') then
                  read(input_line(ileft:iright),*) time_step
                  write(*,*) 'Time Step: ',time_step
              elseif (keyword(1:8).eq.'STOPTIME') then
                  read(input_line(ileft:iright),*) stop_time
                  write(*,*) 'Stop Time: ', stop_time
              elseif (keyword(1:8).eq.'TMONITOR') then
                  read(input_line(ileft:iright),*) ndt_monitor
                  write(*,*) 'Monitor Time Step :',ndt_monitor
              elseif (keyword(1:8).eq.'SNAPSHOT') then
                  is = is +1
                  read(input_line(ileft:iright),*) t_snapshot(is)
                  write(*,*) 'Snapshots: ',t_snapshot(is)
              elseif (keyword(1:7).eq.'OPTIOUT') then
                  read(input_line(ileft:iright),*) option_out_var(1),option_out_var(2),option_out_var(3),&
                      option_out_var(4),option_out_var(5),option_out_var(6)
                  write(*,*) 'Output Var: ',option_out_var(:)
              elseif (keyword(1:4) .eq. 'MLST') then                        
                  n_lst = 1                                
                  read(input_line(ileft:iright),*) depth_search_mon_lst, monfile_lst
              elseif (keyword(1:8) .eq. 'TESTMODE') then                        
                  test = 1                                
              elseif (keyword(1:8) .eq. 'DGMETHOD') then
                  read(input_line(ileft:iright),*) dg_c
                  write(*,*) 'DG Method: ',dg_c
              elseif (keyword(1:8) .eq. 'PENALIZC') then
                  read(input_line(ileft:iright),*) pen_c
                  write(*,*) 'Penalize: ',pen_c
              elseif (keyword(1:7).eq.'TIMEDEG') then
                  read(input_line(ileft:iright),*) time_deg
                  write(*,*) 'Time Deg: ', time_deg
              elseif (keyword(1:7).eq.'NONLINE') then
                  read(input_line(ileft:iright),*) temp
                  if (temp==1) then
                      NLFLAG=.true.
                  else
                      NLFLAG=.false.
                  endif
                  write(*,*) 'Nonlinear: ',NLFLAG
              endif
          enddo

          ! Ordered snapshots

              if (nsnapshots.gt.1) then
                  do i = 1,nsnapshots-1
                  do j = i+1, nsnapshots
                  if (t_snapshot(i).gt.t_snapshot(j)) then
                      val = t_snapshot(i)
                      t_snapshot(i) = t_snapshot(j)
                      t_snapshot(j) = val
                  endif
                  enddo
                  enddo
              else
                  write(*,'(A)')'ERROR: Not snapshot list!' 
                  write(*,*)
            endif
        close(20)

        return
    end subroutine READ_HEADER
