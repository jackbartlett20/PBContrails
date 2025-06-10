!**********************************************************************************************

subroutine psr_pbe()

!**********************************************************************************************
!
! Perfectly stirred reactor for the homogeneous PBE
! Stelios Rigopoulos (06/05/2017)
! Modified 25/06/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision int_time,tin,current_time,meansize,dt

integer i,i_step,n_steps,iflag,flowflag,nin,i_write,n_write,total_writes
logical first_write
integer agg_kernel_update

!**********************************************************************************************

! Initialisation

! Initialise PBE
call pbe_read()
call pbe_grid()
call pbe_init()

! Read PSR input data
open(30,file='psr/psr.in')
do i=1,2
  read(30,*)
end do
read(30,*) int_time
read(30,*) dt
read(30,*) agg_kernel_update
read(30,*) total_writes
close(30)

! Initialise PSR integration
n_steps = int_time/dt
n_write = n_steps/total_writes
current_time= 0.D0
i_write = 1
first_write = .true.

! Temporary:
call pbe_output_general_psd(current_time, first_write)
first_write = .false.

!----------------------------------------------------------------------------------------------

! Integration

do i_step = 1,n_steps

  current_time = current_time + dt

  ! Update temperature and vapour pressure
  call pbe_set_environment(current_time)

  ! The following should be done if the kernel should be updated at each time step due to e.g. 
  ! temperature dependency
  if (agg_kernel_update==1) then
    ! Insert here the expression for updating the kernel
    ! agg_kernel_const = 
    call PBE_agg_beta(2)
  end if

  ! Update freezing (before soot so a soot particle cannot immediately turn to ice)
  ! Currently assumes whole volume is water
  call pbe_freezing(dt)

  ! Update soot activation
  !call pbe_activation()

  ! Integrate - for processes which change over time (not just those which alter ni)
  ! Pvap is reduced by growth here
  call pbe_integ(dt)


  ! Write PSD
  if ((i_write==n_write).or.(i_step==n_steps)) then

    ! General particles
    call pbe_output_general_psd(current_time, first_write)

    ! Crystals
    call pbe_output_crystal_psd(current_time, first_write)

    ! Environment variables
    call pbe_output_env(current_time, first_write)

    first_write = .false.
    i_write = 0
  end if
  i_write = i_write + 1

end do

!----------------------------------------------------------------------------------------------

! Deallocate arrays
call PBE_deallocate()

end subroutine psr_pbe

!**********************************************************************************************