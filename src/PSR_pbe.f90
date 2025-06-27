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

double precision int_time,tin,current_time,meansize,dt_req,dt

integer i,i_step,n_steps,iflag,flowflag,nin,i_write,n_write,total_writes
double precision n_steps_to_double, write_interval, next_write_time
logical first_write, stop_flag
integer variable_dt,agg_kernel_update

character(len=30) filename

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
read(30,*) dt_req
read(30,*) variable_dt
read(30,*) n_steps_to_double
read(30,*) total_writes
read(30,*) agg_kernel_update
close(30)

! Initialise PSR integration
dt = dt_req
!n_steps = int_time/dt
!n_write = n_steps/total_writes
write_interval = int_time / total_writes
next_write_time = write_interval
current_time = 0.D0
!i_write = 1
first_write = .true.

!----------------------------------------------------------------------------------------------

! Integration

!do i_step = 1,n_steps
do
  if (current_time.ge.int_time) then
    exit
  end if

  if (variable_dt==1) then
    call calc_dt(dt_req, n_steps_to_double, dt)
  end if

  write(*,*) "Using dt = ",dt

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

  ! Update volume of water in droplets

  ! Update freezing
  ! Currently assumes whole volume is water
  call pbe_freezing(dt)

  ! Integrate
  call pbe_integ(dt, stop_flag)

  if (stop_flag) then
    write(*,*) "Stopping at t = ",current_time
    stop
  end if

  ! Write outputs
  !if ((i_write==n_write).or.(i_step==n_steps)) then
  if (current_time.ge.next_write_time) then

    ! Droplets
    filename = 'output/psd_droplet.out'
    call pbe_output_psd(ni_droplet, filename, current_time, first_write)

    ! Crystals
    filename = 'output/psd_crystal.out'
    call pbe_output_psd(ni_crystal, filename, current_time, first_write)

    ! Properties
    filename = 'output/hygroscopicity.out'
    call pbe_output_property(kappa, filename, current_time, first_write)
    filename = 'output/density.out'
    call pbe_output_property(rho, filename, current_time, first_write)
    filename = 'output/dry_fraction.out'
    call pbe_output_property(f_dry, filename, current_time, first_write)

    ! Environment variables
    call pbe_output_env(current_time, first_write)

    first_write = .false.
    !i_write = 0
    next_write_time = next_write_time + write_interval
  end if
  !i_write = i_write + 1

end do

!----------------------------------------------------------------------------------------------

! Deallocate arrays
call PBE_deallocate()

end subroutine psr_pbe

!**********************************************************************************************