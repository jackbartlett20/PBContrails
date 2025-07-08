!**********************************************************************************************
!
! PBE solver subroutines
! Stelios Rigopoulos
!
!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ(dt, integ_success)

!**********************************************************************************************
!
! Temporal integration
!
! Stelios Rigopoulos (02/06/2019)
!
! Modified 14/06/06
! Modified 19/12/2017
! Modified 13/05/2018
! Modified 02/06/2019
! Modified 25/06/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: dt
logical, intent(out) :: integ_success

double precision ni_droplet_prime(m), ni_crystal_prime(m)
double precision kappa_prime(m), rho_prime(m), f_dry_prime(m)
double precision g_term, sum_gn
integer index

!----------------------------------------------------------------------------------------------

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_integ_euler(dt, ni_droplet_prime, ni_crystal_prime, kappa_prime, rho_prime, f_dry_prime)

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_integ_RK2(dt, ni_droplet_prime, ni_crystal_prime, kappa_prime, rho_prime, f_dry_prime)

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_integ_RK4(dt, ni_droplet_prime, ni_crystal_prime, kappa_prime, rho_prime, f_dry_prime)

end if


call check_valid_properties(dt,kappa_prime,rho_prime,f_dry_prime,integ_success)


if (integ_success) then

  ! Update Pvap due to growth
  sum_gn = 0.D0
  do index=1,m
    ! Droplets
    g_term = 0.5D0 * (g_droplet(index-1) + g_droplet(index))
    sum_gn = sum_gn + g_term*ni_droplet(index)*dv(index)
    ! Crystals
    g_term = 0.5D0 * (g_crystal(index-1) + g_crystal(index))
    sum_gn = sum_gn + g_term*ni_crystal(index)*dv(index)
  end do
  Pvap = Pvap + (boltzmann_constant * temperature * (-sum_gn / water_molecular_vol)) * dt
  
  ! Update arrays
  ni_droplet = ni_droplet + ni_droplet_prime * dt
  ni_crystal = ni_crystal + ni_crystal_prime * dt
  kappa = kappa + kappa_prime * dt
  !rho = rho + rho_prime * dt
  f_dry = f_dry + f_dry_prime * dt
  
  ! Adjust f_dry within tolerance
  do index=1,m
    if ((f_dry(index)>1.D0).and.(f_dry(index)<1.D0+f_dry_tolerance)) then
      f_dry(index) = 1.D0
    else if ((f_dry(index)<0.D0).and.(f_dry(index)>-f_dry_tolerance)) then
      f_dry(index) = 0.D0
    end if
  end do

end if


end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ_euler(dt, ni_droplet_prime, ni_crystal_prime, kappa_prime, rho_prime, f_dry_prime)

!**********************************************************************************************
!
! Euler explicit temporal integration
!
! Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: dt
double precision, dimension(m), intent(out) :: ni_droplet_prime
double precision, dimension(m), intent(out) :: ni_crystal_prime
double precision, dimension(m), intent(out) :: kappa_prime
double precision, dimension(m), intent(out) :: rho_prime
double precision, dimension(m), intent(out) :: f_dry_prime

integer index

!----------------------------------------------------------------------------------------------

call pbe_ydot_droplet(ni_droplet,ni_droplet_prime)

call pbe_ydot_crystal(ni_crystal,ni_crystal_prime)

call pbe_ydot_kappa(kappa,ni_droplet_prime,kappa_prime)

!call pbe_ydot_rho(rho,ni_droplet_prime,rho_prime)

call pbe_ydot_f_dry(f_dry,ni_droplet_prime,f_dry_prime)


end subroutine pbe_integ_euler

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ_RK2(dt, ni_droplet_prime, ni_crystal_prime, kappa_prime, rho_prime, f_dry_prime)

!**********************************************************************************************
!
! Runge-Kutta 2nd order temporal integration
!
! Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: dt
double precision, dimension(m), intent(out) :: ni_droplet_prime
double precision, dimension(m), intent(out) :: ni_crystal_prime
double precision, dimension(m), intent(out) :: kappa_prime
double precision, dimension(m), intent(out) :: rho_prime
double precision, dimension(m), intent(out) :: f_dry_prime

double precision temp(m)

!----------------------------------------------------------------------------------------------

call pbe_ydot_droplet(ni_droplet,ni_droplet_prime)
temp = ni_droplet + 0.5D0 * ni_droplet_prime * dt
call pbe_ydot_droplet(temp,ni_droplet_prime)

call pbe_ydot_crystal(ni_crystal,ni_crystal_prime)
temp = ni_droplet + 0.5D0 * ni_crystal_prime * dt
call pbe_ydot_crystal(temp,ni_crystal_prime)

call pbe_ydot_kappa(kappa,ni_droplet_prime,kappa_prime)
temp = kappa + 0.5D0 * kappa_prime * dt
call pbe_ydot_kappa(temp,ni_droplet_prime,kappa_prime)

call pbe_ydot_rho(rho,ni_droplet_prime,rho_prime)
temp = rho + 0.5D0 * rho_prime * dt
call pbe_ydot_rho(temp,ni_droplet_prime,rho_prime)

call pbe_ydot_f_dry(f_dry,ni_droplet_prime,f_dry_prime)
temp = f_dry + 0.5D0 * f_dry_prime * dt
call pbe_ydot_f_dry(temp,ni_droplet_prime,f_dry_prime)


end subroutine pbe_integ_RK2

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ_RK4(dt, ni_droplet_prime, ni_crystal_prime, kappa_prime, rho_prime, f_dry_prime)

!**********************************************************************************************
!
! Runge-Kutta 4th order temporal integration
!
! Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: dt
double precision, dimension(m), intent(out) :: ni_droplet_prime
double precision, dimension(m), intent(out) :: ni_crystal_prime
double precision, dimension(m), intent(out) :: kappa_prime
double precision, dimension(m), intent(out) :: rho_prime
double precision, dimension(m), intent(out) :: f_dry_prime

double precision niprime(m),nitemp(m)
double precision k1(m),k2(m),k3(m),k4(m)

!----------------------------------------------------------------------------------------------

call pbe_ydot_droplet(ni_droplet,niprime)
k1 = niprime * dt
nitemp = ni_droplet + 0.5D0 * k1
call pbe_ydot_droplet(nitemp,niprime)
k2 = niprime * dt
nitemp = ni_droplet + 0.5D0 * k2
call pbe_ydot_droplet(nitemp,niprime)
k3 = niprime * dt
nitemp = ni_droplet + k3
call pbe_ydot_droplet(nitemp,niprime)
k4 = niprime * dt
ni_droplet = ni_droplet + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4


end subroutine pbe_integ_RK4

!**********************************************************************************************



!**********************************************************************************************

subroutine check_valid_properties(dt,kappa_prime,rho_prime,f_dry_prime,integ_success)

!**********************************************************************************************
!
! Checks if property values would be valid if updated
!
! Jack Bartlett (27/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in)               :: dt
double precision, dimension(m), intent(in) :: kappa_prime
double precision, dimension(m), intent(in) :: rho_prime
double precision, dimension(m), intent(in) :: f_dry_prime
logical, intent(out)                       :: integ_success

integer i

!----------------------------------------------------------------------------------------------

integ_success = .true.

do i=1,m
  if ((kappa(i)+kappa_prime(i)*dt)<0.D0) then
    write(*,*) "ERROR: Found kappa = ",(kappa(i)+kappa_prime(i)*dt)," at index ",i
    integ_success = .false.
    exit
  else if ((kappa(i)+kappa_prime(i)*dt)>1.D0) then
    write(*,*) "ERROR: Found kappa = ",(kappa(i)+kappa_prime(i)*dt)," at index ",i
    integ_success = .false.
    exit
  else if ((rho(i)+rho_prime(i)*dt)<0.D0) then
    write(*,*) "ERROR: Found rho = ",(rho(i)+rho_prime(i)*dt)," at index ",i
    integ_success = .false.
    exit
  else if ((f_dry(i)+f_dry_prime(i)*dt)<-f_dry_tolerance) then
    write(*,*) "ERROR: Found f_dry = ",(f_dry(i)+f_dry_prime(i)*dt)," at index ",i
    integ_success = .false.
    exit
  else if ((f_dry(i)+f_dry_prime(i)*dt)>1.D0+f_dry_tolerance) then
    write(*,*) "ERROR: Found f_dry = ",(f_dry(i)+f_dry_prime(i)*dt)," at index ",i
    integ_success = .false.
    exit
  end if
end do

end subroutine check_valid_properties

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_droplet(ni_droplet_temp,ni_droplet_prime)

!**********************************************************************************************
!
! Calculates the right hand side of the droplet ODE to be integrated
!
! By Stelios Rigopoulos
! Modified by Jack Bartlett (30/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in)  :: ni_droplet_temp
double precision, dimension(m), intent(out) :: ni_droplet_prime

double precision growth_source

integer index

!----------------------------------------------------------------------------------------------

ni_droplet_prime = 0.D0


! Nucleation
! None unless homogeneous nucleation at high RH is added


! Particle growth
do index=1,m
  call growth_tvd(ni_droplet_temp, index, 0, growth_source)
  ni_droplet_prime(index) = ni_droplet_prime(index) + growth_source
end do


!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ni_droplet_temp,ni_droplet_prime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(ni_droplet_prime,ni_droplet_temp)
end if

end subroutine pbe_ydot_droplet

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_crystal(ni_crystal_temp,ni_crystal_prime)

!**********************************************************************************************
!
! Calculates the right hand side of the crystal ODE to be integrated
!
! By Stelios Rigopoulos
! Modified by Jack Bartlett (30/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in)  :: ni_crystal_temp
double precision, dimension(m), intent(out) :: ni_crystal_prime

double precision growth_source

integer index

!----------------------------------------------------------------------------------------------

ni_crystal_prime = 0.D0

! Nucleation
! None unless homogeneous nucleation at high RH is added


! Particle growth
do index=1,m
  call growth_tvd(ni_crystal_temp, index, 1, growth_source)
  ni_crystal_prime(index) = ni_crystal_prime(index) + growth_source
end do


!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ni_crystal_temp,ni_crystal_prime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(ni_crystal_temp,ni_crystal_prime)
end if

end subroutine pbe_ydot_crystal

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_kappa(kappa_temp,ni_droplet_prime,kappa_prime)

!**********************************************************************************************
!
! Calculates the right hand side of the hygroscopicity ODE to be integrated
!
! By Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in)  :: kappa_temp
double precision, dimension(m), intent(in)  :: ni_droplet_prime
double precision, dimension(m), intent(out) :: kappa_prime

double precision nkappa(m) ! total quantity of kappa
double precision growth_source

integer index

!----------------------------------------------------------------------------------------------

kappa_prime = 0.D0

nkappa = kappa_temp*ni_droplet

! Particle growth
do index=1,m
  call growth_tvd(nkappa, index, 0, growth_source)
  kappa_prime(index) = kappa_prime(index) + growth_source
end do

!Aggregation - make include correct birth/death rates of kappa
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  !call pbe_agg_cfv(dv,v_m,nkappa,kappa_prime)
end if

!Fragmentation
if (break_const>0.) then
  !call pbe_breakage_cfv(kappa_prime,nkappa)
end if

! Change in ni_droplet
kappa_prime = kappa_prime - kappa * ni_droplet_prime

! Scaling
kappa_prime = kappa_prime/ni_droplet

end subroutine pbe_ydot_kappa

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_rho(rho_temp,ni_droplet_prime,rho_prime)

!**********************************************************************************************
!
! Calculates the right hand side of the mass density ODE to be integrated
!
! By Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in)  :: rho_temp
double precision, dimension(m), intent(in)  :: ni_droplet_prime
double precision, dimension(m), intent(out) :: rho_prime

double precision growth_source

integer index

!----------------------------------------------------------------------------------------------

rho_prime = 0.D0

! Particle growth
do index=1,m
  call growth_tvd_rho(ni_droplet, rho_temp, index, growth_source)
  rho_prime(index) = rho_prime(index) + growth_source
end do

!Aggregation - make include correct birth/death rates of rho

!Fragmentation

! Change in ni_droplet
rho_prime = rho_prime - rho * ni_droplet_prime

! Scaling
rho_prime = rho_prime/ni_droplet

end subroutine pbe_ydot_rho

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_f_dry(f_dry_temp,ni_droplet_prime,f_dry_prime)

!**********************************************************************************************
!
! Calculates the right hand side of the dry volume fraction ODE to be integrated
!
! By Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in)  :: f_dry_temp
double precision, dimension(m), intent(in)  :: ni_droplet_prime
double precision, dimension(m), intent(out) :: f_dry_prime

double precision growth_source

integer index

!----------------------------------------------------------------------------------------------

f_dry_prime = 0.D0

! Particle growth
!do index=1,m
!  call growth_tvd_f_dry(ni_droplet, f_dry_temp, index, growth_source)
!  f_dry_prime(index) = f_dry_prime(index) + growth_source
!end do

! Add contribution to growth to interval 1 from outside spectrum
!if (g_droplet(0) > 0.D0) then
!  f_dry_prime(1) = f_dry_prime(1) + ni_droplet(1)*(-f_dry(1) * 0.5D0*g_droplet(0) / v_m(1))
!end if

do index=1,m
  f_dry_prime(index) = -f_dry(index) * 0.5D0*(g_droplet(index-1)+g_droplet(index)) / v_m(index)
end do

!Aggregation - make include correct birth/death rates of f_dry

!Fragmentation

! Change in ni_droplet
!f_dry_prime = f_dry_prime - f_dry * ni_droplet_prime

! Scaling
!f_dry_prime = f_dry_prime/ni_droplet


end subroutine pbe_ydot_f_dry

!**********************************************************************************************