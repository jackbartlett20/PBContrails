!**********************************************************************************************
!
! PBE solver subroutines
! Stelios Rigopoulos
!
!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ(dt)

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

double precision g_term, sum_gn
integer index

!----------------------------------------------------------------------------------------------

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_integ_euler(dt)

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_integ_RK2(dt)

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_integ_RK4(dt)

end if

do index=1,m
  if (kappa(index)<0.D0) then
    write(*,*) "Found kappa = ",kappa(index)," at index ",index,". Stopping."
    stop
  else if (kappa(index)>1.D0) then
    write(*,*) "Found kappa = ",kappa(index)," at index ",index,". Stopping."
    stop
  else if (rho(index)<0.D0) then
    write(*,*) "Found rho = ",rho(index)," at index ",index,". Stopping."
    stop
  else if (f_dry(index)<0.D0) then
    write(*,*) "Found f_dry = ",f_dry(index)," at index ",index,". Stopping."
    stop
  else if ((f_dry(index)>1.D0).and.(f_dry(index)<1.1D0)) then
    write(*,*) "Found f_dry = ",f_dry(index)," at index ",index,". Continuing."
    f_dry(index) = 1.D0
  else if (f_dry(index)>1.1D0) then
    write(*,*) "Found f_dry = ",f_dry(index)," at index ",index,". Stopping."
    stop
  end if
end do


! Change Pvap due to growth
sum_gn = 0.D0
do index=1,m
  call calc_growth_rate_liquid(index, .false., g_term)
  sum_gn = sum_gn + g_term*ni_droplet(index)*dv(index)
end do
Pvap = Pvap + (boltzmann_constant * temperature * (sum_gn / water_molecular_vol)) * dt


! Cap at zero after growth
!where (ni_droplet < 0.D0)
!  ni_droplet = 0.D0
!end where

!sum_VG = 0.

! Adjust water vapour due to droplet growth
!Pvap = Pvap + ((sum_VG / water_molecular_vol) * boltzmann_constant * temperature) * dt
!
!if (Pvap<0.D0) then
!  Pvap = 0.D0 ! just in case
!end if

end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ_euler(dt)

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

double precision ni_droplet_prime(m),kappa_prime(m),rho_prime(m),f_dry_prime(m)

!----------------------------------------------------------------------------------------------

call pbe_ydot_droplet(ni_droplet,ni_droplet_prime,dt)

call pbe_ydot_kappa(kappa,ni_droplet_prime,kappa_prime,dt)

call pbe_ydot_rho(rho,ni_droplet_prime,rho_prime,dt)

call pbe_ydot_f_dry(f_dry,ni_droplet_prime,f_dry_prime,dt)

ni_droplet = ni_droplet + ni_droplet_prime * dt
kappa = kappa + kappa_prime * dt
rho = rho + rho_prime * dt
f_dry = f_dry + f_dry_prime * dt


end subroutine pbe_integ_euler

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ_RK2(dt)

!**********************************************************************************************
!
! !Runge-Kutta 2nd order temporal integration
!
! Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: dt

double precision ni_droplet_prime(m),kappa_prime(m),rho_prime(m),f_dry_prime(m)
double precision temp(m)

!----------------------------------------------------------------------------------------------

call pbe_ydot_droplet(ni_droplet,ni_droplet_prime,dt)
temp = ni_droplet + 0.5D0 * ni_droplet_prime * dt
call pbe_ydot_droplet(temp,ni_droplet_prime,dt)
ni_droplet = ni_droplet + ni_droplet_prime * dt

call pbe_ydot_kappa(kappa,ni_droplet_prime,kappa_prime,dt)
temp = kappa + 0.5D0 * kappa_prime * dt
call pbe_ydot_kappa(temp,ni_droplet_prime,kappa_prime,dt)
kappa = kappa + kappa_prime * dt

call pbe_ydot_rho(rho,ni_droplet_prime,rho_prime,dt)
temp = rho + 0.5D0 * rho_prime * dt
call pbe_ydot_rho(temp,ni_droplet_prime,rho_prime,dt)
rho = rho + rho_prime * dt

call pbe_ydot_f_dry(f_dry,ni_droplet_prime,f_dry_prime,dt)
temp = f_dry + 0.5D0 * f_dry_prime * dt
call pbe_ydot_f_dry(temp,ni_droplet_prime,f_dry_prime,dt)
f_dry = f_dry + f_dry_prime * dt


end subroutine pbe_integ_RK2

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ_RK4(dt)

!**********************************************************************************************
!
! !Runge-Kutta 4th order temporal integration
!
! Jack Bartlett (18/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: dt

double precision niprime(m),nitemp(m)
double precision k1(m),k2(m),k3(m),k4(m)

!----------------------------------------------------------------------------------------------

call pbe_ydot_droplet(ni_droplet,niprime,dt)
k1 = niprime * dt
nitemp = ni_droplet + 0.5D0 * k1
call pbe_ydot_droplet(nitemp,niprime,dt)
k2 = niprime * dt
nitemp = ni_droplet + 0.5D0 * k2
call pbe_ydot_droplet(nitemp,niprime,dt)
k3 = niprime * dt
nitemp = ni_droplet + k3
call pbe_ydot_droplet(nitemp,niprime,dt)
k4 = niprime * dt
ni_droplet = ni_droplet + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4


end subroutine pbe_integ_RK4

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_droplet(ni_droplet_temp,ni_droplet_prime,dt)

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
double precision, intent(in) :: dt

double precision growth_source,growth_mass_source

integer index

!----------------------------------------------------------------------------------------------

ni_droplet_prime = 0. ! d(ni)/dt


! Nucleation
! Is there any homogeneous ice nucleation?


! Particle growth
do index=1,m
  call growth_tvd(ni_droplet_temp, index, dt, growth_source)
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

subroutine pbe_ydot_crystal(ni_crystal_temp,ni_crystal_prime,dt)

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
double precision, intent(in) :: dt

double precision growth_source,growth_mass_source

integer index

!----------------------------------------------------------------------------------------------

ni_crystal_prime = 0. ! d(ni)/dt



end subroutine pbe_ydot_crystal

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_kappa(kappa_temp,ni_droplet_prime,kappa_prime,dt)

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
double precision, intent(in) :: dt

double precision nkappa(m) ! total quantity of kappa
double precision growth_source,growth_mass_source

integer index

!----------------------------------------------------------------------------------------------

kappa_prime= 0.

nkappa = kappa_temp * ni_droplet


! Particle growth
do index=1,m
  call growth_tvd(nkappa, index, dt, growth_source)
  kappa_prime(index) = kappa_prime(index) + growth_source
end do

! No change in hygroscopicity due to condensation

!Aggregation - make include correct birth/death rates of kappa
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,nkappa,kappa_prime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(kappa_prime,nkappa)
end if

! Change in ni_droplet
kappa_prime = kappa_prime - kappa * ni_droplet_prime

! Scaling
kappa_prime = kappa_prime/ni_droplet

end subroutine pbe_ydot_kappa

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_rho(rho_temp,ni_droplet_prime,rho_prime,dt)

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
double precision, intent(in) :: dt

double precision nrho(m) ! total quantity of rho
double precision growth_source,g_term

integer index

!----------------------------------------------------------------------------------------------

rho_prime= 0.

nrho = rho_temp * ni_droplet


! Particle growth
do index=1,m
  call growth_tvd(nrho, index, dt, growth_source)
  rho_prime(index) = rho_prime(index) + growth_source
end do

! Condensation
do index=1,m
  call calc_growth_rate_liquid(index, .false., g_term)
  rho_prime(index) = rho_prime(index) + ni_droplet(index) * (water_density - rho(index)) * g_term / v_m(index)
end do

!Aggregation - make include correct birth/death rates of rho
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,nrho,rho_prime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(rho_prime,nrho)
end if

! Change in ni_droplet
rho_prime = rho_prime - rho * ni_droplet_prime

! Scaling
rho_prime = rho_prime/ni_droplet

end subroutine pbe_ydot_rho

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_f_dry(f_dry_temp,ni_droplet_prime,f_dry_prime,dt)

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
double precision, intent(in) :: dt

double precision nf_dry(m) ! total quantity of f_dry
double precision growth_source,g_term

integer index

!----------------------------------------------------------------------------------------------

f_dry_prime= 0.

nf_dry = f_dry_temp * ni_droplet


! Particle growth
do index=1,m
  call growth_tvd(nf_dry, index, dt, growth_source)
  f_dry_prime(index) = f_dry_prime(index) + growth_source
end do

! Condensation
do index=1,m
  call calc_growth_rate_liquid(index, .false., g_term)
  f_dry_prime(index) = f_dry_prime(index) - ni_droplet(index) * f_dry(index) * g_term / v_m(index)
end do

!Aggregation - make include correct birth/death rates of f_dry
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,nf_dry,f_dry_prime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(f_dry_prime,nf_dry)
end if

! Change in ni_droplet
f_dry_prime = f_dry_prime - f_dry * ni_droplet_prime

! Scaling
f_dry_prime = f_dry_prime/ni_droplet

end subroutine pbe_ydot_f_dry

!**********************************************************************************************