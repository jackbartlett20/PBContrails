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

double precision niprime(m),nitemp(m) ! used for both droplets and crystals
double precision k1(m),k2(m),k3(m),k4(m) ! used for both droplets and crystals
double precision r_m, J, sum_Jn, delta_supersaturation_l, delta_supersaturation_i

integer index

!----------------------------------------------------------------------------------------------

! DROPLETS

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_ydot_droplet(ni_droplet,niprime,dt)
  ni_droplet = ni_droplet + niprime * dt

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_ydot_droplet(ni_droplet,niprime,dt)
  nitemp = ni_droplet + 0.5D0 * niprime * dt
  call pbe_ydot_droplet(nitemp,niprime,dt)
  ni_droplet = ni_droplet + niprime * dt

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
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

end if

! Cap at zero after growth
do index = 1,m
  if (ni_droplet(index) < 0.D0) then
    ni_droplet(index) = 0.D0
  end if
end do


! Deplete supersaturation due to droplet growth
sum_Jn = 0.D0
if (supersaturation_l>0.D0) then
  
  do index=1,m
    ! Calculate H2O flux to particles of size r_m
    r_m = ((3.D0*v_m(index))/(4.D0*pi))**(1.D0/3.D0)
    call calc_J(r_m, supersaturation_l, J)
    sum_Jn = sum_Jn + J * ni_droplet(index) * dv(index) ! Last part to change to absolute number density
  end do

  delta_supersaturation_l = (- sum_Jn / n_sat) * dt ! in brackets is ds/dt
  Pvap = Pvap + delta_supersaturation_l * Psat_l
  if (Pvap<0.D0) then
    Pvap = 0.D0 ! just in case
  end if
end if

!----------------------------------------------------------------------------------------------

! CRYSTALS

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_ydot_crystal(ni_crystal,niprime,dt)
  ni_crystal = ni_crystal + niprime * dt

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_ydot_crystal(ni_crystal,niprime,dt)
  nitemp = ni_crystal + 0.5D0 * niprime * dt
  call pbe_ydot_crystal(nitemp,niprime,dt)
  ni_crystal = ni_crystal + niprime * dt

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_ydot_crystal(ni_crystal,niprime,dt)
  k1 = niprime * dt
  nitemp = ni_crystal + 0.5D0 * k1
  call pbe_ydot_crystal(nitemp,niprime,dt)
  k2 = niprime * dt
  nitemp = ni_crystal + 0.5D0 * k2
  call pbe_ydot_crystal(nitemp,niprime,dt)
  k3 = niprime * dt
  nitemp = ni_crystal + k3
  call pbe_ydot_crystal(nitemp,niprime,dt)
  k4 = niprime * dt
  ni_crystal = ni_crystal + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4

end if

! Cap at zero after growth
do index = 1,m
  if (ni_crystal(index) < 0.D0) then
    ni_crystal(index) = 0.D0
  end if
end do


! Deplete supersaturation due to crystal growth
sum_Jn = 0.D0
if (supersaturation_i>0.D0) then
  
  do index=1,m
    ! Calculate H2O flux to particles of size r_m
    r_m = ((3.D0*v_m(index))/(4.D0*pi))**(1.D0/3.D0)
    call calc_J(r_m, supersaturation_i, J)
    sum_Jn = sum_Jn + J * ni_crystal(index) * dv(index) ! Last part to change to absolute number density
  end do

  delta_supersaturation_i = (- sum_Jn / n_sat) * dt ! in brackets is ds/dt
  Pvap = Pvap + delta_supersaturation_i * Psat_i
  if (Pvap<0.D0) then
    Pvap = 0.D0 ! just in case
  end if
end if

end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_droplet(ni,niprime,dt)

!**********************************************************************************************
!
! Calculates the right hand side of the ODEs to be integrated
!
! By Stelios Rigopoulos
! Modified by Jack Bartlett (30/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
double precision, dimension(m), intent(out) :: niprime

double precision, intent(in) :: dt

double precision dn(m)

double precision growth_source,growth_mass_source,params(1)

integer index

!----------------------------------------------------------------------------------------------

niprime = 0. ! d(ni)/dt
params(1) = 0.


! Nucleation
! Add homogeneous nucleation at high supersaturation? Otherwise none


! Droplet growth
if (supersaturation_l>0) then

  do index = 1,m
    call growth_tvd(ni, index, supersaturation_l, dt, growth_source)
    niprime(index) = niprime(index) + growth_source
  end do

  ! Else niprime(index) does not change
end if


!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ni,niprime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(niprime,ni)
end if

end subroutine pbe_ydot_droplet

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_crystal(ni,niprime,dt)

!**********************************************************************************************
!
! Calculates the right hand side of the ODEs to be integrated
!
! By Stelios Rigopoulos
! Modified by Jack Bartlett (30/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
double precision, dimension(m), intent(out) :: niprime

double precision, intent(in) :: dt

double precision dn(m)

double precision growth_source,growth_mass_source,params(1)

integer index

!----------------------------------------------------------------------------------------------

niprime = 0. ! d(ni)/dt
params(1) = 0.


! Nucleation
! Is there any homogeneous ice nucleation?


! Crystal growth
if (supersaturation_i>0) then

  do index = 1,m
    call growth_tvd(ni, index, supersaturation_i, dt, growth_source)
    niprime(index) = niprime(index) + growth_source
  end do

  ! Else niprime(index) does not change
end if


!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ni,niprime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(niprime,ni)
end if

end subroutine pbe_ydot_crystal

!**********************************************************************************************