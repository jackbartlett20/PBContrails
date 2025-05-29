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

double precision niprime(m),nitemp(m)
double precision k1(m),k2(m),k3(m),k4(m)

!----------------------------------------------------------------------------------------------

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_ydot(ni_droplet,niprime,dt)
  ni_droplet = ni_droplet + niprime * dt

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_ydot(ni_droplet,niprime,dt)
  nitemp = ni_droplet + 0.5D0 * niprime * dt
  call pbe_ydot(nitemp,niprime,dt)
  ni_droplet = ni_droplet + niprime * dt

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_ydot(ni_droplet,niprime,dt)
  k1 = niprime * dt
  nitemp = ni_droplet + 0.5D0 * k1
  call pbe_ydot(nitemp,niprime,dt)
  k2 = niprime * dt
  nitemp = ni_droplet + 0.5D0 * k2
  call pbe_ydot(nitemp,niprime,dt)
  k3 = niprime * dt
  nitemp = ni_droplet + k3
  call pbe_ydot(nitemp,niprime,dt)
  k4 = niprime * dt
  ni_droplet = ni_droplet + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4

end if

end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot(ani,niprime,dt)

!**********************************************************************************************
!
! Calculates the right hand side of the ODEs to be integrated
!
! By Stelios Rigopoulos
! 14/01/2002
! Modified 04/05/2017
! Modified 23/06/2020
! Modified 25/06/2020
! Modified 05/07/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ani
double precision, dimension(m), intent(out) :: niprime

double precision, intent(in) :: dt

double precision dn(m)

double precision growth_source,growth_mass_source,params(1)

double precision n_sat, supersaturation_l, thermal_speed, diff_coeff, r_m, J, sum_Jn

integer index

!----------------------------------------------------------------------------------------------

niprime = 0. ! d(ni)/dt
params(1) = 0.
sum_Jn = 0.D0


! Nucleation
! Add homogeneous nucleation at high supersaturation? Otherwise none


! Droplet growth
supersaturation_l = Pvap/Psat_l - 1.D0
if (supersaturation_l>0) then

  ! H2O number concentration at water saturation - not sure this is correct
  n_sat = avogadro_constant * Pvap / (ideal_gas_constant * temperature)

  thermal_speed = sqrt(3 * boltzmann_constant * temperature / water_molecular_mass)

  diff_coeff = 2.11D-5 * (temperature/273.15D0)**(1.94D0) * (101325D0 / P_ambient)

  do index = 1,m
    call growth_tvd(ani, index, dt, n_sat, supersaturation_l, thermal_speed, diff_coeff,&
                    &growth_source)
    niprime(index) = niprime(index) + growth_source

    ! Calculate H2O flux to particles of size r_m - THIS SHOULD BE IN pbe_integ FOR OTHER SOLVERS TO WORK
    r_m = ((3.D0*v_m(index))/(4.D0*pi))**(1.D0/3.D0)
    call calc_J(r_m, n_sat, supersaturation_l, thermal_speed, diff_coeff, J)
    sum_Jn = sum_Jn + J * ani(index) * dv(index) ! Last part to change to absolute number density
  end do

  ! Reduce supersaturation
  supersaturation_l = supersaturation_l + (- sum_Jn / n_sat) * dt ! in brackets is ds/dt
  Pvap = (supersaturation_l + 1.D0) * Psat_l
  if (Pvap<0.D0) then
    Pvap = 0.D0 ! just in case
  end if

  ! Else niprime(index) does not change
end if


!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ani,niprime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(niprime,ani)
end if

end subroutine pbe_ydot

!**********************************************************************************************