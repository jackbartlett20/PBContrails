!**********************************************************************************************
!
! PBE solver subroutines
! Stelios Rigopoulos
!
!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ(ni,dt)

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

double precision, dimension(m), intent(inout) :: ni
double precision, intent(in)                  :: dt

double precision niprime(m),nitemp(m)
double precision k1(m),k2(m),k3(m),k4(m)

!----------------------------------------------------------------------------------------------

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_ydot(ni,niprime)
  ni = ni + niprime * dt

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_ydot(ni,niprime)
  nitemp = ni + 0.5D0 * niprime * dt
  call pbe_ydot(nitemp,niprime)
  ni = ni + niprime * dt

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_ydot(ni,niprime)
  k1 = niprime * dt
  nitemp = ni + 0.5D0 * k1
  call pbe_ydot(nitemp,niprime)
  k2 = niprime * dt
  nitemp = ni + 0.5D0 * k2
  call pbe_ydot(nitemp,niprime)
  k3 = niprime * dt
  nitemp = ni + k3
  call pbe_ydot(nitemp,niprime)
  k4 = niprime * dt
  ni = ni + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4

end if

end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot(ni,niprime)

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

double precision, dimension(m), intent(in)  :: ni
double precision, dimension(m), intent(out)    :: niprime

double precision dn(m)

double precision growth_source,growth_mass_source,params(1)

double precision n_sat, supersaturation_l, thermal_speed, diff_coeff

integer index

!----------------------------------------------------------------------------------------------

niprime = 0. ! d(ni)/dt
params(1) = 0.

! Nucleation
! Add homogeneous nucleation at high supersaturation? Otherwise none

! Growth

! H2O number concentration at water saturation - not sure this is correct
n_sat = avogadro_constant * Pvap / (ideal_gas_constant * temperature)

supersaturation_l = Pvap/Psat_l - 1.D0

thermal_speed = sqrt(3 * boltzmann_constant * temperature / water_molecular_mass)

diff_coeff = 2.11D-5 * (temperature/273.15D0)**(1.94D0) * (101325D0 / P_ambient)

do index = 1,m
  call growth_tvd(ni,index,n_sat,supersaturation_l,thermal_speed,diff_coeff,growth_source)
  niprime(index) = niprime(index) + growth_source
end do


!if (growth_function>0) then
!  do index = 1,m
!    call growth_tvd(ni,index,growth_source)
!    niprime(index) = niprime(index) + growth_source
!  end do
!  if (i_gm==1) then
!    ! For mass-conservative growth scheme, apply growth source term after the first interval
!    do index=2,m
!      niprime(index) = niprime(index) + (ni(index)*g_coeff1 &
!      /((g_coeff2+1.)*0.5*(v(index)**2-v(index-1)**2))) & 
!      * (v(index)**(g_coeff2+1.)-v(index-1)**(g_coeff2+1.))
!    end do
!  end if
!end if

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

end subroutine pbe_ydot

!**********************************************************************************************