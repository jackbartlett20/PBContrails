!**********************************************************************************************
!
! PBE finite volume discretisation for growth
!
!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd(ni, index, type, growth_source)

!**********************************************************************************************
!
! Growth for finite volume method; used for droplets, crystals, and quantity of kappa
!
! Stelios Rigopoulos, Fabian Sewerin, Binxuan Sun
! Modified by Jack Bartlett
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in)                        :: index
integer, intent(in)                        :: type ! 0=droplet, 1=crystal
double precision, intent(out)              :: growth_source

double precision :: g_terml,g_termr,phi
double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Number density in left cell
double precision :: nll               !< Number density in left-left cell
double precision :: nc                !< Number density in this cell
double precision :: nr                !< Number density in right cell
double precision :: nrr               !< Number density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio (avoids div by zero)
double precision :: rl,rr             !< r+ at left and right surface
double precision :: r                 !< Radius at boundary

parameter(eps = 1.D1*epsilon(1.D0))

!**********************************************************************************************

if (type==0) then
  ! Growth rate at right boundary
  g_termr = g_droplet(index)
  ! Growth rate at left boundary
  g_terml = g_droplet(index-1)
else if (type==1) then
  ! Growth rate at right boundary
  g_termr = g_crystal(index)
  ! Growth rate at left boundary
  g_terml = g_crystal(index-1)
end if

!----------------------------------------------------------------------------------------------
!TVD scheme ref: S.Qamar et al 2006: A comparative study of high resolution schemes for solving
!                population balances in crystallization
!----------------------------------------------------------------------------------------------

if (g_terml>0.D0) then
  ! growth rate at left boundary is in positive direction
  
  if (index==1) then
    gnl = 0.0D0
  else if (index==m) then
    rl = (ni(m) - ni(m-1) + eps) / (ni(m-1) - ni(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (ni(m-1) + 0.5 * phi * (ni(m-1) - ni(m-2)))
  else if (index==2 ) then
    gnl = g_terml * 0.5 * (ni(1)+ni(2))
  else
    nl = ni(index-2)
    nc = ni(index-1)
    nr = ni(index)
    rl = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc + 0.5 * phi * (nc - nl))
  end if

else
  ! growth rate at left boundary is in negative direction
  
  if (index==1) then
    gnl = g_terml * (ni(1) + 0.5 * (ni(1) - ni(2)))
  else if (index==m) then
    gnl = g_terml * 0.5 * (ni(m)+ni(m-1))
  else
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rl = (nl - nc + eps) / (nc - nr + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc + 0.5 * phi * (nc - nr))
  end if
end if


if (g_termr>0.D0) then
  ! growth rate at right boundary is in positive direction

  if (index==1) then
    gnr = g_termr * 0.5 * (ni(1)+ni(2))
  else if (index==m) then
    gnr = g_termr * (ni(m) + 0.5*(ni(m) - ni(m-1)))
  else
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc + 0.5 * phi * (nc - nl))
  end if

else
  ! growth rate at right boundary is in negative direction

  if (index==1) then
    rr = (ni(1) - ni(2) + eps) / (ni(2) - ni(3) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (ni(2) + 0.5 * phi * (ni(2) - ni(3)))
  else if (index==m) then
    gnr = 0.D0
  else if (index==m-1) then
    gnr = g_termr * 0.5 * (ni(m)+ni(m-1))
  else
    nl = ni(index)
    nc = ni(index+1)
    nr = ni(index+2)
    rr = (nl - nc + eps) / (nc - nr + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc + 0.5 * phi * (nc - nr))
  end if
end if

! Calculate growth source
growth_source = (gnl - gnr) / dv(index)

end subroutine growth_tvd

!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd_rho(ni, rho_temp, index, growth_source)

!**********************************************************************************************
!
! Growth of quantity of density. Takes into account change in density due to growth.
!
! Jack Bartlett (20/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
double precision, dimension(m), intent(in) :: rho_temp
integer, intent(in)                        :: index
double precision, intent(out)              :: growth_source

double precision :: g_terml,g_termr,phi
double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Quantity density in left cell
double precision :: nll               !< Quantity density in left-left cell
double precision :: nc                !< Quantity density in this cell
double precision :: nr                !< Quantity density in right cell
double precision :: nrr               !< Quantity density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio (avoids div by zero)
double precision :: rl,rr             !< r+ at left and right surface
double precision :: r                 !< Radius at boundary
double precision :: rho_after_a,rho_after_b !< Temporary variables for holding density

parameter(eps = 1.D1*epsilon(1.D0))

!**********************************************************************************************

! Growth rate at right boundary calculation
g_termr = g_droplet(index)
! Growth rate at left boundary calculation
g_terml = g_droplet(index-1)

!----------------------------------------------------------------------------------------------
!TVD scheme ref: S.Qamar et al 2006: A comparative study of high resolution schemes for solving
!                population balances in crystallization
!----------------------------------------------------------------------------------------------

if (g_terml>0.D0) then
  ! growth rate at left boundary is in positive direction
  
  if (index==1) then
    gnl = 0.0D0
  else if (index==m) then
    call rho_growth(m-1, m, rho_temp(m-1), rho_after_a)
    call rho_growth(m-2, m, rho_temp(m-2), rho_after_b)
    rl = (ni(m)*rho_temp(m) - ni(m-1)*rho_after_a + eps) / (ni(m-1)*rho_after_a - ni(m-2)*rho_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (ni(m-1)*rho_after_a + 0.5 * phi * (ni(m-1)*rho_after_a - ni(m-2)*rho_after_b))
  else if (index==2 ) then
    call rho_growth(1, 2, rho_temp(1), rho_after_a)
    gnl = g_terml * 0.5 * (ni(1)*rho_after_a + ni(2)*rho_temp(2))
  else
    nl = ni(index-2)
    nc = ni(index-1)
    nr = ni(index)
    call rho_growth(index-1, index, rho_temp(index-1), rho_after_a)
    call rho_growth(index-2, index, rho_temp(index-2), rho_after_b)
    rl = (nr*rho_temp(index) - nc*rho_after_a + eps) / (nc*rho_after_a - nl*rho_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc*rho_after_a + 0.5 * phi * (nc*rho_after_a - nl*rho_after_b))
  end if

else
  ! growth rate at left boundary is in negative direction
  
  if (index==1) then
    gnl = g_terml * (ni(1)*rho_temp(1) + 0.5 * (ni(1)*rho_temp(1) - ni(2)*rho_temp(2)))
  else if (index==m) then
    gnl = g_terml * 0.5 * (ni(m)*rho_temp(m) + ni(m-1)*rho_temp(m-1))
  else
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rl = (nl*rho_temp(index-1) - nc*rho_temp(index) + eps) / (nc*rho_temp(index) - nr*rho_temp(index+1) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc*rho_temp(index) + 0.5 * phi * (nc*rho_temp(index) - nr*rho_temp(index+1)))
  end if
end if


if (g_termr>0.D0) then
  ! growth rate at right boundary is in positive direction

  if (index==1) then
    gnr = g_termr * 0.5 * (ni(1)*rho_temp(1) + ni(2)*rho_temp(2))
  else if (index==m) then
    gnr = g_termr * (ni(m)*rho_temp(m) + 0.5*(ni(m)*rho_temp(m) - ni(m-1)*rho_temp(m-1)))
  else
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr*rho_temp(index+1) - nc*rho_temp(index) + eps) / (nc*rho_temp(index) - nl*rho_temp(index-1) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc*rho_temp(index) + 0.5 * phi * (nc*rho_temp(index) - nl*rho_temp(index-1)))
  end if

else
  ! growth rate at right boundary is in negative direction

  if (index==1) then
    call rho_growth(2, 1, rho_temp(2), rho_after_a)
    call rho_growth(3, 1, rho_temp(3), rho_after_b)
    rr = (ni(1)*rho_temp(1) - ni(2)*rho_after_a + eps) / (ni(2)*rho_after_a - ni(3)*rho_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (ni(2)*rho_after_a + 0.5 * phi * (ni(2)*rho_after_a - ni(3)*rho_after_b))
  else if (index==m) then
    gnr = 0.D0
  else if (index==m-1) then
    call rho_growth(m, m-1, rho_temp(m), rho_after_a)
    gnr = g_termr * 0.5 * (ni(m)*rho_after_a + ni(m-1)*rho_temp(m-1))
  else
    nl = ni(index)
    nc = ni(index+1)
    nr = ni(index+2)
    call rho_growth(index+1, index, rho_temp(index+1), rho_after_a)
    call rho_growth(index+2, index, rho_temp(index+2), rho_after_b)
    rr = (nl*rho_temp(index) - nc*rho_after_a + eps) / (nc*rho_after_a - nr*rho_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc*rho_after_a + 0.5 * phi * (nc*rho_after_a - nr*rho_after_b))
  end if
end if

! Calculate growth source
growth_source = (gnl - gnr) / dv(index)

end subroutine growth_tvd_rho

!**********************************************************************************************



!**********************************************************************************************

subroutine rho_growth(i, j, rho_before, rho_after)

!**********************************************************************************************
!
! Calculates change in density due to growth (moving from interval i to interval j)
!
! By Jack Bartlett (20/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

integer, intent(in)                        :: i
integer, intent(in)                        :: j
double precision, intent(in)               :: rho_before
double precision, intent(out)              :: rho_after

!----------------------------------------------------------------------------------------------

rho_after = (rho_before*v_m(i) + water_density*(v_m(j) - v_m(i))) / v_m(j)

end subroutine rho_growth

!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd_f_dry(ni, f_dry_temp, index, growth_source)

!**********************************************************************************************
!
! Growth of quantity of dry fraction. Takes into account change in dry fraction due to growth.
!
! Jack Bartlett (20/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
double precision, dimension(m), intent(in) :: f_dry_temp
integer, intent(in)                        :: index
double precision, intent(out)              :: growth_source

double precision :: g_terml,g_termr,phi
double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Quantity density in left cell
double precision :: nll               !< Quantity density in left-left cell
double precision :: nc                !< Quantity density in this cell
double precision :: nr                !< Quantity density in right cell
double precision :: nrr               !< Quantity density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio (avoids div by zero)
double precision :: rl,rr             !< r+ at left and right surface
double precision :: r                 !< Radius at boundary
double precision :: f_dry_after_a,f_dry_after_b !< Temporary variables for holding f_dry

parameter(eps = 1.D1*epsilon(1.D0))

!**********************************************************************************************

! Growth rate at right boundary calculation
g_termr = g_droplet(index)
! Growth rate at left boundary calculation
g_terml = g_droplet(index-1)

!----------------------------------------------------------------------------------------------
!TVD scheme ref: S.Qamar et al 2006: A comparative study of high resolution schemes for solving
!                population balances in crystallization
!----------------------------------------------------------------------------------------------

if (g_terml>0.D0) then
  ! growth rate at left boundary is in positive direction
  
  if (index==1) then
    gnl = 0.0D0
  else if (index==m) then
    call f_dry_growth(m-1, m, f_dry_temp(m-1), f_dry_after_a)
    call f_dry_growth(m-2, m, f_dry_temp(m-2), f_dry_after_b)
    rl = (ni(m)*f_dry_temp(m) - ni(m-1)*f_dry_after_a + eps) / (ni(m-1)*f_dry_after_a - ni(m-2)*f_dry_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (ni(m-1)*f_dry_after_a + 0.5 * phi * (ni(m-1)*f_dry_after_a - ni(m-2)*f_dry_after_b))
  else if (index==2 ) then
    call f_dry_growth(1, 2, f_dry_temp(1), f_dry_after_a)
    gnl = g_terml * 0.5 * (ni(1)*f_dry_after_a + ni(2)*f_dry_temp(2))
  else
    nl = ni(index-2)
    nc = ni(index-1)
    nr = ni(index)
    call f_dry_growth(index-1, index, f_dry_temp(index-1), f_dry_after_a)
    call f_dry_growth(index-2, index, f_dry_temp(index-2), f_dry_after_b)
    rl = (nr*f_dry_temp(index) - nc*f_dry_after_a + eps) / (nc*f_dry_after_a - nl*f_dry_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc*f_dry_after_a + 0.5 * phi * (nc*f_dry_after_a - nl*f_dry_after_b))
  end if

else
  ! growth rate at left boundary is in negative direction
  
  if (index==1) then
    gnl = g_terml * (ni(1)*f_dry_temp(1) + 0.5 * (ni(1)*f_dry_temp(1) - ni(2)*f_dry_temp(2)))
  else if (index==m) then
    gnl = g_terml * 0.5 * (ni(m)*f_dry_temp(m) + ni(m-1)*f_dry_temp(m-1))
  else
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rl = (nl*f_dry_temp(index-1) - nc*f_dry_temp(index) + eps) / (nc*f_dry_temp(index) - nr*f_dry_temp(index+1) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc*f_dry_temp(index) + 0.5 * phi * (nc*f_dry_temp(index) - nr*f_dry_temp(index+1)))
  end if
end if


if (g_termr>0.D0) then
  ! growth rate at right boundary is in positive direction

  if (index==1) then
    gnr = g_termr * 0.5 * (ni(1)*f_dry_temp(1) + ni(2)*f_dry_temp(2))
  else if (index==m) then
    gnr = g_termr * (ni(m)*f_dry_temp(m) + 0.5*(ni(m)*f_dry_temp(m) - ni(m-1)*f_dry_temp(m-1)))
  else
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr*f_dry_temp(index+1) - nc*f_dry_temp(index) + eps) / (nc*f_dry_temp(index) - nl*f_dry_temp(index-1) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc*f_dry_temp(index) + 0.5 * phi * (nc*f_dry_temp(index) - nl*f_dry_temp(index-1)))
  end if

else
  ! growth rate at right boundary is in negative direction

  if (index==1) then
    call f_dry_growth(2, 1, f_dry_temp(2), f_dry_after_a)
    call f_dry_growth(3, 1, f_dry_temp(3), f_dry_after_b)
    rr = (ni(1)*f_dry_temp(1) - ni(2)*f_dry_after_a + eps) / (ni(2)*f_dry_after_a - ni(3)*f_dry_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (ni(2)*f_dry_after_a + 0.5 * phi * (ni(2)*f_dry_after_a - ni(3)*f_dry_after_b))
  else if (index==m) then
    gnr = 0.D0
  else if (index==m-1) then
    call f_dry_growth(m, m-1, f_dry_temp(m), f_dry_after_a)
    gnr = g_termr * 0.5 * (ni(m)*f_dry_after_a + ni(m-1)*f_dry_temp(m-1))
  else
    nl = ni(index)
    nc = ni(index+1)
    nr = ni(index+2)
    call f_dry_growth(index+1, index, f_dry_temp(index+1), f_dry_after_a)
    call f_dry_growth(index+2, index, f_dry_temp(index+2), f_dry_after_b)
    rr = (nl*f_dry_temp(index) - nc*f_dry_after_a + eps) / (nc*f_dry_after_a - nr*f_dry_after_b + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc*f_dry_after_a + 0.5 * phi * (nc*f_dry_after_a - nr*f_dry_after_b))
  end if
end if

! Calculate growth source
growth_source = (gnl - gnr) / dv(index)

end subroutine growth_tvd_f_dry

!**********************************************************************************************



!**********************************************************************************************

subroutine f_dry_growth(i, j, f_dry_before, f_dry_after)

!**********************************************************************************************
!
! Calculates change in dry fraction due to growth (moving from interval i to interval j)
!
! By Jack Bartlett (20/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

integer, intent(in)                        :: i
integer, intent(in)                        :: j
double precision, intent(in)               :: f_dry_before
double precision, intent(out)              :: f_dry_after

!----------------------------------------------------------------------------------------------

f_dry_after = f_dry_before*v_m(i) / v_m(j)

end subroutine f_dry_growth

!**********************************************************************************************



!**********************************************************************************************

subroutine calc_growth_rate_liquid(r, kappa_i, rho_i, f_dry_i, g_term)

!**********************************************************************************************
!
! Calculates droplet growth rate with given parameters
!
! By Jack Bartlett (10/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in)               :: r
double precision, intent(in)               :: kappa_i
double precision, intent(in)               :: rho_i
double precision, intent(in)               :: f_dry_i
double precision, intent(out)              :: g_term

double precision accom_coeff
double precision raoult_term, kelvin_term, S_droplet
double precision diffusivity_mod, k_air_mod
double precision F_d, F_k

!----------------------------------------------------------------------------------------------

accom_coeff = 1.D0

if (kappa_i.eq.0.D0) then
  raoult_term = 1.D0
else
  raoult_term = (1.D0 - f_dry_i)/(1.D0 - (1.D0 - kappa_i)*f_dry_i)
end if

kelvin_term = exp((2.D0 * sigma_water * water_molar_mass)/(ideal_gas_constant*temperature*water_density*r))

S_droplet = raoult_term * kelvin_term

diffusivity_mod = diffusivity / ( r/(r + 0.7*mfp_air) + &
                                & diffusivity/(r*accom_coeff) * &
                                & sqrt(2*pi*water_molar_mass/(ideal_gas_constant*temperature)) )

k_air_mod = k_air / ( r/(r + 0.7*mfp_air) + &
                    & k_air/(r*accom_coeff*air_density*cp_air) * &
                    & sqrt(2*pi*air_molar_mass/(ideal_gas_constant*temperature)) )

F_d = (water_density * ideal_gas_constant * temperature) / (Psat_l * diffusivity_mod * water_molar_mass)

F_k = (l_v * water_density)/(k_air_mod * temperature) * (l_v*water_molar_mass/(ideal_gas_constant*temperature) - 1)

g_term = (4*pi*r) * 1/(F_d + F_k) * (saturation_l - S_droplet)

end subroutine calc_growth_rate_liquid

!**********************************************************************************************



!**********************************************************************************************

subroutine calc_growth_rate_crystal(index, g_term)

!**********************************************************************************************
!
! Calculates crystal growth rate at indexed boundary
!
! By Jack Bartlett (27/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

integer, intent(in) :: index
double precision, intent(out) :: g_term

double precision r, J

!----------------------------------------------------------------------------------------------

r = ((3.D0*v(index))/(4.D0*pi))**(1.D0/3.D0) ! Find radius of indexed boundary

call calc_J(r, J)

g_term = water_molecular_vol * J

end subroutine calc_growth_rate_crystal

!**********************************************************************************************



!**********************************************************************************************

subroutine calc_J(r, J)

!**********************************************************************************************
!
! Calculates J (flux of water molecules to crystals of size r) from Kärcher 2015
!
! By Jack Bartlett (28/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: r
double precision, intent(out) :: J

double precision accom_coeff, correction_factor

!----------------------------------------------------------------------------------------------

accom_coeff = 1.D0

correction_factor = 1 + accom_coeff * vapour_thermal_speed * r / (4.D0 * diffusivity)

J = (pi*r**2 * accom_coeff * vapour_thermal_speed * (saturation_i-1.D0) * n_sat) / correction_factor

end subroutine calc_J

!**********************************************************************************************



!**********************************************************************************************

subroutine calc_f_dry_equ(r, kappa_i, f_dry_eq)

!**********************************************************************************************
!
! Calculates dry fraction in equilibrium
!
! By Jack Bartlett (02/07/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in)               :: r
double precision, intent(in)               :: kappa_i
double precision, intent(out)              :: f_dry_eq

double precision kelvin_term

!----------------------------------------------------------------------------------------------

kelvin_term = exp((2.D0 * sigma_water * water_molar_mass)/(ideal_gas_constant*temperature*water_density*r))

f_dry_eq = (1 - kelvin_term/saturation_l) / (1 - kappa_i - kelvin_term/saturation_l)

end subroutine calc_f_dry_equ

!**********************************************************************************************



!**********************************************************************************************

subroutine courant_check(dt, courant_success, dt_sugg)

!**********************************************************************************************
!
! Checks if the Courant-Friedrichs-Lewy (CFL) condition for stability (C <= desired number)
! holds.
!
! By Jack Bartlett (24/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in)  :: dt
logical, intent(out)          :: courant_success
double precision, intent(out) :: dt_sugg

double precision courant
integer i

!----------------------------------------------------------------------------------------------

courant_success = .true.
dt_sugg = dt

if (maxval(f_dry) > 0.999D0) then
  do i=1,m
    if (ni_droplet(i).eq.0.D0) then ! Don't bother if there are no droplets there
      cycle
    else
      courant = max(abs(g_droplet(i-1)), abs(g_droplet(i))) * dt / dv(i)
      if (courant>courant_condition_tight) then
        !write(*,*) "Courant number of ",courant," detected at index ",i,"."
        courant_success = .false.
        dt_sugg = dt / (1.01D0 * courant / courant_condition_tight)
        exit
      end if
    end if
  end do
else
  do i=1,m
    if (ni_droplet(i).eq.0.D0) then ! Don't bother if there are no droplets there
      cycle
    else
      courant = max(abs(g_droplet(i-1)), abs(g_droplet(i))) * dt / dv(i)
      if (courant>courant_condition) then
        !write(*,*) "Courant number of ",courant," detected at index ",i,"."
        courant_success = .false.
        dt_sugg = dt / (1.01D0 * courant / courant_condition)
        exit
      end if
    end if
  end do
end if

end subroutine courant_check

!**********************************************************************************************