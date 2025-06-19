!**********************************************************************************************
!
! PBE finite volume discretisation for growth
!
!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd(ni, index, dt, growth_source)

!**********************************************************************************************
!
! Growth for finite volume method
! Size can be diameter, volume etc.
!
! Stelios Rigopoulos, Fabian Sewerin, Binxuan Sun
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in)                        :: index
double precision, intent(in)               :: dt ! only for Courant number check
double precision, intent(out)              :: growth_source

double precision :: g_terml,g_termr,phi
double precision :: courant
double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Number density in left cell
double precision :: nll               !< Number density in left-left cell
double precision :: nc                !< Number density in this cell
double precision :: nr                !< Number density in right cell
double precision :: nrr               !< Number density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio (avoids div by zero)
double precision :: rl,rr             !< r+ at left and right surface

parameter(eps = 1.D1*epsilon(1.D0))

!**********************************************************************************************

! Growth rate at right boundary calculation
call calc_growth_rate_liquid(index, .true., g_termr)
! Growth rate at left boundary calculation
call calc_growth_rate_liquid(index-1, .true., g_terml)

! Courant-Friedrichs-Lewy (CFL) condition (C <= 1 for PBE)
courant = g_termr * dt / dv(index)
if (courant>1) then
  write(*,*) "Courant number of ",courant," detected at index ",index,"."
  write(*,*) "Courant number should be <= 1 for growth function to work."
  write(*,*) "Adjust dt proportionately."
  stop
end if

!----------------------------------------------------------------------------------------------
!TVD scheme ref: S.Qamar et al 2006: A comparative study of high resolution schemes for solving
!                population balances in crystallization
!----------------------------------------------------------------------------------------------

if (g_termr>0.D0) then

  ! growth rate is along positive direction
  if (index==1) then

    gnl = 0.0D0
    gnr = g_termr * 0.5 * (ni(1)+ni(2))

  else if (index==m) then

    rl = (ni(m) - ni(m-1) + eps) / (ni(m-1) - ni(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (ni(m-1) + 0.5 * phi * (ni(m-1) - ni(m-2)))
    gnr = g_termr * (ni(m) + 0.5*(ni(m) - ni(m-1)))

  else

    ! Fluxes at cell right surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc + 0.5 * phi * (nc - nl))

    ! Fluxes at cell left surface
    if (index==2 ) then
      gnl = g_terml * 0.5 * (ni(1)+ni(2))
    else
      nl = ni(index-2)
      nc = ni(index-1)
      nr = ni(index)
      rl = (nr - nc + eps) / (nc - nl + eps)
      phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
      gnl = g_terml * (nc + 0.5 * phi * (nc - nl))
    end if
  end if

else

  ! growth rate is along negative direction
  if (index==1) then

    gnl = g_terml * (ni(1) + 0.5 * (ni(1) - ni(2)))
    rr = (ni(1) - ni(2) + eps) / (ni(2) - ni(3) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (ni(2) + 0.5 * phi * (ni(2) - ni(3)))

  else if (index==m) then

    gnr = 0
    gnl = g_terml * 0.5 * (ni(m)+ni(m-1))

  else

    ! Fluxes at cell right surface
    if (index==m-1) then
      gnr = g_termr * 0.5 * (ni(m)+ni(m-1))
    else
      nl = ni(index)
      nc = ni(index+1)
      nr = ni(index+2)
      rr = (nl - nc + eps) / (nc - nr + eps)
      phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
      gnr = g_termr * (nc + 0.5 * phi * (nc - nr))
    end if

    ! Fluxes at cell left surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rl = (nl - nc + eps) / (nc - nr + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc + 0.5 * phi * (nc - nr))

  end if

end if

! Calculate growth source
if (i_gm==1) then
  ! For mass-conservative growth scheme, apply it after the first interval
  if (index>1) then
    growth_source = (v(index-1)*gnl - v(index)*gnr) / (0.5*(v(index)**2-v(index-1)**2))
  else
    growth_source = (gnl - gnr) / dv(index)
  end if
else
  growth_source = (gnl - gnr) / dv(index)
end if

end subroutine growth_tvd

!**********************************************************************************************



!**********************************************************************************************

subroutine calc_growth_rate_liquid(index, boundary, g_term)

!**********************************************************************************************
!
! Calculates droplet growth rate at indexed boundary
!
! By Jack Bartlett (10/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

integer, intent(in)                        :: index
logical, intent(in)                        :: boundary
double precision, intent(out)              :: g_term

double precision surf_tens, accom_coeff
double precision r, f_dry_i, kappa_i, rho_i
double precision raoult_term, kelvin_term, S_droplet
double precision diff_coeff_mod

!----------------------------------------------------------------------------------------------

surf_tens = 72.8D-3 ! Update!
accom_coeff = 1.D0

if (boundary) then
  ! Find properties of interval boundary
  r = ((3.D0*v(index))/(4.D0*pi))**(1.D0/3.D0)

  if (index==0) then
    f_dry_i = f_dry(1)
    kappa_i = kappa(1)
    rho_i = rho(1)
  else if (index==m) then
    f_dry_i = f_dry(m)
    kappa_i = kappa(m)
    rho_i = rho(m)
  else
    f_dry_i = 0.5D0 * (f_dry(index) + f_dry(index+1))
    kappa_i = 0.5D0 * (kappa(index) + kappa(index+1))
    rho_i = 0.5D0 * (rho(index) + rho(index+1))
  end if

else
  ! Find properties of interval midpoint
  r = d_m(index)/2
  f_dry_i = f_dry(index)
  kappa_i = kappa(index)
  rho_i = rho(index)
end if

if (kappa_i.eq.0.D0) then
  raoult_term = 1.D0
else
  raoult_term = (1.D0 - f_dry_i)/(1.D0 - (1.D0 - kappa_i)*f_dry_i)
end if

kelvin_term = exp((2.D0 * surf_tens * water_molar_mass)/(ideal_gas_constant*temperature*water_density*r))

S_droplet = raoult_term * kelvin_term

diff_coeff_mod = diff_coeff / ( r/(r + 0.7*mfp_air) + &
                              & diff_coeff/(r*accom_coeff) * &
                              & sqrt(2*pi*water_molar_mass/(ideal_gas_constant*temperature)) )

g_term = 1.D-6 * 4.D0*pi*r * (diff_coeff_mod * water_molar_mass)/(rho_i * ideal_gas_constant * temperature) * (Pvap - S_droplet*Psat_l)

end subroutine calc_growth_rate_liquid

!**********************************************************************************************



!**********************************************************************************************

subroutine calc_growth_rate_crystal(index, supersaturation, g_term)

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
double precision, intent(in) :: supersaturation ! either wrt liquid or wrt ice
double precision, intent(out) :: g_term

double precision r, J

!----------------------------------------------------------------------------------------------

r = ((3.D0*v(index))/(4.D0*pi))**(1.D0/3.D0) ! Find radius of indexed boundary

call calc_J(r, supersaturation, J)

g_term = water_molecular_vol * J

end subroutine calc_growth_rate_crystal

!**********************************************************************************************



!**********************************************************************************************

subroutine calc_J(r, supersaturation, J)

!**********************************************************************************************
!
! Calculates J (flux of water molecules to droplets/crystals of size r)
!
! By Jack Bartlett (28/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: r
double precision, intent(in) :: supersaturation ! either wrt liquid or wrt ice
double precision, intent(out) :: J

double precision accom_coeff, correction_factor

!----------------------------------------------------------------------------------------------

accom_coeff = 1.D0

correction_factor = 1 + accom_coeff * vapour_thermal_speed * r / (4.D0 * diff_coeff)

J = (pi * r**2 * accom_coeff * vapour_thermal_speed * supersaturation * n_sat) / correction_factor

end subroutine calc_J

!**********************************************************************************************