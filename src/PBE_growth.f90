!**********************************************************************************************
!
! PBE finite volume discretisation for growth
!
!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd_general(indices, dt, growth_source)

!**********************************************************************************************
!
! Growth of general particles for finite volume method
! Size can be diameter, volume etc.
!
! Stelios Rigopoulos, Fabian Sewerin, Binxuan Sun
! Modified by Jack Bartlett (10/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

integer, dimension(4), intent(in)          :: indices
double precision, intent(in)               :: dt ! only for Courant number check
double precision, intent(out)              :: growth_source

integer, dimension(4) :: indices_l ! for the left boundary
double precision, dimension(m) :: nislice
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
call calc_growth_rate_liquid(indices, g_termr)
! Growth rate at left boundary calculation
indices_l = indices
indices_l(1) = indices_l(1) - 1
call calc_growth_rate_liquid(indices_l, g_terml)

! Courant-Friedrichs-Lewy (CFL) condition (C <= 1 for PBE)
courant = g_termr * dt / dv(indices(1))
if (courant>1) then
  write(*,*) "Courant number of ",courant," detected."
  write(*,*) "Courant number should be <= 1 for growth function to work."
  write(*,*) "Adjust dt proportionately."
  stop
end if

!----------------------------------------------------------------------------------------------
!TVD scheme ref: S.Qamar et al 2006: A comparative study of high resolution schemes for solving
!                population balances in crystallization
!----------------------------------------------------------------------------------------------

nislice = ni(:,indices(2),indices(3),indices(4))

if (g_termr>0.D0) then

  ! growth rate is along positive direction
  if (indices(1)==1) then

    gnl = 0.0D0
    gnr = g_termr * 0.5 * (nislice(1)+nislice(2))

  else if (indices(1)==m) then

    rl = (nislice(m) - nislice(m-1) + eps) / (nislice(m-1) - nislice(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nislice(m-1) + 0.5 * phi * (nislice(m-1) - nislice(m-2)))
    gnr = g_termr * (nislice(m) + 0.5*(nislice(m) - nislice(m-1)))

  else

    ! Fluxes at cell right surface
    nl = nislice(indices(1)-1)
    nc = nislice(indices(1))
    nr = nislice(indices(1)+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc + 0.5 * phi * (nc - nl))

    ! Fluxes at cell left surface
    if (indices(1)==2) then
      gnl = g_terml * 0.5 * (nislice(1)+nislice(2))
    else
      nl = nislice(indices(1)-2)
      nc = nislice(indices(1)-1)
      nr = nislice(indices(1))
      rl = (nr - nc + eps) / (nc - nl + eps)
      phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
      gnl = g_terml * (nc + 0.5 * phi * (nc - nl))
    end if
  end if

else

  ! growth rate is along negative direction
  if (indices(1)==1) then

    gnl = g_terml * (nislice(1) + 0.5 * (nislice(1) - nislice(2)))
    rr = (nislice(1) - nislice(2) + eps) / (nislice(2) - nislice(3) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nislice(2) + 0.5 * phi * (nislice(2) - nislice(3)))

  else if (indices(1)==m) then

    gnr = 0
    gnl = g_terml * 0.5 * (nislice(m)+nislice(m-1))

  else

    ! Fluxes at cell right surface
    if (indices(1)==m-1) then
      gnr = g_termr * 0.5 * (nislice(m)+nislice(m-1))
    else
      nl = nislice(indices(1))
      nc = nislice(indices(1)+1)
      nr = nislice(indices(1)+2)
      rr = (nl - nc + eps) / (nc - nr + eps)
      phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
      gnr = g_termr * (nc + 0.5 * phi * (nc - nr))
    end if

    ! Fluxes at cell left surface
    nl = nislice(indices(1)-1)
    nc = nislice(indices(1))
    nr = nislice(indices(1)+1)
    rl = (nl - nc + eps) / (nc - nr + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (nc + 0.5 * phi * (nc - nr))

  end if

end if

! Calculate growth source
if (i_gm==1) then
  ! For mass-conservative growth scheme, apply it after the first interval
  if (indices(1)>1) then
    growth_source = (v(indices(1)-1)*gnl - v(indices(1))*gnr) / (0.5*(v(indices(1))**2-v(indices(1)-1)**2))
  else
    growth_source = (gnl - gnr) / dv(indices(1))
  end if
else
  growth_source = (gnl - gnr) / dv(indices(1))
end if

end subroutine growth_tvd_general

!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd(ani, index, supersaturation, dt, growth_source)

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

double precision, dimension(m), intent(in) :: ani

integer, intent(in)                        :: index
double precision, intent(in)               :: supersaturation ! either wrt liquid or wrt ice
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
call calc_growth_rate_crystal(index, supersaturation, g_termr)
! Growth rate at left boundary calculation
call calc_growth_rate_crystal(index-1, supersaturation, g_terml)

! Courant-Friedrichs-Lewy (CFL) condition (C <= 1 for PBE)
courant = g_termr * dt / dv(index)
if (courant>1) then
  write(*,*) "Courant number of ",courant," detected."
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
    gnr = g_termr * 0.5 * (ani(1)+ani(2))

  else if (index==m) then

    rl = (ani(m) - ani(m-1) + eps) / (ani(m-1) - ani(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (ani(m-1) + 0.5 * phi * (ani(m-1) - ani(m-2)))
    gnr = g_termr * (ani(m) + 0.5*(ani(m) - ani(m-1)))

  else

    ! Fluxes at cell right surface
    nl = ani(index-1)
    nc = ani(index)
    nr = ani(index+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc + 0.5 * phi * (nc - nl))

    ! Fluxes at cell left surface
    if (index==2 ) then
      gnl = g_terml * 0.5 * (ani(1)+ani(2))
    else
      nl = ani(index-2)
      nc = ani(index-1)
      nr = ani(index)
      rl = (nr - nc + eps) / (nc - nl + eps)
      phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
      gnl = g_terml * (nc + 0.5 * phi * (nc - nl))
    end if
  end if

else

  ! growth rate is along negative direction
  if (index==1) then

    gnl = g_terml * (ani(1) + 0.5 * (ani(1) - ani(2)))
    rr = (ani(1) - ani(2) + eps) / (ani(2) - ani(3) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (ani(2) + 0.5 * phi * (ani(2) - ani(3)))

  else if (index==m) then

    gnr = 0
    gnl = g_terml * 0.5 * (ani(m)+ani(m-1))

  else

    ! Fluxes at cell right surface
    if (index==m-1) then
      gnr = g_termr * 0.5 * (ani(m)+ani(m-1))
    else
      nl = ani(index)
      nc = ani(index+1)
      nr = ani(index+2)
      rr = (nl - nc + eps) / (nc - nr + eps)
      phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
      gnr = g_termr * (nc + 0.5 * phi * (nc - nr))
    end if

    ! Fluxes at cell left surface
    nl = ani(index-1)
    nc = ani(index)
    nr = ani(index+1)
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

subroutine calc_growth_rate_liquid(indices, g_term)

!**********************************************************************************************
!
! Calculates droplet growth rate at indexed boundary
!
! By Jack Bartlett (10/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

integer, dimension(4), intent(in) :: indices
double precision, intent(out) :: g_term

double precision surf_tens,r,dry_frac,kappa,raoult_term,kelvin_term,S_droplet,particle_density

!----------------------------------------------------------------------------------------------

surf_tens = 72.8D-3 ! Update!

r = ((3.D0*v_m(indices(1)))/(4.D0*pi))**(1.D0/3.D0) ! Find radius of indexed boundary

dry_frac = max(vf_m(indices(1)) + vf_m(indices(2)) + vf_m(indices(3)), 1.D0)

kappa = (vf_m(indices(1))*comp_kappas(1) + vf_m(indices(2))*comp_kappas(2) + vf_m(indices(3))*comp_kappas(3)) / dry_frac

if (dry_frac.eq.1.D0) then
  raoult_term = (1.D0 - dry_frac)/(1.D0 - (1.D0 - kappa)*dry_frac)
else
  raoult_term = 1.D0
end if

kelvin_term = exp((2.D0 * surf_tens * water_molar_mass)/(ideal_gas_constant*temperature*comp_densities(4)*r))

S_droplet = raoult_term * kelvin_term

particle_density = vf_m(indices(1))*comp_densities(1) + vf_m(indices(2))*comp_densities(2) + vf_m(indices(3))*comp_densities(3) + vf_m(indices(4))*comp_densities(4)

g_term = 4*pi*r * (diff_coeff * water_molar_mass)/(particle_density * ideal_gas_constant * temperature) * (Pvap - S_droplet*Psat_l)

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