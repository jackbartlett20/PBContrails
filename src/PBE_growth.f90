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

subroutine growth_tvd_rho(ni, rho_temp, index, dt, growth_source)

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
double precision, intent(in)               :: dt ! only for Courant number check
double precision, intent(out)              :: growth_source

double precision :: g_terml,g_termr,phi
double precision :: courant
double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Quantity density in left cell
double precision :: nll               !< Quantity density in left-left cell
double precision :: nc                !< Quantity density in this cell
double precision :: nr                !< Quantity density in right cell
double precision :: nrr               !< Quantity density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio (avoids div by zero)
double precision :: rl,rr             !< r+ at left and right surface
double precision :: rho_after_a,rho_after_b !< Temporary variables for holding density

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
    call rho_growth(2, 1, rho_temp(2), rho_after_a)
    gnr = g_termr * 0.5 * (ni(1)*rho_temp(1) + ni(2)*rho_after_a)

  else if (index==m) then

    rl = (ni(m) - ni(m-1) + eps) / (ni(m-1) - ni(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    call rho_growth(m-1, m, rho_temp(m-1), rho_after_a)
    call rho_growth(m-2, m, rho_temp(m-2), rho_after_b)
    gnl = g_terml * (ni(m-1)*rho_after_a + 0.5 * phi * (ni(m-1)*rho_after_a - ni(m-2)*rho_after_b))
    gnr = g_termr * (ni(m)*rho_temp(m) + 0.5*(ni(m)*rho_temp(m) - ni(m-1)*rho_after_a))

  else

    ! Fluxes at cell right surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    call rho_growth(index-1, index, rho_temp(index-1), rho_after_a)
    gnr = g_termr * (nc*rho_temp(index) + 0.5 * phi * (nc*rho_temp(index) - nl*rho_after_a))

    ! Fluxes at cell left surface
    if (index==2 ) then
      call rho_growth(1, 2, rho_temp(1), rho_after_a)
      gnl = g_terml * 0.5 * (ni(1)*rho_after_a + ni(2)*rho_temp(2))
    else
      nl = ni(index-2)
      nc = ni(index-1)
      nr = ni(index)
      rl = (nr - nc + eps) / (nc - nl + eps)
      phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
      call rho_growth(index-1, index, rho_temp(index-1), rho_after_a)
      call rho_growth(index-2, index, rho_temp(index-2), rho_after_b)
      gnl = g_terml * (nc*rho_after_a + 0.5 * phi * (nc*rho_after_a - nl*rho_after_b))
    end if
  end if

else

  ! growth rate is along negative direction
  if (index==1) then

    call rho_growth(2, 1, rho_temp(2), rho_after_a)
    gnl = g_terml * (ni(1)*rho_temp(1) + 0.5 * (ni(1)*rho_temp(1) - ni(2)*rho_after_a))
    rr = (ni(1) - ni(2) + eps) / (ni(2) - ni(3) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    call rho_growth(3, 1, rho_temp(3), rho_after_b)
    gnr = g_termr * (ni(2)*rho_after_a + 0.5 * phi * (ni(2)*rho_after_a - ni(3)*rho_after_b))

  else if (index==m) then

    gnr = 0
    call rho_growth(m-1, m, rho_temp(m-1), rho_after_a)
    gnl = g_terml * 0.5 * (ni(m)*rho_temp(m) + ni(m-1)*rho_after_a)

  else

    ! Fluxes at cell right surface
    if (index==m-1) then
      call rho_growth(m, m-1, rho_temp(m), rho_after_a)
      gnr = g_termr * 0.5 * (ni(m)*rho_after_a + ni(m-1)*rho_temp(m-1))
    else
      nl = ni(index)
      nc = ni(index+1)
      nr = ni(index+2)
      rr = (nl - nc + eps) / (nc - nr + eps)
      phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
      call rho_growth(index+1, index, rho_temp(index+1), rho_after_a)
      call rho_growth(index+2, index, rho_temp(index+2), rho_after_b)
      gnr = g_termr * (nc*rho_after_a + 0.5 * phi * (nc*rho_after_a - nr*rho_after_b))
    end if

    ! Fluxes at cell left surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rl = (nl - nc + eps) / (nc - nr + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    call rho_growth(index+1, index, rho_temp(index+1), rho_after_a)
    gnl = g_terml * (nc*rho_temp(index) + 0.5 * phi * (nc*rho_temp(index) - nr*rho_after_a))

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

subroutine growth_tvd_f_dry(ni, f_dry_temp, index, dt, growth_source)

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
double precision, intent(in)               :: dt ! only for Courant number check
double precision, intent(out)              :: growth_source

double precision :: g_terml,g_termr,phi
double precision :: courant
double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Quantity density in left cell
double precision :: nll               !< Quantity density in left-left cell
double precision :: nc                !< Quantity density in this cell
double precision :: nr                !< Quantity density in right cell
double precision :: nrr               !< Quantity density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio (avoids div by zero)
double precision :: rl,rr             !< r+ at left and right surface
double precision :: f_dry_after_a,f_dry_after_b !< Temporary variables for holding f_dry

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
    call f_dry_growth(2, 1, f_dry_temp(2), f_dry_after_a)
    gnr = g_termr * 0.5 * (ni(1)*f_dry_temp(1) + ni(2)*f_dry_after_a)

  else if (index==m) then

    rl = (ni(m) - ni(m-1) + eps) / (ni(m-1) - ni(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    call f_dry_growth(m-1, m, f_dry_temp(m-1), f_dry_after_a)
    call f_dry_growth(m-2, m, f_dry_temp(m-2), f_dry_after_b)
    gnl = g_terml * (ni(m-1)*f_dry_after_a + 0.5 * phi * (ni(m-1)*f_dry_after_a - ni(m-2)*f_dry_after_b))
    gnr = g_termr * (ni(m)*f_dry_temp(m) + 0.5*(ni(m)*f_dry_temp(m) - ni(m-1)*f_dry_after_a))

  else

    ! Fluxes at cell right surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    call f_dry_growth(index-1, index, f_dry_temp(index-1), f_dry_after_a)
    gnr = g_termr * (nc*f_dry_temp(index) + 0.5 * phi * (nc*f_dry_temp(index) - nl*f_dry_after_a))

    ! Fluxes at cell left surface
    if (index==2 ) then
      call f_dry_growth(1, 2, f_dry_temp(1), f_dry_after_a)
      gnl = g_terml * 0.5 * (ni(1)*f_dry_after_a + ni(2)*f_dry_temp(2))
    else
      nl = ni(index-2)
      nc = ni(index-1)
      nr = ni(index)
      rl = (nr - nc + eps) / (nc - nl + eps)
      phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
      call f_dry_growth(index-1, index, f_dry_temp(index-1), f_dry_after_a)
      call f_dry_growth(index-2, index, f_dry_temp(index-2), f_dry_after_b)
      gnl = g_terml * (nc*f_dry_after_a + 0.5 * phi * (nc*f_dry_after_a - nl*f_dry_after_b))
    end if
  end if

else

  ! growth rate is along negative direction
  if (index==1) then

    call f_dry_growth(2, 1, f_dry_temp(2), f_dry_after_a)
    gnl = g_terml * (ni(1)*f_dry_temp(1) + 0.5 * (ni(1)*f_dry_temp(1) - ni(2)*f_dry_after_a))
    rr = (ni(1) - ni(2) + eps) / (ni(2) - ni(3) + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    call f_dry_growth(3, 1, f_dry_temp(3), f_dry_after_b)
    gnr = g_termr * (ni(2)*f_dry_after_a + 0.5 * phi * (ni(2)*f_dry_after_a - ni(3)*f_dry_after_b))

  else if (index==m) then

    gnr = 0
    call f_dry_growth(m-1, m, f_dry_temp(m-1), f_dry_after_a)
    gnl = g_terml * 0.5 * (ni(m)*f_dry_temp(m) + ni(m-1)*f_dry_after_a)

  else

    ! Fluxes at cell right surface
    if (index==m-1) then
      call f_dry_growth(m, m-1, f_dry_temp(m), f_dry_after_a)
      gnr = g_termr * 0.5 * (ni(m)*f_dry_after_a + ni(m-1)*f_dry_temp(m-1))
    else
      nl = ni(index)
      nc = ni(index+1)
      nr = ni(index+2)
      rr = (nl - nc + eps) / (nc - nr + eps)
      phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
      call f_dry_growth(index+1, index, f_dry_temp(index+1), f_dry_after_a)
      call f_dry_growth(index+2, index, f_dry_temp(index+2), f_dry_after_b)
      gnr = g_termr * (nc*f_dry_after_a + 0.5 * phi * (nc*f_dry_after_a - nr*f_dry_after_b))
    end if

    ! Fluxes at cell left surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rl = (nl - nc + eps) / (nc - nr + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    call f_dry_growth(index+1, index, f_dry_temp(index+1), f_dry_after_a)
    gnl = g_terml * (nc*f_dry_temp(index) + 0.5 * phi * (nc*f_dry_temp(index) - nr*f_dry_after_a))

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
double precision diffusivity_mod, k_air_mod
double precision F_d, F_k

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

diffusivity_mod = diffusivity / ( r/(r + 0.7*mfp_air) + &
                                & diffusivity/(r*accom_coeff) * &
                                & sqrt(2*pi*water_molar_mass/(ideal_gas_constant*temperature)) )

k_air_mod = k_air / ( r/(r + 0.7*mfp_air) + &
                    & k_air/(r*accom_coeff*air_density*cp_air) * &
                    & sqrt(2*pi*air_molar_mass/(ideal_gas_constant*temperature)) )

F_d = (water_density * ideal_gas_constant * temperature) / (Psat_l * diffusivity_mod * water_molar_mass)

F_k = (l_v * water_density)/(k_air_mod * temperature) * (l_v*water_molar_mass/(ideal_gas_constant*temperature) - 1)

g_term = (4*pi*r) * 1/(F_d + F_k) * (Pvap/Psat_l - S_droplet)

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

correction_factor = 1 + accom_coeff * vapour_thermal_speed * r / (4.D0 * diffusivity)

J = (pi * r**2 * accom_coeff * vapour_thermal_speed * supersaturation * n_sat) / correction_factor

end subroutine calc_J

!**********************************************************************************************