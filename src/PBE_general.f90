!**********************************************************************************************
!
! PBE general subroutines
! Stelios Rigopoulos, Anxiong Liu, Binxuan Sun, Daniel O'Sullivan
! Modified 19/12/2017
! Modified 15/01/2018
! Modified 09/05/2018
!
!**********************************************************************************************



!**********************************************************************************************

module pbe_mod

!**********************************************************************************************
!
! Declaration of common variables related to grid and kernels
!
! by Stelios Rigopoulos
! 17/12/2001
! 19/04/2002 new version
!
! Modified 07/05/2017
! Modified 31/05/2017
! Modified 14/07/2020
!
!**********************************************************************************************

implicit none

save

double precision, allocatable, dimension(:) :: v
double precision, allocatable, dimension(:) :: dv
double precision, allocatable, dimension(:) :: v_m
double precision, allocatable, dimension(:) :: d_m
double precision, allocatable, dimension(:) :: nuc
double precision, allocatable, dimension(:) :: g_droplet
double precision, allocatable, dimension(:) :: g_crystal
double precision, allocatable, dimension(:) :: f_droplet

double precision, allocatable, dimension(:) :: ni_droplet
double precision, allocatable, dimension(:) :: ni_crystal
double precision, allocatable, dimension(:) :: kappa
double precision, allocatable, dimension(:) :: rho
double precision, allocatable, dimension(:) :: f_dry

double precision, allocatable, dimension(:) :: savgol_coeffs


double precision v0,grid_lb,grid_rb
double precision agg_kernel_const
double precision break_const
double precision f_dry_tolerance
double precision courant_condition,courant_condition_tight

integer m,grid_type
integer i_gm
integer agg_kernel
integer growth_function
integer order_of_gq
integer solver_pbe
integer do_smoothing, smoothing_window, half_sm_win

! Environment variables
double precision temperature,T_exhaust,T_ambient
double precision P_ambient
double precision Pvap,Pvap_exhaust,Pvap_ambient
double precision Psat_l,Psat_i
double precision saturation_l,saturation_i
double precision air_density
double precision n_sat
double precision vapour_thermal_speed
double precision diffusivity, k_air
double precision mfp_air
double precision l_v
double precision sigma_water

! Plume diffusion constants
double precision r_0,u_0,eps_diffusivity
double precision x_m ! unentrained length (m)
double precision tau_m ! mixing timescale (s)


! Double precision kind
integer, parameter :: dp = selected_real_kind(15, 307)

! Mathematical constants
real(kind=dp), parameter :: pi = 3.14159265358979D0
real(kind=dp), parameter, dimension(-2:2) :: &
    & savgol_coeffs_5 = (/ -3.D0/35.D0, 12.D0/35.D0, 17.D0/35.D0, 12.D0/35.D0, -3.D0/35.D0 /)
real(kind=dp), parameter, dimension(-3:3) :: &
    & savgol_coeffs_7 = 1.D0/21.D0 * (/ -2.D0, 3.D0, 6.D0, 7.D0, 6.D0, 3.D0, -2.D0 /)
real(kind=dp), parameter, dimension(-4:4) :: &
    & savgol_coeffs_9 = 1.D0/231.D0 * (/ -21.D0, 14.D0, 39.D0, 54.D0, 59.D0, 54.D0, 39.D0, 14.D0, -21.D0 /)
real(kind=dp), parameter, dimension(-5:5) :: &
    & savgol_coeffs_11 = 1.D0/429.D0 * (/ -36.D0, 9.D0, 44.D0, 69.D0, 84.D0, 89.D0, 84.D0, 69.D0, 44.D0, 9.D0, -36.D0 /)

! Physical constants
real(kind=dp), parameter :: ideal_gas_constant = 8.314D0 ! (J mol-1 K-1)
real(kind=dp), parameter :: boltzmann_constant = 1.380649D-23 ! (J K-1)
real(kind=dp), parameter :: water_molar_mass = 18.D-3 ! (kg mol-1)
real(kind=dp), parameter :: air_molar_mass = 28.9D-3 ! (kg mol-1)
real(kind=dp), parameter :: avogadro_constant = 6.02214D23 ! (mol-1)
real(kind=dp), parameter :: water_molecular_vol = 2.97D-29 ! (m3)
real(kind=dp), parameter :: water_molecular_mass = 2.99D-26 ! (kg)
real(kind=dp), parameter :: water_density = 1.E3 ! (kgm-3)
real(kind=dp), parameter :: cp_air = 1004 ! (J kg-1 K-1)

end module pbe_mod

!**********************************************************************************************



!**********************************************************************************************

module agg_cfv_mod

implicit none

save

integer :: ID_Ajk, ID_xj, ID_xk, ID_dgj, ID_dgk
integer :: ID_j, ID_k, ID_i
integer :: NDou_AggSec, NInt_AggSec
integer :: Nbeta_Agg
integer :: N_AggSec
integer :: index_Sec

double precision, allocatable, dimension(:) :: Dou_AggSec
integer, allocatable, dimension(:)          :: Int_AggSec
double precision, allocatable, dimension(:) :: beta_AggSec
double precision, allocatable, dimension(:) :: beta_AggSecOrig

double precision :: Pbe_AggrInfo(3)

parameter (ID_Ajk = 1, ID_xj = 2, ID_xk = 3, ID_dgj=4, ID_dgk=5)
parameter (ID_j = 1, ID_k = 2, ID_i = 3)
parameter (NDou_AggSec = 5, NInt_AggSec = 3)

end module agg_cfv_mod

!**********************************************************************************************



!**********************************************************************************************

module frag_cfv_mod

implicit none

save

double precision :: break_a, break_b, break_c, h, delta, lambda ! break dist parameters
integer :: break_section_count, break_global_cntr
double precision, allocatable :: break_kernel_store(:,:)

end module frag_cfv_mod

!**********************************************************************************************



!**********************************************************************************************

module gauss_mod

! Data for Gaussian integration

implicit none

save

double precision, allocatable, dimension(:,:) :: gw, gp

! Gauss data 2-point
double precision, parameter, dimension(4,2) :: gp_2 &
& = reshape((/-1.D0/(3.D0**0.5D0),-1.D0/(3.D0**0.5D0), &
&1.D0/(3.D0**0.5D0),1.D0/(3.D0**0.5D0),-1.D0/(3.D0**0.5D0),1.D0/3.D0**0.5D0,-1.D0/(3.D0**0.5D0),1.D0/(3.D0**0.5D0)/),shape(gp_2))
double precision, parameter :: gw_2 = 1.D0

! Gauss data 3 point
double precision, parameter, dimension(9,2) :: gp_3 &
& = reshape((/ -(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,&
&0.D0,0.D0,0.D0,(6.D0/10.D0)**0.5D0,(6.D0/10.D0)**0.5D0,(6.D0/10.D0)**0.5D0,&
&-(6.D0/10.D0)**0.5D0,0.D0,(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,&
& 0.D0,(6.D0/10.D0)**0.5D0,-(6.D0/10.D0)**0.5D0,0.D0,(6.D0/10.D0)**0.5D0/),shape(gp_3))
double precision, parameter, dimension(9,2) :: gw_3 &
& = reshape((/ 5.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 8.D0/9.D0 , 8.D0/9.D0 , 8.D0/9.D0,&
&5.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 8.D0/9.D0 , 5.D0/9.D0,&
& 5.D0/9.D0 , 8.D0/9.D0 , 5.D0/9.D0 , 5.D0/9.D0 , 8.D0/9.D0 ,5.D0/9.D0/),shape(gw_3))

end module gauss_mod

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_read()

!**********************************************************************************************
!
! Reads PBE data
!
! Stelios Rigopoulos 03/12/2014
! Modified 06/05/2017
! Modified 25/06/2020
! Modified 14/07/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

integer i

!----------------------------------------------------------------------------------------------

! Read PBE parameters

open(30,file='pbe/pbe.in')
do i=1,2
  read(30,*)
end do
read(30,*) agg_kernel
read(30,*) agg_kernel_const
read(30,*) r_0
read(30,*) u_0
read(30,*) eps_diffusivity
read(30,*) T_exhaust
read(30,*) T_ambient
read(30,*) P_ambient
read(30,*) Pvap_exhaust
read(30,*) Pvap_ambient
read(30,*) i_gm
read(30,*) break_const
read(30,*) order_of_gq
read(30,*) grid_type
read(30,*) m
read(30,*) grid_lb
read(30,*) grid_rb
read(30,*) v0
read(30,*) solver_pbe
read(30,*) f_dry_tolerance
read(30,*) courant_condition
read(30,*) courant_condition_tight
read(30,*) do_smoothing
read(30,*) smoothing_window
close(30)


end subroutine pbe_read

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_init()

!**********************************************************************************************
!
! Initialises PBE data
!
! Stelios Rigopoulos 03/12/2014
! Modified 06/05/2017
! Modified 25/06/2020
! Modified 14/07/2020
!
!**********************************************************************************************

use pbe_mod
use frag_cfv_mod
use gauss_mod

implicit none

!----------------------------------------------------------------------------------------------

! Initialise temperature, vapour pressure, and other environment variables
call pbe_set_environment(0.D0)

! Initialise particle distributions

ni_droplet = 0.D0
ni_crystal = 0.D0

call pbe_read_species()


! Initialise nucleation
nuc = 0.D0

! Initialise aggregation
if (agg_kernel>0) then
  call PBE_agg_fvLocate(v,1)
  call PBE_agg_fvLocate(v,2)
  call PBE_agg_beta(1)
  call PBE_agg_beta(2)
end if

if (break_const>0.) then
  call gauss_init(order_of_gq)
  call pbe_breakage_calc(1)
  allocate(break_kernel_store(break_global_cntr,3))
  call pbe_breakage_calc(2)
end if

half_sm_win = (smoothing_window-1)/2
allocate(savgol_coeffs(-half_sm_win:half_sm_win))
if (smoothing_window == 5) then
  savgol_coeffs = savgol_coeffs_5
else if (smoothing_window == 7) then
  savgol_coeffs = savgol_coeffs_7
else if (smoothing_window == 9) then
  savgol_coeffs = savgol_coeffs_7
else if (smoothing_window == 11) then
  savgol_coeffs = savgol_coeffs_11
else
  write(*,*) "Smoothing window of size ",smoothing_window," not accepted. Stopping."
  stop
end if

end subroutine pbe_init

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_read_species()

!**********************************************************************************************
!
! Initialises species data
!
! Jack Bartlett (19/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision r_m, dlogr, n
double precision n_tot, r_mean, sigma, kappa_i, rho_i, f_dry_i
integer i, num_species

character(len=256) :: line

!----------------------------------------------------------------------------------------------

num_species = 0

open(30,file='pbe/species.in')
do i=1,4
  read(30,*)
end do

do
  read(30, fmt='(A)') line ! fmt ensures the whole line is read as a string

  if (trim(line) == "END") then
    exit
  end if

  read(line, fmt=*) n_tot, r_mean, sigma, kappa_i, rho_i, f_dry_i

  if (kappa_i < 0.D0) then
    write(*,*) "Hygroscopicity of ",kappa_i," found in input."
    write(*,*) "Hygroscopicity must be more than 0. Stopping."
    stop
  else if (rho_i < 0.D0) then
    write(*,*) "Density of ",rho_i," found in input."
    write(*,*) "Density must be more than 0. Stopping."
    stop
  else if ((f_dry_i < 0.D0).or.(f_dry_i > 1.D0)) then
    write(*,*) "Dry volume fraction of ",f_dry_i," found in input."
    write(*,*) "Dry volume fraction must be between 0 and 1. Stopping."
    stop
  end if

  ! Add species to number density array
  do i=1,m
    r_m = d_m(i)/2.D0
    dlogr = 1.D0/3.D0 * log(v(i)/v(i-1))
    ! Particle number density per interval width to be added
    n = 1.D0/(dv(i)) * n_tot * dlogr * (1.D0/(sqrt(2.D0*pi)*sigma)) * exp(-(log(r_m/r_mean))**2.D0/(2.D0*sigma**2.D0))

    ! Update average properties in interval
    kappa(i) = (kappa(i)*ni_droplet(i) + kappa_i*n) / (ni_droplet(i) + n)
    rho(i) = (rho(i)*ni_droplet(i) + rho_i*n) / (ni_droplet(i) + n)
    f_dry(i) = (f_dry(i)*ni_droplet(i) + f_dry_i*n) / (ni_droplet(i) + n)

    ! Update number density per interval width
    ni_droplet(i) = ni_droplet(i) + n
  end do

  num_species = num_species + 1
end do
close(30)

write(*,*) "Read in ",num_species," species."

end subroutine pbe_read_species

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_grid()

!**********************************************************************************************
!
! Subroutine grid
! Generates/reads the grid
!
! By Stelios Rigopoulos
! 31/10/2001 - initial version
! Modified 01/06/2017
! Modified 12/05/2018
! Modified 14/07/2020
!
! m                     number of points
! ni_droplet            population density of droplets
! v                     independent variable - volume
! dv                    lengths of intervals
! v_m                   interval mid-points
! d_m                   diameter corresponding to mid-point
! v0                    size of nuclei
! grid_lb               left boundary
! grid_rb               right boundary
!
!**********************************************************************************************e

use pbe_mod

implicit none

double precision alpha,v1,v2

integer allocated_size_bytes
integer i

!----------------------------------------------------------------------------------------------

! Allocate arrays
allocated_size_bytes = 11*8*m
if (allocated_size_bytes > 1000**3) then
  write(*,*) "Allocating around",(allocated_size_bytes/(1000**3))," GB to arrays."
else if (allocated_size_bytes > 1000**2) then
  write(*,*) "Allocating around",(allocated_size_bytes/(1000**2))," MB to arrays."
else if (allocated_size_bytes > 1000) then
  write(*,*) "Allocating around",(allocated_size_bytes/(1000))," kB to arrays."
else
  write(*,*) "Allocating around",(allocated_size_bytes)," bytes to arrays."
end if
allocate(v(0:m),dv(m),v_m(m),d_m(m),nuc(m),g_droplet(0:m),g_crystal(0:m),f_droplet(0:m))
allocate(ni_droplet(m),ni_crystal(m),kappa(m),rho(m),f_dry(m))

! Calculate grid
if (grid_type==1) then

  !Option 1: geometric grid
  if (grid_lb==0) then
    write(*,*) "Left grid boundary must be greater than zero for geometric grid."
    stop
  end if
  alpha = exp(log(grid_rb/grid_lb)/m) ! Calculates geometric ratio
  v(0) = grid_lb
  do i=1,m
    v(i) = v(i-1) * alpha
  end do

  !Previous solution:
  !v(0) = grid_lb
  !v(1) = v0 + (v0 - grid_lb)
  !v1 = v(1) - grid_lb
  !v2 = grid_rb - grid_lb
  !call inc_ratio(v1,v2,m,alpha)
  !do i=2,m
  !  v(i) = v(0)+(v(1)-v(0))*(1-alpha**real(i))/(1-alpha)
  !end do

  write(*,*) 'left boundary: ',v(0)
  write(*,*) 'first node: ',v(1)
  write(*,*) 'right boundary: ',v(m)
  write(*,*) 'ratio: ',alpha

else if (grid_type==2) then

  !Option 2: uniform grid
  alpha = (grid_rb - grid_lb)/m
  v(0) = grid_lb
  do i=1,m
    v(i) = v(i-1) + alpha
  end do
  write(*,*) 'left boundary: ',v(0)
  write(*,*) 'first node: ',v(1)
  write(*,*) 'right boundary: ',v(m)
  write(*,*) 'increment: ',alpha

end if

!Interval length
do i=1,m
  dv(i) = v(i)-v(i-1)
end do

!Determine mid-points
do i=1,m
  v_m(i) = v(i-1)+0.5D0*dv(i)
end do

! Determine diameters of mid-points
do i=1,m
  d_m(i) = (6.D0/pi*v_m(i))**(1.D0/3.D0)
end do

end subroutine pbe_grid

!**********************************************************************************************



!**********************************************************************************************

subroutine inc_ratio(lb,rb,k,q)

!**********************************************************************************************
!
! Calculation of common ratio in exponential grid (obsolete)
!
! By Binxuan Sun
! 28/02/2018
! lb: left boundary, which is also the first term a1 in geometric progression
! rb: right boundary
! k:  number of total elements
! q:  common ratio
!
! Sum a(1,2,3,....,m) = a1 * ( 1 - q^m )/ ( 1 - q )
! Sum = the length of the interval [0,rb]
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in)    :: lb
double precision, intent(in)    :: rb
integer, intent(in)             :: k
double precision, intent(inout) :: q

double precision s,r1,r2,rm,rt
double precision t1,t2,a1,t
double precision abs_error

integer i,IMAX

!**********************************************************************************************

IMAX = 1000 !maximum iteration times
r1 = 1.0001
r2 = 20
abs_error= 1.0D-6

! a1 is the first number in the series
a1 = lb
!s = ( rb - lb )/ a1
s = rb / a1
t1 = r1**k - s * r1 + s -1
t2 = r2**k - s * r2 + s -1

if ((t1.GE.0.0).OR.(t2.LE.0.0)) then
  write(*,*) "Error in the right grid boundary and the number of nodes, re-select them"
  stop
end if
rt = 0

if ((t1.LT.0.0).AND.(t2.GT.0.0)) then
  do i=1,IMAX
    rm = ( r1+ r2 )*0.5
    if (abs(rm-rt).LE.abs_error) then
      go to 666
    end if
    rt = rm
    q = rm
    t = q**k - s * q + s - 1
    if (t.GT.0) then
      r2 = rm
    else
      r1 = rm
    end if
  end do
end if

666 q = rm

write(*,*) "alpha = ", q

end subroutine inc_ratio

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_set_environment(current_time)

!**********************************************************************************************
!
! Set (or update) environment variables (temperature and pressure)
!
! By Jack Bartlett (22/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: current_time

double precision dilution_factor, mixing_grad, temperature_new

!----------------------------------------------------------------------------------------------

! Initialisation
if (current_time.eq.0.D0) then
  temperature = T_exhaust
  Pvap = Pvap_exhaust
  x_m = r_0 * sqrt(2.D0/eps_diffusivity)
  tau_m = x_m / u_0

! Updating
else

  if (current_time.le.tau_m) then
    dilution_factor = 1.D0
  else
    dilution_factor = (tau_m / current_time)**(0.9D0)
  end if

  mixing_grad = (Pvap - Pvap_ambient) / (temperature - T_ambient) ! Schumann's G

  temperature_new = T_ambient + (T_exhaust - T_ambient) * dilution_factor
  ! Pvap can't be calculated the same way because it is also updated by growth
  ! However, it does decrease like temperature due to mixing
  Pvap = Pvap + mixing_grad * (temperature_new - temperature)

  temperature = temperature_new

end if

! Saturation vapour pressures - Buck
Psat_l = 6.1121D2*exp((18.678D0 - (temperature-273.15D0)/234.5D0) * ((temperature-273.15D0)/(temperature-16.01D0)))
Psat_i = 6.1115D2*exp((23.036D0 - (temperature-273.15D0)/333.7D0) * ((temperature-273.15D0)/(temperature+6.67D0)))

! Saturation ratios
saturation_l = Pvap/Psat_l
saturation_i = Pvap/Psat_i

! Density of air (kg m-3)
air_density = P_ambient * air_molar_mass / (ideal_gas_constant * temperature)

! Mean free path of air (m) - assumes effective collision diameter 0.37 nm
mfp_air = (boltzmann_constant*temperature)/(sqrt(2*pi)*(0.37D-9)**2*P_ambient)

! H2O density
!comp_densities(4) = 1.D3

! H2O number concentration at water saturation (m-3)
n_sat = avogadro_constant * Pvap / (ideal_gas_constant * temperature)

! Water surface tension (N m-1) - IAPWS
sigma_water = 235.8D-3 * (1-temperature/647.096D0)**(1.256D0) * (1.D0 - 0.625D0*(1-temperature/647.096D0))

! Thermal speed of water vapour (m s-1)
vapour_thermal_speed = sqrt(8 * boltzmann_constant * temperature / (pi*water_molecular_mass))

! Diffusivity (m2 s-1) - Pruppacher and Klett eq. 13-3
diffusivity = 2.11D-5 * (temperature/273.15D0)**(1.94D0) * (101325.D0 / P_ambient)

! Thermal conductivity of air (J m-1 s-1 K-1) - Seinfeld and Pandis eq. 17.54
k_air = 1.D-3 * (4.39D0 + 0.071D0 * temperature)

! Specific latent heat of vaporisation of water (J kg-1)
l_v = 2.501D6 - 2.37D3 * (temperature - 273.15D0)


end subroutine pbe_set_environment

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_update_g_arrays()

!**********************************************************************************************
!
! Update arrays of growth terms
!
! By Jack Bartlett (01/07/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision r
integer index

!----------------------------------------------------------------------------------------------

do index=0,m
  ! Droplets
  r = ((3.D0*v(index))/(4.D0*pi))**(1.D0/3.D0)
  if (index==m) then
    call calc_growth_rate_liquid(r, kappa(index), rho(index), f_dry(index), g_droplet(index))
  else
    call calc_growth_rate_liquid(r, kappa(index+1), rho(index+1), f_dry(index+1), g_droplet(index))
  end if

  ! Crystals
  call calc_growth_rate_crystal(index, g_crystal(index))
end do


end subroutine pbe_update_g_arrays

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_smooth_props()

!**********************************************************************************************
!
! Smooth property arrays to avoid numerical errors
!
! By Jack Bartlett (03/07/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision kappa_smooth(m),rho_smooth(m),f_dry_smooth(m)
double precision sum_savgol
integer i, j

!----------------------------------------------------------------------------------------------

! Hygroscopicity

!do i=1,half_sm_win
!  kappa_smooth(i) = kappa(i)
!  kappa_smooth(m+1-i) = kappa(m+1-i)
!end do

!do i=half_sm_win+1,m-half_sm_win
!  sum_savgol = 0.D0
!  do j=-half_sm_win,half_sm_win
!    sum_savgol = sum_savgol + savgol_coeffs(j) * kappa(i+j)
!  end do
!  kappa_smooth(i) = sum_savgol
!end do

!kappa = kappa_smooth


! Density

do i=1,half_sm_win
  rho_smooth(i) = rho(i)
  rho_smooth(m+1-i) = rho(m+1-i)
end do

do i=half_sm_win+1,m-half_sm_win
  sum_savgol = 0.D0
  do j=-half_sm_win,half_sm_win
    sum_savgol = sum_savgol + savgol_coeffs(j) * rho(i+j)
  end do
  rho_smooth(i) = sum_savgol
end do

rho = rho_smooth


! Dry fraction

do i=1,half_sm_win
  f_dry_smooth(i) = f_dry(i)
  f_dry_smooth(m+1-i) = f_dry(m+1-i)
end do

do i=half_sm_win+1,m-half_sm_win
  sum_savgol = 0.D0
  do j=-half_sm_win,half_sm_win
    sum_savgol = sum_savgol + savgol_coeffs(j) * f_dry(i+j)
  end do
  f_dry_smooth(i) = sum_savgol
end do

f_dry = f_dry_smooth


end subroutine pbe_smooth_props

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_freezing(dt)

!**********************************************************************************************
!
! Update droplet freezing each time step
!
! By Jack Bartlett (30/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: dt

double precision ice_germ_rate

integer index

!----------------------------------------------------------------------------------------------

! Add Koop 2000 freezing rate based on water activity

! Ice germ rate is number of ice germs formed per droplet vol per second
ice_germ_rate = 1.D6 * exp(-3.5714D0 * temperature + 858.719D0)

do index=1,m
  ! Condition to have >= 1 ice germ per droplet
  if ((ice_germ_rate * v_m(index) * (1.D0 - f_dry(index)) * dt).ge.(1.D0)) then
    ni_crystal(index) = ni_crystal(index) + ni_droplet(index)
    ni_droplet(index) = 0.D0
  end if
end do

end subroutine pbe_freezing

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_moments(ni,moment,meansize)

!**********************************************************************************************
!
! Calculation of zeroth and first moments
!
! By Stelios Rigopoulos
! Modified 06/05/2017
! Modified 04/07/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in)    :: ni
double precision, dimension(0:1), intent(out) :: moment
double precision, intent(out)                 :: meansize

double precision M1_lp,lp

integer i

!----------------------------------------------------------------------------------------------

moment(0) = 0.D0
moment(1) = 0.D0

do i=1,m
  moment(0) = moment(0) + ni(i)*dv(i)
  moment(1) = moment(1) + 0.5D0*(v(i-1)+v(i))*ni(i)*dv(i)
end do

M1_lp = 0.D0
do i=m-5,m
  M1_lp = M1_lp + 0.5D0*(v(i-1)+v(i))*ni(i)*dv(i)
end do

!lp = M1_lp/moment(1)
!if (lp.gt.0.001) then
!  write(*,*) 'warning, more than 0.1% of mass in the last five nodes'
!end if

if (moment(0).gt.1.D-10) then
  meansize = moment(1)/moment(0)
else
  meansize = 0.D0
end if

end subroutine pbe_moments

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_output_psd(ni,filename,current_time,first_write)

!**********************************************************************************************
!
! Writes PSD
!
! By Jack Bartlett (29/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ni
character(len=30), intent(in) :: filename
double precision, intent(in) :: current_time
logical, intent(in) :: first_write

double precision :: nitemp(m)
double precision, dimension(0:1) :: moment

double precision meansize

integer i

!----------------------------------------------------------------------------------------------

call pbe_moments(ni,moment,meansize)

do i=1,m
  if (abs(ni(i))<1.D-16) then
    nitemp(i) = 0.D0
  else
    nitemp(i) = ni(i)
  end if
end do


if (first_write) then
  open(99,file=filename,status='replace')
else
  open(99,file=filename,status='old',position='append')
end if
do i=1,m
  write(99,1001) current_time, v_m(i), dv(i), d_m(i), nitemp(i), &
  & nitemp(i)*dv(i)/moment(0), v_m(i)*nitemp(i), v_m(i)*nitemp(i)*dv(i)/moment(1)
end do
close(99)


1001 format(8E20.10)

end subroutine pbe_output_psd

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_output_property(property,filename,current_time,first_write)

!**********************************************************************************************
!
! Writes property
!
! By Jack Bartlett (20/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: property
character(len=30), intent(in) :: filename
double precision, intent(in) :: current_time
logical, intent(in) :: first_write

integer i

!----------------------------------------------------------------------------------------------

if (first_write) then
  open(99,file=filename,status='replace')
else
  open(99,file=filename,status='old',position='append')
end if
do i=1,m
  write(99,1003) current_time, v_m(i), property(i)
end do
close(99)


1003 format(3E20.10)

end subroutine pbe_output_property

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_output_env(current_time,first_write)

!**********************************************************************************************
!
! Writes environment variables
!
! By Jack Bartlett (29/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: current_time
logical, intent(in) :: first_write

integer i

!----------------------------------------------------------------------------------------------


! Write environment variables to end of environment_variables.out
if (first_write) then
  open(99,file='output/environment_variables.out',status='replace')
else
  open(99,file='output/environment_variables.out',status='old',position='append')
end if
write(99,1002) current_time, temperature, Pvap, Psat_l, Psat_i
close(99)


1002 format(5E20.10)

end subroutine pbe_output_env

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_output_growth(current_time,first_write)

!**********************************************************************************************
!
! Writes growth rate at interval boundaries
!
! By Jack Bartlett (01/07/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: current_time
logical, intent(in) :: first_write

integer i

!----------------------------------------------------------------------------------------------


! Write environment variables to end of growth_rate.out
if (first_write) then
  open(99,file='output/growth_rate.out',status='replace')
else
  open(99,file='output/growth_rate.out',status='old',position='append')
end if
do i=0,m
  write(99,1004) current_time, v(i), g_droplet(i), g_crystal(i)
end do
close(99)

1004 format(4E20.10)

end subroutine pbe_output_growth

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_deallocate()

!**********************************************************************************************
!
! Deallocate PBE variables
!
! By Stelios Rigopoulos
! Modified 06/05/2017
! Modified 04/07/2020
!
!**********************************************************************************************

use pbe_mod

deallocate(v,dv,v_m,d_m,nuc,g_droplet,g_crystal,f_droplet,ni_droplet,ni_crystal,kappa,rho,f_dry)
deallocate(savgol_coeffs)

end subroutine pbe_deallocate

!**********************************************************************************************