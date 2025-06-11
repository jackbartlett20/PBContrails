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
double precision, allocatable, dimension(:) :: vf
double precision, allocatable, dimension(:) :: vf_m
double precision, allocatable, dimension(:) :: nuc

double precision, allocatable, dimension(:) :: ni_droplet
double precision, allocatable, dimension(:) :: ni_crystal
double precision, allocatable, dimension(:) :: ni_soot
double precision, allocatable, dimension(:,:,:,:) :: ni
double precision, allocatable, dimension(:,:,:,:) :: niprime
double precision, allocatable, dimension(:) :: nislice_m
double precision, allocatable, dimension(:) :: nislice_vf

double precision v0,grid_lb,grid_rb
double precision agg_kernel_const
double precision break_const

integer m,grid_type
integer i_gm
integer agg_kernel
integer growth_function
integer order_of_gq

! Environment variables
double precision temperature,T_exhaust,T_ambient
double precision P_ambient
double precision Pvap,Pvap_exhaust,Pvap_ambient
double precision Psat_l,Psat_i
double precision supersaturation_l,supersaturation_i
double precision mixing_grad
double precision n_sat
double precision vapour_thermal_speed
double precision diff_coeff

! Soot initial distribution
double precision n_soot,r_mean_soot,sigma_soot

! Particle composition parameters
integer n_vf
double precision vf_width
integer, parameter :: n_components = 4
double precision, dimension(n_components) :: comp_densities
double precision, dimension(n_components) :: comp_kappas


! Double precision kind
integer, parameter :: dp = selected_real_kind(15, 307)

! Mathematical constants
real(kind=dp), parameter :: pi = 3.141592654D0

! Physical constants
real(kind=dp), parameter :: ideal_gas_constant = 8.314D0 ! (J mol-1 K-1)
real(kind=dp), parameter :: boltzmann_constant = 1.380649D-23 ! (J K-1)
real(kind=dp), parameter :: water_molar_mass = 0.01801528D0 ! (kg mol-1)
real(kind=dp), parameter :: avogadro_constant = 6.02214D23 ! (mol-1)
real(kind=dp), parameter :: water_molecular_vol = 2.97D-29 ! (m3)
real(kind=dp), parameter :: water_molecular_mass = 2.99D-26 ! (kg)

! Plume diffusion constants
real(kind=dp), parameter :: eps_diffusivity = 0.0285D0 ! turbulent diffusivity ()
real(kind=dp), parameter :: r_0 = 0.5D0 ! jet radius at exhaust (m)
real(kind=dp), parameter :: x_m = r_0 * sqrt(2.D0/eps_diffusivity) ! unentrained length (m)
real(kind=dp), parameter :: u_0 = 4.D2 ! exhaust velocity (m/s)
real(kind=dp), parameter :: tau_m = x_m / u_0 ! mixing timescale (s)

! Soot activation constants
real(kind=dp), parameter :: r_k = 1.D-9 ! Approx Kelvin radius for dry particle activation (m)
real(kind=dp), parameter :: soot_solubility = 0.1D0 ! Solubility parameter for soot, a guess ()

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
read(30,*) T_exhaust
read(30,*) T_ambient
read(30,*) P_ambient
read(30,*) Pvap_exhaust
read(30,*) Pvap_ambient
read(30,*) n_soot
read(30,*) r_mean_soot
read(30,*) sigma_soot
read(30,*) i_gm
read(30,*) break_const
read(30,*) order_of_gq
read(30,*) grid_type
read(30,*) m
read(30,*) grid_lb
read(30,*) grid_rb
read(30,*) n_vf
read(30,*) v0
close(30)

! Read component parameters

open(30,file='pbe/components.in')
do i=1,4
  read(30,*)
end do
do i=1,(n_components-1)
  read(30,*) comp_densities(i), comp_kappas(i)
end do
close(30)

comp_kappas(4) = 0.D0


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

double precision r_m
double precision n_tot, r_mean, sigma, vf1, vf2, vf3, vf4

character(len=256) :: line

integer i, i_vf1, i_vf2, i_vf3

!----------------------------------------------------------------------------------------------

! Initialise temperature, vapour pressure, and other environment variables
call pbe_set_environment(0.D0)

! Initialise particle composition distribution

ni = 0.D0

open(30,file='pbe/species.in')
do i=1,4
  read(30,*)
end do

do
  read(30, fmt='(A)') line ! fmt ensures the whole line is read as a string

  if (trim(line) == "END") then
    exit
  end if

  read(line, fmt=*) n_tot, r_mean, sigma, vf1, vf2, vf3, vf4

  if ((vf1+vf2+vf3+vf4 < 0.9999).or.(vf1+vf2+vf3+vf4 > 1.0001)) then
    write(*,*) "Volume fractions for a species must sum to 1. Stopping."
    stop
  end if

  ! Determine relevant volume fraction intervals
  ! Add in something clever to optimise volume fractions for intervals
  do i=1,n_vf
    if ((vf1.ge.vf(i-1)).and.(vf1.le.vf(i))) then
      i_vf1 = i
    end if
    if ((vf2.ge.vf(i-1)).and.(vf2.le.vf(i))) then
      i_vf2 = i
    end if
    if ((vf3.ge.vf(i-1)).and.(vf3.le.vf(i))) then
      i_vf3 = i
    end if
  end do
  write(*,*) vf_m(i_vf1), vf_m(i_vf2), vf_m(i_vf3)

  ! Add species to number density array
  do i=1,m
    r_m = d_m(i)/2.D0
    ni(i,i_vf1,i_vf2,i_vf3) = ni(i,i_vf1,i_vf2,i_vf3) + (1.D0/(dv(i)*vf_width**3)) * &
    & n_tot * (1.D0/(sqrt(2.D0*pi)*sigma)) * exp(-(log(r_m/r_mean))**2.D0/(2.D0*sigma**2.D0))
    ! scaling to per interval width * total number density * log-normal distribution
  end do

end do
close(30)

if (sum(ni) /= 0.D0) then
  write(*,*) "not all zero"
end if


! Initialise ice crystal distribution
ni_crystal = 0.D0

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


end subroutine pbe_init

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

integer array_size_bytes
integer i

!----------------------------------------------------------------------------------------------

! Allocate arrays
allocate(v(0:m),dv(m),v_m(m),d_m(m),nuc(m),ni_droplet(m),ni_crystal(m),ni_soot(m))
allocate(vf(0:n_vf), vf_m(n_vf))
allocate(nislice_m(m), nislice_vf(n_vf))

array_size_bytes = 8*m*n_vf**3
if (array_size_bytes > 1000**3) then
  write(*,*) "Creating 2 arrays of size ",(8*m*n_vf**3/(1000**3))," GB."
else if (array_size_bytes > 1000**2) then
  write(*,*) "Creating 2 arrays of size ",(8*m*n_vf**3/(1000**2))," MB."
else
  write(*,*) "Creating 2 arrays of size ",(8*m*n_vf**3/(1000))," kB."
end if
allocate(ni(m,n_vf,n_vf,n_vf))
allocate(niprime(m,n_vf,n_vf,n_vf))

if (grid_type==1) then

  !Option 1: geometric grid
  if (grid_lb==0) then
    write(*,*) "Left grid boundary must be greater than zero for geometric grid"
    stop
  end if
  alpha = exp(1.D0/m * log(grid_rb/grid_lb)) ! Calculates geometric ratio
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

! Calculate volume fraction dimension
vf_width = 1./n_vf
vf(0) = 0.
do i=1,n_vf
  vf(i) = vf(i-1) + vf_width
end do

! Determine volume fraction mid-points
vf_m(1) = vf_width/2.
do i=2,n_vf
  vf_m(i) = vf_m(i-1) + vf_width
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

double precision dilution_factor, temperature_new

!----------------------------------------------------------------------------------------------

! Initialisation
if (current_time.eq.0.D0) then
  temperature = T_exhaust
  Pvap = Pvap_exhaust

! Updating
else

  if (current_time.LE.tau_m) then
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

Psat_l = 6.108D2*exp(17.27D0 * (temperature - 273.15D0)/(temperature - 35.86D0))
Psat_i = 6.108D2*exp(21.87D0 * (temperature - 273.15D0)/(temperature - 7.66D0))

supersaturation_l = Pvap/Psat_l - 1.D0
supersaturation_i = Pvap/Psat_i - 1.D0

! Update H2O density
comp_densities(4) = 1.D3

! H2O number concentration at water saturation - not totally sure this is correct
n_sat = avogadro_constant * Pvap / (ideal_gas_constant * temperature)

vapour_thermal_speed = sqrt(3 * boltzmann_constant * temperature / water_molecular_mass)

diff_coeff = 2.11D-5 * (temperature/273.15D0)**(1.94D0) * (101325D0 / P_ambient)


end subroutine pbe_set_environment

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

! Ice germ rate is number of ice germs formed per droplet vol per second
ice_germ_rate = 1.D6 * exp(-3.5714D0 * temperature + 858.719D0)

do index=1,m
  ! Condition to have >= 1 ice germ per droplet
  if ((ice_germ_rate * v_m(index) * dt).ge.(1.D0)) then
    ni_crystal(index) = ni_crystal(index) + vf_width**3 * sum(ni(index,:,:,:))
    ni(index,:,:,:) = 0.D0
  end if
end do

end subroutine pbe_freezing

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_activation()

!**********************************************************************************************
!
! Update (soot) activation each time step
!
! By Jack Bartlett (27/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision r_act, v_act, sat_ratio_l

integer index

!----------------------------------------------------------------------------------------------

if (Pvap>Psat_l) then

  ! Calculate minimum soot radius to have activated
  sat_ratio_l = Pvap/Psat_l
  r_act = r_k / (54.D0*soot_solubility*(log(sat_ratio_l))**2)**(1.D0/3.D0)

  ! Convert radius to volume
  v_act = 4.D0/3.D0 * pi * r_act**3.D0

  ! Move ni_soot to ni_droplet for all volumes larger than critical
  do index=1,m
    if (v_m(index) > v_act) then
      ni_droplet(index) = ni_droplet(index) + ni_soot(index)
      ni_soot(index) = 0.D0
    end if
  end do
end if

end subroutine pbe_activation

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_moments(ani,moment,meansize)

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

double precision, dimension(m), intent(in)    :: ani
double precision, dimension(0:1), intent(out) :: moment
double precision, intent(out)                 :: meansize

double precision M1_lp,lp

integer i

!----------------------------------------------------------------------------------------------

moment(0) = 0.0
moment(1) = 0.0

do i=1,m
  moment(0) = moment(0) + ani(i)*dv(i)
  moment(1) = moment(1) + 0.5D0*(v(i-1)+v(i))*ani(i)*dv(i)
end do

M1_lp = 0.D0
do i=m-5,m
  M1_lp = M1_lp + 0.5D0*(v(i-1)+v(i))*ani(i)*dv(i)
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

subroutine pbe_output_psd(ani,filename,current_time,first_write)

!**********************************************************************************************
!
! Writes PSD each time step
!
! By Jack Bartlett (29/05/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, dimension(m), intent(in) :: ani
character(len=30), intent(in) :: filename
double precision, intent(in) :: current_time
logical, intent(in) :: first_write

double precision :: nitemp(m)
double precision, dimension(0:1) :: moment

double precision meansize

integer i

!----------------------------------------------------------------------------------------------

call pbe_moments(nitemp,moment,meansize)

do i=1,m
  if (abs(ani(i))<1.D-16) then
    nitemp(i) = 0.D0
  else
    nitemp(i) = ani(i)
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

subroutine pbe_output_general_psd(current_time,first_write)

!**********************************************************************************************
!
! Writes general particle PSD
!
! By Jack Bartlett (10/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: current_time
logical, intent(in) :: first_write

double precision :: nitemp(m)
double precision, dimension(0:1) :: moment

double precision meansize

integer i

!----------------------------------------------------------------------------------------------

! Integrate over component dimensions
do i=1,m
  nitemp(i) = vf_width**3 * sum(ni(i,:,:,:))
end do

call pbe_moments(nitemp,moment,meansize)

do i=1,m
  if (abs(nitemp(i))<1.D-16) then
    nitemp(i) = 0.D0
  end if
end do


if (first_write) then
  open(99,file='output/psd_general.out',status='replace')
else
  open(99,file='output/psd_general.out',status='old',position='append')
end if
do i=1,m
  write(99,1001) current_time, v_m(i), dv(i), d_m(i), nitemp(i), &
  & nitemp(i)*dv(i)/moment(0), v_m(i)*nitemp(i), v_m(i)*nitemp(i)*dv(i)/moment(1)
end do
close(99)


1001 format(8E20.10)

end subroutine pbe_output_general_psd

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_output_crystal_psd(current_time,first_write)

!**********************************************************************************************
!
! Writes ice crystal PSD
!
! By Jack Bartlett (10/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

double precision, intent(in) :: current_time
logical, intent(in) :: first_write

double precision :: nitemp(m)
double precision, dimension(0:1) :: moment

double precision meansize

integer i

!----------------------------------------------------------------------------------------------

nitemp = ni_crystal

call pbe_moments(nitemp,moment,meansize)

do i=1,m
  if (abs(nitemp(i))<1.D-16) then
    nitemp(i) = 0.D0
  end if
end do


if (first_write) then
  open(99,file='output/psd_crystal.out',status='replace')
else
  open(99,file='output/psd_crystal.out',status='old',position='append')
end if
do i=1,m
  write(99,1001) current_time, v_m(i), dv(i), d_m(i), nitemp(i), &
  & nitemp(i)*dv(i)/moment(0), v_m(i)*nitemp(i), v_m(i)*nitemp(i)*dv(i)/moment(1)
end do
close(99)


1001 format(8E20.10)

end subroutine pbe_output_crystal_psd

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_output_env(current_time,first_write)

!**********************************************************************************************
!
! Writes environment variables each time step
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
 
deallocate(v,dv,v_m,d_m,nuc,ni_droplet,ni_crystal,ni_soot)
deallocate(vf, vf_m)
deallocate(ni)
deallocate(niprime)
deallocate(nislice_m, nislice_vf)

end subroutine pbe_deallocate

!**********************************************************************************************