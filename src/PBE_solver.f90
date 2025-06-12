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

double precision sum_VG ! sum of volume times growth source
double precision r_m, J, sum_Jn, delta_supersaturation_l, delta_supersaturation_i ! now only used for crystals

integer index

!----------------------------------------------------------------------------------------------

! GENERAL PARTICLES

!Euler explicit
call pbe_ydot_general(dt)
ni = ni + niprime * dt ! niprime is now stored in pbe_mod

! Cap at zero after growth
where (ni < 0.D0)
  ni = 0.D0
end where

sum_VG = 0.

! Adjust water vapour due to droplet growth
Pvap = Pvap + ((sum_VG / water_molecular_vol) * boltzmann_constant * temperature) * dt

if (Pvap<0.D0) then
  Pvap = 0.D0 ! just in case
end if

!----------------------------------------------------------------------------------------------

! CRYSTALS

!Euler explicit
!call pbe_ydot_crystal(niprime_crystal,dt)
!ni_crystal = ni_crystal + niprime_crystal * dt

! Cap at zero after growth
!do index = 1,m
!  if (ni_crystal(index) < 0.D0) then
!    ni_crystal(index) = 0.D0
!  end if
!end do


! Deplete supersaturation due to crystal growth
!sum_Jn = 0.D0
!if (supersaturation_i>0.D0) then
!  
!  do index=1,m
!    ! Calculate H2O flux to particles of size r_m
!    r_m = ((3.D0*v_m(index))/(4.D0*pi))**(1.D0/3.D0)
!    call calc_J(r_m, supersaturation_i, J)
!    sum_Jn = sum_Jn + J * ni_crystal(index) * dv(index) ! Last part to change to absolute number density
!  end do
!
!  delta_supersaturation_i = (- sum_Jn / n_sat) * dt ! in brackets is ds/dt
!  Pvap = Pvap + delta_supersaturation_i * Psat_i
!  if (Pvap<0.D0) then
!    Pvap = 0.D0 ! just in case
!  end if
!end if

end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_general(dt)

!**********************************************************************************************
!
! Calculates the right hand side of the ODEs to be integrated
!
! By Stelios Rigopoulos
! Modified by Jack Bartlett (10/06/2025)
!
!**********************************************************************************************

use pbe_mod

implicit none

! Interface to get nislice array of assumed dimension to work
interface
  subroutine growth_tvd_general(nislice, i, i_left, rb_index, lb_index, max_index, interval_width, dt, growth_source)
    double precision, dimension(:), intent(in) :: nislice
    integer, dimension(4), intent(in)          :: i
    integer, dimension(4), intent(in)          :: i_left
    integer, intent(in)                        :: rb_index
    integer, intent(in)                        :: lb_index
    integer, intent(in)                        :: max_index
    double precision, intent(in)               :: interval_width
    double precision, intent(in)               :: dt
    double precision, intent(out)              :: growth_source
  end subroutine growth_tvd_general
end interface

double precision, intent(in) :: dt

double precision growth_source,growth_mass_source,params(1),growth_rate,growth_source_tot
double precision new_V,vf1,vf2,vf3
double precision interval_width

integer, dimension(4) :: i, i_left
integer i1,i2,i3,i4,i3_max,i4_max,rb_index,lb_index,max_index

!----------------------------------------------------------------------------------------------

niprime = 0. ! d(ni)/dt
params(1) = 0.


! Nucleation
! Add homogeneous nucleation at high supersaturation? Otherwise none


! Droplet growth

do i1 = 1,m
  do i2 = 1,n_vf
    i3_max = n_vf + 1 - i2
    do i3 = 1,i3_max
      i4_max = n_vf + 1 - i2 - i3
      do i4 = 1,i4_max
        i = (/i1,i2,i3,i4/)

        growth_source_tot = 0.D0

        ! Volume growth source
        nislice_m = ni(:, i(2), i(3), i(4))
        rb_index = i(1)
        lb_index = i(1) - 1
        max_index = m
        i_left = (/lb_index, i(2), i(3), i(4)/)
        interval_width = dv(i(1))
        call growth_tvd_general(nislice_m, i, i_left, rb_index, lb_index, max_index, interval_width, dt, growth_source)
        growth_source_tot = growth_source_tot + growth_source

        ! Component 1 growth source
        nislice_vf = ni(i(1), :, i(3), i(4))
        rb_index = i(2)
        lb_index = i(2) - 1
        max_index = n_vf
        i_left = (/i(1), lb_index, i(3), i(4)/)
        interval_width = vf_width
        call growth_tvd_general(nislice_vf, i, i_left, rb_index, lb_index, max_index, interval_width, dt, growth_source)
        growth_source_tot = growth_source_tot + growth_source

        ! Component 2 growth source
        nislice_vf = ni(i(1), i(2), :, i(4))
        rb_index = i(3)
        lb_index = i(3) - 1
        max_index = n_vf
        i_left = (/i(1), i(2), lb_index, i(4)/)
        interval_width = vf_width
        call growth_tvd_general(nislice_vf, i, i_left, rb_index, lb_index, max_index, interval_width, dt, growth_source)
        growth_source_tot = growth_source_tot + growth_source

        ! Component 3 growth source
        nislice_vf = ni(i(1), i(2), i(3), :)
        rb_index = i(4)
        lb_index = i(4) - 1
        max_index = n_vf
        i_left = (/i(1), i(2), i(3), lb_index/)
        interval_width = vf_width
        call growth_tvd_general(nislice_vf, i, i_left, rb_index, lb_index, max_index, interval_width, dt, growth_source)
        growth_source_tot = growth_source_tot + growth_source
        
        niprime(i(1), i(2), i(3), i(4)) = niprime(i(1), i(2), i(3), i(4)) + growth_source_tot



        !call calc_growth_rate_liquid(i, growth_rate)
        !new_V = v_m(i(1)) + (growth_rate * dt)
        !vf1 = vf_m(i(2)) * (v_m(i(1)) / new_V) ! nvPM
        !vf2 = vf_m(i(3)) * (v_m(i(1)) / new_V) ! H2SO4
        !vf3 = vf_m(i(4)) * (v_m(i(1)) / new_V) ! organics

        ! Determine new intervals
        !do index=1,m
        !  if ((new_V.ge.v(index-1)).and.(new_V.le.v(index))) then
        !    new_i(1) = index
        !  end if
        !end do
        !do index=1,n_vf
        !  if ((vf1.ge.vf(index-1)).and.(vf1.le.vf(index))) then
        !    new_i(2) = index
        !  end if
        !  if ((vf2.ge.vf(index-1)).and.(vf2.le.vf(index))) then
        !    new_i(3) = index
        !  end if
        !  if ((vf3.ge.vf(index-1)).and.(vf3.le.vf(index))) then
        !    new_i(4) = index
        !  end if
        !end do

        ! Take everything in current indices and move to new indices
        !if ((new_i(1) /= i(1)).or.(new_i(2) /= i(2)).or.(new_i(3) /= i(3)).or.(new_i(4) /= i(4))) then
        !  sum_VG = sum_VG + growth_rate * ni(i(1), i(2), i(3), i(4))
        !  ni(new_i(1), new_i(2), new_i(3), new_i(4)) = &
        !  & ni(new_i(1), new_i(2), new_i(3), new_i(4)) + ni(i(1), i(2), i(3), i(4))
        !  ni(i(1), i(2), i(3), i(4)) = 0.D0
        !end if
      end do
    end do
  end do
end do


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

end subroutine pbe_ydot_general

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot_crystal(niprime_crystal,dt)

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

double precision, dimension(m), intent(out) :: niprime_crystal

double precision, intent(in) :: dt

double precision growth_source,growth_mass_source,params(1)

integer index

!----------------------------------------------------------------------------------------------

niprime_crystal = 0. ! d(ni)/dt
params(1) = 0.


! Nucleation
! Is there any homogeneous ice nucleation?


! Crystal growth
if (supersaturation_i>0) then

  do index = 1,m
    ! growth_tvd has been modified to only work with crystals to allow compile
    call growth_tvd(ni_crystal, index, supersaturation_i, dt, growth_source)
    niprime_crystal(index) = niprime_crystal(index) + growth_source
  end do

  ! Else niprime_crystal(index) does not change
end if


!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ni_crystal,niprime_crystal)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(niprime_crystal,ni_crystal)
end if

end subroutine pbe_ydot_crystal

!**********************************************************************************************