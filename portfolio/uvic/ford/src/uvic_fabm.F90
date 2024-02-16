! Copyright (C) 2024 Bolding & Bruggeman

!> Some native UVic_ESCM variables will have to be converted to be compatible with FABM.
!> This can either be because of different dimensionality or different units.
!> This is done by creating module level private variables that will be calculated/updated
!> based on the original UVic_ECSM variables. Some will only need to be calculated once - like
!> layer heights - and some will have to be updated every time step - like density.
!> @warning
!> This module is still under development.
!> API and functioning might change without notice.
!> @endwarning
!>
!> @history
!> A list of important UVic (MOM2) variables used by FABM:  
!>
!> - \(t(imt,km,jmt,nt,-1:1)\) - all tracers [source/mom/mw.h]
!> - \(src(imt,km,jsmw:jemw,nsrc)\): tracers with sources [source/mom/tracer.f]
!> - \(sbc(imt,jmt,numsbc)\) surftace boundary conditions [source/common/csbc.h]
!> - \(rho(imt,km,jsmw:jmw)\):  density [source/mom/mw.h]
!>
!> with
!>
!> - \((imt,km,jmt) = (102,19,102)\) [commom/size.h]
!> - \(nt = 2+??\) number of tracers [commom/size.h]
!> - \(jsmw:jemw = (2, jmw=jmt) or (2, jmw=(3,4,5)\) [commom/size.h]
!> - \(nsrc\): number of tracers with source terms [common/csbc.h]
!> - \(numsbc\): total number of surface boundary conditions - 
!>   list in [common/csbc.h] set in [common/UVic_ESCM.F].
!>
!> Updating FABM is done i the tracer subroutine
!>
!> \(call tracer (joff, jstrac, jetrac, is, ie)\)

!> with
!>
!> - \((joff, jstrac, jetrac, is, ie) = (0,2,101,2,101)\)
!>
!> @endhistory
!>
!> @note
!> The FABM calculation domain in UVic reference is \( (2:imt-1,km,2:jmt-1) \)
!>
!> Is it best to allocate arrays correspondingly or do the mapping 
!> in the do-loops?
!> @endnote
!>
!> @note
!> Dimension of z, dz and pressure in a z-coordinate model?  
!> Option to be 1D?
!>
!> @endnote

module uvic_fabm

#ifdef O_fabm

use uvic_common_blocks
use fabm
use fabm_types
#define DEBUG
#ifdef DEBUG
use fabm_debug
#endif

IMPLICIT NONE

private

class (type_fabm_model), pointer :: model
   !! This variable will contain all FABM configuration and
   !! give access to FABM routines

#if 1
! The following is a list of all FABM standard variables - the definitions can be updated
! during the implementation. The order of the variables is mainained from:
! https://github.com/fabm-model/fabm/wiki/List-of-standard-variables
! For mandatory variables - e.g. cell thinkness it is not necessary to obtain an id.

! Interior variables
!! id_alkalinity_expressed_as_mole_equivalent
!! id_attenuation_coefficient_of_photosynthetic_radiative_flux
!! id_attenuation_coefficient_of_shortwave_flux
!type (type_fabm_interior_variable_id) :: id_cell_thickness
real(rke), allocatable, target :: dz(:,:,:)
type (type_fabm_interior_variable_id) :: id_density
real(rke), allocatable, target :: rho_fabm(:,:,:)
!type (type_fabm_interior_variable_id) :: id_depth
real(rke), allocatable, target :: depth(:,:,:)
!! id_downwelling_photosynthetic_radiative_flux
!! id_downwelling_shortwave_flux
!! id_fractional_saturation_of_oxygen
!! id_mass_concentration_of_suspended_matter
!! id_mole_concentration_of_ammonium
!! id_mole_concentration_of_carbonate_expressed_as_carbon
!! id_mole_concentration_of_dissolved_inorganic_carbon
!! id_mole_concentration_of_dissolved_iron
!! id_mole_concentration_of_nitrate
!! id_mole_concentration_of_phosphate
!! id_mole_concentration_of_silicate
!! id_net_rate_of_absorption_of_shortwave_energy_in_layer
!! id_ph_reported_on_total_scale
type (type_fabm_interior_variable_id) :: id_practical_salinity
real(rke), allocatable, target :: salt(:,:,:)
type (type_fabm_interior_variable_id) :: id_pressure
real(rke), allocatable, target :: pressure(:,:,:)
!! id_secchi_depth
!type (type_fabm_interior_variable_id) :: id_temperature

! Surface variables
!! id_cloud_area_fraction
!! id_ice_area_fraction
type (type_fabm_horizontal_variable_id) :: id_mole_fraction_of_carbon_dioxide_in_air
real(rke), allocatable, target :: mole_fraction_of_carbon_dioxide_in_air(:,:)
!! id_surface_air_pressure
!! id_surface_albedo
!! id_surface_downwelling_photosynthetic_radiative_flux
!! id_surface_downwelling_photosynthetic_radiative_flux_in_air
type (type_fabm_horizontal_variable_id) :: id_surface_downwelling_shortwave_flux
real(rke), allocatable, target :: surface_downwelling_shortwave_flux(:,:)
!! id_surface_downwelling_shortwave_flux_in_air
!! id_surface_drag_coefficient_in_air
!! id_surface_specific_humidity
!! id_surface_temperature
type (type_fabm_horizontal_variable_id) :: id_windspeed
real(rke), allocatable, target :: windspeed(:,:)

! Bottom variables
!! id_bottom_depth
!! id_bottom_depth_below_geoid
!! id_bottom_roughness_length
type (type_fabm_horizontal_variable_id) :: id_bottom_stress
real(rke), allocatable, target :: bottom_stress(:,:)

! Global variables
!! id_number_of_days_since_start_of_the_year

! Universal variables
!! id_total_carbon
!! id_total_iron
!! id_total_nitrogen
!! id_total_phosphorus
!! id_total_silicate
! horizontal FABM ids
! interior FABM variables - calculated from UVic_ESCM variables
! horizontal FABM variables
#endif

! public available routines
public fabm_configure, fabm_link_data, fabm_update, fabm_list, fabm_clean

! module level variables - static or allocatable? Same goes with e.g. windspeed and rho_fabm


#define _DOMAIN_  2:imt-1,km,2:jmt-1
#define _I_  2:imt-1
#define _J_  2:jmt-1
#define _K_  km

real(rke) :: surface_flux(_I_,_J_,nt)
  !! surface fluxes
real(rke) :: surface_sms(_I_,_J_,nt)
  !! surface source-sink terms
real(rke) :: bottom_flux(_I_,_J_,nt)
  !! bottom fluxes
real(rke) :: bottom_sms(_I_,_J_,nt)
  !! bottom source-sink terms
!real(rke) :: w(imt,km,jsmw:jemw,nt)
real(rke) :: w(_I_,_K_,_J_,nt)
  !! vertical velocity in m/s

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine fabm_configure(dt)
   !! read fabm.yaml and call FABM configuration subroutines

   real(rke), intent(in) :: dt
      !! bio-geochemical time step as set by MOM2 [s]

   print*, '== Initializing FABM component with nt=',nt

   model => fabm_create_model('fabm.yaml')

   !parameter (jsmw=2, jemw=jmw-1) - parameter (jmw=jmt)
   ! joff,js,je,is,ie 0, 2, 101, 2, 101 - fabm_update()
   print*, imt,jmt,jmw
   print*, 'zt: ',shape(zt)
   print*, 't: ',shape(t)
   print*, 'sbc: ',shape(sbc)
   print*, 'src: ',shape(src) ! imt,km,jsmw:jemw,nsrc
   print*, 'source: ',shape(source) ! imt,km,jsmw:jemw
   print*, 'rho: ',shape(rho)
   !stop 112

   if (nt-2 .ne. size(model%interior_state_variables)) then
      print*, nt-2,size(model%interior_state_variables)
!KB      stop 'aa'
   end if

!   call model%set_domain(imt,km,jmt,dt)
   call model%set_domain(imt-2,km,jmt-2,dt) ! jmt or jmw?
   call model%set_domain_start(2,1,2)
   call model%set_domain_stop(imt-1,km,jmt-1)

!   call model%set_mask(tmask(2:jmt-1,:,2:jmt-1),tmask(2:jmt-1,1,2:jmt-1)) 
   call model%set_mask(tmask(_I_,:,_J_),tmask(_I_,1,_J_)) 
   !> @note
   !! seems tmask is not initialised until called in mom() - 
   !! i.e. after initialization - so all values are 0 here
   !! @endnote
   !KBprint*, tmask(53,:,53)

   call model%set_bottom_index(kmt(_I_,_J_))
end subroutine fabm_configure

!-----------------------------------------------------------------------

subroutine fabm_link_data()
   !! link all FABM configured external dependencies - and call
   !! model%start() to assure proper configuration
   integer :: n

   ! link to time in-dependent data that do require transformation
   call link_grid()

   ! link to time dependent data that do NOT require transformation
   call model%link_interior_data(fabm_standard_variables%temperature,t(_I_,:,_J_,itemp,0))

   ! link to time dependent data that do require transformation
   ! initialize and update time changing environmental variables
#if 1
   call link_wind()
   call link_mole_fraction_of_carbon_dioxide_in_air()
   call link_surface_downwelling_shortwave_flux()
   call link_bottom_stress()
   call link_salinity()
   call link_density()
#endif

   ! link to FABM's interior state variables
   do n = 1, size(model%interior_state_variables)
      call model%link_interior_state_data(n, t(_I_,:,_J_,2+n,0))
      mapt(2+n) = trim(model%interior_state_variables(n)%name)
      !KBmapst(2+n) = 's'//trim(mapt(2+n))
   end do

   call model%start()

end subroutine fabm_link_data

!-----------------------------------------------------------------------

subroutine fabm_list()
   !! lists all FABM configured variables
   integer :: n

   print*, 'FABM interior state variables:'
   do n = 1,size(model%interior_state_variables)
      print*, n, &
             trim(model%interior_state_variables(n)%name), '  ', &
             trim(model%interior_state_variables(n)%units),'  ',&
             trim(model%interior_state_variables(n)%long_name)
   end do

   print*, 'FABM surface-bound state variables:'
   do n=1,size(model%surface_state_variables)
      print*, n, &
             trim(model%surface_state_variables(n)%name), '  ', &
             trim(model%surface_state_variables(n)%units),'  ', &
             trim(model%surface_state_variables(n)%long_name)
   end do

   print*, 'FABM bottom-bound state variables:'
   do n=1,size(model%bottom_state_variables)
      print*, n, &
             trim(model%bottom_state_variables(n)%name), '  ', &
             trim(model%bottom_state_variables(n)%units),'  ', &
             trim(model%bottom_state_variables(n)%long_name)
   end do

#if 0
   print*, 'FABM diagnostic variables defined on the full model domain:'
   do n=1,size(model%interior_diagnostic_variables)
      print*, n, &
             trim(model%interior_diagnostic_variables(n)%name), '  ', &
             trim(model%interior_diagnostic_variables(n)%units),'  ', &
             trim(model%interior_diagnostic_variables(n)%long_name)
   end do

   print*, 'FABM diagnostic variables defined on a horizontal slice of the model domain:'
   do n=1,size(model%horizontal_diagnostic_variables)
      print*, n, &
             trim(model%horizontal_diagnostic_variables(n)%name), '  ', &
             trim(model%horizontal_diagnostic_variables(n)%units),'  ', &
             trim(model%horizontal_diagnostic_variables(n)%long_name)
   end do
#endif
end subroutine fabm_list

!-----------------------------------------------------------------------

subroutine fabm_update(joff, js, je, is, ie)
   !! update the environment and calculate the source/sink terms
   !! call with the same argument list as mom() calls tracer() i.e. 
   !! the specification on the active UVic window - typically the 
   !! full domain on modern hardware

   integer, intent(in) :: joff
     !! offset row in global window
   integer, intent(in) :: js
     !! start row
   integer, intent(in) :: je
     !! end row
   integer, intent(in) :: is
     !! start column
   integer, intent(in) :: ie
     !! end column

   integer :: i,j,k,n
     ! local loop counters

   print*, 'fabm_update',joff,js,je,is,ie
   ! t(:,:,:,var,0) is updated in loadmw() in mom()
   call update_data(joff)
#if 0
   print*, 'wind ',sbc(53,53,iws)/100._rke
   print*, 'salt ',35.+ 1000*t(53,:,53,isalt,0)
   stop 'egon'
#endif

!KB   call model%prepare_inputs(t=real(imal,rk))
   call model%prepare_inputs()
   stop 111

   ! here the surface is updated
   surface_flux = 0._rke
   surface_sms = 0._rke
#if 0
   do j=js,je
      call model%get_surface_sources(is,ie,j, &
                 surface_flux(is:ie,j,:),surface_sms(is:ie,j,:))
   end do
#endif

   ! here the pelagic is updated
   do j=js,je
      do k=1,km
         call model%get_interior_sources(is,ie,k,j,src(is:ie,k,j,:))
      end do
   end do

   ! here the bottom is updated
   bottom_flux = 0._rke
   bottom_sms = 0._rke
#if 0
   do j=js,je
      call model%get_bottom_sources(is,ie,j, &
                 bottom_flux(is:ie,j,:),bottom_sms(is:ie,j,:))
   end do
#endif

   print*, shape(bottom_flux)
   print*, shape(bottom_sms)
   print*, shape(src)
   stop 110

#if 1
   ! fold the surface and bottom flux terms - src keeps track on 
   ! which variables actually have sources - itrc(n) and the size
   ! of src reflects this - tracer.F90 line 1122
   do n=1,nt
   do j=js,je
      do i=is,ie
         if (kmt(i,j) > 0) then
            ! surface
            k=1
            src(i,k,j,:)=src(i,k,j,n)+surface_flux(i,j,n)/dz(i,k,j)
            ! bottom
            k=kmt(i,j)
            src(i,k,j,:)=src(i,k,j,n)+bottom_flux(i,j,n)/dz(i,k,j)
         end if
      end do
   end do
   end do
#endif

   ! vertical velocities
   do j=js,je
      do k=1,km
!         call model%get_vertical_movement(is,ie,k,j,w)
      end do
   end do

   call model%finalize_outputs()

   !stop 120

end subroutine fabm_update

!-----------------------------------------------------------------------

subroutine fabm_clean()
   !! de-allocate all allocated arrays
   if (allocated(windspeed)) deallocate(windspeed)
   !if (allocated(windspeed)) deallocate(windspeed)
   !if (allocated(windspeed)) deallocate(windspeed)
   if (allocated(salt)) deallocate(salt)
   if (allocated(rho_fabm)) deallocate(rho_fabm)
   ! :
   ! :
end subroutine fabm_clean

!-----------------------------------------------------------------------

subroutine update_data(joff)
   integer, intent(in) :: joff
     !! offset row in global window
   !! update all time varying FABM configured external dependencies
   !! by calling individual update routines - tests done in routines
   call update_wind()
   call update_mole_fraction_of_carbon_dioxide_in_air()
   call update_surface_downwelling_shortwave_flux()
   call update_bottom_stress(joff)
   call update_salinity()
   call update_density()
end subroutine update_data

!-----------------------------------------------------------------------

subroutine link_grid()
   !! Allocate and link grid related FABM standard variables that 
   !! are being transformed from UVic native variables [cm -> m].

   integer :: rc
   integer :: i,j,k

#if 1
   allocate(depth(_I_,_K_,_J_),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (depth)'
   depth = 0._rke
   allocate(pressure(_I_,_K_,_J_),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (pressure)'
   pressure = 0._rke
   allocate(dz(_I_,_K_,_J_),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (dz)'
   dz = 0._rke
#if 0
   do j=1,jmt
      do i=1,imt
#else
   do j=2,jmt-1
      do i=2,imt-1
#endif
         if (kmt(i,j) > 0) then
            depth(i,:,j) = zt/100._rke
            pressure(i,:,j) = depth(i,:,j)/10._rke
            dz(i,:,j) = dzt/100._rke
         end if
      end do
   end do
#else
   allocate(depth(imt,km,jmt),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (depth)'
   depth = 0._rke
   allocate(pressure(imt,km,jmt),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (pressure)'
   pressure = 0._rke
   allocate(dz(imt,km,jmt),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (dz)'
   dz = 0._rke
   do j=1,jmt
      do i=1,imt
         if (kmt(i,j) > 0) then
            depth(i,:,j) = zt/100._rke
            pressure(i,:,j) = depth(i,:,j)/10._rke
            dz(i,:,j) = dzt/100._rke
         end if
      end do
   end do
#endif
   call model%link_interior_data(fabm_standard_variables%depth,depth)
   call model%link_interior_data(fabm_standard_variables%pressure,pressure)
   call model%link_interior_data(fabm_standard_variables%cell_thickness,dz)
#if 0
   print*, depth(53,:,53),pressure(:53,:,53),dz(53,:,53)
   stop 'egon'
#endif
end subroutine link_grid

!-----------------------------------------------------------------------

subroutine link_wind()
   !! get wind speed FABM standard variable and if needed by FABM 
   !! allocate memory
   integer rc
   id_windspeed = model%get_horizontal_variable_id(standard_variables%wind_speed)
   if (model%variable_needs_values(id_windspeed)) then
      allocate(windspeed(_I_,_J_),stat=rc)
      if (rc /= 0) stop 'link_wind(): Error allocating (windspeed)'
      windspeed = 0._rke
      call model%link_horizontal_data(id_windspeed,windspeed)
   end if
end subroutine link_wind

!-----------------------------------------------------------------------

subroutine update_wind()
   !! calculate wind speed in m/s according to 
   !! $$w = w_{UVic}/100$$
   ! calculate the wind speed in m/s - iws, iaws

   integer i,j
   if (model%variable_needs_values(id_windspeed)) then
#if 0
   do j=1,jmt
      do i=1,imt
#else
   do j=2,jmt-1
      do i=2,imt-1
#endif
            if (kmt(i,j) > 0) windspeed(i,j) = sbc(i,j,iws)/100._rke
         end do
      end do
   end if
end subroutine update_wind

!-----------------------------------------------------------------------

subroutine link_mole_fraction_of_carbon_dioxide_in_air()
   integer i,j,rc
   id_mole_fraction_of_carbon_dioxide_in_air = model% &
       get_horizontal_variable_id(standard_variables% &
       mole_fraction_of_carbon_dioxide_in_air)
   if (model%variable_needs_values(id_mole_fraction_of_carbon_dioxide_in_air)) then
      allocate(mole_fraction_of_carbon_dioxide_in_air(_I_,_J_),stat=rc)
      if (rc /= 0) stop 'link_mole_fraction_of_carbon_dioxide_in_air(): Error allocating (mole_fraction_of_carbon_dioxide_in_air)'
      mole_fraction_of_carbon_dioxide_in_air = 0._rke
      call model%link_horizontal_data(id_mole_fraction_of_carbon_dioxide_in_air,mole_fraction_of_carbon_dioxide_in_air)
   end if
end subroutine link_mole_fraction_of_carbon_dioxide_in_air

!-----------------------------------------------------------------------

subroutine update_mole_fraction_of_carbon_dioxide_in_air()
   !! calculate the ?????? in W/m^2
   integer i,j
   if (model%variable_needs_values(id_mole_fraction_of_carbon_dioxide_in_air)) then
#if 0
      do j=1,jmt
         do i=1,imt
#else
      do j=2,jmt-1
         do i=2,imt-1
#endif
            if (kmt(i,j) > 0) mole_fraction_of_carbon_dioxide_in_air(i,j) =  10._rke !KB
         end do
      end do
   end if
end subroutine update_mole_fraction_of_carbon_dioxide_in_air

!-----------------------------------------------------------------------

subroutine link_surface_downwelling_shortwave_flux()
   integer rc
   id_surface_downwelling_shortwave_flux = model% &
        get_horizontal_variable_id(standard_variables% &
        surface_downwelling_shortwave_flux)
   if (model%variable_needs_values(id_surface_downwelling_shortwave_flux)) then
      allocate(surface_downwelling_shortwave_flux(_I_,_J_),stat=rc)
      if (rc /= 0) stop 'link_surface_downwelling_shortwave_flux(): Error allocating (surface_downwelling_shortwave_flux)'
      surface_downwelling_shortwave_flux = 0._rke
      call model%link_horizontal_data(id_surface_downwelling_shortwave_flux,surface_downwelling_shortwave_flux)
   end if
end subroutine link_surface_downwelling_shortwave_flux

!-----------------------------------------------------------------------

subroutine update_surface_downwelling_shortwave_flux()
   !! calculate the short wave flux in W/m^2

   integer i,j
   if (model%variable_needs_values(id_surface_downwelling_shortwave_flux)) then
#if 0
      do j=1,jmt
         do i=1,imt
#else
      do j=2,jmt-1
         do i=2,imt-1
#endif
            if (kmt(i,j) > 0) surface_downwelling_shortwave_flux(i,j) = 200._rke !KB
         end do
      end do
   end if
end subroutine update_surface_downwelling_shortwave_flux

!-----------------------------------------------------------------------

subroutine link_bottom_stress()
   !! get bottom stress FABM standard variable and if needed by FABM 
   !! allocate memory
   integer rc
   id_bottom_stress = model%get_horizontal_variable_id(standard_variables%bottom_stress)
   if (model%variable_needs_values(id_bottom_stress)) then
      allocate(bottom_stress(_I_,_J_),stat=rc)
      if (rc /= 0) stop 'link_bottom_stress(): Error allocating (bottom_stress)'
      bottom_stress = 0._rke
      call model%link_horizontal_data(id_bottom_stress,bottom_stress)
   end if
end subroutine link_bottom_stress

!-----------------------------------------------------------------------

subroutine update_bottom_stress(joff)
   !! calculate the bottom stress in Pa
   integer, intent(in) :: joff
     !! offset row in global window
   real(rke), parameter :: x=10._rke ! dynes/cm2 --> Pa
   integer i,j,jrow
   if (model%variable_needs_values(id_bottom_stress)) then
#if 0
      do j=1,jmt
         jrow = j+joff
         do i=1,imt
#else
      do j=2,jmt-1
         jrow = j+joff
         do i=2,imt-1
#endif
            if (kmt(i,jrow) > 0) bottom_stress(i,j) = x*sqrt( &
                bmf(i,jrow,1)**2 + bmf(i,jrow,2)**2)
!            if (kmt(i,j) > 0) bottom_stress(i,j) =  0.001_rke !KB
         end do
      end do
   end if
end subroutine update_bottom_stress

!-----------------------------------------------------------------------

subroutine link_salinity()
   !! get salinity FABM standard variable and if needed by FABM allocate
   !! memory
   integer rc
   id_practical_salinity = model%get_interior_variable_id(fabm_standard_variables%practical_salinity)
   if (model%variable_needs_values(id_practical_salinity)) then
      allocate(salt(_I_,_K_,_J_),stat=rc)
      if (rc /= 0) stop 'link_salinity(): Error allocating (salt)'
      salt = 0._rke
      call model%link_interior_data(id_practical_salinity,salt)
   end if
end subroutine link_salinity

!-----------------------------------------------------------------------

subroutine update_salinity()
   !! calculate salinity in PSU according to $$S = 35 + 1000*S_{UVic}$$
   integer i,j,k
   if (model%variable_needs_values(id_practical_salinity)) then
#if 0
      do j=1,jmt
         do k=1,km
            do i=1,imt
#else
      do j=2,jmt-1
         do k=1,km
            do i=2,imt-1
#endif
               if (kmt(i,j) > 0) salt(i,k,j) = 35._rk+1000._rke*t(i,k,j,isalt,0)
            end do
         end do
      end do
   end if
end subroutine update_salinity

!-----------------------------------------------------------------------

subroutine link_density()
   !! get density FABM standard variable and if needed by FABM allocate
   !! memory
   integer rc
   id_density = model%get_interior_variable_id(fabm_standard_variables%density)
   if (model%variable_needs_values(id_density)) then
      allocate(rho_fabm(_I_,_K_,_J_),stat=rc)
      if (rc /= 0) stop 'link_density(): Error allocating (rho_fabm)'
      rho_fabm = 0._rke
      call model%link_interior_data(id_density,rho_fabm)
   end if
end subroutine link_density

!-----------------------------------------------------------------------

subroutine update_density()
   !! calculate density in kg/m³ according to 
   !! $$\rho = 1000*(\rho_0 + \rho_{UVic})$$
   !! with \(\rho_0 = 1.035\). MUST match rho0 from UVic_ESCM.F90

   !! @note
   !! loadmw.F: l 154
   !!
   !! imt=102, km=19, jsmw=2, jmw=jmt --- jemw=jmw-1  
   !!
   !! declared: rho(imt,km,jsmw:jmw) calculated rho(1:102,2:102) rho
   !! @endnote

   !! @note
   !! KB check depth dependent reference density
   !! @endnote

   integer i,j,k
   real(rke), parameter :: rho0=1.035
   if (model%variable_needs_values(id_density)) then
      !do j=1,jmt
      do j=jsmw,jmw ! must be 2:102
         do k=1,km
#if 0
            do i=1,imt
#else
            do i=2,imt-1
#endif
               if (kmt(i,j) > 0) rho_fabm(i,k,j) = 1000._rke*(rho0+rho(i,k,j))
            end do
         end do
      end do
   end if
#if 0
   print*, rho_fabm(53,:,53)
   stop 'kurt'
#endif
end subroutine update_density

!-----------------------------------------------------------------------

#endif

end module uvic_fabm

!-----------------------------------------------------------------------
!MVV: associate( MNV => model%variable_needs_values )
!MHG: associate( MHGV => model%get_horizontal_variable_id )
!MIG: associate( MIGV => model%get_horizontal_variable_id )

