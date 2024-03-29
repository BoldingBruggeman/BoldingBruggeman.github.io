!> Copyright (C) 2024 Bolding & Bruggeman
!>
!> @warning
!> This module is still under development.
!> API and functioning might change without notice.
!> @endwarning
!>
!> @history
!> A list of important UVic (MOM2) variables used by FABM:  
!>
!> - \(t(imt,km,jmt,nt,-1:1)\) - all tracers [source/mom/mw.h] -
!> [T,S] = [1,2]
!> - \(src(imt,km,jsmw:jemw,nsrc)\): tracers with sources [source/mom/tracer.f]
!> - \(sbc(imt,jmt,numsbc)\) surface boundary conditions [source/common/csbc.h]
!> - \(rho(imt,km,jsmw:jmw)\):  density [source/mom/mw.h]
!>
!> with
!>
!> - \((imt,km,jmt) = (102,19,102)\) [commom/size.h]
!> - \(nt = 2+??\) number of tracers [commom/size.h]
!> - \(jsmw:jemw = (2, jmw=jmt) or (2, jmw=(3,4,5)\) [commom/size.h]
!> - \(nsrc\): number of tracers with source terms [common/size.h]
!> - \(numsbc\): total number of surface boundary conditions - 
!>   list in [common/csbc.h] set in [common/UVic_ESCM.F].
!>
!> Updating FABM is done in the tracer() subroutine called like:
!>
!> - \(call\ tracer (joff, jstrac, jetrac, is, ie)\)
!> with
!> \((joff, jstrac, jetrac, is, ie) = (0,2,101,2,101)\)
!>
!> The arguments to tracer() are passed directly to fabm_update(). Note
!> that potentially jstrac:jetrac does not cover the entire domain.  
!> Focus - initially - will be on getting it to work with !O_min_window
!> i.e. the entire domain calculated in one go.
!>
!> @endhistory
!>
!> @note
!> The FABM calculation domain in UVic reference is 
!> \( t(2:imt-1,km,2:jmt-1,3:nt) \)
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
!>
!> @note
!> What are the proper links to FABM modules?  
!> @endnote
!>
!> @note
!> Is there a way to check if a FABM pelagic variable has source terms?
!>
!> With UVic model all tracer variables - except T, S and CFC gases.
!> All respecting #ifdef's.
!> @endnote
!>
!> Some native UVic_ESCM variables will have to be converted to be compatible with FABM.
!> This can either be because of different dimensionality or different units.
!> This is done by creating module level private variables that will be calculated/updated
!> based on the original UVic_ECSM variables. Some will only need to be calculated once - like
!> layer heights - and some will have to be updated every time step - like density.

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

! The following is a list of all FABM standard variables - the definitions can be updated
! during the implementation. The order of the variables is mainained from:
! https://github.com/fabm-model/fabm/wiki/List-of-standard-variables
! For mandatory variables - e.g. cell thinkness it is not necessary to obtain an id.

! Interior variables

! id_alkalinity_expressed_as_mole_equivalent
! id_attenuation_coefficient_of_photosynthetic_radiative_flux
! id_attenuation_coefficient_of_shortwave_flux
!type (type_fabm_interior_variable_id) :: id_cell_thickness
real(rke), allocatable, target :: dz(:,:,:)
type (type_fabm_interior_variable_id) :: id_density
real(rke), allocatable, target :: rho_fabm(:,:,:)
!type (type_fabm_interior_variable_id) :: id_depth
real(rke), allocatable, target :: depth(:,:,:)
type (type_fabm_interior_variable_id) id_downwelling_photosynthetic_radiative_flux
real(rke), allocatable, target :: downwelling_photosynthetic_radiative_flux(:,:,:)
! id_downwelling_shortwave_flux
! id_fractional_saturation_of_oxygen
! id_mass_concentration_of_suspended_matter
! id_mole_concentration_of_ammonium
! id_mole_concentration_of_carbonate_expressed_as_carbon
! id_mole_concentration_of_dissolved_inorganic_carbon
! id_mole_concentration_of_dissolved_iron
! id_mole_concentration_of_nitrate
! id_mole_concentration_of_phosphate
! id_mole_concentration_of_silicate
! id_net_rate_of_absorption_of_shortwave_energy_in_layer
! id_ph_reported_on_total_scale
type (type_fabm_interior_variable_id) :: id_practical_salinity
real(rke), allocatable, target :: salt(:,:,:)
type (type_fabm_interior_variable_id) :: id_pressure
real(rke), allocatable, target :: pressure(:,:,:)
! id_secchi_depth
!type (type_fabm_interior_variable_id) :: id_temperature

! Surface variables
! id_cloud_area_fraction
! id_ice_area_fraction
type (type_fabm_horizontal_variable_id) :: id_mole_fraction_of_carbon_dioxide_in_air
real(rke), allocatable, target :: mole_fraction_of_carbon_dioxide_in_air(:,:)
! id_surface_air_pressure
! id_surface_albedo
!KB ipsw - 1 cal cm-2 s-1 = 41868 W/m2
type (type_fabm_horizontal_variable_id) :: id_surface_downwelling_photosynthetic_radiative_flux
real(rke), allocatable, target :: surface_downwelling_photosynthetic_radiative_flux(:,:)
! id_surface_downwelling_photosynthetic_radiative_flux_in_air
type (type_fabm_horizontal_variable_id) :: id_surface_downwelling_shortwave_flux
real(rke), allocatable, target :: surface_downwelling_shortwave_flux(:,:)
! id_surface_downwelling_shortwave_flux_in_air
! id_surface_drag_coefficient_in_air
! id_surface_specific_humidity
! id_surface_temperature
type (type_fabm_horizontal_variable_id) :: id_windspeed
real(rke), allocatable, target :: windspeed(:,:)

! Bottom variables
! id_bottom_depth
! id_bottom_depth_below_geoid
! id_bottom_roughness_length
type (type_fabm_horizontal_variable_id) :: id_bottom_stress
real(rke), allocatable, target :: bottom_stress(:,:)

! Global variables
! id_number_of_days_since_start_of_the_year

! Universal variables
! id_total_carbon
! id_total_iron
! id_total_nitrogen
! id_total_phosphorus
! id_total_silicate
! horizontal FABM ids
! interior FABM variables - calculated from UVic_ESCM variables
! horizontal FABM variables

! public available routines
public fabm_configure
public fabm_link_data
public fabm_update
public fabm_list
public fabm_clean

! module level variables - static or allocatable? Same goes with e.g. windspeed and rho_fabm

!integer, parameter :: npel=nt-2

real(rke) :: pelagic_sms(imt,km,1,nt-2)
  !! pelagic source-sink terms in one j-stride
real(rke) :: surface_flux(imt,jmt,nt-2)
  !! surface fluxes
real(rke) :: surface_sms(imt,jmt,nt-2)
  !! surface source-sink terms
real(rke) :: bottom_flux(imt,jmt,nt-2)
  !! bottom fluxes
real(rke) :: bottom_sms(imt,jmt,nt-2)
  !! bottom source-sink terms
!real(rke) :: w(imt,km,jsmw:jemw,nt)
real(rke) :: w(imt,km,jmt,nt-2)
  !! vertical velocity in m/s

integer :: nsurface
integer :: npelagic
integer :: nsediment

!> @note
!> The variables to hold surface, pelagic and bottom state variables
!> comes from mw.h - included in a O_fabm compilation clause
!> @endnote
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine fabm_configure(dt,yaml_file)
   !! read fabm.yaml and call FABM configuration subroutines

   real(rke), intent(in) :: dt
      !! bio-geochemical time step as set by MOM2 [s]
   character(len=*), intent(in), optional :: yaml_file
      !! name of alternativ FABM configuration file

   print*, '==> Initializing FABM component with nt=',nt

   if (present(yaml_file)) then
      model => fabm_create_model(trim(yaml_file))
   else
      model => fabm_create_model('fabm.yaml')
   end if

   nsurface = size(model%surface_state_variables)
   npelagic = size(model%interior_state_variables)
   nsediment = size(model%bottom_state_variables)

   if (nt-2 .ne. npelagic) then
      print*, 'nt (UVic) = ',nt-2
      print*, 'nsurface  = ',nsurface
      print*, 'npelagic  = ',npelagic
      print*, 'nsediment = ',nsediment
      stop 'fabm_configure()'
   end if

   !parameter (jsmw=2, jemw=jmw-1) - parameter (jmw=jmt)
#ifdef DEBUG
   print*, imt,jmt,jmw
   !print*, jrow,js,je,is,ie
   print*, 'zt: ',shape(zt)
   print*, 't: ',shape(t)
   print*, 'sbc: ',shape(sbc)
   print*, 'src: ',shape(src) ! imt,km,jsmw:jemw,nsrc
   print*, 'source: ',shape(source) ! imt,km,jsmw:jemw
   print*, 'rho: ',shape(rho)
   !stop 112
#endif   

   call model%set_domain(imt,km,jmt,dt)
   call model%set_domain_start(2,1,2)
   call model%set_domain_stop(imt-1,km,jmt-1)
   call model%set_mask(tmask,tmask(:,1,:))
   call model%set_bottom_index(kmt)

   !! @note
   !! seems tmask is not initialised until called in mom() - 
   !! i.e. after initialization - so all values are 0 here
   !! @endnote
end subroutine fabm_configure

!-----------------------------------------------------------------------

subroutine fabm_link_data()
   !! link all FABM configured external dependencies - and call
   !! model%start() to assure proper configuration
   integer :: j,k,n

   ! link to time in-dependent data that do require transformation
   call link_grid()

   ! link to time dependent data that do NOT require transformation
   call model%link_interior_data(fabm_standard_variables%temperature,t(:,:,:,itemp,0))

   ! link to time dependent data that do require transformation
   ! initialize and update time changing environmental variables
   call link_wind()
   call link_mole_fraction_of_carbon_dioxide_in_air()
   call link_surface_downwelling_photosynthetic_radiative_flux()
   call link_surface_downwelling_shortwave_flux()
   call link_bottom_stress()
   call link_downwelling_photosynthetic_radiative_flux()
   call link_salinity()
   call link_density()

   ! link to FABM's surface state variables
   do n = 1,nsurface
      !call model%link_surface_state_data(n, sed(:,:,n))
   end do

   ! link to FABM's interior state variables
   do n = 1, size(model%interior_state_variables)
      call model%link_interior_state_data(n, t(:,:,:,2+n,0))
      mapt(2+n) = trim(model%interior_state_variables(n)%name)
      itrc(n+2) = n+2
      !KBmapst(2+n) = 's'//trim(mapt(2+n))
   end do

   ! link to FABM's bottom state variables
   do n = 1,nsediment
      call model%link_bottom_state_data(n, sed(:,:,n))
   end do

   do j = 2, jmt-1
      do k = 1, km
!         call model%initialize_interior_state(2, imt-1, k, j)
      end do
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
   !! update the environment and calculate the source/sink terms -
   !! is called with the same argument list as mom() calls tracer()
   !! i.e. the specification of the active UVic window - typically
   !! the full domain on modern hardware
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

#ifdef DEBUG
   print*, 'fabm_update:',joff,js,je,is,ie
#endif

   ! t(:,:,:,var,0) is updated in loadmw() in mom()
   ! this is done before the call to tracer() - and thus
   ! data are ready here
   call update_data(joff)

   call model%prepare_inputs()

   ! update the surface
   surface_flux = 0._rke
   surface_sms = 0._rke
   if (nsurface > 0) then
      do j=js,je
         call model%get_surface_sources(is,ie,j, &
                    surface_flux(is:ie,j,:),surface_sms(is:ie,j,:))
      end do
   end if

   ! update the pelagic
   do j=js,je
      do k=1,km
         call model%get_interior_sources(is,ie,k,j,pelagic_sms(is:ie,k,1,:))
      end do
      src(is:ie,:,j,3:) = pelagic_sms(is:ie,:,1,:)
   end do

   ! update the bottom
   bottom_flux = 0._rke
   bottom_sms = 0._rke
   if (nsediment > 0) then
      do j=js,je
         call model%get_bottom_sources(is,ie,j, &
                    bottom_flux(is:ie,j,:),bottom_sms(is:ie,j,:))
      end do
   end if

   ! fold the surface and bottom flux terms - src keeps track on 
   ! which variables actually have sources - itrc(n) and the size
   ! of src reflects this - tracer.F90 line 1122
   do n=1,nt-2
      do j=js,je
         do i=is,ie
            if (kmt(i,j) > 0) then
               ! surface
               k=1
               src(i,k,j,3:)=src(i,k,j,n)+surface_flux(i,j,n)/dz(i,k,j)
               ! bottom
               k=kmt(i,j)
               src(i,k,j,3:)=src(i,k,j,n)+bottom_flux(i,j,n)/dz(i,k,j)
            end if
         end do
      end do
   end do

   ! vertical velocities
   do j=js,je
      do k=1,km
!         call model%get_vertical_movement(is,ie,k,j,w)
      end do
   end do

   call model%finalize_outputs()
   stop "fabm_update"
end subroutine fabm_update

!-----------------------------------------------------------------------

subroutine fabm_clean()
   !! de-allocate all allocated arrays
   if (allocated(windspeed)) deallocate(windspeed)
   if (allocated(mole_fraction_of_carbon_dioxide_in_air)) &
           deallocate(mole_fraction_of_carbon_dioxide_in_air)
   if (allocated(surface_downwelling_photosynthetic_radiative_flux)) &
           deallocate(surface_downwelling_photosynthetic_radiative_flux)
   if (allocated(surface_downwelling_shortwave_flux)) &
           deallocate(surface_downwelling_shortwave_flux)
   if (allocated(bottom_stress)) deallocate(bottom_stress)
   if (allocated(salt)) deallocate(salt)
   if (allocated(downwelling_photosynthetic_radiative_flux)) &
           deallocate(downwelling_photosynthetic_radiative_flux)
   if (allocated(rho_fabm)) deallocate(rho_fabm)
   ! :
   ! :
end subroutine fabm_clean

!-----------------------------------------------------------------------

subroutine update_data(joff)
   !! update all time varying FABM configured external dependencies
   !! by calling individual update routines - tests done in routines
   integer, intent(in) :: joff
     !! offset row in global window

   call update_wind()
   call update_mole_fraction_of_carbon_dioxide_in_air()
   call update_surface_downwelling_photosynthetic_radiative_flux()
   call update_surface_downwelling_shortwave_flux()
   call update_bottom_stress(joff)
   call update_downwelling_photosynthetic_radiative_flux()
   call update_salinity()
   call update_density()
end subroutine update_data

!-----------------------------------------------------------------------

subroutine link_grid()
   !! Allocate and link grid related FABM standard variables that 
   !! are being transformed from UVic native variables [cm -> m].

   integer :: rc
   integer :: i,j,k

   allocate(depth(imt,km,jmt),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (depth)'
   depth = 0._rke
   allocate(pressure(imt,km,jmt),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (pressure)'
   pressure = 0._rke
   allocate(dz(imt,km,jmt),stat=rc)
   if (rc /= 0) stop 'link_grid(): Error allocating (dz)'
   dz = 0._rke
   do j=2,jmt-1
      do i=2,imt-1
         if (kmt(i,j) > 0) then
            depth(i,:,j) = zt/100._rke
            pressure(i,:,j) = depth(i,:,j)/10._rke
            dz(i,:,j) = dzt/100._rke
         end if
      end do
   end do
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
      allocate(windspeed(imt,jmt),stat=rc)
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
   do j=2,jmt-1
      do i=2,imt-1
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
      allocate(mole_fraction_of_carbon_dioxide_in_air(imt,jmt),stat=rc)
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
      do j=2,jmt-1
         do i=2,imt-1
            if (kmt(i,j) > 0) mole_fraction_of_carbon_dioxide_in_air(i,j) =  10._rke !KB
         end do
      end do
   end if
end subroutine update_mole_fraction_of_carbon_dioxide_in_air

!-----------------------------------------------------------------------

subroutine link_surface_downwelling_photosynthetic_radiative_flux()
   integer rc
   id_surface_downwelling_photosynthetic_radiative_flux = model% &
        get_horizontal_variable_id(standard_variables% &
        surface_downwelling_photosynthetic_radiative_flux)
   if (model%variable_needs_values(id_surface_downwelling_photosynthetic_radiative_flux)) then
      allocate(surface_downwelling_photosynthetic_radiative_flux(imt,jmt),stat=rc)
      if (rc /= 0) stop 'link_surface_downwelling_photosynthetic_radiative_flux(): &
              Error allocating (surface_downwelling_photosynthetic_radiative_flux)'
      surface_downwelling_photosynthetic_radiative_flux = 0._rke
      call model%link_horizontal_data( &
              id_surface_downwelling_photosynthetic_radiative_flux, &
              surface_downwelling_photosynthetic_radiative_flux)
   end if
end subroutine link_surface_downwelling_photosynthetic_radiative_flux

!-----------------------------------------------------------------------

subroutine update_surface_downwelling_photosynthetic_radiative_flux()
   !! calculate the ?????????? flux in W/m^2

   integer i,j
   if (model%variable_needs_values(id_surface_downwelling_photosynthetic_radiative_flux)) then
      do j=2,jmt-1
         do i=2,imt-1
            if (kmt(i,j) > 0) surface_downwelling_photosynthetic_radiative_flux(i,j) = 200._rke !KB
         end do
      end do
   end if
end subroutine update_surface_downwelling_photosynthetic_radiative_flux

!-----------------------------------------------------------------------

subroutine link_surface_downwelling_shortwave_flux()
   integer rc
   id_surface_downwelling_shortwave_flux = model% &
        get_horizontal_variable_id(standard_variables% &
        surface_downwelling_shortwave_flux)
   if (model%variable_needs_values(id_surface_downwelling_shortwave_flux)) then
      allocate(surface_downwelling_shortwave_flux(imt,jmt),stat=rc)
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
      do j=2,jmt-1
         do i=2,imt-1
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
      allocate(bottom_stress(imt,jmt),stat=rc)
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
      do j=2,jmt-1
         jrow = j+joff
         do i=2,imt-1
            if (kmt(i,jrow) > 0) bottom_stress(i,j) = x*sqrt( &
                bmf(i,jrow,1)**2 + bmf(i,jrow,2)**2)
!            if (kmt(i,j) > 0) bottom_stress(i,j) =  0.001_rke !KB
         end do
      end do
   end if
end subroutine update_bottom_stress

!-----------------------------------------------------------------------

subroutine link_downwelling_photosynthetic_radiative_flux()
   !! get salinity FABM standard variable and if needed by FABM allocate
   !! memory
   integer rc
   id_downwelling_photosynthetic_radiative_flux = &
   model%get_interior_variable_id(fabm_standard_variables%downwelling_photosynthetic_radiative_flux)
   if (model%variable_needs_values(id_downwelling_photosynthetic_radiative_flux)) then
      allocate(downwelling_photosynthetic_radiative_flux(imt,km,jmt),stat=rc)
      if (rc /= 0) stop 'link_salinity(): Error allocating (downwelling_photosynthetic_radiative_flux)'
      downwelling_photosynthetic_radiative_flux = 0._rke
      call model%link_interior_data( &
              id_downwelling_photosynthetic_radiative_flux, &
              downwelling_photosynthetic_radiative_flux)
   end if
end subroutine link_downwelling_photosynthetic_radiative_flux

!-----------------------------------------------------------------------

subroutine update_downwelling_photosynthetic_radiative_flux()
   !! calculate salinity in PSU according to $$S = 35 + 1000*S_{UVic}$$
   integer i,j,k
   if (model%variable_needs_values(id_downwelling_photosynthetic_radiative_flux)) then
      do j=2,jmt-1
         do k=1,km
            do i=2,imt-1
               if (kmt(i,j) > 0) downwelling_photosynthetic_radiative_flux(i,k,j) = &
                            35._rk+1000._rke*t(i,k,j,isalt,0)
            end do
         end do
      end do
   end if
end subroutine update_downwelling_photosynthetic_radiative_flux
!-----------------------------------------------------------------------

subroutine link_salinity()
   !! get salinity FABM standard variable and if needed by FABM allocate
   !! memory
   integer rc
   id_practical_salinity = model%get_interior_variable_id(fabm_standard_variables%practical_salinity)
   if (model%variable_needs_values(id_practical_salinity)) then
      allocate(salt(imt,km,jmt),stat=rc)
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
      do j=2,jmt-1
         do k=1,km
            do i=2,imt-1
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
      allocate(rho_fabm(imt,km,jmt),stat=rc)
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
            do i=2,imt-1
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
