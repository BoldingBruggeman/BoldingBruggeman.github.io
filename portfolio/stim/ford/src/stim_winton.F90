!> The winton ice model
!>
!> After Michael Winton
!>
!> authors: Karsten Bolding, Adolf Stips and Jesper Larsen

   MODULE stim_winton

!>
!>  The model consists of a zero heat capacity snow layer overlying two equally
!>  thick sea ice layers. The upper ice layer has a variable heat capacity to
!>  represent brine pockets. The lower ice layer has a fixed heat capacity.
!>  The prognostic variables are hs (snow layer thickness), hi (ice layer
!>  thickness), T1 and T2, the upper and lower ice layer temperatures located
!>  at the midpoints of the layers. The ice model performs two functions, the
!>  first is to calculate the ice temperature and the second is to calculate
!>  changes in the thickness of ice and snow.
!>
!>------------------------------------------------------------------------------
!>                                                                              
!>                       THREE-LAYER VERTICAL THERMODYNAMICS                    
!>                                                                              
!> Reference:  M. Winton, 2000: "A reformulated three-layer sea ice model",     
!>            Journal of Atmospheric and Oceanic Technology, 17, 525-531.       
!>                                                                              
!>------------------------------------------------------------------------------
!>        -> +---------+ <- Ts - diagnostic surface temperature ( <= 0C )       
!>       /   |         |                                                        
!>     hs    |  snow   | <- 0-heat capacity snow layer                          
!>       \   |         |                                                       
!>        -> +---------+                                                        
!>       /   |         |                                                        
!>      /    |         | <- T1 - upper 1/2 ice temperature; this layer has      
!>     /     |         |         a variable (T/S dependent) heat capacity       
!>     \     |         |                                                       
!>      \    |         | <- T2 - lower 1/2 ice temp. (fixed heat capacity)      
!>       \   |         |                                                        
!>        -> +---------+ <- Tf - base of ice fixed at seawater freezing temp.   
!>                                                                              
!>                                                     Mike Winton (mw@gfdl.gov)
!>------------------------------------------------------------------------------

!>  Note: in this implementation the equations are multiplied by hi to improve
!>  thin ice accuracy
!>
!>  The code is based on the open source sea ice model included in the Modular
!>  Ocean Model.
!>

!   hi      |   ice   | <- non 0-heat capacity ice layer                           

   use stim_variables
   use ice_thm_mod
   IMPLICIT NONE

   private

   public :: init_stim_winton
   public :: do_stim_winton

   real(rk), pointer :: Ts,T1,T2
   real(rk), pointer :: hi,hs,dh1,dh2
   real(rk), pointer :: trn
   real(rk)          :: pen
   real(rk), pointer :: tmelt,bmelt
   real(rk), pointer :: fb    ! heat flux from ocean to ice bottom (W/m^2)

   contains

   SUBROUTINE init_stim_winton(Ta)

   real(rk), intent(in) :: Ta
     !! Air temperature [C]
!-----------------------------------------------------------------------
   trn => transmissivity
   Ts => Tice_surface
   Ts = Ta
   T1 => Tice(1)
   T1 = Ta
   T2 => Tice(2)
   T2 = Tf
   hs => Hsnow
   hi => Hice
   dh1 => dHis
   dh2 => dHib
   tmelt => surface_ice_energy
   bmelt => bottom_ice_energy
   fb => ocean_ice_flux
END SUBROUTINE init_stim_winton

!-----------------------------------------------------------------------

   SUBROUTINE do_stim_winton(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
     !! This SUBROUTINE updates the sea ice prognostic variables. The updated
     !! variables are upper ice layer temperature (T1), lower ice layer temperature
     !! (T2), snow thickness (hs), and ice thickness (hi).
     !!
     !! The ice model performs this in two steps. First the temperatures are updated
     !! and secondly the changes in ice and snow thickness are calculated.
     !!
     !! Any surplus energy that is not used for melting is returned in tmelt and
     !! bmelt.
     !!
     !! Evaporation and bottom ablation formation are not included in
     !! this version of the model. Furthermore we do not keep an explicit water
     !! and salt budget for the sea ice and how that affects the water and salt
     !! budgets in the ocean.

   real(rk), intent(in)    :: dz,dt,Ta,S,precip,Qsw
   integer, intent(inout)  :: ice_cover
   real(rk), intent(inout) :: Tw
   interface
      SUBROUTINE Qfluxes(T,qh,qe,qb)
         integer, parameter                   :: rk = kind(1.d0)
         real(rk), intent(in)                 :: T
         real(rk), intent(out)                :: qh,qe,qb
      END SUBROUTINE
   end interface

   real(rk)        :: I     ! solar absorbed by upper ice (W/m^2)
   real(rk)        :: evap  ! evaporation of ice (m/s)
   real(rk)        :: snow
   real(rk)        :: A, B, dts=0.01
   real(rk)        :: qe,qh,qb
   real(rk)        :: h1,h2
   real(rk)        :: ts_new
   real(rk)        :: frazil
   real(rk)        :: heat_to_ocn, h2o_to_ocn, h2o_from_ocn, snow_to_ice
!-----------------------------------------------------------------------
   tmelt = 0._rk
   bmelt = 0._rk

   ! Calculate seawater freezing temperature
   Tf = -0.0575_rk*S

   if (ice_cover .gt. 0) then
      call ice_optics(albedo_ice, pen, trn, hs, hi, ts, Tf)
      I = Qsw*(1._rk-trn)
      h1 = hi/2._rk
      h2 = h1

      ! check this out
      call Qfluxes(Ts,qe,qh,qb)
      A = -(qe+qh+qb) ! (7-)
      call Qfluxes(Ts+dts,qe,qh,qb)
      B = -(qe+qh+qb)
      B = (B-A)/dts ! (8)
      A = A-I - Ts*B ! (-7)

      !https://github.com/mom-ocean/MOM5/blob/08266af73b04d2334be4a52d0c45c174f447cee4/src/ice_sis/ice_model.F90
      call ice3lay_temp(hs,hi,t1,t2,ts_new,A,B,pen*I,Tf,fb,dt,tmelt,bmelt)
      ts = ts_new
 !     frazil = 0._rk
      Hfrazil = 0._rk
   else
      frazil = -(Tw-Tf)*dz*Cw
      if (frazil .gt. 0._rk) Hfrazil = frazil/(rho_ice*L_ice)
   end if

   if (ice_cover .gt. 0 .or. frazil .gt. 0._rk) then
      call ice3lay_resize(hs, hi, t1, t2, snow, frazil, evap, tmelt, bmelt, &
                          tf, heat_to_ocn, h2o_to_ocn, h2o_from_ocn,       &
                          snow_to_ice)
   !                       snow_to_ice, bablt)
!write(*,*) 'CC ',heat_to_ocn, h2o_to_ocn, h2o_from_ocn
   end if

hs = 0._rk
   if (hi .gt. 0._rk) then
      ice_cover = 2
   else
      ice_cover = 0
   end if
END SUBROUTINE do_stim_winton

!-----------------------------------------------------------------------

   END MODULE stim_winton

!-----------------------------------------------------------------------
! Copyright by the GETM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
