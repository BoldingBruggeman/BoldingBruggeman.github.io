!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 'obsice' ice module
!
! !INTERFACE:
   module stim_obsice
!
! !DESCRIPTION:
!  https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016JC012199
!
! !USES:
   use stim_variables, only: transmissivity, albedo_ice
   use stim_variables, only: rk, Hice, dHis, dHib, Tf, fdd
   use stim_variables, only: rho_ice, L_ice, ocean_ice_flux
   IMPLICIT NONE
!  Default all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_stim_obsice, do_stim_obsice, clean_stim_obsice
!
! !PUBLIC DATA MEMBERS:
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialisation of the ice variables
!
! !INTERFACE:
   subroutine init_stim_obsice(ice_cover)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)  :: ice_cover
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-----------------------------------------------------------------------
!BOC
   end subroutine init_stim_obsice
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do the 'obsice' ice calculations
!
! !INTERFACE:
   subroutine do_stim_obsice(ice_cover,Tw,S,Ta)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in)    :: S,Ta
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)  :: ice_cover
   real(rk), intent(inout) :: Tw
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   ! calculate ice thickness according to obsice 1938 
   Tf = -0.0575_rk*S
write(*,*) Tf
   if (ice_cover .eq. 2) then
      Tw = Tf
!KB      Tice = (alpha*Tf+Ta)/(1._rk+alpha) ! From mylake - but why?
      albedo_ice = 0.6_rk
   else
      Hice = 0._rk
      albedo_ice = 0._rk
      transmissivity = 1._rk
   end if
   end subroutine do_stim_obsice
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleaning up the 'obsice' ice variables
!
! !INTERFACE:
   subroutine clean_stim_obsice()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !LEVEL2 'clean_stim_obsice'

   end subroutine clean_stim_obsice
!EOC

!-----------------------------------------------------------------------

   end module stim_obsice

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
