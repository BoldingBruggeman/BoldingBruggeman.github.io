! Copyright (C) 2024 Bolding & Bruggeman

!> This program is used to assure that the sizes of arrays
!> specified in size.h and the size deducted from fabm.yaml
!> are in aggreement.
!>
!> @note The way the include files are formatted does not
!> allow to be included in F90 free format source code files.
!> Still to be solved is the issue of getting the variables in mw.h to
!> double precission - for now mw.h in the base UVic source code folder
!> has a kind=8 on all real arguments @endnote

program UVic_ESCM_FABM

use uvic_common_blocks
use uvic_fabm, only: fabm_configure, fabm_list

implicit none

#if 0
print*, size(u)
print*, shape(u)
print*, size(t)
print*, shape(t)
#endif

! from UVic_ESCM.f: tracer_init
itemp=1
isalt=2

print*, '== UNIVERSITY OF VICTORIA EARTH SYSTEM CLIMATE MODEL with FABM =='
print*, '== Configuration from: ',_DIR_UVIC_CODE_
! tracer_init 
call configure_fabm()
call list_fabm()
end

