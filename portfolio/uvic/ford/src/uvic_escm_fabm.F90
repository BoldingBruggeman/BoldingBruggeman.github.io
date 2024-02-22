!> Copyright (C) 2024 Bolding & Bruggeman
!>

!> This program is used to assure that the sizes of arrays
!> specified in size.h and the size deducted from fabm.yaml
!> are in agreement.
!>
!> @note The way the include files are formatted does not
!> allow to be included in F90 free format source code files.
!>
!> Still to be solved is the issue of getting the variables in mw.h to
!> double precission - for now mw.h in the base UVic source code folder
!> has a kind=8 on all real arguments
!>
!> An alternative method is to provide CMake with options to compile
!> in 8 (4) byte. This is used for now - and no changes to mw.h
!> are necessary
!> @endnote

program UVic_ESCM_FABM

use uvic_common_blocks
use uvic_fabm, only: fabm_configure, fabm_list

implicit none

integer :: n
character(len=128) :: fname = 'fabm.yaml'
logical :: exists

n=command_argument_count()
if (n .eq. 1) then
   call get_command_argument(n,fname)
end if

INQUIRE(FILE=trim(fname), EXIST=exists)
if ( .not. exists) then
   write(*,*) 'usage: uvic_escm_fabm.x <yaml_file>'
   stop 0
endif

itemp=1
isalt=2

print*, '== UNIVERSITY OF VICTORIA EARTH SYSTEM CLIMATE MODEL with FABM =='
print*, '== Source configuration from: ',_DIR_UVIC_CODE_
print*, '== FABM configuration from:   ',trim(fname)
call fabm_configure(dble(0.),trim(fname))
call fabm_list()
end
