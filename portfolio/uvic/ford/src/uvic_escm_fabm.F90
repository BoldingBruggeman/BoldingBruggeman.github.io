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

integer :: n
character(len=128) :: fname = 'fabm.yaml'
logical :: file_exists

n=command_argument_count()
if (n .eq. 1) then
   call get_command_argument(n,fname)
else
   write(*,*) 'using default fabm.yaml'
end if

inquire(FILE=trim(fname), EXIST=file_exists)

if (.not. file_exists) then
   write(*,*) 'usage: uvic_escm_fabm.x <yaml_file> - default fabm.yaml'
   stop 1
end if


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
print*, '== nt:   ',nt
print*, '== nsrc: ',nsrc
call fabm_configure(dble(3600.),yaml_file=trim(fname))
#if 0
if ((size(model%interior_state_variables) .ne. nt) then
   stop 'kaj ...'
end if
#endif
call fabm_list()
end
