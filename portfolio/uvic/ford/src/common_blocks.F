c Copyright (C) 2024 Bolding & Bruggeman

c> @note
c> We need this F77 file to get the information of array sizes and 
c> a number of static allocated MOM2 variables. 
c>
c> The way the UVic include files are formatted does not
c> allow to be included in F90 free format source code files.
c>
c> Still to be solved is the issue of getting the variables in mw.h to
c> double precission - for now updates/fabm/source/mom/mw.h has a 
c! kind=8 on all real arguments
c> @endnote

      module uvic_common_blocks
      include "size.h"
      include "coord.h"
      include "csbc.h"
      include "mw.h"
      include "levind.h"
      end module uvic_common_blocks

