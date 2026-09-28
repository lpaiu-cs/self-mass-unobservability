! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton, Ed Brown
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
!
! ***********************************************************************
 
      module chem_isos_io
      use chem_def
      
      implicit none


      contains
      

      subroutine do_read_chem_isos(isotopes_filename, ierr)
      	use utils_lib
      	use crlibm_lib, only: str_to_vector
      	use const_def, only: mesa_data_dir
         character (len=*), intent(in) :: isotopes_filename
      	integer, intent(out) :: ierr
      	integer :: i, j, k, iounit, pass, nvec, Z, N
      	character (len=256) :: filename, buf
      	real(dp), target :: vec_ary(256)
      	real(dp), pointer :: vec(:)
      	
      	ierr = 0
      	vec => vec_ary
	
      	iounit = alloc_iounit(ierr)
      	if (ierr /= 0) then
      	   write(*,*) 'failed to alloc iounit'
      	   return
      	end if
      	
      	filename = trim(mesa_data_dir) // '/chem_data/' // trim(isotopes_filename)
         num_chem_isos = 0
      	
      	do pass = 1, 2
      	
         	open(unit=iounit, file=trim(filename), iostat=ierr, status='old',action='read')
         	if ( ierr /= 0 ) then
         		write(*,*) 'unable to open '// trim(filename)
         		return
         	end if
         	read(iounit,*) ! skip line 1
         	
         	if (pass == 1) then
         	
         	   do ! 4 lines per nuclide
            		read(iounit, *, iostat=ierr) 
            		if (ierr /= 0) exit
            		read(iounit, *, iostat=ierr) 
            		if (ierr /= 0) exit
            		read(iounit, *, iostat=ierr) 
            		if (ierr /= 0) exit
            		read(iounit, *, iostat=ierr) 
            		if (ierr /= 0) exit
            		num_chem_isos = num_chem_isos+1
         	   end do
         	   if (num_chem_isos == 0) then
							write (*,*) 'unable to retrieve isotopes from '//trim(filename)
							return
						end if
         	else
         	
            	call allocate_nuclide_data(chem_isos, num_chem_isos, ierr)
            	if (ierr /= 0) then
            	   write(*,*) 'unable to allocate nuclide data'
            	   return
            	end if
	
            	do i = 1, num_chem_isos
            	
            		read(iounit, *, iostat=ierr) &
            		   chem_isos% name(i), chem_isos% W(i), chem_isos% Z(i), chem_isos% N(i), &
            		 	chem_isos% spin(i), chem_isos% mass_excess(i)
            		if (ierr /= 0) exit
            		
            		do k=1,3
               		read(iounit,'(a)',iostat=ierr) buf
               		if (ierr == 0) then
               		   call str_to_vector(buf, vec, nvec, ierr)
               		   if (nvec /= 8) ierr = -1
               		end if
               		if (ierr /= 0) exit
               		do j=1,8
               		   chem_isos% pfcn(j+(k-1)*8, i) = vec(j)
               		end do
               	end do
            		if (ierr /= 0) exit
            		
            		chem_isos% chem_id(i) = i
            		chem_isos% nuclide(i) = i
            		chem_isos% isomeric_state(i) = get_isomeric_state(chem_isos% name(i), ierr)
            		
            	end do
            	if (ierr /= 0) then
								write (*,*) 'something went wrong in read of '//trim(filename)
								return
							end if

         	end if
         	
         	close(iounit)
         	
         end do
         
      	call free_iounit(iounit)

      	if (ierr /= 0) return
      	
         call do_create_nuclides_dict(chem_isos, chem_isos_dict, ierr)
      	if (ierr /= 0) return
	
         chem_isos% Z_plus_N = chem_isos% Z + chem_isos% N
         
      	chem_isos% binding_energy = chem_isos% Z*del_Mp + chem_isos% N*del_Mn - chem_isos% mass_excess
      	where (chem_isos% Z_plus_N <= 1) chem_isos% binding_energy = 0.0 ! neutrons and protons by various names
      	
      	element_min_N = 99999
      	element_max_N = -1
      	do i = 1, num_chem_isos
      	   Z = chem_isos% Z(i)
      	   N = chem_isos% N(i)
      	   if (N < element_min_N(Z)) element_min_N(Z) = N
      	   if (N > element_max_N(Z)) element_max_N(Z) = N
      	end do      	
      	
      	contains
      	
      	integer function get_isomeric_state(name, ierr)
      	   character (len=*), intent(in) :: name
      	   integer, intent(out) :: ierr
      	   integer :: i, len
      	   ierr = 0
      	   get_isomeric_state = 0
      	   len = len_trim(name)
      	   do i=1,len
      	      if (name(i:i) == '-') then
      	         read(name(i+1:len),*,iostat=ierr) get_isomeric_state
      	         if (ierr /= 0 .or. get_isomeric_state < 0 .or. get_isomeric_state > 99) then
      	            write(*,*) 'ERROR: invalid name for iso ' // trim(name) // ' in ' // trim(filename)
      	            return
      	         end if
      	         return
      	      end if
      	   end do
      	end function get_isomeric_state
      	

      end subroutine do_read_chem_isos
            
      
      subroutine do_create_nuclides_dict(nuclides, nuclides_dict, ierr)
         use utils_lib, only: integer_dict_define, integer_dict_create_hash, integer_dict_lookup
      	type(nuclide_data), intent(in) :: nuclides
         type (integer_dict), pointer :: nuclides_dict ! will be allocated
      	integer, intent(out) :: ierr
         integer :: i
         
         ierr = 0
         nullify(nuclides_dict)
         do i=1,nuclides% nnuclides
            call integer_dict_define(nuclides_dict, nuclides% name(i), i, ierr)
            if (ierr /= 0) return
         end do

         call integer_dict_create_hash(nuclides_dict, ierr)
         if (ierr /= 0) return
         
      end subroutine do_create_nuclides_dict
      

      end module chem_isos_io