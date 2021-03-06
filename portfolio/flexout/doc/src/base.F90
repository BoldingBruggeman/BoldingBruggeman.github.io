module output_operators_base

   use output_manager_core
   use field_manager
   use yaml_types

   implicit none

   private

   public type_operator_result, type_universal_operator_result

   type, extends(type_base_output_field) :: type_operator_result
      class (type_base_operator), pointer :: operator => null()
      class (type_base_output_field), pointer :: source => null()
      real(rk)              :: result_0d
      real(rk), allocatable :: result_1d(:)
      real(rk), allocatable :: result_2d(:,:)
      real(rk), allocatable :: result_3d(:,:,:)
      type (type_dimension_pointer), allocatable :: dimensions(:)
   contains
      procedure :: new_data
      procedure :: before_save
      procedure :: get_metadata
      procedure :: flag_as_required
      procedure :: fill
      procedure :: get_field
   end type

   type, extends(type_operator_result) :: type_universal_operator_result
   contains
      procedure :: get_field => universal_get_field
   end type

   contains

   recursive subroutine new_data(self)
      class (type_operator_result), intent(inout) :: self
      call self%source%new_data()
   end subroutine

   recursive subroutine before_save(self)
      class (type_operator_result), intent(inout) :: self
      call self%source%before_save()
   end subroutine

   recursive subroutine get_metadata(self, long_name, units, dimensions, minimum, maximum, fill_value, standard_name, path, attributes)
      class (type_operator_result), intent(in) :: self
      character(len=:), allocatable, intent(out), optional :: long_name, units, standard_name, path
      type (type_dimension_pointer), allocatable, intent(out), optional :: dimensions(:)
      real(rk), intent(out), optional :: minimum, maximum, fill_value
      type (type_attributes), intent(out), optional :: attributes

      ! Workaround for gfortran BUG 88511 - passing optional allocatable deferred length character arguments to the next routine causes memory corruption
      character(len=:), allocatable :: long_name2, units2, standard_name2, path2

      call self%source%get_metadata(long_name2, units2, dimensions, minimum, maximum, fill_value, standard_name2, path2, attributes)

      ! Workaround for gfortran
      if (present(long_name) .and. allocated(long_name2)) long_name = long_name2
      if (present(units) .and. allocated(units2)) units = units2
      if (present(standard_name) .and. allocated(standard_name2)) standard_name = standard_name2
      if (present(path) .and. allocated(path2)) path = path2

      if (present(dimensions) .and. allocated(self%dimensions)) dimensions(:) = self%dimensions(:)
   end subroutine

   recursive subroutine flag_as_required(self, required)
      class (type_operator_result), intent(inout) :: self
      logical,                      intent(in)    :: required

      call self%source%flag_as_required(required)
   end subroutine

   subroutine fill(self, value)
      class (type_operator_result), intent(inout) :: self
      real(rk),                     intent(in)    :: value

      if (allocated(self%result_3d)) then
         self%result_3d(:,:,:) = value
      elseif (allocated(self%result_2d)) then
         self%result_2d(:,:) = value
      elseif (allocated(self%result_1d)) then
         self%result_1d(:) = value
      else
         self%result_0d = value
      end if
   end subroutine

   recursive function get_field(self, field) result(output_field)
      class (type_operator_result), intent(in) :: self
      type (type_field), target                :: field
      class (type_base_output_field), pointer  :: output_field
      output_field => self%source%get_field(field)
   end function

   recursive function universal_get_field(self, field) result(output_field)
      class (type_universal_operator_result), intent(in) :: self
      type (type_field), target                          :: field
      class (type_base_output_field), pointer            :: output_field
      output_field => self%type_operator_result%get_field(field)
      if (associated(output_field)) output_field => self%operator%apply(output_field)
   end function

end module