module iso_c_utilities

use iso_c_binding
implicit none

contains

function char_array_to_string(char_array)
    character(c_char) :: char_array(:)
    character(len=size(char_array)) :: char_array_to_string
    integer :: i
    do i = 1, size(char_array)
        char_array_to_string(i:i) = char_array(i)
    enddo
end function

end module
