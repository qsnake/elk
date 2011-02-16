module cdirac

use iso_c_utilities
use iso_c_binding
!use elk

implicit none

contains

subroutine c_rdirac(sol, n, l, k, np, nr, r, vr, eval, g0, f0) bind(c)
real(c_double), intent(in) :: sol
integer(c_int), intent(in) :: n
integer(c_int), intent(in) :: l
integer(c_int), intent(in) :: k
integer(c_int), intent(in) :: np
integer(c_int), intent(in) :: nr
real(c_double), intent(in) :: r(nr)
real(c_double), intent(in) :: vr(nr)
real(c_double), intent(inout) :: eval
real(c_double), intent(out) :: g0(nr)
real(c_double), intent(out) :: f0(nr)

call rdirac(sol, n, l, k, np, nr, r, vr, eval, g0, f0)
! !INPUT/OUTPUT PARAMETERS:
!   (sol,n,l,k,np,nr,r,vr,eval,g0,f0)
!   sol  : speed of light in atomic units (in,real)
!   n    : principal quantum number (in,integer)
!   l    : quantum number l (in,integer)
!   k    : quantum number k (l or l+1) (in,integer)
!   np   : order of predictor-corrector polynomial (in,integer)
!   nr   : number of radial mesh points (in,integer)
!   r    : radial mesh (in,real(nr))
!   vr   : potential on radial mesh (in,real(nr))
!   eval : eigenvalue without rest-mass energy (inout,real)
!   g0   : major component of the radial wavefunction (out,real(nr))
!   f0   : minor component of the radial wavefunction (out,real(nr))
end subroutine

end module
