
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmat(vmt,vir,vmat)
! generates potential matrix elements for all states and k-points
use modmain
use modmpi
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: vir(ngrtot)
complex(8), intent(out) :: vmat(nstsv,nstsv,nkpt)
! local variables
integer ld,is,ia,ias
integer ik,ispn,n,lp
! local arrays
real(8), allocatable :: rfmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
! allocate local arrays
allocate(rfmt(lmmaxvr,nrcmtmax,natmtot))
! convert muffin-tin potential to spherical coordinates
ld=lmmaxvr*lradstp
do is=1,nspecies
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ias)
!$OMP DO
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
     vmt(:,:,ias),ld,0.d0,rfmt(:,:,ias),lmmaxvr)
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv) &
!$OMP PRIVATE(wfmt,wfir,ispn)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngrtot,nspinor,nstsv))
!$OMP CRITICAL
  write(*,'("Info(genvmat): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! calculate the second-variational wavefunctions for all states
  call genwfsv(.false.,.false.,.false.,ngk(1,ik),igkig(:,1,ik),evalsv(:,ik), &
   apwalm,evecfv,evecsv,wfmt,ngrtot,wfir)
  call genvmatk(rfmt,vir,wfmt,wfir,vmat(:,:,ik))
  deallocate(apwalm,evecfv,evecsv,wfmt,wfir)
end do
!$OMP END DO
!$OMP END PARALLEL
! broadcast matrix elements to every process
n=nstsv*nstsv
do ik=1,nkpt
  lp=mod(ik-1,np_mpi)
  call mpi_bcast(vmat(:,:,ik),n,mpi_double_complex,lp,mpi_comm_world,ierror)
end do
deallocate(rfmt)
return
end subroutine

