
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqnhf(ikp,evecsvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(inout) :: evecsvp(nstsv,nstsv)
! local variables
integer is,ia,ias
integer ngknr,ik,jk,iq
integer ig,iv(3),igq0,igk
integer ist1,ist2,ist3
integer lmax,ir,irc,idm
integer lwork,info,ld
real(8) cfq,v(3),t1
complex(8) zrho01,zrho02,zt1,zt2
! automatic arrays
complex(8) sfacgq0(natmtot)
! allocatable arrays
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:),vgkcnr(:,:)
real(8), allocatable :: gkcnr(:),tpgkcnr(:,:)
real(8), allocatable :: vgqc(:,:),tpgqc(:,:),gqc(:)
real(8), allocatable :: jlgqr(:,:,:),jlgq0r(:,:,:)
real(8), allocatable :: evalsvp(:),evalsvnr(:)
real(8), allocatable :: vmt(:,:,:),vir(:)
real(8), allocatable :: bmt(:,:,:,:),bir(:,:)
real(8), allocatable :: rfmt(:,:),rwork(:)
complex(8), allocatable :: h(:,:),vmat(:,:)
complex(8), allocatable :: sfacgknr(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
complex(8), allocatable :: zfmt(:,:),work(:)
! external functions
complex(8) zfinp
external zfinp
!$OMP CRITICAL
write(*,'("Info(seceqnhf): ",I6," of ",I6," k-points")') ikp,nkpt
!$OMP END CRITICAL
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax),vgkcnr(3,ngkmax))
allocate(gkcnr(ngkmax),tpgkcnr(2,ngkmax))
allocate(vgqc(3,ngvec),tpgqc(2,ngvec),gqc(ngvec))
allocate(jlgqr(0:lmaxvr+npsden+1,ngvec,nspecies))
allocate(jlgq0r(0:lmaxvr,nrcmtmax,nspecies))
allocate(evalsvp(nstsv),evalsvnr(nstsv))
allocate(rwork(3*nstsv))
allocate(h(nstsv,nstsv),vmat(nstsv,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(sfacgknr(ngkmax,natmtot))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zfmt(lmmaxvr,nrcmtmax))
lwork=2*nstsv
allocate(work(lwork))
! coefficient of long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(:,ikp),evalsvp)
call getevecfv(vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for all states for the input k-point
call genwfsv(.false.,.false.,.false.,ngk(1,ikp),igkig(:,1,ikp),evalsvp,apwalm, &
 evecfv,evecsvp,wfmt1,ngrtot,wfir1)
!---------------------------------!
!     kinetic matrix elements     !
!---------------------------------!
call zgemm('N','N',nstsv,nstsv,nstsv,zone,kinmatc(:,:,ikp),nstsv,evecsvp, &
 nstsv,zzero,vmat,nstsv)
call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsvp,nstsv,vmat,nstsv,zzero,h, &
 nstsv)
!-----------------------------------------!
!     local potential matrix elements     !
!-----------------------------------------!
allocate(vmt(lmmaxvr,nrcmtmax,natmtot))
ld=lmmaxvr*lradstp
if (hybmix.lt.1.d0) then
! hybrid functional case
  allocate(rfmt(lmmaxvr,nrcmtmax))
  allocate(vir(ngrtot))
  t1=1.d0-hybmix
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        rfmt(:,irc)=vclmt(:,ir,ias)+t1*vxcmt(:,ir,ias)
      end do
      call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr,rfmt, &
       lmmaxvr,0.d0,vmt(:,:,ias),lmmaxvr)
    end do
  end do
  vir(:)=vclir(:)+t1*vxcir(:)
  call genvmatk(vmt,vir,wfmt1,wfir1,vmat)
  h(:,:)=h(:,:)+vmat(:,:)
  deallocate(rfmt,vir)
  if (spinpol) then
    allocate(bmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
    allocate(bir(ngrtot,ndmag))
    do idm=1,ndmag
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,t1,rbshtvr,lmmaxvr, &
           bxcmt(:,:,ias,idm),ld,0.d0,bmt(:,:,ias,idm),lmmaxvr)
        end do
      end do
      bir(:,idm)=t1*bxcir(:,idm)
    end do
    call genbmatk(bmt,bir,wfmt1,wfir1,vmat)
    h(:,:)=h(:,:)+vmat(:,:)
    deallocate(bmt,bir)
  end if
else
! normal Hartree-Fock
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
       vclmt(:,:,ias),ld,0.d0,vmt(:,:,ias),lmmaxvr)
    end do
  end do
  call genvmatk(vmt,vclir,wfmt1,wfir1,vmat)
  h(:,:)=h(:,:)+vmat(:,:)
end if
deallocate(vmt)
!-------------------------------------------!
!     non-local Coulomb matrix elements     !
!-------------------------------------------!
vmat(:,:)=0.d0
! start loop over non-reduced k-point set
do ik=1,nkptnr
! find the equivalent reduced k-point
  iv(:)=ivk(:,ik)
  jk=ikmap(iv(1),iv(2),iv(3))
! generate the G+k vectors
  call gengpvec(vkl(:,ik),vkc(:,ik),ngknr,igkignr,vgklnr,vgkcnr)
! generate the spherical coordinates of the G+k vectors
  do igk=1,ngknr
    call sphcrd(vgkcnr(:,igk),gkcnr(igk),tpgkcnr(:,igk))
  end do
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! get the eigenvalues/vectors from file for non-reduced k-point
  call getevalsv(vkl(:,ik),evalsvnr)
  call getevecfv(vkl(:,ik),vgklnr,evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! determine q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ngvec
! determine G+q vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vector
    call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  lmax=lmaxvr+npsden+1
  call genjlgpr(lmax,gqc,jlgqr)
  call genjlgq0r(gqc(igq0),jlgq0r)
! calculate the wavefunctions for all states
  call genwfsv(.false.,.false.,.false.,ngknr,igkignr,evalsvnr,apwalm,evecfv, &
   evecsv,wfmt2,ngrtot,wfir2)
  do ist3=1,nstsv
    if (occsv(ist3,jk).gt.epsocc) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zvclmt,zvclir) &
!$OMP PRIVATE(zrho01,zrho02,ist1,zt1,zt2,t1)
!$OMP DO
      do ist2=1,nstsv
        allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngrtot))
        allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot),zvclir(ngrtot))
! calculate the complex overlap density
        call genzrho(.true.,wfmt2(:,:,:,:,ist3),wfmt1(:,:,:,:,ist2), &
         wfir2(:,:,ist3),wfir1(:,:,ist2),zrhomt,zrhoir)
! calculate the Coulomb potential
        call genzvclmt(nrcmt,nrcmtmax,rcmt,nrcmtmax,zrhomt,zvclmt)
        call zpotcoul(nrcmt,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq,zrhoir, &
         nrcmtmax,zvclmt,zvclir,zrho02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
        do ist1=1,ist2
! calculate the complex overlap density
          call genzrho(.true.,wfmt2(:,:,:,:,ist3),wfmt1(:,:,:,:,ist1), &
           wfir2(:,:,ist3),wfir1(:,:,ist1),zrhomt,zrhoir)
          zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
! compute the density coefficient of the smallest G+q-vector
          call zrhogp(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt,zrhoir,zrho01)
          zt2=cfq*wiq2(iq)*(conjg(zrho01)*zrho02)
          t1=occsv(ist3,jk)/occmax
!$OMP CRITICAL
          vmat(ist1,ist2)=vmat(ist1,ist2)-t1*(wkptnr*zt1+zt2)
!$OMP END CRITICAL
        end do
        deallocate(zrhomt,zrhoir,zvclmt,zvclir)
      end do
!$OMP END DO
!$OMP END PARALLEL
    end if
  end do
! end loop over non-reduced k-point set
end do
! scale the non-local matrix elements in the case of a hybrid functional
if (hybmix.lt.1.d0) then
  vmat(:,:)=hybmix*vmat(:,:)
end if
! add the non-local matrix elements to Hamiltonian
h(:,:)=h(:,:)+vmat(:,:)
! diagonalise the Hartree-Fock Hamiltonian (eigenvalues in global array)
call zheev('V','U',nstsv,h,nstsv,evalsv(:,ikp),work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(seceqnhf): diagonalisation of the Hartree-Fock Hamiltonian &
   &failed")')
  write(*,'(" for k-point ",I8)') ikp
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
! apply unitary transformation to second-variational states
evecsv(:,:)=evecsvp(:,:)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,h,nstsv,zzero,evecsvp, &
 nstsv)
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,tpgqc,gqc,jlgqr,jlgq0r)
deallocate(evalsvp,evalsvnr,evecfv,evecsv,rwork)
deallocate(h,vmat,apwalm,sfacgknr,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,zfmt,work)
return
end subroutine

