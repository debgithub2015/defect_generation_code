!PROGRAMMER ::Dilpuneet S Aidhy
!revised version in f90: 12/16/11, IBM India
!Sitebysite analysis to locate point defects in fluroite and rocksalt lattice
! the code locates vacancies, interstitials of host and doped ions

Module hell_variables

IMPLICIT NONE
INTEGER::idcr,ipbc,lsurf,i_v_save,ltype1,ltype2,ipot0,ieam,i3b,in
REAL::amass1,amass2,rsmax0,real_time,PotEnrg
INTEGER::nntypes,islab,ipoly,ipolysize,ndoub1,ndoub2,noub3,iatoms
REAL,dimension(3,3)::cbeta
REAL,dimension(2,2)::rnn
REAL::hx,hy,hz,shx,shy,shz
INTEGER,dimension(2,2)::lnn
INTEGER,dimension(64000)::igrain,igo,ithermo
REAL,dimension(64000)::x5,y5,z5,sx5,sy5,sz5

END Module


Module global_variables

IMPLICIT NONE
REAL,dimension(64000)::rx,ry,rz,charge,sx,sy,sz,scharge
INTEGER,dimension(64000)::imass,ident,ic,simass,itick_vac,itick_int
INTEGER,dimension(64000)::isubs
INTEGER::natoms,nnatoms,i
REAL::acharge, bcharge, rrcut
REAL :: lp
INTEGER::iU,iO,iDy,iUn,iDyn,iOn,intU,intY,intO,ivacU,ivacO,ivacDy

END Module

Program compare
USE global_variables
USE hell_variables
IMPLICIT NONE
!Write(*,*)'enter latice paramter'
!READ(*,*)lp
lp=3.52
WRITE(*,*)'Lattice parameter is=',lp

! Calling the files with no defect

!Call ReadStructure_lammps_atomformat
!Call ReadStructure_lammps_atomformat_tempfile
Call ReadStructure_lammps_dumpformat_tempfile

! Calling the files with defects

!Call ReadStructure_lammps_atomformat_strfile
Call ReadStructure_lammps_dumpformat_strfile

!Call ReadStructure_hell

Call InputParameters
Call CountVacancies
Call CountInterstitials
!Call CountSubstitution
Call OutputFiles

End program compare


!**********************

Subroutine ReadStructure_lammps_atomformat_tempfile
!Initial file with no defects

Use hell_variables
Use global_variables
Implicit None
REAL::j,hxo,hyo,hzo
CHARACTER::a,d
OPEN (UNIT=1,STATUS='OLD',FILE='quenched.temp.noFP')
READ(1,*)a
READ(1,*)natoms,a
READ(1,*)j,a
READ(1,*)hxo,hx,a,a
READ(1,*)hyo,hy,a,a
READ(1,*)hzo,hz,a,a
READ(1,*)d
hx = hx - hxo
hy = hy - hyo
hz = hy - hzo
WRITE(*,*)hx,hy,hz
hx=hx/lp
hy=hy/lp
hz=hz/lp
WRITE(*,*)natoms
Do i=1,natoms
 READ(1,*)ident(i),imass(i),charge(i),rx(i),ry(i),rz(i)
rx(i) = rx(i)/lp
ry(i) = ry(i)/lp
rz(i) = rz(i)/lp
ENDDO

End Subroutine ReadStructure_lammps_atomformat_tempfile
!**********************

Subroutine ReadStructure_lammps_dumpformat_tempfile
!Initial file with no defects
!Files in dump format

Use hell_variables
Use global_variables
Implicit None
REAL::j,hxo,hyo,hzo
CHARACTER::a,d
OPEN (UNIT=1,STATUS='OLD',FILE='quenched.temp.noFP')
READ(1,*)a
READ(1,*)j
READ(1,*)a
READ(1,*)natoms
READ(1,*)a
READ(1,*)hxo,hx
READ(1,*)hyo,hy
READ(1,*)hzo,hz
READ(1,*)d
hx = hx - hxo
hy = hy - hyo
hz = hz - hzo
hx=hx/lp
hy=hy/lp
hz=hz/lp
WRITE(*,*)natoms
Do i=1,natoms
!READ(1,*)ident(i),imass(i),charge(i),rx(i),ry(i),rz(i)
READ(1,*)ident(i),imass(i),rx(i),ry(i),rz(i)

!rx(i) = rx(i)/lp
!ry(i) = ry(i)/lp
!rz(i) = rz(i)/lp

rx(i) = rx(i)*hx
ry(i)=ry(i)*hy
rz(i)=rz(i)*hz

ENDDO
ENd Subroutine ReadStructure_lammps_dumpformat_tempfile

!**********************

Subroutine ReadStructure_lammps_atomformat_strfile
!MD time progressed file

Use hell_variables
Use global_variables
Implicit None
REAL::j,shxo,shyo,shzo
CHARACTER::a,d

OPEN (UNIT=2,STATUS='OLD',FILE='quenched.str')
READ(2,*)a
READ(2,*)nnatoms, a
READ(2,*)j,a
READ(2,*)shxo,shx,a,a
READ(2,*)shyo,shy,a,a
READ(2,*)shzo,shz,a,a
READ(2,*)d
shx = shx - shxo
shy = shy - shyo
shz = shz - shzo
shx=shx/lp
shy=shy/lp
shz=shz/lp
Do i=1,nnatoms
READ(2,*)ic(i),simass(i),scharge(i),sx(i),sy(i),sz(i)
sx(i) = sx(i)/lp
sy(i) = sy(i)/lp
sz(i) = sz(i)/lp
ENDDO

Do i=1,nnatoms
sx(i) = sx(i)*hx/shx
sy(i) = sy(i)*hy/shy
sz(i) = sz(i)*hz/shz
ENDDO

WRITE(*,*)'structure reading done'
WRITE(*,*)' '
END Subroutine ReadStructure_lammps_atomformat_strfile

!***************

Subroutine ReadStructure_lammps_dumpformat_strfile
!Files in dump format
Use hell_variables
Use global_variables
Implicit None
REAL::j,shxo,shyo,shzo
CHARACTER::a,d

OPEN (UNIT=2,STATUS='OLD',FILE='quenched.str')
READ(2,*)a
READ(2,*)j
READ(2,*)a
READ(2,*)nnatoms
READ(2,*)a
READ(2,*)shxo,shx
READ(2,*)shyo,shy
READ(2,*)shzo,shz
READ(2,*)a
shx = shx - shxo
shy = shy - shyo
shz = shz - shzo
WRITE(*,*)shx,shy,shz

shx=shx/lp
shy=shy/lp
shz=shz/lp
Do i=1,nnatoms
!READ(2,*)ic(i),simass(i),scharge(i),sx(i),sy(i),sz(i)
READ(2,*)ic(i),simass(i),sx(i),sy(i),sz(i)
ENDDO

Do i=1,nnatoms
!sx(i) = sx(i)/lp
!sy(i) = sy(i)/lp
!sz(i) = sz(i)/lp

sx(i) = sx(i)*shx
sy(i) = sy(i)*shy
sz(i)=sz(i)*shz
sx(i) = sx(i)*hx/shx
sy(i) = sy(i)*hy/shy
sz(i) = sz(i)*hz/shz
ENDDO

WRITE(*,*)'structure reading done'
WRITE(*,*)' '

END Subroutine ReadStructure_lammps_dumpformat_strfile

!***************

subroutine InputParameters
Use global_variables
Implicit None

WRITE(*,*)'enter rcut - cutoff for identifying vac and int.'
!WRITE(*,*)'for fluorite = 0.19'
!write(*,*)'for MgO = 0.25'
write(*,*) ' '
read(*,*)rrcut
!rrcut = 0.23
!rrcut  = rrcut /lp   !Lammps output coordinates are in fractional coordinates.
                     !so the rrcut is divided by lp

DO i = 1, natoms
 IF (imass(i).eq.1) iU = iU + 1
 IF (imass(i).eq.2) iO = iO + 1
 IF (imass(i).eq.3) iDy = iDy + 1
ENDDO
DO i = 1, nnatoms
 IF (simass(i).eq.1) iUn = iUn + 1
 IF (simass(i).eq.3) iDyn = iDyn + 1
 IF (simass(i).eq.2) iOn = iOn + 1
ENDDO
WRITE(*,*)'In the MD defected file, there are'
WRITE(*,*)'         cations    subs     anions'
WRITE(*,*)iUn,iDyn,iOn

END Subroutine InputParameters

!*******************

Subroutine CountVacancies 
Use global_variables
Use hell_variables
Implicit None
INTEGER::icounter,ki,j,k
REAL::dx,dy,dz,rrr
REAL,dimension(64000)::repeatvac
icounter=0

Do ki=1,natoms
itick_vac(ki) = 0
ENDDO
write(*,*)'counting vacancies'
DO ki=1,natoms
 IF (mod(ki,20000).eq.0) WRITE(*,*) 'Atoms counted',ki
  DO j= 1,nnatoms
!   DO k=1,icounter
!    If (repeatvac(k).eq.j) THEN
!     go to 103
!    ENDIF
!   ENDDO
   dx  =  0.0d0 ; dy  =  0.0 ;   dz = 0.0d0 ;
   dx  = rx(ki) - sx(j)
   dy  = ry(ki) - sy(j)
   dz  = rz(ki) - sz(j)
   IF  (dx.GE.hx/2)      dx = dx - hx
   IF  (dx.LT.-hx/2)     dx = dx + hx
   IF  (dy.GE.hy/2)      dy = dy - hy
   IF  (dy.LT.-hy/2)     dy = dy + hy
   IF  (dz.GE.hz/2)      dz = dz - hz
   IF  (dz.LT.-hz/2)     dz = dz + hz
   rrr = SQRT(dx*dx+dy*dy+dz*dz)
    IF (rrr.le.rrcut) THEN
     itick_vac(ki) = 1
     icounter = icounter + 1
     repeatvac (icounter) = j
     IF (imass(ki).eq.1) ivacU = ivacU+1
     IF (imass(ki).eq.2) ivacO = ivacO+1
    IF (imass(ki).eq.3) ivacDy = ivacDy+1
     go to 100
    ENDIF
103 ENDDO
100 ENDDO
WRITE(*,*)'Number of Vacancies found'
WRITE(*,*)'Cation vacancies:',iU-ivacU
WRITE(*,*)'Oxygen vacancies:',iO-ivacO
WRITE(*,*)'Imass 3 vacancies', iDy-ivacDy
!WRITE(*,*)'Lattice atoms : ', 'cation =',ivacU,'Oxygen',ivacO
WRITE(*,*)' '

End Subroutine CountVacancies
!********************

Subroutine CountInterstitials
Use global_variables 
Use hell_variables
IMPLICIT NONE
INTEGER::j,ki,k,icounter
REAL::dx,dy,dz,rrr
INTEGER,dimension(64000)::repeatint

Do ki=1,natoms
itick_int(ki) = 0
ENDDO

DO j=1,nnatoms
 IF (mod(j,20000).eq.0) WRITE(*,*) 'Atoms counted', j
 DO ki=1,natoms
!  DO k=1,icounter
!   If (repeatint(k).eq.ki) THEN
!       go to 104
!  ENDIF
!  ENDDO
  dx  =  0.0d0 ; dy  =  0.0 ;   dz = 0.0d0 ;
  dx  = sx(j) - rx(ki)
  dy  = sy(j) - ry(ki)
  dz  = sz(j) - rz(ki)
  IF  (dx.GE.shx/2)      dx = dx - shx
  IF  (dx.LT.-shx/2)     dx = dx + shx
  IF  (dy.GE.shy/2)      dy = dy - shy
  IF  (dy.LT.-shy/2)     dy = dy + shy
  IF  (dz.GE.shz/2)      dz = dz - shz
  IF  (dz.LT.-shz/2)     dz = dz + shz
  rrr = SQRT(dx*dx + dy*dy + dz*dz)
  IF (rrr.le.rrcut) THEN
   itick_int(j)=1
   icounter = icounter + 1
   repeatint (icounter) = ki
   IF (simass(j).eq.1) intU = intU + 1
   IF (simass(j).eq.3) intY = intY + 1
   IF (simass(j).eq.2) intO = intO + 1
  
   go to 200

  ENDIF
104 ENDDO
200 ENDDO

WRITE(*,*)'Number of Interstitials found '
WRITE(*,*)'Cation4+ Interstitials:',iUn-intU
WRITE(*,*)'Cation3+ Interstitials:',iDyn-intY
WRITE(*,*)'Oxygen Interstitials:',iOn-intO
WRITE(*,*)' '
END subroutine CountInterstitials

!*****************************


Subroutine Countsubstitution
Use global_variables 
Use hell_variables
IMPLICIT NONE
INTEGER::ki,j,subs,isubU,icounter,k
REAL::rrr,dx,dy,dz
INTEGER,dimension(64000)::repeatsubs
DO ki=1,natoms
 IF (mod(ki,20000).eq.0) WRITE(*,*) 'Atoms counted',ki
 DO j=1,nnatoms
!  DO k=1,icounter
!   If (repeatsubs(k).eq.j) THEN
!    go to 301
!   ENDIF
!  ENDDO
  dx  =  0.0d0 ; dy  =  0.0 ;   dz = 0.0d0 ;
  dx  = rx(ki) - sx(j)
  dy  = ry(ki) - sy(j)
  dz  = rz(ki) - sz(j)
   IF  (dx.GE.hx/2)      dx = dx - hx
   IF  (dx.LT.-hx/2)     dx = dx + hx
   IF  (dy.GE.hy/2)      dy = dy - hy
   IF  (dy.LT.-hy/2)     dy = dy + hy
   IF  (dz.GE.hz/2)      dz = dz - hz
   IF  (dz.LT.-hz/2)     dz = dz + hz
  rrr = SQRT(dx*dx+dy*dy+dz*dz)
  IF(rrr.le.rrcut) THEN
    IF(imass(ki).ne.simass(j))THEN 
     isubs(ki) = 1
     subs = subs + 1
     icounter = icounter + 1
     repeatsubs (icounter) = j
!     WRITE(*,*)ki,charge(ki),scharge(j)
!     IF(imass(ki).ne.acharge) isubU = isubU + 1
     go to 300
    ENDIF
   ENDIF
301  ENDDO
300 ENDDO
WRITE(*,*)' '
WRITE(*,*)'Number of Substitution ions found '
WRITE(*,*)'Cation substitutions:',subs
WRITE(*,*)' '
End Subroutine Countsubstitution

!***********************

Subroutine OutputFiles
Use global_variables
Use hell_variables 
IMPLICIT NONE
INTEGER::ki,j,icount
OPEN (UNIT=3,STATUS='unknown',FILE='defects.xyz')
WRITE(*,*)'lattice parameter', lp
DO ki=1,natoms
  IF (itick_vac(ki).eq.0) THEN
    IF (imass (ki).eq.1) WRITE(3,101)'Ce',rx(ki)*lp,ry(ki)*lp,rz(ki)*lp   ! U vacancies
    IF (imass (ki).eq.2) WRITE(3,101)'Bi',rx(ki)*lp,ry(ki)*lp,rz(ki)*lp    ! O vacancies
    IF (imass (ki).eq.3) WRITE(3,101)'Ir',rx(ki)*lp,ry(ki)*lp,rz(ki)*lp    ! imass 3 vacancies
 ENDIF

  IF (itick_int(ki).eq.0)THEN
!    IF(simass(ki).eq.1) WRITE(*,*)'identi',ident(ki)
    IF (simass(ki).eq.1) WRITE(3,101)'Pb',sx(ki)*lp,sy(ki)*lp,sz(ki)*lp    ! U interstitials
    IF (simass(ki).eq.3) WRITE(3,101)'Dy',sx(ki)*lp,sy(ki)*lp,sz(ki)*lp    ! Ce3+ interstitials
    IF (simass(ki).eq.2) WRITE(3,101)'Se',sx(ki)*lp,sy(ki)*lp,sz(ki)*lp     ! O interstitials
  ENDIF

  IF (itick_int(ki).eq.1)THEN
!   IF(isubs(ki).eq.0)THEN
    IF (simass (ki).eq.1) WRITE(3,101)'Si',sx(ki)*lp,sy(ki)*lp,sz(ki)*lp       ! U lattice atoms
    IF (simass (ki).eq.2) WRITE(3,101)'O',sx(ki)*lp,sy(ki)*lp,sz(ki)*lp      ! O lattice atoms
!  ELSEIF(isubs(ki).eq.1)THEN
 !    WRITE(3,101)'Yb',rx(ki)*lp,ry(ki)*lp,rz(ki)*lp   ! Ce3+ is substitution for Ce4+
   ENDIF
! ENDIF

IF (simass(ki).eq.3) WRITE(3,101) 'Yb',sx(ki)*lp,sy(ki)*lp,sz(ki)*lp

! IF(isubs(ki).eq.1) THEN
!   icount = icount+1
!   IF (imass(ki).eq.simass(ki)) WRITE(3,101)'Y',rx(ki),ry(ki),rz(ki)   ! Ce3+ is Y
! ENDIF
ENDDO

101   FORMAT(a,3(f12.5))

END subroutine OutputFiles


!**********************

Subroutine ReadStructure_hell
Use hell_variables
Use global_variables
Implicit None
INTEGER::dum
REAL::dum1
OPEN (UNIT=1,STATUS='OLD',FILE='quenched.temp.noFP')
READ(1,*)dum,dum,dum,dum,dum
READ(1,*)dum,dum,dum,dum,dum
READ(1,*)natoms
READ(1,*)dum1, dum1
READ(1,*)dum1,dum1
READ(1,*)dum,dum,dum
READ(1,*)dum,dum
READ(1,*)dum1,dum1,dum1
READ(1,*)dum1,dum1,dum1
READ(1,*)dum1,dum1,dum1
READ(1,*)hx,dum1,dum1
READ(1,*)dum1,hy,dum1
READ(1,*)dum1,dum1,hz
READ(1,*)dum,dum,dum
READ(1,*)dum1
READ(1,*)dum1,dum1
READ(1,*)dum1,dum1
READ(1,*)dum, dum
READ(1,*)dum,dum
WRITE(*,*)natoms
Do i=1,natoms
 READ(1,*)ident(i),rx(i),ry(i),rz(i),x5(i),y5(i),z5(i),charge(i),imass(i),dum,dum,dum
ENDDO

OPEN (UNIT=2,STATUS='OLD',FILE='quenched.str')
READ(2,*)dum,dum,dum,dum,dum
READ(2,*)dum,dum,dum,dum,dum
READ(2,*)nnatoms
READ(2,*)dum1, dum1
READ(2,*)dum1,dum1
READ(2,*)dum,dum,dum
READ(2,*)dum,dum
READ(2,*)dum1,dum1,dum1
READ(2,*)dum1,dum1,dum1
READ(2,*)dum1,dum1,dum1
READ(2,*)shx,dum1,dum1
READ(2,*)dum1,shy,dum1
READ(2,*)dum1,dum1,shz
READ(2,*)dum,dum,dum
READ(2,*)dum1
READ(2,*)dum1,dum1
READ(2,*)dum1,dum1
READ(2,*)dum, dum
READ(2,*)dum,dum

WRITE(*,*)nnatoms
Do i=1,nnatoms
 READ(2,*)ic(i),sx(i),sy(i),sz(i),sx5(i),sy5(i),sz5(i),scharge(i),simass(i),dum,dum,dum
ENDDO

Do i=1,nnatoms
sx(i) = sx(i)*hx/shx
sy(i) = sy(i)*hy/shy
sz(i) = sz(i)*hz/shz
ENDDO

WRITE(*,*)'structur reading done'

END Subroutine ReadStructure_hell

!***************
