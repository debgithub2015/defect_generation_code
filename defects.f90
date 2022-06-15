!Modified March 2013
!Modified April 1, 2014 - for adding Frenkel pairs in SrTiO3 for Eva
!Modified May 5, 2014 for creating Gd2Ti2O7 from fluorite

MODULE variables

IMPLICIT NONE
REAL,dimension(7000000)::rx,ry,rz,charge,x5,y5,z5,x6,y6,z6
INTEGER,dimension(7000000)::imass,ident,igrain,igo,ithermo
INTEGER,dimension(7000000)::ivac,ireplace
INTEGER::natoms,i,zatom, coortype,typeatoms
REAL,dimension(3,3)::h
REAL::lp1,lp2,lp3,hx,hxo,hy,hyo,hz,hzo,xnum,ynum,znum

CHARACTER*2:: atomname (3)

!coortype :for dlpoly - indicates if velocites and forces are present in the structure file
!atomname used on Dl-poly structure read-write
!Check cation and anion charge in Create_void subroutine for your system
! Check Lattice parameter, lp for your system


END Module

Module str

IMPLICIT NONE
INTEGER::idcr,ipbc,lsurf,i_v_save,ltype1,ltype2,ipot0,ieam,i3b,ion
REAL::amass1,amass2,rsmax0,real_time,PotEnrg, amass3
INTEGER::iatoms,nntypes,islab,ipoly,ipolysize,ndoub1,ndoub2,noub3
REAL,dimension(3,3)::cbeta
REAL,dimension(2,2)::rnn
INTEGER,dimension(2,2)::lnn

END Module



Program Defects
USE variables
IMPLICIT NONE
INTEGER::option

WRITE(*,*) 'Give the size of the box at all 3 dimensions'
READ(*,*) xnum,ynum,znum

!WRITE(*,*)'enter lattice parameter in ANGSTROMS'
!READ(*,*)lp
!lp = 5.56                ! Tho2 BD08 Behera potential
!lp = 5.411               ! CeO2 lattice parameter (Gotter07 potentiial)
!lp = 3.918               ! SrTiO3 lattice parameter (Sekiguchi potenital)
!lp=5.41
!lp=10.185
!lp = 3.52    ! Ni FCC
!WRITE(*,*)' '
!WRITE(*,*)'CHECK IF LP IS CORRECT', lp
!WRITE(*,*)' '
!Call ReadStructure_hell
Call ReadStructure_lammps_atomformat
!CAll ReadStructure_lammps_dumpformat
 
!Call ReadStructure_Dlpoly_NoVelFors

!CAll ReadStructure_Dlpoly_VelFors

Call random_seed

Write(*,*)'Choose what would you like to do'
WRITE(*,*)'Substitute CeO2 with Ce2O3       - type 1'
WRITE(*,*)'Add charge-neutral interstitials - type 2'
WRITE(*,*)'Create Frenkel Pairs in Rocksalt - type 3'
WRITE(*,*)'Create charge-neutral voids      - type 4'
WRITE(*,*)'Create defects in SrTiO3         - type 5'
!STO subroutine can be used for Pyrochlore as well
WRITE(*,*)'Create Schottky defects in CeO2  - type 6'
WRITE(*,*)'Create Frenkel pairs in SrTiO3   - type 7 '
WRITE(*,*)'Convert ZrO2-fluorite to Ge2Ti2O7 -type 8'
WRITE(*,*)'Convert Pyrochlore to Defect-fluorite - type 9'
WRITE(*,*)'Convert Metal to Alloy           - type 10'
WRITE(*,*)'Create vacancies in metal or alloy-type 11'
WRITE(*,*)'Create interstitials in metal or alloy-type 12'
WRITE(*,*)'Create void in FCC Ni or metal      - type 13' 
WRITE(*,*)'Create ZrO2 box inside CeO2 box     - type 14'
 
READ(*,*)option
!option=9
If (option == 1) Call  substitute
IF (option == 2) Call  addinterstitials
IF (option == 3) Call  FProcksalt
IF (option == 4) Call  Create_void
IF (option == 5) Call  SrTiO3_defects
IF (option == 6) Call  CeO2_Schottky
IF (option == 7) Call  STO_FPs
IF (option == 8) Call  ZrO2_Gd2Ti2O7
IF (option == 9) Call  PyrotoDefectfluo
IF (option == 10) Call random_alloy
IF (option == 11) Call random_alloy_vacancies
IF (option == 12) Call random_alloy_interstitials
IF (option == 13) Call Create_metalvoid
IF (option == 14) Call Create_zro2box

!Call Overlapping_atoms
Call WriteStructure

END PROGRAM Defects

!***************************


Subroutine ReadStructure_lammps_atomformat
! Atoms structure format
!based on LAmmps input file for SrTiO3, Alex-UF format
Use variables
Use str
Implicit None
REAL::j,dum1 
INTEGER::atom_type
CHARACTER::a,d
!OPEN (UNIT=1,STATUS='OLD',FILE='lammpsstrfile')
OPEN (UNIT=1,STATUS='OLD',FILE='/project/RD-HEA/dchakrab/initial_code/create-defects-code-debajit/lammpsstrfile')
READ(1,*)a
READ(1,*)natoms,a
READ(1,*)atom_type,a
READ(1,*)hxo,hx, a
READ(1,*)hyo,hy, a
READ(1,*)hzo,hz, a
READ(1,*)d
lp1=(hx-hxo)/xnum
lp2=(hy-hyo)/ynum
lp3=(hz-hzo)/znum
write(*,*) 'the lattice parameters are', lp1,lp2,lp3
h(1,1)=(hx-hxo)/lp1
h(2,2)=(hy-hyo)/lp2
h(3,3)=(hz-hzo)/lp3
Do i=1,natoms

 READ(1,*)ident(i),imass(i),charge(i),rx(i),ry(i),rz(i)

rx(i) = (rx(i)-hxo)/(hx-hxo)
ry(i) = (ry(i)-hyo)/(hy-hyo)
rz(i) = (rz(i)-hzo)/(hz-hzo)

ENDDO

WRITE(*,*)hxo,hx,hyo,hy,hzo,hz
WRITE(*,*)h(1,1),h(2,2),h(3,3)
END Subroutine ReadStructure_lammps_atomformat

!***************************

Subroutine ReadStructure_lammps_dumpformat
!based on Dump file format used initially for Sigma5 GB simulations to restart sim from output file
Use variables
Use str
Implicit None
REAL::j,dum1
INTEGER::atom_type,ab
CHARACTER::a,d
OPEN (UNIT=1,STATUS='OLD',FILE='lammpsstrfile')
READ(1,*)a
READ(1,*)ab
READ(1,*)a
READ(1,*)natoms
READ(1,*)a
READ(1,*)hxo,hx
READ(1,*)hyo,hy
READ(1,*)hzo,hz
READ(1,*)a
lp1=(hx-hxo)/xnum
lp2=(hy-hyo)/ynum
lp3=(hz-hzo)/znum
h(1,1)=(hx-hxo)/lp1
h(2,2)=(hy-hyo)/lp2
h(3,3)=(hz-hzo)/lp3

Do i=1,natoms

 READ(1,*)ident(i),imass(i),charge(i),rx(i),ry(i),rz(i)

rx(i) = rx(i)
ry(i) = ry(i)
rz(i) = rz(i)

ENDDO

END Subroutine ReadStructure_lammps_dumpformat



!***************************

Subroutine ReadStructure_hell

Use variables
Use str
IMPLICIT NONE
INTEGER:: dum
REAL::dum1

OPEN (UNIT=1,STATUS='OLD',FILE='hellstrfile')
READ(1,*) idcr, ipbc, lsurf, i_v_save, ltype1
REAd(1,*) ltype2, ipot0, ieam, i3b, ion
READ(1,*) natoms
READ(1,*) amass1, amass2
READ(1,*) rsmax0,  real_time
READ(1,*) iatoms, nntypes, islab
READ(1,*) ipoly, ipolysize
READ(1,*) cbeta(1,1), cbeta(2,1), cbeta(3,1)
READ(1,*) cbeta(1,2), cbeta(2,2), cbeta(3,2)
READ(1,*) cbeta(1,3), cbeta(2,3), cbeta(3,3)
READ(1,*) h(1,1), h(2,1), h(3,1)
READ(1,*) h(1,2), h(2,2), h(3,2)
READ(1,*) h(1,3), h(2,3), h(3,3)
READ(1,*) ndoub1,ndoub2,noub3
READ(1,*) PotEnrg
READ(1,*) rnn(1,1), rnn(2,1)
READ(1,*) rnn(1,2), rnn(2,2)
READ(1,*) lnn(1,1), lnn(2,1)
READ(1,*) lnn(1,2), lnn(2,2)

WRITE(*,*)natoms
Do i=1,natoms
 READ(1,*)ident(i),rx(i),ry(i),rz(i),x5(i),y5(i),z5(i),&
 charge(i),imass(i),igrain(i),igo(i),ithermo(i)
ENDDO

zatom = natoms
END Subroutine ReadStructure_hell

!******************************

Subroutine ReadStructure_Dlpoly_NoVelFors
!based on DL_Poly format, the structure that has no Velocities and forces
Use variables
Use str
Implicit None
REAL::j,dum1
!INTEGER::typeatoms
CHARACTER::a
!CHARACTER*2:: atomname (3)
OPEN (UNIT=1,STATUS='OLD',FILE='dl-polystr')
READ(1,*)a
READ(1,*)coortype,typeatoms
READ(1,*) h(1,1), h(2,1), h(3,1)
READ(1,*) h(1,2), h(2,2), h(3,2)
READ(1,*) h(1,3), h(2,3), h(3,3)

Write(*,*)'enter total number of atoms'
READ(*,*)natoms

Do i=1,natoms
  READ(1,*)atomname(i),ident(i)
  READ(1,*)rx(i),ry(i),rz(i)
  If(atomname(i) == 'Sr') imass(i) = 1
  If(atomname(i) == 'Ti') imass(i) = 3
  If(atomname(i) == 'O') imass(i) = 2
!  WRITE(*,*)imass(i)
rx(i) = rx(i)/lp1
ry(i) = ry(i)/lp2
rz(i) = rz(i)/lp3

ENDDO

Do i=1,natoms
!WRITE(*,*)atomname(i),ident(i)
!WRITE(*,*)rx(i)*lp,ry(i)*lp,rz(i)*lp
ENDDO

END Subroutine ReadStructure_Dlpoly_NoVelFors

!******************************
Subroutine ReadStructure_Dlpoly_VelFors
!based on DL_Poly format, the structure that has no Velocities and forces
Use variables
Use str
Implicit None
REAL::dum
INTEGER::mdsteps
CHARACTER::a
!CHARACTER*2:: atomname (3)
OPEN (UNIT=1,STATUS='OLD',FILE='dl-polystr')
READ(1,*)a
READ(1,*)coortype,typeatoms,mdsteps,dum,dum
READ(1,*) h(1,1), h(2,1), h(3,1)
READ(1,*) h(1,2), h(2,2), h(3,2)
READ(1,*) h(1,3), h(2,3), h(3,3)

Write(*,*)'enter total number of atoms'
READ(*,*)natoms

Do i=1,natoms
  READ(1,*)atomname(i),ident(i)
  READ(1,*)rx(i),ry(i),rz(i)
  READ(1,*)x5(i),y5(i),z5(i)
  READ(1,*)x6(i),y6(i),z6(i)
  If(atomname(i) == 'Sr') imass(i) = 1
  If(atomname(i) == 'Ti') imass(i) = 3
  If(atomname(i) == 'O') imass(i) = 2

rx(i) = rx(i)/lp1
ry(i) = ry(i)/lp2
rz(i) = rz(i)/lp3

ENDDO

!Do i=1,natoms
!WRITE(*,*)atomname(i),ident(i)
!WRITE(*,*)rx(i)*lp,ry(i)*lp,rz(i)*lp

!ENDDO

END Subroutine ReadStructure_Dlpoly_VelFors
!********************************************

Subroutine WriteStructure
USE variables
USE str
IMPLICIT NONE
Integer::ik

!OPEN (UNIT=2,STATUS='unknown',FILE='defect.str')
OPEN (UNIT=3,STATUS='unknown',FILE='defect-lammps-atomformat.str')
!OPEN (UNIT=33,STATUS='unknown',FILE='defect-lammps-dumpformat.str')

!OPEN (UNIT=34,STATUS='unknown',FILE='defect-dlpoly.str')
!OPEN (UNIT=35,STATUS='unknown',FILE='gulp.str')
OPEN (UNIT=4,STATUS='unknown',FILE='defect.xyz')



!WRITE(2,*) idcr, ipbc, lsurf, i_v_save, ltype1
!WRITE(2,*) ltype2, ipot0, ieam, i3b, ion
!WRITE(2,*) zatom
!WRITE(2,*) amass1, amass2
!WRITE(2,*) rsmax0,  real_time
!WRITE(2,*) iatoms, nntypes, islab
!WRITE(2,*) ipoly, ipolysize
!WRITE(2,*) cbeta(1,1), cbeta(2,1), cbeta(3,1)
!WRITE(2,*) cbeta(1,2), cbeta(2,2), cbeta(3,2)
!WRITE(2,*) cbeta(1,3), cbeta(2,3), cbeta(3,3)
!WRITE(2,*) h(1,1), h(2,1), h(3,1)
!WRITE(2,*) h(1,2), h(2,2), h(3,2)
!WRITE(2,*) h(1,3), h(2,3), h(3,3)
!WRITE(2,*) ndoub1,ndoub2,noub3
!WRITE(2,*) PotEnrg
!WRITE(2,*) rnn(1,1), rnn(2,1)
!WRITE(2,*) rnn(1,2), rnn(2,2)
!WRITE(2,*) lnn(1,1), lnn(2,1)
!WRITE(2,*) lnn(1,2), lnn(2,2)

WRITE(4,*) natoms
WRITE(4,*) ' '
Do i=1,natoms
 IF (ivac(i) == 0) THEN
! WRITE(2,101)ident(i),rx(i),ry(i),rz(i),x5(i),y5(i),z5(i),&
! charge(i),imass(i),igrain(i),igo(i),ithermo(i)

!DC modifies according to the system of interest 07/12/2016

  If(imass(i)==1) Write(4,*)'Ni', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo  !Host Lattice(Ni)
!  If(imass(i)==1) Write(4,*)'Al', rx(i)*lp,ry(i)*lp,rz(i)*lp  !Cation lattice atom

  If(imass(i)==2) Write(4,*)'Fe', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo   ! 1st Alyloy component(Fe)
!  If(imass(i)==2) Write(4,*)'O', rx(i)*lp,ry(i)*lp,rz(i)*lp   ! Oxygen lattice atom

  IF(imass(i)==3) Write(4,*)'Cr', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo  ! 2nd Alloy component(Cr)
!  IF(imass(i)==3) Write(4,*)'Pt', rx(i)*lp,ry(i)*lp,rz(i)*lp  ! Cation interstitial

  IF(imass(i)==4) Write(4,*)'Pd', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo  ! 3rd Alloy Component(Pd)
!  IF(imass(i)==4) Write(4,*)'Si', rx(i)*lp,ry(i)*lp,rz(i)*lp  ! for Oxygen interstitial
 
  IF(imass(i)==5) Write(4,*)'Co', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo  ! 4th Alloy Component(Co)

  IF(imass(i)==6) Write(4,*)'Mn', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo  ! 5th Alloy Component(Mn)

 ELSEIF (ivac(i) == 1) THEN
  If(imass(i)==1) Write(4,*)'Sr', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo ! Ni-vacancy
!  If(imass(i)==1) Write(4,*)'Sr', rx(i)*lp,ry(i)*lp,rz(i)*lp ! imass = 1 cation vacancy

  If(imass(i)==2) Write(4,*)'Er', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo  ! Fe-vacancy
!  If(imass(i)==2) Write(4,*)'Er', rx(i)*lp,ry(i)*lp,rz(i)*lp  ! O vacancy

  IF(imass(i)==3) Write(4,*)'Ti', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo ! Cr-vacancy
!  IF(imass(i)==3) Write(4,*)'Ti', rx(i)*lp,ry(i)*lp,rz(i)*lp ! imass = 3 vacancy

  IF(imass(i)==4) Write(4,*)'Th', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo ! Pd-vacancy
!  IF(imass(i)==4) Write(4,*)'Th', rx(i)*lp,ry(i)*lp,rz(i)*lp ! for replacing Ti in SrTiO3 for TC calculations

  IF(imass(i)==5) Write(4,*)'Th', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo ! Co-vacancy

  IF(imass(i)==6) Write(4,*)'Th', (rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo ! Mn-vacancy

 ENDIF 
ENDDO

!DC ends

101 format (i7,1x,6(f12.5,1x),f6.2,1x,4(i4,1x))

WRITE(3,*)'LAMMPS data file'
WRITE(3,*)' '
WRITE(3,*) zatom, 'atoms'
WRITE(3,*)'6 atom types'
WRITE(3,*) hxo,hx,'xlo xhi'
WRITE(3,*) hyo,hy,'ylo yhi'
WRITE(3,*) hzo,hz,'zlo zhi'
!WRITE(3,*) '0.0', h(1,1)*lp1,'xlo xhi'
!WRITE(3,*)'0.0', h(2,2)*lp2,'ylo yhi'
!WRITE(3,*)'0.0', h(3,3)*lp3,'zlo zhi'
WRITE(3,*)' '
WRITE(3,*) 'Atoms'
WRITE(3,*) ' '

!WRITE(33,*)'ITEM: TIMESTEP'
!WRITE(33,*)'00 '
!WRITE(33,*)'ITEM: NUMBER OF ATOMS'
!WRITE(33,*)zatom
!WRITE(33,*)'ITEM:BOX BOUNDS'
!WRITE(33,*) '0.0', h(1,1)*lp1
!WRITE(33,*)'0.0', h(2,2)*lp2
!WRITE(33,*)'0.0', h(3,3)*lp3
!WRITE(33,*)'ITEM: Atoms if type q x y z '

!dl-poly structure
!WRITE(34,*) 'title'
!If( coortype == 0) THEN
!  WRITE(34,*)coortype,typeatoms
!ELSEIF(coortype ==2) THEN
!  WRITE(34,*)coortype,typeatoms,natoms,0,0,0
!ENDIF
!
!WRITE(34,*) h(1,1), h(1,2), h(1,3)
!WRITE(34,*) h(2,1), h(2,2), h(2,3)
!WRITE(34,*) h(3,1), h(3,2), h(3,3)
!    
!
!WRITE(35,'(a)')'opti conp'
!WRITE(35,'(a)')'title'
!WRITE(35,'(a)')'pyro'
!WRITE(35,'(a)')'end'
!WRITE(35,'(a)')'cell 10'
!WRITE(35,'(3(f12.5,1x),3(i5))')h(1,1)*lp1, h(2,2)*lp2, h(3,3)*lp3, 90, 90, 90
!WRITE(35,'(a)')'cartesian'
!
!DO i=1,natoms
!  IF (ivac(i) == 0) THEN
!   IF(imass(i) == 1) Write(35,103)'Gd ', 'c ', rx(i)*lp1,ry(i)*lp2,rz(i)*lp3 
!   IF(imass(i) == 2) Write(35,103)'O ', 'c ', rx(i)*lp1,ry(i)*lp2,rz(i)*lp3 
!   IF(imass(i) == 3) Write(35,103)'Zr ', 'c ', rx(i)*lp1,ry(i)*lp2,rz(i)*lp3 
!  ENDIF
!ENDDO
!
103  format ((a),1x,a,3(f12.5,1x))
104  format (3(f12.5,1x),3(i5))
ik=0
DO i=1,natoms
  IF (ivac(i) == 0) THEN
    ik=ik+1
    WRITE(3,105)ik,imass(i),charge(i),(rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo
  ELSEIF (ireplace(i) == 1) THEN
    WRITE(3,105)ik,imass(i),charge(i),(rx(i)*(hx-hxo))+hxo,(ry(i)*(hy-hyo))+hyo,(rz(i)*(hz-hzo))+hzo
  ENDIF
ENDDO
!
!ik=0
!
!DO i=1,natoms
!  IF (ivac(i) == 0) THEN
!     ik=ik+1
!     WRITE(33,102)ik,imass(i),charge(i),rx(i)*lp1,ry(i)*lp2,rz(i)*lp3
!     WRITE(34,*)atomname(i),ik
!     WRITE(34,*)rx(i)*lp1,ry(i)*lp2,rz(i)*lp3
!    IF(coortype == 2)THEN
!      WRITE(34,*)x5(i),y5(i),z5(i)
!      WRITE(34,*)x6(i),y6(i),z6(i)
!   ENDIF 
!  ELSEIF (ireplace(i) == 1) THEN
!     WRITE(33,102)ik,imass(i),charge(i),rx(i)*lp1,ry(i)*lp2,rz(i)*lp3
!     WRITE(34,*)atomname(i),ik
!     WRITE(34,*)rx(i)*lp1,ry(i)*lp2,rz(i)*lp3
!    IF(coortype == 2)THEN
!       WRITE(34,*)x5(i),y5(i),z5(i)
!       WRITE(34,*)x6(i),y6(i),z6(i)
!    ENDIF 
!  
!  ENDIF
!ENDDO

102  format (i7,1x,i3,1x,4(f13.6,1x))
105  format (i7,1x,i3,1x,4(f13.6,1x))

END Subroutine WriteStructure

!***************************

Subroutine Substitute
USE variables
IMPLICIT NONE
INTEGER::nmols,ions,k,j,count_ce,count_o
REAL::rn1,rn2,cat_charge,x1,x2,y1,y2
INTEGER::index_atom,vac_atom
INTEGER,dimension(7000000)::inter_atom,inter_vac

Write(*,*)'Enter the number of Ce2O3 molecules to be substituted'
READ(*,*)nmols
ions = 2*nmols  ! 2Ce3+ ions replace 2 Ce4+ ion, and 3 Oxy replace 4 Oxy creating 1 Ovac
WRITE(*,*)'Enter charge of Ce3+'
READ(*,*) cat_charge

Call random_seed
WRITE(*,*)'substituting'
count_ce = 0; count_o = 0;

x1=0.33*h(1,1)
x2=0.66*h(1,1)
y1=0.33*h(2,2)
y2=0.66*h(2,2)

GO TO 10001

Do k=1,ions
 
   1001  Call random_number(rn1)
   index_atom = int(natoms)*rn1 + 1             !Randomly choosing an atom to be replaced

   IF(imass(index_atom) /= 1) Go to 1001

    DO j = 1,k
     IF(inter_atom(j) == index_atom) Go to 1001 ! Not to repeat in selecting an atom for intersitial
    ENDDO 
!   WRITE(*,*)index_atom,rz(index_atom),charge(index_atom)

!  IF((rx(index_atom)>=x1.and.rx(index_atom)<=x2).and.(ry(index_atom)>=y1.and.(ry(index_atom)<=y2))) THEN
!   IF (rz(index_atom)>9.and.rz(index_atom)<10) THEN
     inter_atom(k) = index_atom
     charge(index_atom) = cat_charge                     ! Ce4+ is made 3+
     imass(index_atom) = 3
     count_ce = count_ce + 1
!  ELSEIF((rx(index_atom)<x1.or.rx(index_atom)>x2).or.(ry(index_atom)<y1.or.(ry(index_atom)>y2))) THEN
!   ELSEIF (rz(index_atom)<9.or.rz(index_atom)>10) THEN
!     go to 1001
!   ENDIF
  
ENDDO 

10001 CONTINUE 

DO k=1,natoms
 ivac(k) = 0
ENDDO 
 

Do k=1,nmols
 
    2001 Call random_number(rn2)
    vac_atom = int(natoms)*rn2 + 1

    IF(imass(vac_atom) /= 2) Go to 2001

    DO j = 1,k
     IF(inter_vac(j) == vac_atom) Go to 2001 ! Not to repeat in selecting an atom for vacancy
    ENDDO

  IF((rx(vac_atom)>=x1.and.rx(vac_atom)<=x2).and.(ry(vac_atom)>=y1.and.(ry(vac_atom)<=y2))) THEN
!    IF (rz(vac_atom)>3.and.rz(vac_atom)<8.) THEN
      inter_vac(k) = vac_atom
      ident(vac_atom) = 0
      count_o = count_o + 1
      ivac(vac_atom) = 1
     WRITE(*,*)count_o
  ELSEIF((rx(vac_atom)<x1.or.rx(vac_atom)>x2).or.(ry(vac_atom)<y1.or.(ry(vac_atom)>y2))) THEN
!    ELSEIF (rz(vac_atom)<3.or.rz(vac_atom)>8) THEN
      go to 2001
    WRITE(*,*)ident(vac_atom)
  ENDIF

ENDDO
zatom = natoms - count_o
WRITE(*,*)'# of Ce3+ added = ', count_ce
WRITE(*,*)'# of O vacs created =',count_o

End Subroutine Substitute

!*******************************

Subroutine addinterstitials
USE variables
USE str
IMPLICIT NONE
INTEGER::nmols,k,j,index_atom, oxyions, ceions,icount,type_ints,add_mols
REAL::rn1,cat_charge,ani_charge
INTEGER,dimension(7000000)::inter_atom
ceions = 0; oxyions = 0

WRITE(*,*)'For Ce2O3 interstitials, type - 1'
WRITE(*,*)'For CeO2 interstitials, type - 2'
READ(*,*)type_ints

IF (type_ints == 1) WRITE (*,*)'Enter the number of Ce2O3 molecules to be added'
IF (type_ints == 2) WRITE (*,*)'Enter the number of CeO2 molecules to be added'
READ(*,*)nmols

WRITE(*,*)'Enter cation charge'
READ (*,*)cat_charge
WRITE(*,*)'Enter anion charge'
READ(*,*) ani_charge


IF (type_ints == 1) add_mols = 5*nmols
IF (type_ints == 2) add_mols = 3*nmols

Call random_seed

Do k=1,add_mols


 100 Call random_number(rn1)
 index_atom  = int (natoms) * rn1 + 1
! Randomly choosing an atom near which an added interstitial is placed
 
 IF(imass(index_atom) /= 1) Go to 100
 
 !IF((rx(index_atom)<-3.14.or.rx(index_atom)>3.4).or.(ry(index_atom)< -2.56.or.(ry(index_atom)>3.0))) THEN
 !  go to 100
 !ENDIF
 
 DO j = 1,k
  IF(inter_atom(j) == index_atom) Go to 100 
! Not to repeat in selecting an atom near which interstital is placed
 ENDDO 
! WRITE(*,*)natoms
 natoms = natoms + 1
 rx(natoms)  =  rx(index_atom) + 0.5
 ry(natoms)  =  ry(index_atom) + 0.5
 rz(natoms)  =  rz(index_atom) + 0.5

 IF (type_ints == 1 ) THEN

   IF (k <= (nmols*2)) imass(natoms)    =  3
   IF (k <= (nmols*2)) charge(natoms)   =  cat_charge
   IF (k > (nmols*2))  imass(natoms)    =  2
   IF (k > (nmols*2))  charge(natoms)   = ani_charge
   IF (k > (nmols*2))  oxyions = oxyions + 1
   IF (k <= (nmols*2)) ceions = ceions + 1
   icount = icount + 1

 ELSEIF (type_ints == 2) THEN

   IF (k <= (nmols*1)) imass(natoms)    =  1
   IF (k <= (nmols*1)) charge(natoms)   =  cat_charge
   IF (k > (nmols*1))  imass(natoms)    =  2
   IF (k > (nmols*1))  charge(natoms)   = ani_charge
   IF (k > (nmols*1))  oxyions = oxyions + 1
   IF (k <= (nmols*1)) ceions = ceions + 1
   icount = icount + 1

 ENDIF


 IF(rx(natoms) > h(1,1)) rx (natoms) =rx(natoms)-h(1,1)
 IF(ry(natoms) > h(2,2)) ry (natoms) =ry(natoms)-h(2,2)
 IF(rz(natoms) > h(3,3)) rz (natoms) =rz(natoms)-h(3,3)
 
 ident(natoms)   = natoms 
 x5(natoms)      =  0.0
 y5(natoms)      =  0.0
 z5(natoms)      =  0.0
 igrain(natoms)  =  igrain(1)
 igo(natoms)     =  igo(1)
 ithermo(natoms) =  ithermo(1)
 
 inter_atom(k) = index_atom
ENDDO
 zatom = natoms
 WRITE(*,*)zatom, ceions, oxyions 

End Subroutine addinterstitials


!****************

Subroutine FProcksalt

USE variables
USE str
IMPLICIT NONE
INTEGER::fps, k, index_atom,index_posi,ki,vac_num,int_num
INTEGER,dimension(7000000)::inter_atom,inter_posi,itick
REAL::rn1,rn2
REAL,dimension(7000000)::vac_x,vac_y,vac_z
OPEN (UNIT=4,STATUS='unknown',FILE='defect.xyz')

WRITE(*,*)'Enter the equal number of cation and anion FPs to be created'
READ(*,*)fps
fps = fps*2

vac_num=0;int_num=0

Do k=1,natoms+fps
 itick(k) = 0
ENDDO 

WRITE(4,*) natoms+fps   !because vacancies are also given an atom identity to be able to visualize them 
WRITE(4,*)' '


! Loop # 999
Do k=1,fps
Call random_seed
 100 Call random_number(rn1)
 200 Call random_number(rn2)

 index_atom  = int (natoms) * rn1 + 1   
! Randomly choosing an atom to make it an interstitial
 index_posi  = int (natoms) * rn2 + 1
! Randomly choosing an atom near which interstitial is placed 
  DO ki = 1, k
   IF (inter_atom(ki) == index_atom) go to 100 !Not to repeat in selecting an atom for interstitial 
   If (inter_atom(ki) == index_posi) go to 100 !Here, interstitial is prevented from being chosen as the atom near which interstitial os placed
   If (inter_posi(ki) == index_posi) go to 200 !Not to repeat the atom near which an interstitial has already been placed
  ENDDO
  inter_atom(k) = index_atom                   !Atom that is made interstitial
  inter_posi(k) = index_posi                   !Atom near which interstitial is placed 
  
  IF (imass(inter_atom(k)) == imass(inter_atom(k-1))) GO TO 100 !Alternating cation and anion

IF(((rx(index_posi)<-2.or.rx(index_posi)>2)).or.(ry(index_posi)<0.or.(ry(index_posi)>1))) GO TO 100
! Keeping within a certain volume
  
  vac_x(k) = rx(inter_atom(k))
  vac_y(k) = ry(inter_atom(k))
  vac_z(k) = rz(inter_atom(k))
  vac_num = vac_num + 1   
! Postion of the vacancy
  
  rx (inter_atom(k)) = rx (inter_posi(k)) + 0.25
  ry (inter_atom(k)) = ry (inter_posi(k)) + 0.25
  rz (inter_atom(k)) = rz (inter_posi(k)) + 0.25  
! Atom made interstitial

  IF (rx(inter_atom(k)) > h(1,1)) rx (inter_atom(k)) = rx (inter_atom(k)) - h (1,1)
  IF (ry(inter_atom(k)) > h(2,2)) ry (inter_atom(k)) = ry (inter_atom(k)) - h (2,2)
  IF (rz(inter_atom(k)) > h(3,3)) rz (inter_atom(k)) = rz (inter_atom(k)) - h (3,3)
  
  itick (ident(inter_atom(k))) = 1
  
  int_num = int_num + 1
   
  IF (imass(inter_atom(k)) == 1) THEN
    WRITE(4,*)'Ca', vac_x(k)*lp1,vac_y(k)*lp2,vac_z(k)*lp3                            !Ca is MG vacancy 
    WRITE(4,*)'Fe', rx(inter_atom(k))*lp1,ry(inter_atom(k))*lp2,rz(inter_atom(k))*lp3 ! Fe is Mg interstitial
  ELSEIF (imass(inter_atom(k)) == 2) THEN
    WRITE(4,*)'N', vac_x(k)*lp1,vac_y(k)*lp2,vac_z(k)*lp3                             !N is O vacancy 
    WRITE(4,*)'C', rx(inter_atom(k))*lp1,ry(inter_atom(k))*lp2,rz(inter_atom(k))*lp3 !C is O interstitial
  ENDIF
  
ENDDO
! Loop # 999 closed       
  
  WRITE(*,*) vac_num, int_num
  
  DO k=1,natoms
   IF (itick(ident(k)).ne.1) THEN
    IF  (imass(k) == 1) WRITE(4,*)'Ce',rx(k)*lp1,ry(k)*lp2,rz(k)*lp3
    If  (imass(k) == 2) WRITE(4,*)'O', rx(k)*lp1,ry(k)*lp2,rz(k)*lp3
   ENDIF
  ENDDO  
zatom = natoms 

END Subroutine FProcksalt



!****************


Subroutine STO_FPs

USE variables
USE str
IMPLICIT NONE
INTEGER::fps, k, index_atom,index_posi,ki,vac_num,int_num,stofps
INTEGER,dimension(7000000)::inter_atom,inter_posi,itick
REAL::rn1,rn2
REAL,dimension(7000000)::vac_x,vac_y,vac_z

INTEGER :: clock,ik,n
INTEGER, DIMENSION(:), ALLOCATABLE :: seed

OPEN (UNIT=51,STATUS='unknown',FILE='defect-2.xyz')
stofps = 0



WRITE(*,*)'For Sr FPs     - type 1 '
WRITE(*,*)'For O FPs      - type 2 '
WRITE(*,*)'For Ti FPs     - type 3 '
WRITE(*,*)'For SrO FPs    - type 4 '
WRITE(*,*)'For TiO2 FPs   - type 5 '
WRITE(*,*)'For SrTiO3 FPs - type 6 '
READ(*,*) stofps

WRITE(*,*)'Enter the number of FPs'
WRITE(*,*)'for SrO, TiO2 and SrTiO3, stoichiometric units of FPs will be created'
READ(*,*)fps

WRITE(*,*) fps
IF (stofps == 1) fps = fps*1
IF (stofps == 2) fps = fps*1
IF (stofps == 3) fps = fps*1
IF (stofps == 4) fps = fps*2
!IF (stofps == 5) fps = fps*3
IF (stofps == 5) fps = fps*2
!IF (stofps == 6) fps = fps*5
IF (stofps == 6) fps = fps*3

vac_num=0;int_num=0


Do k=1,natoms+fps
itick(k) = 0
ENDDO

WRITE(51,*) natoms+fps   !because vacancies are also given an atom identity to be able to visualize them
WRITE(51,*)' '

WRITE(*,*) 'entering loop'


! Loop # 999


CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))


DO k=1,fps

   CALL SYSTEM_CLOCK(COUNT=clock)
   seed = clock + 37 * (/ (ik - 1, ik = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)

!  Call random_seed
100 Call random_number(rn1)
200 Call random_number(rn2)
  
  index_atom  = int (natoms) * rn1 + 1
!    Randomly choosing an atom to make it an interstitial
  index_posi  = int (natoms) * rn2 + 1
!    Randomly choosing an atom near which interstitial is placed

   DO ki = 1, k
     IF (inter_atom(ki) == index_atom) go to 100 !Not to repeat in selecting an atom for interstitial
     If (inter_atom(ki) == index_posi) go to 100 !Here, interstitial is prevented from being chosen as the atom near which interstitial os placed
     If (inter_posi(ki) == index_posi) go to 200 !Not to repeat the atom near which an interstitial has already been placed
   ENDDO
   inter_atom(k) = index_atom                   !Atom that is made interstitial
   inter_posi(k) = index_posi                   !Atom near which interstitial is placed
  
   If(imass(inter_posi(k)).ne.1) go to 200    !Select only Sr near which interstitial is placed
   
   IF (stofps == 1) THEN
     IF  (imass(inter_atom(k)).ne.1) go to 100  ! This willl let it select only Sr for FP/
   ENDIF

   IF (stofps == 2) THEN
     IF  (imass(inter_atom(k)).ne.2) go to 100  ! This willl let it select only O for FP/
   ENDIF

   IF (stofps == 3) THEN
     IF  (imass(inter_atom(k)).ne.3) go to 100  ! This willl let it select only Ti for FP/
   ENDIF

   IF (stofps == 4) THEN
     IF (k<=(fps/2)) THEN
       IF (imass(inter_atom(k)).ne.1)   GO TO 100 !This willl let it select only Ti for FP from TIO2
   WRITE(*,*)ident(index_atom)
     ELSEIF (k>(fps/2)) THEN
       IF (imass(inter_atom(k)).ne.2)   GO TO 100 !This willl let it select only Ti for FP from TIO2
   WRITE(*,*)ident(index_atom)
     ENDIF

   ENDIF

   IF (stofps == 5) THEN
!     IF (k<=(fps/3)) THEN
      IF (k<=(fps/2)) THEN     ! for radiation damage, equal O and ce3+ are craeted
       IF (imass(inter_atom(k)).ne.3)   GO TO 100 !This willl let it select only Ti for FP from TIO2
!     ELSEIF (k>(fps/3)) THEN
      ELSEIF (k>(fps/2)) THEN
       IF (imass(inter_atom(k)).ne.2)   GO TO 100 !This willl let it select only Ti for FP from TIO2
     ENDIF
   ENDIF

   IF (stofps == 6) THEN
   
!     IF (k<=(fps/5)) THEN
!       IF (imass(inter_atom(k)).ne.1)  GO TO 100 !This willl let it select only Sr for FP from SrTiO2
!     ELSEIF (k>(fps/5).and. k<=(2*fps/5)) THEN
!       IF (imass(inter_atom(k)).ne.3)  GO TO 100 !This willl let it select only Ti for FP from SrTiO2
!     ELSEIF (k>(2*fps/5)) THEN
!       IF (imass(inter_atom(k)).ne.2)   GO TO 100 !This willl let it select only O for FP from SrTiO2
!     ENDIF

        IF (k<=(fps/3)) THEN
         IF (imass(inter_atom(k)).ne.1)  GO TO 100 !This willl let it select only Sr for FP from SrTiO2
          ELSEIF (k>(fps/3).and. k<=(2*fps/3)) THEN
         IF (imass(inter_atom(k)).ne.3)  GO TO 100 !This willl let it select only Ti for FP from SrTiO2
          ELSEIF (k>(2*fps/3)) THEN
         IF (imass(inter_atom(k)).ne.2)   GO TO 100 !This willl let it select only O for FP from SrTiO2
        ENDIF

   ENDIF

!IF(((rx(index_posi)<-2.or.rx(index_posi)>2)).or.(ry(index_posi)<0.or.(ry(index_posi)>1))) GO TO 100
! Keeping within a certain volume

vac_x(k) = rx(inter_atom(k))
vac_y(k) = ry(inter_atom(k))
vac_z(k) = rz(inter_atom(k))
vac_num = vac_num + 1
! Postion of the vacancy

!for fluorite
rx (inter_atom(k)) = rx (inter_posi(k)) + 0.25
ry (inter_atom(k)) = ry (inter_posi(k)) + 0.25
rz (inter_atom(k)) = rz (inter_posi(k)) + 0.25
!WRITE(*,*)ident(inter_atom(k)),vac_x(k)*lp,rx(inter_atom(k))*lp

!for STO
!rx (inter_atom(k)) = rx (inter_posi(k)) + 0.5
!ry (inter_atom(k)) = ry (inter_posi(k)) + 0.0
!rz (inter_atom(k)) = rz (inter_posi(k)) + 0.0
! Atom made interstitial

IF (rx(inter_atom(k)) > h(1,1)) rx (inter_atom(k)) = rx (inter_atom(k)) - h (1,1)
IF (ry(inter_atom(k)) > h(2,2)) ry (inter_atom(k)) = ry (inter_atom(k)) - h (2,2)
IF (rz(inter_atom(k)) > h(3,3)) rz (inter_atom(k)) = rz (inter_atom(k)) - h (3,3)

itick (ident(inter_atom(k))) = 1

int_num = int_num + 1

IF (imass(inter_atom(k)) == 1) THEN 
WRITE(51,*)'Ca', vac_x(k)*lp1,vac_y(k)*lp2,vac_z(k)*lp3                            ! Ca is Sr vacancy
WRITE(51,*)'Mg', rx(inter_atom(k))*lp1,ry(inter_atom(k))*lp2,rz(inter_atom(k))*lp3 ! Mg is Sr interstitial
ELSEIF (imass(inter_atom(k)) == 2) THEN
WRITE(51,*)'Na', vac_x(k)*lp1,vac_y(k)*lp2,vac_z(k)*lp3                             ! Na is O vacancy
WRITE(51,*)'C', rx(inter_atom(k))*lp1,ry(inter_atom(k))*lp2,rz(inter_atom(k))*lp3  ! C is O interstitial
ELSEIF (imass(inter_atom(k)) == 3) THEN
WRITE(51,*)'Zr', vac_x(k)*lp1,vac_y(k)*lp2,vac_z(k)*lp3                             ! Zr is Ti vacancy
WRITE(51,*)'Sc', rx(inter_atom(k))*lp1,ry(inter_atom(k))*lp2,rz(inter_atom(k))*lp3  !  Sc is Ti interstitial
ENDIF

ENDDO
! Loop # 999 closed

WRITE(*,*) '# of vacancies  ',vac_num
WRITE(*,*)'# of interstitials  ', int_num

DO k=1,natoms
IF (itick(ident(k)).ne.1) THEN
IF  (imass(k) == 1) WRITE(51,*)'Sr',rx(k)*lp1,ry(k)*lp2,rz(k)*lp3
If  (imass(k) == 2) WRITE(51,*)'O', rx(k)*lp1,ry(k)*lp2,rz(k)*lp3
If  (imass(k) == 3) WRITE(51,*)'Ti', rx(k)*lp1,ry(k)*lp2,rz(k)*lp3
ENDIF
ENDDO
zatom = natoms


! Write structure(*****
WRITE(3,*)'LAMMPS data file'
WRITE(3,*)' '
WRITE(3,*) zatom, 'atoms'
WRITE(3,*)'3 atom types'
WRITE(3,*) '0.0', h(1,1)*lp1,'xlo xhi'
WRITE(3,*)'0.0', h(2,2)*lp2,'ylo yhi'
WRITE(3,*)'0.0', h(3,3)*lp3,'zlo zhi'
WRITE(3,*)' '
WRITE(3,*)'Masses'
WRITE(3,*)' '
WRITE(3,*)'1 58.69'
WRITE(3,*)'2 55.84'
WRITE(3,*)'3 51.99'
WRITE(3,*)' '
WRITE(3,*) 'Atoms'
WRITE(3,*) ' '

DO i=1,natoms
     WRITE(3,102)ident(i),imass(i),charge(i),rx(i)*lp1,ry(i)*lp2,rz(i)*lp3
ENDDO
102  format (i7,1x,i3,1x,4(f13.6,1x))

END Subroutine STO_FPs



!****************

Subroutine ZrO2_Gd2Ti2O7

! Converts 1x1x1 ZrO2 fluorite to Pyrochlore
! By creating 1 O vac per unit cell, and by converting half of Zr atoms to Gd (or +3 atoms)
USE variables
USE str
IMPLICIT NONE
INTEGER:: k,ki,novacs,nzrtogd,index_atom, count, count_zr, index_posi
INTEGER,dimension (7000000)::itick, inter_atom, inter_posi
REAL::rn1,rn2

      INTEGER :: clock,ik,n
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

OPEN (UNIT=4,STATUS='unknown',FILE='defect.xyz')

Do i=1,natoms
ivac(i) = 0
ENDDO

WRITE(*,*)'enter number of Zr atoms to be converted'
READ(*,*)nzrtogd

novacs = 172      

!novacs = natoms/12        ! # of O vacancies
!nzrtogd =  2*natoms/12     ! # of Zr (or +4) atoms converted to Gd (+3)
count=0; count_zr=0

!create O vacs first

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

DO k=1,natoms

   CALL SYSTEM_CLOCK(COUNT=clock)
   seed = clock + 37 * (/ (ik - 1, ik = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)

!!!  Call random_seed
   100 Call random_number(rn1)
  
    index_atom  = int (natoms) * rn1 + 1    ! Randonly choosing an O atom to make it a vacancy
   
    IF (imass(index_atom) /= 1) GO to 100     !CHECK IF OXYGEN IMASS IS 1 OR 2 IN THE STRUCTURE FILE
   
    DO ki = 1, k
       IF (inter_atom(ki) == index_atom) go to 100 !Not to repeat in selecting an atom 
    ENDDO
    
    inter_atom(k) = index_atom                   !Atom that is made vacancy
    ivac(index_atom) = 1
    count = count + 1
    IF (count == novacs) go to 101
    writE(*,*)count,novacs,'rn1',rn1
 
ENDDO
101 continue

DO k=1,natoms

   CALL SYSTEM_CLOCK(COUNT=clock)
   seed = clock + 37 * (/ (ik - 1, ik = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)
!  Call random_seed
   200 Call random_number(rn2)   

    index_posi  = int (natoms) * rn2 + 1    ! Randomly choosing Zr atom to convert to +3 charge 

    IF (imass(index_posi) /= 2) GO to 200

    DO ki = 1, k
      IF (inter_posi(ki) == index_posi) go to 200 ! Not to repeat in selecting an atom 
    ENDDO

    inter_posi(k) = index_posi                      
    charge (index_posi) = 4.0    ! CHANGE CHARGE AS REQUIRED
    imass (index_posi) = 4       ! CHANGE IMASS AS REQUITED 
    count_zr = count_zr + 1

    write(*,*)count_zr, 'rn2',rn2
    IF (count_zr == nzrtogd) Go TO 201

ENDDO    
201 continue

zatom = natoms - count
write(*,*)zatom
!DEALLOCATE(seed)
END Subroutine ZrO2_Gd2Ti2O7

!***************

Subroutine PyrotoDefectfluo
!randomizes cations by swaping +4 and +3 positions and creates antisite cations
! Input structure is a pyrochlore structure
USE variables
USE str
IMPLICIT NONE
REAL::rn1,rn2
REAL,dimension(7000000)::temp_x,temp_y,temp_z
INTEGER:: k, index_atom,index_posi,ki, swap, swappings
INTEGER,dimension(7000000)::inter_atom,inter_posi

      INTEGER :: clock,ik,n
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

WRITE(*,*)'enter number of cation swappings'
read(*,*)swappings
Write(*,*)swappings
swap=0

Do i=1,natoms
ivac(i) = 0
ENDDO

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))


DO k=1,natoms

   CALL SYSTEM_CLOCK(COUNT=clock)
   seed = clock + 37 * (/ (ik - 1, ik = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)

!  Call random_seed
  100 Call random_number(rn1)
  200 Call random_number (rn2)

    index_atom  = int (natoms) * rn1 + 1  
    index_posi  = int (natoms) * rn2 + 1     

    IF (imass(index_atom) /= 1) GO to 100  ! Randomly choosing Gd atom to convert to Ti charge
    IF (imass(index_posi) /= 3) GO to 200  ! Randomly choosing ti atom to convert to Gd charge
   
     DO ki = 1, k
       IF (inter_atom(ki) == index_atom) go to 100 !Not to repeat in selecting an atom 
       IF (inter_posi(ki) == index_posi) go to 200 ! Not to repeat in selecting an atom
     ENDDO
       inter_atom(k) = index_atom 
       inter_posi(k) = index_posi   

       temp_x(k) = rx(inter_atom(k))
       temp_y(k) = ry(inter_atom(k))
       temp_z(k) = rz(inter_atom(k))

         rx(index_atom) = rx(index_posi)
         ry(index_atom) = ry(index_posi)
         rz(index_atom) = rz(index_posi)

         rx(index_posi) = temp_x(k)
         ry(index_posi) = temp_y(k)
         rz(index_posi) = temp_z(k)

         swap = swap + 1

WRITE(*,*)ident(index_atom), rz(index_atom)
WRITE(*,*)ident(index_posi), rz(index_posi)
         
         IF (swap == swappings) Go TO 10001



ENDDO         
10001 Continue

zatom = natoms
end Subroutine PyrotoDefectfluo
!***************

Subroutine Overlapping_atoms

USE variables
USE str
IMPLICIT NONE
INTEGER::j,k,icount
REAL::dx,dy,dz,dist
WRITE(*,*)'Entering atom-overlapping check loop' 
Do j=1,zatom
 Do k=1,zatom
  IF(j.ne.k) THEN
    dx = abs(rx(j) - rx(k))
    dy = abs(ry(j) - ry(k))
    dz = abs(rz(j) - rz(k))
    IF ( dx > h(1,1)/2 ) dx = dx-h(1,1)/2
    IF ( dx.le.h(1,1)/2 ) dx = dx+h(1,1)/2
    IF ( dy > h(2,2)/2 ) dy = dy-h(2,2)/2  
    IF ( dy.le.h(2,2)/2 ) dy = dy+h(2,2)/2
    IF ( dz > h(3,3)/2 ) dz = dz-h(3,3)/2
    IF ( dz.le.h(3,3)/2 ) dz = dz+h(3,3)/2
   dist = SQRT (dx*dx + dy*dy + dz*dz)
   IF (dist < 0.2) THEN
   WRITE(*,*)ident(j),ident(k),dx,dy,dz
   WRITE(*,*) dist
   icount = icount + 1
   WRITE(*,*)'# of overlapping atoms =', icount
  ENDIF
  ENDIF
 ENDDO
ENDDO
 
END Subroutine Overlapping_atoms

!****************

Subroutine Create_void

USE variables
USE str
IMPLICIT NONE
INTEGER::j,k,ik,num_voids, extra_atom, count,int_charge,new_count
REAL,dimension(7000000)::cx,cy
REAL::d_void, rn1, charge_void,dist,dx,dy,cation,anion,new_charge_void


!WRITE(*,*)'Enter the number of voids wanted'
!READ(*,*)num_voids
num_voids = 1

WRITE(*,*)'Enter the center of each void'
 Do j = 1,num_voids
  READ(*,*)cx(j),cy(j)     ! void is in two dimension, x and y for the Columnar grain structure. For three dimension, add z coordinate
 ENDDO
 
WRITE(*,*)'Enter the diameter of each void in NANOMETERS'
READ(*,*)d_void

!WRITE(*,*)'Enter cation and anion charge'   ! cation and anion charges
!READ (*,*)cation, anion
cation= 4.0; anion = -2.0
charge_void = 0.0
! DC 11/15/2016
d_void = d_void / (0.1*lp1)     ! diameter of void converted to ao units

Call random_seed

Do j=1,num_voids
 Do k=1,natoms
    dx = abs (rx(k) - cx(j))
    dy = abs (ry(k) - cy(j))
    dist = SQRT (dx*dx + dy*dy)
    IF (dist <= d_void) THEN
     ivac (k) = 1
     charge_void = charge_void + charge(k)
     count = count + 1
    ENDIF
    IF (k == natoms) THEN
!     WRITE(*,*)count
      WRITE(*,*)'charge_void before nullifying', charge_void    
      IF (charge_void > 0) THEN
        int_charge = (charge_void/cation)
!       WRITE(*,*)int_charge 
     
         DO ik = 1,int_charge
           100 Call random_number(rn1)
           extra_atom = int(natoms)*rn1 + 1
           IF (ivac(extra_atom) == 1) Go to 100
           IF (charge(extra_atom).lt.0) Go TO 100
             ivac (extra_atom) = 1
             new_charge_void = new_charge_void + charge(extra_atom)
             new_count = new_count + 1
         ENDDO
         WRITE(*,*)'charge_void after nullifying- cation', charge_void       

      ELSEIF (charge_void < 0) THEN
        int_charge = (charge_void/anion)
!       WRITE(*,*)int_charge      

         DO ik = 1,int_charge
           200 Call random_number(rn1)
           extra_atom = int(natoms)*rn1 + 1
            IF (ivac(extra_atom) == 1) Go to 200
            IF (charge(extra_atom).gt.0) Go TO 200
            ivac (extra_atom) = 1
            new_charge_void = new_charge_void + charge(extra_atom)
            new_count = new_count + 1
         ENDDO
         WRITE(*,*)'charge_void after nullifying-anion', new_charge_void

      ENDIF     
    ENDIF
 ENDDO
 WRITE(*,*)'Total atoms removed = ',count+new_count
 WRITE(*,*)'Total charge =', charge_void - new_charge_void
ENDDO
zatom = natoms - (count+new_count)
WRITE(*,*)zatom

End Subroutine Create_void    

!***************************

Subroutine SrTiO3_defects
Use variables
Use str
IMPLICIT NONE
INTEGER:: number,defect_kind, ij,sro_vac, random_vac, vacs, ti_replace, count_ti_replace, replace 
INTEGER:: sr_vac, o_vac, ti_vac, count_ti, count_o, count_sr
INTEGER,dimension(7000000)::tag_vac
REAL:: rn3     

sr_vac = 0; o_vac = 0; ti_vac = 0
count_ti = 0; count_o = 0; count_sr = 0
count_ti_replace = 0
WRITE(*,*)'For charge-neutral SrO vacancies,    type - 1'
WRITE(*,*)'For charge-neutral TiO2 vacancies,   type - 2'
WRITE(*,*)'For charge-neutral SrTiO3 vacancies, type - 3'
WRITE(*,*)'Replace Ti with another +4 ion,      type - 4'

READ(*,*)defect_kind

Do ij = 1, natoms
 ivac(ij) = 0
 ireplace (ij) = 0
ENDDO

IF (defect_kind == 1) THEN
  WRITE(*,*)'Enter the number of SrO vacancies required'
  READ(*,*)vacs
  sr_vac = vacs
  ti_vac = 0
  o_vac  = vacs
  
ELSEIF (defect_kind == 2) THEN
  WRITE(*,*)'Enter the number of TiO2 vacancies required'
  READ(*,*)vacs
  sr_vac = 0
  ti_vac = vacs
  o_vac  = 2*vacs
  
ELSEIF (defect_kind == 3) THEN
  WRITE(*,*)'Enter the number of SrTiO3 vacancies required'
  READ(*,*)vacs
  sr_vac  = vacs
  ti_vac  = vacs
  o_vac   = 3*vacs  
ENDIF

WRITE(*,*)sr_vac,o_vac,ti_vac

Do ij = 1,sr_vac          ! creating Sr vacancies
   
   400 Call random_number (rn3)        
   random_vac = int (natoms) * rn3 + 1  ! Randomly choosing an atom to make it a vacancy
 
   IF(ivac(random_vac) .ne. 0) Go To 400    ! atom that has not yet been tagged as a vacancy' prevents choosing the same atoms twice

     IF (imass(random_vac).ne.1) GO to 400 !Only Sr-vacs allowed here
       count_sr = count_sr +1
       ivac(random_vac) = 1           ! atom is tagged to be a vacancy

ENDDO


Do ij = 1, o_vac          ! creating O vacancies
   
   500 Call random_number (rn3)        
   random_vac = int (natoms) * rn3 + 1  
 
   IF(ivac(random_vac) .ne. 0) Go to 500    

     IF (imass(random_vac).ne.2) Go to 500 
       count_o = count_o +1
       ivac(random_vac) = 1           

ENDDO


Do ij = 1, ti_vac          ! creating Ti vacancies
   600 Call random_number (rn3)        
   random_vac = int (natoms) * rn3 + 1  

    IF (imass(random_vac).ne.3) Go to 600 
    
      IF(ivac(random_vac) .ne. 0) GO to 600    
       count_ti = count_ti + 1
       ivac(random_vac) = 1          

ENDDO

WRITE(*,*)' '
WRITE(*,*)'Sr-vacs',count_sr, 'O-vacs', count_o, 'Ti-vacs', count_ti

zatom = natoms-(count_sr + count_o + count_ti)

IF (defect_kind == 4) THEN
  WRITE(*,*)'Enter the number of Ti replacements required'
  READ(*,*)ti_replace
ENDIF


DO ij = 1, ti_replace

   700 Call random_number (rn3)
       replace = int (natoms) * rn3 + 1
        IF(ireplace(replace) .ne. 0) GO to 700
        IF (imass(replace).ne.3) Go to 700
        imass(replace) = 4
        ireplace (replace) = 1
        count_ti_replace = count_ti_replace + 1
ENDDO

WRITE(*,*) 'Ti ions replaced', count_ti_replace

END Subroutine SrTiO3_defects



!***************************

Subroutine CeO2_Schottky
Use variables
Use str
IMPLICIT NONE
INTEGER:: random_vac, vacs, ij
INTEGER:: ce_vac, o_vac, count_ce, count_o
INTEGER,dimension(7000000)::tag_vac
REAL:: rn3

o_vac = 0; ce_vac = 0
count_ce = 0; count_o = 0;


Do ij = 1, natoms
ivac(ij) = 0
ENDDO

WRITE(*,*)'Enter the number of CeO2 vacancies required'
READ(*,*)vacs
ce_vac = vacs
o_vac  = 2*vacs

WRITE(*,*)vacs

Do ij = 1, o_vac          ! creating O vacancies

  500 Call random_number (rn3)
  random_vac = int (natoms) * rn3 + 1  ! Randomly choosing an atom to make it a vacancy

  IF(ivac(random_vac) .ne. 0) Go to 500    ! atom that has not yet been tagged as a vacancy' prevents choosing the same atoms twice

  IF (imass(random_vac).ne.2) Go to 500 !Only O-vacs allowed here

!  IF (rz(random_vac)>4.and.rz(random_vac)<8.) THEN
   count_o = count_o +1
   ivac(random_vac) = 1           ! atom is tagged to be a vacancy
!  ELSEIF (rz(random_vac)<4.or.rz(random_vac)>8.) THEN
!   GO tO 500
!  ENDIF

ENDDO



Do ij = 1, ce_vac          ! creating Ce vacancies

  600 Call random_number (rn3)
  random_vac = int (natoms) * rn3 + 1  ! Randomly choosing an atom to make it a vacancy

  IF (imass(random_vac).ne.1) Go to 600 !Only Ce-vacs allowed here

  IF(ivac(random_vac) .ne. 0) GO to 600    ! atom that has not yet been tagged as a vacancy' prevents choosing the same atoms twice
!  IF (rz(random_vac)>4.and.rz(random_vac)<8.) THEN
    count_ce = count_ce + 1
    ivac(random_vac) = 1           ! atom is tagged to be a vacancy
!  ELSEIF (rz(random_vac)<4.or.rz(random_vac)>8.) THEN
!    GO tO 600
!  ENDIF

ENDDO

WRITE(*,*)' '
WRITE(*,*)'O-VACS', count_o, ':', '  Ce-VACS', count_ce

zatom = natoms-(count_o + count_ce)


END Subroutine CeO2_Schottky



!***************************

Subroutine random_alloy
! converts pure metal to alloy by randomly replacing some metal atoms' imass =1 to imass = 2

USE variables
IMPLICIT NONE
INTEGER::nmols,k,count_ce, j
REAL::rn1
! DC 09/27/2017
INTEGER::index_atom,vac_atom, alloycomp, comp_ni, comp_fe,comp_cr,comp_pd,&
         comp_co,comp_mn
INTEGER:: natoms_ni, natoms_fe, natoms_cr, natoms_pd, natoms_co, natoms_mn,&
          count_fe, count_cr, count_pd, count_co, count_mn
          
INTEGER,dimension(7000000)::inter_atom,inter_vac
      
      INTEGER :: clock,ik,n
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      

WRITE(*,*) 'Enter composition in percentages of Ni Fe Cr Pd Co Mn (e.g. 16.67 16.67 16.67 16.67 16.67 16.67)'
READ(*,*) comp_ni, comp_fe, comp_cr, comp_pd, comp_co, comp_mn 

natoms_ni = INT ((comp_ni * natoms)/100)
natoms_fe = INT ((comp_fe * natoms)/100)
natoms_cr = INT ((comp_cr * natoms)/100)
natoms_pd = INT ((comp_pd * natoms)/100)
natoms_co = INT ((comp_co * natoms)/100)
natoms_mn = INT ((comp_mn * natoms)/100)

WRITE(*,*)natoms_ni, natoms_fe, natoms_cr, natoms_pd, natoms_co, natoms_mn
nmols = natoms_fe + natoms_cr + natoms_pd + natoms_co + natoms_mn

WRITE(*,*)nmols

count_ce=0; count_fe=0; count_cr =0; count_pd =0; count_co= 0; count_mn=0;

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

DO k=1,natoms

   CALL SYSTEM_CLOCK(COUNT=clock)
   seed = clock + 37 * (/ (ik - 1, ik = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)
   1001 Call random_number(rn1)
   index_atom = int(natoms)*rn1 + 1             !Randomly choosing an atom to be replaced

! write(*,*)'index_atom', index_atom

   IF(imass(index_atom) /= 1) Go to 1001

    DO j = 1,k
     IF(inter_atom(j) == index_atom) Go to 1001 ! Not to repeat in selecting an atom for intersitial
    ENDDO 

      inter_atom(k) = index_atom
     
       IF (count_ce .lt. natoms_fe) THEN
       imass(index_atom) = 2 
       count_fe = count_fe + 1    
       ELSEIF (count_ce .ge. natoms_fe .and. count_ce.lt.(natoms_fe+natoms_cr))THEN
       imass(index_atom) = 3 
       count_cr = count_cr + 1
       ELSEIF (count_ce .ge.(natoms_fe+natoms_cr).and.&
               count_ce.lt.(natoms_fe+natoms_cr+natoms_pd)) THEN
       imass(index_atom) = 4
       count_pd = count_pd + 1
       ELSEIF (count_ce .ge.(natoms_fe+natoms_cr+natoms_pd).and.&
               count_ce.lt.(natoms_fe+natoms_cr+natoms_pd+natoms_co)) THEN
       imass(index_atom) = 5
       count_co = count_co + 1
       ELSEIF (count_ce .ge.(natoms_fe+natoms_cr+natoms_pd+natoms_co).and.&
               count_ce.lt.(natoms_fe+natoms_cr+natoms_pd+natoms_co+natoms_mn)) THEN
       imass(index_atom) = 6
       count_mn = count_mn + 1
       ENDIF
       count_ce = count_ce + 1
!     write(*,*)'atoms converted',count_ce
     If(count_ce == nmols) go to 1002
ENDDO
1002 continue
write(*,*)'replaced atoms', count_ce, count_fe, count_cr, count_pd, count_co, count_mn

zatom=natoms

End subroutine random_alloy



!***************************

Subroutine random_alloy_vacancies
! creates vacancies in pure metal or in a random alloy

USE variables
IMPLICIT NONE
INTEGER::nmols,k,count, j
REAL::rn1
INTEGER::index_atom,vac_atom
INTEGER,dimension(7000000)::inter_atom,inter_vac
      
      INTEGER :: clock,ik,n
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

 
 Write(*,*)'Enter the # of Ni vacancies to be created'
 READ(*,*)nmols

 Do i=1,natoms
 ivac(i) = 0
 ENDDO

 count=0
 CALL RANDOM_SEED(size = n)
 ALLOCATE(seed(n))

 DO k=1,natoms

   CALL SYSTEM_CLOCK(COUNT=clock)
   seed = clock + 37 * (/ (ik - 1, ik = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)
   1001 Call random_number(rn1)
   index_atom = int(natoms)*rn1 + 1             !Randomly choosing an atom to be replaced

!  write(*,*)'index_atom', index_atom, count

!DC changes 07/12/2016
!   IF(imass(index_atom) /= 1) Go to 1001
!   IF(imass(index_atom) /= 2) Go to 1001
!   IF(imass(index_atom) /= 3) Go to 1001
!   IF(imass(index_atom) /= 4) Go to 1001
!   IF(imass(index_atom) /= 6) Go to 1001
!   IF(imass(index_atom) /= 6) Go to 1001
   

    DO j = 1,k
     IF(inter_atom(j) == index_atom) Go to 1001 ! Not to repeat in selecting an atom for intersitial
    ENDDO 

     inter_atom(k) = index_atom
    ivac(index_atom) = 1
    count = count + 1
!     write(*,*)'atoms converted',count_ce
     If(count == nmols) go to 1002
!     write(*,*) index_atom, count
ENDDO
1002 continue
write(*,*)'vacancies created', count

zatom=natoms-count

End subroutine random_alloy_vacancies

!***************************

Subroutine random_alloy_interstitials
! creates interstitials in pure metal or in a random alloy

USE variables
IMPLICIT NONE
INTEGER::nmols,k,count, j, katoms
REAL::rn1, u
INTEGER::index_atom,vac_atom, type_ints
INTEGER,dimension(7000000)::inter_atom,inter_vac
      
      INTEGER :: clock,ik,n
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      

WRITE(*,*)'For Ni interstitials, type - 1'
WRITE(*,*)'For Fe interstitials, type - 2'
WRITE(*,*) 'For Cr interstitials, type - 3'
WRITE(*,*) 'For Pd interstitials, type - 4'
WRITE(*,*) 'For random Ni-Fe interstitials,type - 5'
WRITE(*,*) 'For random Ni-Cr interstitials,type - 6'
WRITE(*,*) 'For random Ni-Pd interstitials,type - 7'
WRITE(*,*) 'For random Fe-Cr interstitials,type - 8'
WRITE(*,*) 'For random Fe-Pd interstitials,type - 9'
WRITE(*,*) 'For random Cr-Pd interstitials,type - 10'
WRITE(*,*) 'For random Ni-Fe-Cr interstitials,type - 11'
WRITE(*,*) 'For random Ni-Fe-Pd interstitials,type - 12'
WRITE(*,*) 'For random Ni-Cr-Pd interstitials,type - 13'
WRITE(*,*) 'For random Fe-Cr-Pd interstitials,type - 14'
WRITE(*,*) 'For random Ni-Fe-Cr-Pd interstitials,type - 15'
WRITE(*,*) 'For random Ni-Fe-Cr-Pd-Co interstitials,type - 16'
WRITE(*,*) 'For random Ni-Fe-Cr-Pd-Co-Mn interstitials,type - 17'
WRITE(*,*) 'Additional changes can be added later'
READ(*,*)type_ints

Write(*,*)'Enter the # of interstitials to be created'
READ(*,*)nmols


Do i=1,natoms
ivac(i) = 0
ENDDO
katoms = natoms
count=0
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

DO k=1,nmols

   CALL SYSTEM_CLOCK(COUNT=clock)
   seed = clock + 37 * (/ (ik - 1, ik = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)
   1001 Call random_number(rn1)
   index_atom = int(natoms)*rn1 + 1             !Randomly choosing an atom to be replaced

! write(*,*)'index_atom', index_atom

!   IF(imass(index_atom) /= 1) Go to 1001
!   IF(imass(index_atom) /= 2) Go to 1001
!   IF(imass(index_atom) /= 3) Go to 1001

    DO j = 1,k
     IF(inter_atom(j) == index_atom) Go to 1001 ! Not to repeat in selecting an atom for intersitial
    ENDDO 
     
     katoms = katoms + 1
!changed by DC 26th July 2017

      rx(katoms)  =  ((rx(index_atom)*(hx-hxo))+hxo)/lp1
      ry(katoms)  =  ((ry(index_atom)*(hy-hyo))+hyo)/lp2
      rz(katoms)  =  ((rz(index_atom)*(hz-hzo))+hzo)/lp3

      rx(katoms)  =  rx(katoms) + 0.5
      ry(katoms)  =  ry(katoms) + 0.5
      rz(katoms)  =  rz(katoms) + 0.5

 IF(rx(katoms) > h(1,1)) rx (katoms) =rx(katoms)-h(1,1)
 IF(ry(katoms) > h(2,2)) ry (katoms) =ry(katoms)-h(2,2)
 IF(rz(katoms) > h(3,3)) rz (katoms) =rz(katoms)-h(3,3)     
       
     inter_atom(k) = index_atom
     write(*,*) index_atom

      IF (type_ints == 1 )  imass(katoms) = 1
      write(*,*)  imass(katoms), count
      IF (type_ints == 2 )  imass(katoms) = 2
      write(*,*)  imass(katoms), count
      IF (type_ints == 3 )  imass(katoms) = 3 
      write(*,*) imass(katoms), count
      IF (type_ints == 4 )  imass(katoms) = 4
      write(*,*) imass(katoms), count
!      call random_number(u)
!      write(*,*) u
      IF (type_ints == 5)  then
        call random_number(u)
        imass(katoms) = 1+FLOOR(2*u)   
      write(*,*)  imass(katoms), count
      elseif (type_ints == 6 ) then
1200         call random_number(u)
             imass(katoms)=1+FLOOR(3*u)
             if (imass(katoms)==2) go to 1200; continue
      write(*,*) imass(katoms), count
      elseif (type_ints == 7 ) then
1300         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==2.or.imass(katoms)==3) go to 1300; continue
              write(*,*) imass(katoms), count
      elseif (type_ints == 8 ) then
1400         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==1.or.imass(katoms)==4) go to 1400; continue
              write(*,*) imass(katoms), count
     elseif (type_ints == 9 ) then
1500         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==1.or.imass(katoms)==3) go to 1500; continue
              write(*,*) imass(katoms), count
     elseif (type_ints == 10 ) then
1600         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==1.or.imass(katoms)==2) go to 1600; continue
             write(*,*) imass(katoms), count
     elseif (type_ints == 11 ) then
1700         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==4) go to 1700; continue
             write(*,*) imass(katoms), count
     elseif (type_ints == 12 ) then
1800         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==3) go to 1800; continue
             write(*,*) imass(katoms), count
     elseif (type_ints == 13 ) then
1900         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==2) go to 1900; continue
             write(*,*) imass(katoms), count
     elseif (type_ints == 14 ) then
2000         call random_number(u)
             imass(katoms) = 1+FLOOR(4*u)
             if (imass(katoms)==1) go to 2000; continue
             write(*,*) imass(katoms), count
      elseif (type_ints == 15) then
              call random_number(u)
              imass(katoms) = 1+FLOOR(4*u)
              write(*,*)  imass(katoms), count
      elseif (type_ints == 16) then
              call random_number(u)
              imass(katoms) = 1+FLOOR(5*u)
              write(*,*)  imass(katoms), count
      elseif (type_ints == 17) then
              call random_number(u)
              imass(katoms) = 1+FLOOR(6*u)
              write(*,*)  imass(katoms), count
       end if

      charge(katoms) = 0.0
!      ivac(katoms) = 1

    count = count + 1
!     If(count == nmols) go to 1002
     rx(katoms)= ((rx(katoms)*lp1)-hxo)/(hx-hxo)
     ry(katoms)= ((ry(katoms)*lp2)-hyo)/(hy-hyo)
     rz(katoms)= ((rz(katoms)*lp3)-hzo)/(hz-hzo)
   
ENDDO
1002 continue
     write(*,*) index_atom, count
write(*,*)'interstitials created', count
write(*,*)katoms
zatom=katoms 
natoms = zatom
write(*,*)zatom
End subroutine random_alloy_interstitials

!*******************

Subroutine Create_metalvoid

USE variables
USE str
IMPLICIT NONE
INTEGER::j,k,ik,num_voids, extra_atom, count,int_charge,new_count,nvacs
REAL,dimension(7000000)::cx,cy, cz
REAL::d_void, rn1, charge_void,dist,dx,dy,dz, cation,anion,new_charge_void

 Do i=1,natoms
 ivac(i) = 0
 ENDDO
 
num_voids = 1
count = 0
WRITE(*,*)'Enter the center of each void'
 Do j = 1,num_voids
  READ(*,*)cx(j),cy(j),cz(j)     ! void is in two dimension, x and y for the Columnar grain structure. For three dimension, add z coordinate
!center is in l.p units
WRITE(*,*)cx(j),cy(j),cz(j)
 ENDDO
 
WRITE(*,*)'Enter the diameter of each void in L.P. units'
READ(*,*)d_void

WRITE(*,*)'Enter number of vacancies to be created in the void'
READ(*,*)nvacs

!d_void = d_void / (0.1*lp)     ! diameter of void converted to ao units



Do j=1,num_voids
 Do k=1,natoms
     IF (count == nvacs) go to 100

     rx(k)  =  ((rx(k)*(hx-hxo))+hxo)
     ry(k)  =  ((ry(k)*(hy-hyo))+hyo)
     rz(k)  =  ((rz(k)*(hx-hzo))+hzo)
     
    dx = abs (rx(k) - cx(j))
    dy = abs (ry(k) - cy(j))
    dz = abs (rz(k) - cz(j))

    write(*,*)k,rx(k),dx,ry(k),dy,rz(k),dz
    IF (dx <= d_void.and.dy <=d_void.and.dz<=d_void) THEN
     ivac (k) = 1
     count = count + 1
     write(*,*) count
    ENDIF
    
     rx(k)= (rx(k)-hxo)/(hx-hxo)
     ry(k)= (ry(k)-hyo)/(hy-hyo)
     rz(k)= (rz(k)-hzo)/(hz-hzo)

 ENDDO
100 continue
ENDDO
zatom = natoms - (count)
WRITE(*,*)zatom

End Subroutine Create_metalvoid  

!***************************

Subroutine Create_zro2box
! converts a central region-box of CeO2 to ZrO@ as done in Aidhy et al., JPCC, 2014-118,4207

USE variables
USE str
IMPLICIT NONE
INTEGER::k,count_ce, j
REAL:: bsize_x1, bsize_x2, bsize_y1, bsize_y2,bsize_z1, bsize_z2

bsize_x1 = 0.33*h(1,1)
bsize_x2 = 0.66*h(1,1)
bsize_y1 = 0.33*h(2,2)
bsize_y2 = 0.66*h(2,2)
bsize_z1 = 0.33*h(3,3)
bsize_z2 = 0.66*h(3,3)      
      
  count_ce = 0
  write(*,*) natoms    
Do k=1,natoms
 If(imass(k)==1) THEN
   IF (rx(k) > bsize_x1 .and. rx(k) < bsize_x2) THEN
   IF (ry(k) > bsize_y1 .and. ry(k) < bsize_y2) THEN
!   IF (rz(k) > bsize_z1 .and. rz(k) < bsize_z2) THEN   
    imass(k) = 3
    count_ce=count_ce + 1
   ENDIF
   ENDIF
!   ENDIF
 ENDIF

ENDDO   
write(*,*) count_ce
zatom = natoms
End subroutine Create_zro2box



!***************************


