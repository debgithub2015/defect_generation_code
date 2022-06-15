!Dilpuneet Aidhy
! ORNL
!March 14, 2013
!code to create single crystal structure
! It takes coordinates and charges of atoms in a molecule, e.g, cor CeO2, coordinates of 2 O and 1 Ce atom along with -2 and +4 charges
! IT also takes a bravais lattice to translate the molecule
! use dipolemoment.f90 code, in case to create random oxygen vacacnies in a structure

!to make CeO2 in ZrO@ precipitate, use Bravais_lattice_embed_precipitate,Attach_molecule, embed_precipitate subroutines. In Attach_molecule, dont use the  If (tag(i) == 1) THEN .

Module variables
IMPLICIT NONE
INTEGER::i,ii,num_lp,natoms,uc,vacs,nx,ny,nz,ntypes
INTEGER,dimension(5000000):: tag_cell,tag_atom,ident,imass
REAL,dimension(50000000)::ai,aj,ak,rx,ry,rz,charge
REAL::lp_a,lp_b,lp_c,mass_a,mass_b,mass_c
REAL,dimension(1000000)::tag
End Module

Program makestr
USE variables
IMPLICIT NONE


CALL Bravais_lattice
!CALL Bravais_lattice_embed_precipitate
CALL Attach_molecule 
!CALL embed_precipitate

!CALL Attach_molecule_embed_precipitate

CALL Write_HEll2dstr
CALL lammps_strfile
Call dlpoly_strfile

END Program makestr


!****************************************

Subroutine Bravais_lattice

USE variables
IMPLICIT NONE
INTEGER::num_bs,x,y,z,ij
OPEN(UNIT=1,STATUS='OLD',FILE='bravais-basis') ! defines bravais basis set

! bravais-fcc for Bi2O3
! bravais-cubic for Pyrochlore

READ(1,*)num_bs     !Number of bravais basis set points, e.g., fcc has 4, bcc has 2

Do ij=1,num_bs 
 READ(1,*)ai(ij),aj(ij),ak(ij)
ENDDO 

WRITE(*,*)'Enter number of unit cells in each x, y and z directions'
READ(*,*)nx,ny,nz

Do x=1,nx
 Do y=1,ny
  Do z=1,nz
  uc = uc + 1            ! number of unit cells (uc)
   Do ij=1,num_bs
    num_lp = num_lp + 1  ! total number of lattice points, i.e., num_bs + in every unit cell
    ai(num_lp) = ai(ij) + (x-1)  ! x coordinate of a lattice point 
    aj(num_lp) = aj(ij) + (y-1)
    ak(num_lp) = ak(ij) + (z-1)
    tag_cell (num_lp) = uc
!    Write(*,*)ai(num_lp),aj(num_lp),ak(num_lp)
    ENDDO
  ENDDO
 ENDDO
ENDDO

!WRITE(*,*)num_lp

DO ij=1,num_lp
!WRITE(*,101)ij,ai(ij),aj(ij),ak(ij), tag_cell(ij), uc
ENDDO         

101  format (i5,3(f7.3,1x),i5, 1x,i5)
End Subroutine Bravais_lattice

!****************************************


Subroutine Bravais_lattice_embed_precipitate

USE variables
IMPLICIT NONE
INTEGER::num_bs,x,y,z,ij
REAL::pp1_i,pp1_o,pp2_i,pp2_o,pp3_i,pp3_o
OPEN(UNIT=1,STATUS='OLD',FILE='bravais-basis') ! defines bravais basis set

! bravais-fcc for Bi2O3
! bravais-cubic for Pyrochlore

READ(1,*)num_bs     !Number of bravais basis set points, e.g., fcc has 4, bcc has 2

Do ij=1,num_bs 
 READ(1,*)ai(ij),aj(ij),ak(ij)
ENDDO 

WRITE(*,*)'Enter number of unit cells in each x, y and z directions'
READ(*,*)nx,ny,nz

Do x=1,nx
 Do y=1,ny
  Do z=1,nz
  uc = uc + 1            ! number of unit cells (uc)
   Do ij=1,num_bs
    num_lp = num_lp + 1  ! total number of lattice points, i.e., num_bs + in every unit cell
    ai(num_lp) = ai(ij) + (x-1)  ! x coordinate of a lattice point 
    aj(num_lp) = aj(ij) + (y-1)
    ak(num_lp) = ak(ij) + (z-1)
    tag_cell (num_lp) = uc
    tag(num_lp) = 1 
   ENDDO 
  ENDDO
 ENDDO
ENDDO

pp1_i=7; pp1_o=14
pp2_i=27; pp2_o=34
pp3_i=47; pp3_o=54

Do i=1,num_lp
 IF ((ai(i).ge.nx/3.and.ai(i).lt.(2*nx/3)).and.(aj(i).ge.ny/3.and.aj(i).lt.(2*ny/3))) THEN
!for three precipitates,and for 60 unit cells in X
!IF (((ai(i).ge.pp1_i.and.ai(i).lt.pp1_o).or. &
!(ai(i).ge.pp2_i.and.ai(i).lt.pp2_o).or. &
!(ai(i).ge.pp3_i.and.ai(i).lt.pp3_o)).and. &
!(aj(i).ge.ny/3.and.aj(i).lt.(2*ny/3))) THEN
 tag(i) = 0
ENDIF
ENDDO 
 
Do i=1,num_lp
 IF (tag(i) .eq.1 )THEN
! WRITE(*,*)'Ce', ai(i),aj(i),ak(i)
 ENDIF
ENDDO  

!WRITE(*,*)num_lp

DO ij=1,num_lp
!WRITE(*,101)ij,ai(ij),aj(ij),ak(ij), tag_cell(ij), uc
ENDDO         

101  format (i5,3(f7.3,1x),i5, 1x,i5)
End Subroutine Bravais_lattice_embed_precipitate

!****************************************

Subroutine  embed_precipitate
USE variables
IMPLICIT NONE
INTEGER::ij
REAL::pp1_i,pp1_o,pp2_i,pp2_o,pp3_i,pp3_o

pp1_i=7; pp1_o=14
pp2_i=27; pp2_o=34
pp3_i=47; pp3_o=54

Do ij = 1,natoms
 IF (imass(ij).eq.1) THEN
 IF (((rx(ij).ge.(nx/3).and.rx(ij).le.(2*nx/3))).and.(ry(ij).ge.(ny/3).and.ry(ij).le.(2*ny/3))) THEN
!IF (((rx(ij).ge.pp1_i.and.rx(ij).lt.pp1_o).or. &
!(rx(ij).ge.pp2_i.and.rx(ij).lt.pp2_o).or. &
!(rx(ij).ge.pp3_i.and.rx(ij).lt.pp3_o)).and. &
!(ry(ij).ge.ny/3.and.ry(ij).lt.(2*ny/3))) THEN
   imass(ij) = 3
 ENDIF
 ENDIF
 
ENDDO

END Subroutine  embed_precipitate


!****************************************


Subroutine Attach_molecule
USE variables
IMPLICIT NONE
INTEGER:: num_atoms_mol,k,j,ij
REAL,dimension(20000)::mx,my,mz
OPEN(UNIT=2,STATUS='OLD',FILE='molecule') ! defines molecule

READ(2,*)num_atoms_mol,ntypes, lp_a,lp_b,lp_c     !number of atoms in a molecule
WRITE(*,*)lp_a,lp_b,lp_c
IF(ntypes == 1) THEN
  READ(2,*)mass_a
!DC changes
!WRITE(*,*)'num2'
WRITE(*,*) 'num1'
ELSEIF(ntypes == 2) THEN
  READ(2,*)mass_a,mass_b
WRITE(*,*)'num2'
ELSEIF(ntypes == 3) THEN
  READ(2,*)mass_a,mass_b,mass_c
WRITE(*,*)'num3'

ENDIF

write(*,*) mass_a,mass_b

DO i=1,num_atoms_mol
 READ(2,*)mx(i),my(i),mz(i), charge(i),imass(i)
ENDDO


DO i=1,num_lp

! If (tag(i) == 1) THEN   !!!tag USED only when using  Bravais_lattice_embed_precipitate
 
 DO j=1,num_atoms_mol
  natoms=natoms+1
  rx(natoms) = mx(j) + ai(i)
  ry(natoms) = my(j) + aj(i)
  rz(natoms) = mz(j) + ak(i)
  charge(natoms) = charge(j)
  tag_atom (natoms) = tag_cell (i)   ! tagging the atoms to the unit cell
  ident(natoms) = natoms
  imass(natoms) = imass(j)
 ENDDO
 
! ENDIF
ENDDO
                       ! natoms: Total number of atoms in the systems by now
102  format (a,1x,3(f5.2,1x))

END Subroutine Attach_molecule 





!****************************************

Subroutine Read_xyzstructure
USE variables
IMPLICIT NONE
INTEGER::ij
OPEN(UNIT=5,STATUS='OLD',FILE='xyzstructure') ! defines structure in xyz file with charge

READ(5,*)natoms

Do ij=1,natoms
 READ(5,*)rx(ij),ry(ij),rz(ij),charge(ij)
ENDDO

END Subroutine Read_xyzstructure

!****************************************

Subroutine Write_HEll2dstr
USE variables
IMPLICIT NONE
INTEGER::ij
OPEN(UNIT=6,STATUS='unknown',FILE='hell2d_structure') ! hell2d format structure file
OPEN(UNIT=7,STATUS='unknown',FILE='structure.xyz') ! xyz format structure file

WRITE(6,*)'1  3  0  2  1'
WRITE(6,*)'1  2  0  0  1'
WRITE(6,*)natoms
WRITE(6,103)mass_a,mass_b,mass_c
WRITE(6,*)'0.0  0.0'
WRITE(6,*)'0 1 0'
WRITE(6,*)'0 0'
WRITE(6,*)'1.0 0.0 0.0'
WRITE(6,*)'0.0 1.0 0.0'
WRITE(6,*)'0.0 0.0 1.0'
WRITE(6,*)nx, 0.0, 0.0
WRITE(6,*)0.0, ny, 0.0
WRITE(6,*)0.0,0.0, nz
WRITE(6,*)nx,ny,nz
WRITE(6,*)'-33.0'
WRITE(6,*)'0.5 0.5'
WRITE(6,*)'0.7 0.7'
WRITE(6,*)6, 6
WRITE(6,*)12, 12     

WRITE(7,*)natoms
WRITE(7,*)' '

Do ij = 1, natoms
 WRITE(6,105)ident(ij),rx(ij),ry(ij),rz(ij),0.0, 0.0, 0.0, charge(ij),imass(ij),0, 0, 0
  IF (imass(ij) == 1) WRITE(7,106)'Sr',rx(ij)*lp_a,ry(ij)*lp_b,rz(ij)*lp_c
  IF (imass(ij) == 2) WRITE(7,106)'O',rx(ij)*lp_a,ry(ij)*lp_b,rz(ij)*lp_c
  IF (imass(ij) == 3) WRITE(7,106)'Ti',rx(ij)*lp_a,ry(ij)*lp_b,rz(ij)*lp_c
ENDDO
103  format (3(f8.3,1x))
105  format (i8,6(f9.4,1x),f6.3, 4(i3,1x))
106  format (a,3(f9.4,1x))
 
END Subroutine Write_HEll2dstr
 
!****************************************

Subroutine lammps_strfile
USE variables
IMPLICIT NONE
INTEGER::ij

OPEN(UNIT=8,STATUS='unknown',FILE='lammps_structure') ! LAMMPS format structure file

OPEN(UNIT=9,STATUS='unknown',FILE='gulp.str') ! GULP str file

WRITE(8,*)'Title'
WRITE(8,*)' '
WRITE(8,*)natoms,'atoms'
WRITE(8,*)ntypes,'atom types'
WRITE(8,*)0.0,nx*lp_a,'xlo xhi'
WRITE(8,*)0.0,ny*lp_b,'ylo yhi'
WRITE(8,*)0.0,nz*lp_c,'zlo zhi'

WRITE(8,*)' '
WRITE(8,*)'Atoms'
WRITE(8,*)' '

Do ij = 1, natoms
 WRITE(8,107)ident(ij),imass(ij),0.0,rx(ij)*lp_a,ry(ij)*lp_b,rz(ij)*lp_c
! WRITE(8,106)ident(ij),imass(ij),rx(ij)*lp_a,ry(ij)*lp_b,rz(ij)*lp_c

ENDDO

107 format (i7,1x,i3,1x,4(f13.6,1x))
106 format (i7,1x,i3,1x,3(f13.6,1x))

WRITE(9,'(a)')'opti conp'
WRITE(9,'(a)')'title'
WRITE(9,'(a)')'pyro'
WRITE(9,'(a)')'end'
WRITE(9,'(a)')'cell 10'
WRITE(9,'(3(f12.5,1x),3(i5))')nx*lp_a, ny*lp_b, nz*lp_c, 90, 90, 90
WRITE(9,'(a)')'cartesian'

DO i=1,natoms

   IF(imass(i) == 1) Write(9,103)'Gd ', 'c ', rx(i)*lp_a,ry(i)*lp_b,rz(i)*lp_c 
   IF(imass(i) == 2) Write(9,103)'O ', 'c ', rx(i)*lp_a,ry(i)*lp_b,rz(i)*lp_c 
   IF(imass(i) == 3) Write(9,103)'Zr ', 'c ', rx(i)*lp_A,ry(i)*lp_b,rz(i)*lp_c 

ENDDO

103  format ((a),1x,a,3(f12.5,1x))
104  format (3(f12.5,1x),3(i5))



END Subroutine lammps_strfile       
       
       
!****************************************


Subroutine dlpoly_strfile
USE variables
IMPLICIT NONE


OPEN(UNIT=10,STATUS='unknown',FILE='dl_poly') ! LAMMPS format structure file


DO i=1,natoms

   IF(imass(i) == 1) THEN
       Write(10,109)'Gd ', i
       WRITE(10,111) rx(i)*lp_a,ry(i)*lp_b,rz(i)*lp_c 
    ELSEIF(imass(i) == 2) THEN
       Write(10,109)'O ', i
       WRITE(10,111) rx(i)*lp_a,ry(i)*lp_b,rz(i)*lp_c   
    ELSEIF(imass(i) == 3) THEN
       Write(10,109)'Ti ', i
       WRITE(10,111) rx(i)*lp_a,ry(i)*lp_b,rz(i)*lp_c 
    ENDIF
    
 ENDDO        
109  format (a,1x,i10)
111   format (3(f18.5,1x))

END Subroutine dlpoly_strfile     
       
       
!****************************************

Subroutine Attach_molecule_embed_precipitate
USE variables
IMPLICIT NONE
INTEGER:: num_atoms_mol,k,j,ij
REAL::lp_a_p,lp_b_p,lp_c_p
REAL,dimension(20000)::mx,my,mz
!OPEN(UNIT=2,STATUS='OLD',FILE='molecule') ! defines 'molecule'
!OPEN(UNIT=3,STATUS='OLD',FILE='molecule.precipitate') ! defines 'precipitate-molecule' to be embedded in a matrix of 'molecule'
OPEN(UNIT=7,STATUS='unknown',FILE='structure.xyz') ! xyz format structure file

DO i = 1, num_lp
OPEN(UNIT=2,STATUS='OLD',FILE='molecule') ! defines 'molecule'
OPEN(UNIT=3,STATUS='OLD',FILE='molecule.precipitate')

!IF (ai(i).gt.(nx/2)) THEN
IF ((ai(i).ge.(nx/3).and.ai(i).le.(2*nx/3)).and.(ak(i).ge.(nz/3).and.ak(i).le.(2*nz/3))) THEN
  
   READ(3,*)num_atoms_mol,ntypes, lp_a_p,lp_b_p,lp_c_p   !number of atoms in a molecule, lp_a_p: lp of precipitate
     
    IF(ntypes == 2) THEN
     READ(3,*)mass_a,mass_b
    ELSEIF(ntypes == 3) THEN
     READ(3,*)mass_a,mass_b,mass_c
    ENDIF
    DO k=1,num_atoms_mol
      READ(3,*)mx(k),my(k),mz(k), charge(k),imass(k)
    ENDDO
  close (3)  
 ELSEIF ((ai(i).lt.(nx/3).or.ai(i).gt.(2*nx/3)).and.(ak(i).lt.(nz/3).or.ak(i).gt.(2*nz/3))) THEN 
!ELSEIF (ai(i).le.(nx/2)) THEN
  READ(2,*)num_atoms_mol,ntypes,lp_a,lp_b,lp_c     !number of atoms in a molecule
    IF(ntypes == 2) THEN
     READ(2,*)mass_a,mass_b
    ELSEIF(ntypes == 3) THEN
     READ(2,*)mass_a,mass_b,mass_c
    ENDIF
    DO k=1,num_atoms_mol
      READ(2,*)mx(k),my(k),mz(k), charge(k),imass(k)
    ENDDO
close (2)
 ENDIF 

 DO j=1,num_atoms_mol
  natoms=natoms+1
  rx(natoms) = mx(j) + ai(i)
  ry(natoms) = my(j) + aj(i)
  rz(natoms) = mz(j) + ak(i)
  charge(natoms) = charge(j)
  tag_atom (natoms) = tag_cell (i)   ! tagging the atoms to the unit cell
  ident(natoms) = natoms
  imass(natoms) = imass(j)
 ENDDO
ENDDO

WRITE(7,*)natoms
WRITE(7,*)' '
Do ij = 1, natoms
  IF (imass(ij) == 1) WRITE(7,106)'Ce',rx(ij)*lp_a,ry(ij)*lp_b,rz(ij)*lp_c
  IF (imass(ij) == 2) WRITE(7,106)'O',rx(ij)*lp_a,ry(ij)*lp_b,rz(ij)*lp_c
  IF (imass(ij) == 3) WRITE(7,106)'Sr',rx(ij)*lp_a_p,ry(ij)*lp_b_p,rz(ij)*lp_c_p
  IF (imass(ij) == 4) WRITE(7,106)'Ti',rx(ij)*lp_a_p,ry(ij)*lp_b_p,rz(ij)*lp_c_p
   IF (imass(ij) == 5) WRITE(7,106)'N',rx(ij)*lp_a_p,ry(ij)*lp_b_p,rz(ij)*lp_c_p
ENDDO
106  format (a,3(f9.4,1x))

                       ! natoms: Total number of atoms in the systems by now
102  format (a,1x,3(f5.2,1x))

END Subroutine Attach_molecule_embed_precipitate

!****************************************
   

       
       
       
       
       
       
       
       
