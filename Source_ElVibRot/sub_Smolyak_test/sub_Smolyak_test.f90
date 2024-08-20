!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
! MIT License
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!    Copyright 2015 David Lauvergnat [1]
!      with contributions of
!        Josep Maria Luis (optimization) [2]
!        Ahai Chen (MPI) [1,4]
!        Lucien Dupuy (CRP) [5]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Institut de Química Computacional and Departament de Química,
!        Universitat de Girona, Catalonia, Spain
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!
!===========================================================================
!===========================================================================
SUBROUTINE sub_main_Smolyak_test()

  CALL sub_main_testnD_new()
  !CALL sub_main_testnD_DeltaSRep()

  !CALL sub_main_test3D_old()
  !CALL sub_main_test3D()
  !CALL sub_main_testHpsi()
  !!!!!CALL sub_main_testSmat()

  !CALL sub_main_test()
  !CALL sub_main_test3D_BtoG()

END SUBROUTINE sub_main_Smolyak_test

SUBROUTINE sub_main_testHpsi()
USE EVR_system_m
USE mod_Smolyak_DInd
USE mod_Smolyak_RDP
USE mod_Smolyak_ba
USE mod_Smolyak_test
IMPLICIT NONE

integer :: id

TYPE (Type_AllParam) :: AllPara

real(kind=Rkind) :: Norm,coef

TYPE(TypeRDP), allocatable    :: WPB(:)
TYPE(TypeRDP), allocatable    :: WPG(:)
TYPE(TypeRDP), allocatable    :: V_ON_G(:)
TYPE(TypeRDP), allocatable    :: T2_ON_G(:)
TYPE(TypeRDP), allocatable    :: Temp_ON_G(:)

TYPE(TypeRDP), allocatable    :: Hpsi_ON_G(:)
TYPE(TypeRDP), allocatable    :: Hpsi_ON_B(:)

real(kind=Rkind), allocatable :: H(:,:),Vec(:,:),Ene(:)
integer  :: nb_ba,nb_qa


integer :: iG,ibbb,LGmin,iq_d,jq_d
logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.




read(in_unit,*) AllPara%D

read(in_unit,*) AllPara%LB,AllPara%LG


!-- the 1D-basis ---------------------------
CALL Set_tab_ba(AllPara%tab_ba,AllPara%D,AllPara%LB,AllPara%LG)

!-- the nD basis ---------------------------
write(out_unit,*)
write(out_unit,*) 'Basis ind'
CALL Set_nDInd_01order(AllPara%ind_Basis,AllPara%D,0,AllPara%LB,0)

!-- the numbers of nD-grids ... ------------
write(out_unit,*)
write(out_unit,*) 'Grid ind'
LGmin = AllPara%LG-AllPara%D+1
CALL Set_nDInd_10order(AllPara%ind_Grid,AllPara%D,LGmin,AllPara%LG)


!-- weight of the Smolyak grids ------------
CALL Set_SmolyakWeight(AllPara%WSG,AllPara%ind_Grid(0),AllPara%D,AllPara%LG)

!-- Potential on the grid ------------------
CALL V_ON_GRID(V_ON_G,AllPara)
  !write(out_unit,*) 'coucou potG',ibbb
  !CALL Write_TabRDP_pack(V_ON_G)
  !write(out_unit,*) 'end coucou potG',ibbb


CALL Set_BgG_FOR_id(Hpsi_ON_G,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,0)

!-- Initialization WP on the basis ----------
CALL Set_BgG_FOR_id(WPB,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,AllPara%D)

write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
nb_ba = WPB(1)%n3
nb_qa = sum(Hpsi_ON_G(:)%n2)
write(out_unit,*) 'nb_ba,nb_qa',nb_ba,nb_qa
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='


!TEST B=>G=>B
  !--------- WP0 on the basis ----------------
  ibbb = 1
  WPB(1)%RDP(1,1,:)    = ONE
  WPB(1)%RDP(1,1,ibbb) = TWO

  !--------- WP on the grid ------------------
  CALL time_perso('sub_main_testHpsi: B=>G')
  CALL sub_B_TO_G(WPB,WPG,AllPara)
  CALL time_perso('sub_main_testHpsi: B=>G')

  !--------- WP on the basis -----------------
  CALL time_perso('sub_main_testHpsi: G=>B')
  CALL sub_G_TO_B(WPG,WPB,AllPara)
  CALL time_perso('sub_main_testHpsi: G=>B')



STOP
allocate(H(nb_ba,nb_ba))
H(:,:) = ZERO
allocate(Vec(nb_ba,nb_ba))
allocate(Ene(nb_ba))

write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '=======D:',AllPara%D,'================='
CALL time_perso('sub_main_testHpsi')
write(out_unit,*) '====================================='
flush(out_unit)
write(out_unit,*) 'basis set list:'
DO ibbb=1,AllPara%ind_Basis(AllPara%D+1)%MaxnD
  write(out_unit,*) 'ibb,tab_i',ibbb,':',AllPara%ind_Basis(AllPara%D+1)%tab_ind(:,ibbb)
END DO


DO ibbb=1,AllPara%ind_Basis(AllPara%D+1)%MaxnD
!ibbb=2


  !--------- WP0 on the basis ----------------
  WPB(1)%RDP(1,1,:)    = ZERO
  WPB(1)%RDP(1,1,ibbb) = ONE

  !--------- WP on the grid ------------------
  CALL sub_B_TO_G(WPB,WPG,AllPara)

  write(out_unit,*) 'coucou WPG',ibbb
  CALL Write_TabRDP_pack(WPG)
  write(out_unit,*) 'end coucou WPG',ibbb

  IF (debug) CALL Norm_OF_BgG(WPG,AllPara%WSG,AllPara%ind_Grid(0),AllPara%tab_ba,AllPara%D,AllPara%LG)
  IF (debug) CALL time_perso('sub_main_testHpsi: B=>G')

  !--------- Hpsi on the basis ----------------
  DO iG=1,size(WPG)
    Hpsi_ON_G(iG)%RDP = V_ON_G(iG)%RDP * WPG(iG)%RDP
  END DO
  IF (debug) CALL time_perso('sub_main_testHpsi: Vpsi')


!  DO iq_d=1,AllPara%D
!
!    CALL T2Psi_ON_GRID(Temp_ON_G,WPG,iq_d,AllPara)
!
!    DO iG=1,size(WPG)
!      Hpsi_ON_G(iG)%RDP = Hpsi_ON_G(iG)%RDP -HALF * Temp_ON_G(iG)%RDP
!    END DO
!
!  END DO
!  IF (debug) CALL time_perso('sub_main_testHpsi: All T2psi')
!
!  DO iq_d=1,AllPara%D
!
!    CALL T1Psi_ON_GRID(Temp_ON_G,WPG,iq_d,AllPara)
!
!    DO iG=1,size(WPG)
!      Hpsi_ON_G(iG)%RDP = Hpsi_ON_G(iG)%RDP -HALF*ONETENTH**10 * Temp_ON_G(iG)%RDP
!    END DO
!
!  END DO
!  IF (debug) CALL time_perso('sub_main_testHpsi: All T1psi')


!  DO iq_d=1,AllPara%D
!  DO jq_d=iq_d,AllPara%D
!
!    CALL T11Psi_ON_GRID(Temp_ON_G,WPG,iq_d,jq_d,AllPara)
!
!    DO iG=1,size(WPG)
!      Hpsi_ON_G(iG)%RDP = Hpsi_ON_G(iG)%RDP -HALF*ONETENTH**10 * Temp_ON_G(iG)%RDP
!    END DO
!
!  END DO
!  END DO
!  IF (debug) CALL time_perso('sub_main_testHpsi: All T11psi')

  !--------- WP on the basis -----------------
  CALL sub_G_TO_B(Hpsi_ON_G,Hpsi_ON_B,AllPara)
  IF (debug) CALL Analysis_BbG_ON_basis(Hpsi_ON_B,AllPara)
  IF (debug) CALL time_perso('sub_main_testHpsi: G=>B')

  H(:,ibbb) = Hpsi_ON_B(1)%RDP(1,1,:)
  write(out_unit,*) 'ibb',ibbb,Hpsi_ON_B(1)%RDP(1,1,ibbb)
  flush(out_unit)

END DO

write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
CALL time_perso('sub_main_testHpsi')
write(out_unit,*) '====================================='
flush(out_unit)
write(out_unit,*) 'V matrix'
CALL write_mat(H,6,5)
write(out_unit,*) 'END V matrix'


Vec = H-transpose(H)
write(out_unit,*) 'non symmetric?',maxval(abs(Vec))
CALL diagonalization(H,Ene,Vec,nb_ba,2,1,.TRUE.)
write(out_unit,'(a,10f12.6)') 'Ene HL',Ene(1:min(10,nb_ba))

END SUBROUTINE sub_main_testHpsi

SUBROUTINE sub_main_testSmat()
USE EVR_system_m
USE mod_Smolyak_DInd
USE mod_Smolyak_RDP
USE mod_Smolyak_ba
USE mod_Smolyak_test
IMPLICIT NONE

integer :: id

TYPE (Type_AllParam) :: AllPara

real(kind=Rkind) :: Norm,coef

TYPE(TypeRDP), allocatable    :: WPB(:)
TYPE(TypeRDP), allocatable    :: WPG(:)
TYPE(TypeRDP), allocatable    :: V_ON_G(:)
TYPE(TypeRDP), allocatable    :: T2_ON_G(:)
TYPE(TypeRDP), allocatable    :: Temp_ON_G(:)

TYPE(TypeRDP), allocatable    :: Hpsi_ON_G(:)
TYPE(TypeRDP), allocatable    :: Hpsi_ON_B(:)

real(kind=Rkind), allocatable :: H(:,:),Vec(:,:),Ene(:)
integer  :: nb_ba,nb_qa


integer :: iG,ibbb,LGmin,iq_d,jq_d
!logical, parameter :: debug = .FALSE.
logical, parameter :: debug = .TRUE.




read(in_unit,*) AllPara%D

read(in_unit,*) AllPara%LB,AllPara%LG


write(out_unit,*) 'nb of threads',BasisTOGrid_maxth


!-- the 1D-basis ---------------------------
CALL Set_tab_ba(AllPara%tab_ba,AllPara%D,AllPara%LB,AllPara%LG)

!-- the nD basis ---------------------------
write(out_unit,*)
write(out_unit,*) 'Basis ind'
CALL Set_nDInd_01order(AllPara%ind_Basis,AllPara%D,0,AllPara%LB,0)

!-- the numbers of nD-grids ... ------------
write(out_unit,*)
write(out_unit,*) 'Grid ind'
LGmin = AllPara%LG-AllPara%D+1
CALL Set_nDInd_10order(AllPara%ind_Grid,AllPara%D,LGmin,AllPara%LG)


!-- weight of the Smolyak grids ------------
CALL Set_SmolyakWeight(AllPara%WSG,AllPara%ind_Grid(0),AllPara%D,AllPara%LG)

CALL Set_BgG_FOR_id(Hpsi_ON_G,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,0)


!-- Initialization WP on the basis ----------
CALL Set_BgG_FOR_id(WPB,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,AllPara%D)

write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
nb_ba = WPB(1)%n3
nb_qa = sum(Hpsi_ON_G(:)%n2)
write(out_unit,*) 'nb_ba,nb_qa',nb_ba,nb_qa
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
allocate(H(nb_ba,nb_ba))
H(:,:) = ZERO

write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '=======D:',AllPara%D,'================='
CALL time_perso('sub_main_testSmat')
write(out_unit,*) '====================================='
flush(out_unit)

DO ibbb=1,AllPara%ind_Basis(AllPara%D+1)%MaxnD
!ibbb=2

  !--------- WP0 on the basis ----------------
  WPB(1)%RDP(1,1,:)    = ZERO
  WPB(1)%RDP(1,1,ibbb) = ONE

  !--------- WP on the grid ------------------
  IF (debug) CALL time_perso('sub_main_testSmat: B=>G')
  CALL sub_B_TO_G(WPB,WPG,AllPara)
  !IF (debug) CALL Norm_OF_BgG(WPG,AllPara%WSG,AllPara%ind_Grid(0),AllPara%tab_ba,AllPara%D,AllPara%LG)
  IF (debug) CALL time_perso('sub_main_testSmat: B=>G')

  !--------- Hpsi = WPG ----------------
  DO iG=1,size(WPG)
    Hpsi_ON_G(iG)%RDP = WPG(iG)%RDP
  END DO



  !--------- WP on the basis -----------------
  IF (debug) CALL time_perso('sub_main_testSmat: G=>B')
  CALL sub_G_TO_B(Hpsi_ON_G,Hpsi_ON_B,AllPara)
  IF (debug) CALL time_perso('sub_main_testSmat: G=>B')
  H(:,ibbb)  = Hpsi_ON_B(1)%RDP(1,1,:)
  write(out_unit,*) 'ibb',ibbb
  CALL Write_VecMat(Hpsi_ON_B(1)%RDP(1,1,:),out_unit,5)
  !write(out_unit,*) 'ibb',ibbb,Hpsi_ON_B(1)%RDP(1,1,:)
  flush(out_unit)

END DO
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
CALL time_perso('sub_main_testSmat')
write(out_unit,*) '====================================='
flush(out_unit)


write(out_unit,*) 'H partial'
CALL Write_VecMat(H,out_unit,5)

END SUBROUTINE sub_main_testSmat

SUBROUTINE sub_main_test3D()
USE EVR_system_m
USE mod_Smolyak_DInd
USE mod_Smolyak_RDP
USE mod_Smolyak_ba
USE mod_Smolyak_test
IMPLICIT NONE

integer :: id

TYPE (Type_AllParam) :: AllPara

real(kind=Rkind) :: Norm,coef

integer, allocatable          :: tablb0(:)
real(kind=Rkind), allocatable :: RWPG(:)
TYPE(TypeRDP), allocatable    :: NDPBgG(:)
TYPE(TypeRDP), allocatable    :: NDPBbG(:)
TYPE(TypeRDP), allocatable    :: V_ON_G(:)
TYPE(TypeRDP), allocatable    :: T2_ON_G(:)


integer :: ibb,ibbb,nb_mult_id,nb_BG,nb_not_zero,LGmin




read(in_unit,*) AllPara%D
allocate(tablb0(AllPara%D))

read(in_unit,*) AllPara%LB,AllPara%LG


!-------------------------------------------
!-------------------------------------------
! the 1D-basis
CALL Set_tab_ba(AllPara%tab_ba,AllPara%D,AllPara%LB,AllPara%LG)
!-------------------------------------------
!-------------------------------------------


!-------------------------------------------
!-------------------------------------------
! the nD basis
CALL Set_nDInd_01order(AllPara%ind_Basis,AllPara%D,0,AllPara%LB,0)
!-------------------------------------------
!-------------------------------------------

!-------------------------------------------
!-------------------------------------------
! the numbers of nD-grids ...
write(out_unit,*)
write(out_unit,*) 'Grid ind'
LGmin = AllPara%LG-AllPara%D+1
CALL Set_nDInd_10order(AllPara%ind_Grid,AllPara%D,LGmin,AllPara%LG)

!weight of the Smolyak grids
CALL Set_SmolyakWeight(AllPara%WSG,AllPara%ind_Grid(0),AllPara%D,AllPara%LG)
!-------------------------------------------
!-------------------------------------------
DO ibbb=1,AllPara%ind_Basis(AllPara%D+1)%MaxnD
!ibbb=3

  !-------------------------------------------
  ! WP0

  !--------- WP0 on the basis ----------------
  tablb0(:) = AllPara%ind_Basis(AllPara%D+1)%tab_ind(:,ibbb)

  !--------- WP0 on the grid -----------------
  CALL WP0_ON_GRID(RWPG,tablb0,AllPara%WSG,AllPara%ind_Grid(0),AllPara%tab_ba,AllPara%D,AllPara%LG)

  !--------- transfer of WP0  ----------------
  CALL Set_BgG_FOR_id(NDPBgG,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,0)
  CALL Transfer_WP0_TO_BgG(RWPG,NDPBgG)
  !CALL Norm_OF_BgG(NDPBgG,AllPara%WSG,AllPara%ind_Grid(0),AllPara%tab_ba,AllPara%D,AllPara%LG)
  !CALL Norm_OFF_Diff_WP0_BgG(RWPG,NDPBgG)
  !CALL SumSq_TabRDP(NDPBgG)
  CALL Write_TabRDP(NDPBgG)


  !--------- WP on the basis -----------------
  CALL sub_G_TO_B(NDPBgG,NDPBbG,AllPara)
  !CALL Compare_WP0_BbG_ON_basis(tablb0,NDPBbG,AllPara)
  CALL Write_TabRDP(NDPBgG)

  !--------- WP on the grid ------------------
  CALL sub_B_TO_G(NDPBbG,NDPBgG,AllPara)
  CALL Norm_OFF_Diff_WP0_BgG(RWPG,NDPBgG)
  CALL Norm_OF_BgG(NDPBgG,AllPara%WSG,AllPara%ind_Grid(0),AllPara%tab_ba,AllPara%D,AllPara%LG)

  !--------- WP on the basis -----------------
  CALL sub_G_TO_B(NDPBgG,NDPBbG,AllPara)
  CALL Compare_WP0_BbG_ON_basis(tablb0,NDPBbG,AllPara)


END DO


END SUBROUTINE sub_main_test3D

SUBROUTINE sub_main_test3D_old()
USE EVR_system_m
USE mod_Smolyak_DInd
USE mod_Smolyak_RDP
USE mod_Smolyak_ba
USE mod_Smolyak_test
IMPLICIT NONE

integer :: D,id

integer :: LB,LG

TYPE(TypeDInd),  allocatable :: ind_Grid(:)
TYPE(TypeDInd),  allocatable :: ind_Basis(:)
TYPE(TypeBa),    allocatable :: tab_ba(:,:) ! tab_ba(L,D)

real(kind=Rkind), allocatable :: WSG(:) ! WSG(nG123)


real(kind=Rkind) :: Norm,coef

integer, allocatable          :: tablb0(:)
real(kind=Rkind), allocatable :: RWPG(:)
TYPE(TypeRDP), allocatable    :: NDPBgG(:)
TYPE(TypeRDP), allocatable    :: NDPBbG(:)

integer :: nb_mult,ibb,ibbb,nb_mult_id,nb_BG,nb_not_zero



D=4 !! dimension
LG=5
LB=3
nb_mult = 0

read(in_unit,*) D
allocate(tablb0(D))

read(in_unit,*) LB,LG
!read(in_unit,*) tablb0(:)


!-------------------------------------------
!-------------------------------------------
! the 1D-basis
CALL Set_tab_ba(tab_ba,D,LB,LG)
!-------------------------------------------
!-------------------------------------------


!-------------------------------------------
!-------------------------------------------
! the nD basis
CALL Set_nDInd_01order(ind_Basis,D,0,LB,0)
!-------------------------------------------
!-------------------------------------------

!-------------------------------------------
!-------------------------------------------
! the numbers of nD-grids ...
write(out_unit,*)
write(out_unit,*) 'Grid ind'
CALL Set_nDInd_10order(ind_Grid,D,(LG-D+1),LG)

!weight of the Smolyak grids
CALL Set_SmolyakWeight(WSG,ind_Grid(0),D,LG)
!-------------------------------------------
!-------------------------------------------
DO ibbb=1,ind_Basis(D+1)%MaxnD
!ibbb=3
 tablb0(:) = ind_Basis(D+1)%tab_ind(:,ibbb)
!-------------------------------------------
!-------------------------------------------
! Allocation and initalization of WP on the grid (old way)
CALL WP0_ON_GRID(RWPG,tablb0,WSG,ind_Grid(0),tab_ba,D,LG)
CALL Set_BgG_FOR_id(NDPBgG,ind_Grid,ind_Basis,tab_ba,D,LG,0)
CALL Transfer_WP0_TO_BgG(RWPG,NDPBgG)
!CALL Norm_OF_BgG(NDPBgG,WSG,ind_Grid(0),tab_ba,D,LG)
!CALL Norm_OFF_Diff_WP0_BgG(RWPG,NDPBgG)
CALL SumSq_TabRDP(NDPBgG)
!CALL Write_TabRDP(NDPBgG)
!-------------------------------------------
!-------------------------------------------
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
CALL time_perso('sub_main_test3D: G=>B')
write(out_unit,*) '====================================='

DO id=1,D
  CALL Size_TabRDP(NDPBgG,nb_BG)
  write(out_unit,*) 'id, size NDPBgG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024._Rkind**2)),' MB'
  flush(out_unit)

  CALL BgG_TO_BbG(NDPBgG,NDPBbG,WSG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id,nb_mult_id)

  CALL Size_TabRDP(NDPBbG,nb_BG)
  write(out_unit,*) 'id, size NDPBbG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024._Rkind**2)),' MB'
  write(out_unit,*) 'id, nb_mult_id ',id,nb_mult_id
  flush(out_unit)

  CALL Transfer_BbG_TO_BgG(NDPBbG,NDPBgG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)
END DO
write(out_unit,*) '====================================='
CALL time_perso('sub_main_test3D: G=>B')
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
!-------------------------------------------
! the print DPbbb
!write(out_unit,*) 'shape NDPBgG(1)%RDP',shape(NDPBgG(1)%RDP)
!Norm = ZERO
!nb_not_zero = 0
!DO ibb=1,NDPBgG(1)%n3 ! nbb (the others are 1)
!  coef = NDPBgG(1)%RDP(1,1,ibb)
!  Norm = Norm + coef**2
!  IF (abs(coef) > 1.d-10) THEN
!    nb_not_zero = nb_not_zero + 1
!    write(out_unit,*) 'WP',ind_Basis(D+1)%tab_ind(:,ibb),coef
!    IF (sum(abs(tablb0-ind_Basis(D+1)%tab_ind(:,ibb))) == 0) THEN
!      write(out_unit,*) 'WP0 OK',ibb,'/',shape(NDPBgG(1)%RDP)
!    ELSE
!      write(out_unit,*) 'WP0 NOT OK',ibb,'/',shape(NDPBgG(1)%RDP)
!    END IF
!  END IF
!END DO
!IF (nb_not_zero == 0) THEN
!  write(out_unit,*) 'WP0 NOT OK (0)',ibb,'/',shape(NDPBgG(1)%RDP)
!END IF
!write(out_unit,*) 'Norm',Norm

write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
CALL time_perso('sub_main_test3D: B=>G')
write(out_unit,*) '====================================='

DO id=D,1,-1
  CALL Size_TabRDP(NDPBgG,nb_BG)
  write(out_unit,*) 'id, size NDPBgG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024._Rkind**2)),' MB'
  flush(out_unit)

  CALL Transfer_BgG_TO_BbG(NDPBgG,NDPBbG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)


  CALL Size_TabRDP(NDPBbG,nb_BG)
  write(out_unit,*) 'id, size NDPBbG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024._Rkind**2)),' MB'
  write(out_unit,*) 'id, nb_mult_id ',id,nb_mult_id
  flush(out_unit)

  CALL BbG_TO_BgG(NDPBbG,NDPBgG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id,nb_mult_id)

END DO
write(out_unit,*) '====================================='
CALL time_perso('sub_main_test3D: B=>G')
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='

CALL Norm_OFF_Diff_WP0_BgG(RWPG,NDPBgG)
CALL Norm_OF_BgG(NDPBgG,WSG,ind_Grid(0),tab_ba,D,LG)


write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
CALL time_perso('sub_main_test3D: G=>B')
write(out_unit,*) '====================================='

DO id=1,D
  CALL Size_TabRDP(NDPBgG,nb_BG)
  write(out_unit,*) 'id, size NDPBgG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024._Rkind**2)),' MB'
  flush(out_unit)

  CALL BgG_TO_BbG(NDPBgG,NDPBbG,WSG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id,nb_mult_id)

  CALL Size_TabRDP(NDPBbG,nb_BG)
  write(out_unit,*) 'id, size NDPBbG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024._Rkind**2)),' MB'
  write(out_unit,*) 'id, nb_mult_id ',id,nb_mult_id
  flush(out_unit)

  CALL Transfer_BbG_TO_BgG(NDPBbG,NDPBgG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)
END DO
write(out_unit,*) '====================================='
CALL time_perso('sub_main_test3D: G=>B')
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
write(out_unit,*) '====================================='
!-------------------------------------------
! the print DPbbb
write(out_unit,*) 'shape NDPBgG(1)%RDP',shape(NDPBgG(1)%RDP)
Norm = ZERO
nb_not_zero = 0
DO ibb=1,NDPBgG(1)%n3 ! nbb (the others are 1)
  coef = NDPBgG(1)%RDP(1,1,ibb)
  Norm = Norm + coef**2
  IF (abs(coef) > 1.d-10) THEN
    nb_not_zero = nb_not_zero + 1
    write(out_unit,*) 'WP',ind_Basis(D+1)%tab_ind(:,ibb),coef
    IF (sum(abs(tablb0-ind_Basis(D+1)%tab_ind(:,ibb))) == 0) THEN
      write(out_unit,*) 'WP0 OK',ibb,'/',shape(NDPBgG(1)%RDP)
    ELSE
      write(out_unit,*) 'WP0 NOT OK',ibb,'/',shape(NDPBgG(1)%RDP)
    END IF
  END IF
END DO
IF (nb_not_zero == 0) THEN
  write(out_unit,*) 'WP0 NOT OK (0)',ibb,'/',shape(NDPBgG(1)%RDP)
END IF
write(out_unit,*) 'Norm',Norm



END DO


END SUBROUTINE sub_main_test3D_old
SUBROUTINE sub_main_testnD_new()
USE EVR_system_m
USE mod_Smolyak_DInd
USE mod_Smolyak_RDP
USE mod_Smolyak_ba
USE mod_Smolyak_test
IMPLICIT NONE

integer :: D,id

integer :: LB,LG

TYPE(TypeDInd),  allocatable :: ind_Basis(:)
TYPE(TypeDInd) :: B_nDind,Smolyak_nDind
TYPE(TypeBa),    allocatable :: tab_ba(:,:) ! tab_ba(L,D)
TYPE(TypeBa),    allocatable :: tab_DelatBa(:,:) ! tab_DelatBa(L,D)

real(kind=Rkind), allocatable :: WSG(:) ! WSG(nG123)


real(kind=Rkind) :: R,coef

real(kind=Rkind), allocatable :: RPsi(:)
TYPE(Type_SmolyakRep)         :: SRep1,SRep2
TYPE(Type_SmolyakRep)         :: SRepWeight
TYPE(Type_SmolyakRep)         :: SBRep1,VSRep

integer :: i,iG,nb_B,nb_G,ibb,jbb,nb_Bold
integer, allocatable :: tab_i(:)
real (kind=Rkind), allocatable :: S(:,:),H(:,:),Vec(:,:),Ene(:)
logical :: Grid=.TRUE.
logical :: debug=.FALSE.

D=4 !! dimension
LG=3
LB=3

read(in_unit,*) D

read(in_unit,*) LB,LG


!-------------------------------------------
!-------------------------------------------
! the 1D-basis
CALL Set_tab_ba(tab_ba,D,LG,LG)
CALL Set_tab_DelatBa(tab_DelatBa,D,LG,LG)
!-------------------------------------------
!-------------------------------------------

write(out_unit,*)  'D,LB,LG',D,LB,LG
!-------------------------------------------

write(out_unit,*) 'Smolyak Rep (old way): Basis'
CALL Set_nDInd_01order(ind_Basis,D,0,LB,0)
B_nDind = ind_Basis(D+1)
B_nDind%tab_ind = B_nDind%tab_ind -1
!CALL Write_TypeDInd(B_nDind)
flush(out_unit)
CALL alloc_SmolyakRep(SRep1,B_nDind%tab_ind,tab_DelatBa,Grid=.FALSE.,Delta=.TRUE.)
SRep1 = ZERO
nb_Bold = Size_SmolyakRep(SRep1)
IF (debug) CALL Write_SmolyakRep(Srep1)
write(out_unit,*) 'size basis (old)',nb_Bold


write(out_unit,*)
write(out_unit,*) 'Smolyak ind (new way)'
CALL Set_Smolyak_nDInd(Smolyak_nDind,D,max(0,(LG-D+1)),LG)
!CALL Write_TypeDInd(Smolyak_nDind)
flush(out_unit)


write(out_unit,*) 'Weight of the Smolyak grids'
CALL Set_SmolyakWeight(WSG,Smolyak_nDind,D,LG)
!write(out_unit,*) 'WSG',WSG
!-------------------------------------------
!-------------------------------------------
write(out_unit,*) 'the potential'
VSRep = Set_V_TO_SmolyakRep(Smolyak_nDind%tab_ind,tab_ba)
!write(out_unit,*) 'coucou pot'
!CALL Write_SmolyakRep_pack(VSRep)
!write(out_unit,*) 'end coucou pot'

write(out_unit,*) 'END the potential'
!STOP

!-------------------------------------------
!-------------------------------------------
write(out_unit,*) 'Test B=>G=>B'

!write(out_unit,*) 'Alloc Smolyak Rep'
CALL alloc_SmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_ba,Grid=.FALSE.)
nb_B = Size_SmolyakRep(SRep1)
write(out_unit,*) 'size basis (old)',nb_Bold
write(out_unit,*) 'size Smolyak Rep',nb_B
flush(out_unit)

!=======================!TEST B=>G

tab_i = B_nDind%tab_ind(:,1)
DO i=1,D
  tab_i(i) = l_TO_n(tab_i(i),1)
END DO
SRep1 = ONE
CALL R2_TO_SmolyakRep1_with_tab_i(SRep1,TWO,tab_i,tab_ba,tab_ind=Smolyak_nDind%tab_ind)
SRep2 = SRep1
write(out_unit,*) 'Write Smolyak Rep 1: Basis'
IF (debug) CALL Write_SmolyakRep(Srep1)


CALL time_perso('sub_main_testSmat: B=>G') ; flush(out_unit)
CALL BSmolyakRep_TO_GSmolyakRep_01(SRep1,Smolyak_nDind%tab_ind,tab_ba)
!CALL BSmolyakRep_TO_GSmolyakRep_01_v2(SRep1,Smolyak_nDind%tab_ind,tab_ba)

!CALL BSmolyakRep_TO_GSmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_ba)
CALL time_perso('sub_main_testSmat: B=>G') ; flush(out_unit)
RETURN
nb_G = Size_SmolyakRep(SRep1)

write(out_unit,*) 'Write Smolyak Rep 1: Grid'
IF (debug) CALL Write_SmolyakRep(Srep1)

CALL time_perso('sub_main_testSmat: G=>B') ; flush(out_unit)
!CALL GSmolyakRep_TO_BSmolyakRep_01(SRep1,Smolyak_nDind%tab_ind,tab_ba)
CALL GSmolyakRep_TO_BSmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_ba)
CALL time_perso('sub_main_testSmat: G=>B') ; flush(out_unit)

write(out_unit,*) 'Write Smolyak Rep 1: Basis (after B=>G=>B)'
IF (debug) CALL Write_SmolyakRep(Srep1)

write(out_unit,*) 'Write Smolyak Rep 1: Basis max diff',MaxVal_SmolyakRep(Srep1 - Srep2)

write(out_unit,*) 'size basis (old)',nb_Bold
write(out_unit,*) 'size Smolyak Rep (basis + Grid)',nb_B,nb_G

write(out_unit,*)
write(out_unit,*) 'nb_mult_BTOG,nb_mult_GTOB',nb_mult_BTOG,nb_mult_GTOB

write(out_unit,*) 'END Test B=>G=>B'

RETURN

write(out_unit,*) 'Write Smolyak Rep 1 of SRepWeight: Grid'
SRepWeight = Set_weight_TO_SmolyakRep(Smolyak_nDind%tab_ind,tab_ba)
IF (debug) CALL Write_SmolyakRep(SRepWeight)
!STOP

!======================= basis set list
write(out_unit,*) 'basis set list:'
DO ibb=1,B_nDind%MaxnD
  CALL alloc_SmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_ba)

  tab_i = B_nDind%tab_ind(:,ibb)
  DO i=1,D
    tab_i(i) = l_TO_n(tab_i(i),1)
  END DO
  write(out_unit,*) 'ibb,tab_i',ibb,':',tab_i
END DO

write(out_unit,*) 'TEST: V matrix'

allocate(H(B_nDind%MaxnD,B_nDind%MaxnD))

DO ibb=1,B_nDind%MaxnD
!ibb=2
  CALL alloc_SmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_ba)
  CALL alloc_SmolyakRep(SRep2,Smolyak_nDind%tab_ind,tab_ba)

  tab_i = B_nDind%tab_ind(:,ibb)
  DO i=1,D
    tab_i(i) = l_TO_n(tab_i(i),1)
  END DO
  SRep1 = ZERO
  CALL R2_TO_SmolyakRep1_with_tab_i(SRep1,ONE,tab_i,tab_ba,tab_ind=Smolyak_nDind%tab_ind)

  !write(out_unit,*) 'psi (on basis)',ibb,tab_i
  !CALL Write_SmolyakRep(Srep1)

  CALL BSmolyakRep_TO_GSmolyakRep_01(SRep1,Smolyak_nDind%tab_ind,tab_ba)
  write(out_unit,*) 'coucou psi (on grid)',ibb
  CALL Write_SmolyakRep_pack(Srep1)
  write(out_unit,*) 'end coucou psi (on grid)',ibb

  SRep1 = SRep1 * VSRep

  write(out_unit,*) 'coucou V.psi (on grid)',ibb
  CALL Write_SmolyakRep_pack(Srep1)
  write(out_unit,*) 'end coucou V.psi (on grid)',ibb

  CALL GSmolyakRep_TO_BSmolyakRep_01(SRep1,Smolyak_nDind%tab_ind,tab_ba)

  !write(out_unit,*) 'V.psi (on basis)',ibb,tab_i
  !CALL Write_SmolyakRep(Srep1)

  DO jbb=1,B_nDind%MaxnD

    tab_i = B_nDind%tab_ind(:,jbb)
    DO i=1,D
      tab_i(i) = l_TO_n(tab_i(i),1)
    END DO
    SRep2 = ZERO
    CALL R2_TO_SmolyakRep1_with_tab_i(SRep2,ONE,tab_i,tab_ba,tab_ind=Smolyak_nDind%tab_ind)

    H(ibb,jbb)=dot_product_SmolyakRep(SRep1,SRep2,WSG)

  END DO
END DO
write(out_unit,*) 'V matrix'
CALL write_mat(H,6,5)
write(out_unit,*) 'END V matrix'


Vec = H-transpose(H)
Ene = H(:,1) ! for the allocation
write(out_unit,*) 'non symmetric?',maxval(abs(Vec))
CALL diagonalization(H,Ene,Vec,B_nDind%MaxnD,2,1,.TRUE.)
write(out_unit,'(a,10f12.6)') 'Ene HL',Ene(1:min(10,B_nDind%MaxnD))
RETURN

deallocate(H,Vec,Ene)
write(out_unit,*) 'end TEST: V matrix'

STOP

!======================= The overlap Matrix
write(out_unit,*) 'The overlap matrix'
allocate(S(B_nDind%MaxnD,B_nDind%MaxnD))

DO ibb=1,B_nDind%MaxnD
DO jbb=1,B_nDind%MaxnD

  CALL alloc_SmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_ba)
  CALL alloc_SmolyakRep(SRep2,Smolyak_nDind%tab_ind,tab_ba)

  tab_i = B_nDind%tab_ind(:,ibb)
  DO i=1,D
    tab_i(i) = l_TO_n(tab_i(i),1)
  END DO
  SRep1 = ZERO
  CALL R2_TO_SmolyakRep1_with_tab_i(SRep1,ONE,tab_i,tab_ba,tab_ind=Smolyak_nDind%tab_ind)

  IF (debug) THEN
    write(out_unit,*) 'Write Smolyak Rep 1: Basis'
    CALL Write_SmolyakRep(Srep1)
  END IF

  tab_i = B_nDind%tab_ind(:,jbb)
  DO i=1,D
    tab_i(i) = l_TO_n(tab_i(i),1)
  END DO
  SRep2 = ZERO
  CALL R2_TO_SmolyakRep1_with_tab_i(SRep2,ONE,tab_i,tab_ba,tab_ind=Smolyak_nDind%tab_ind)

  IF (debug) THEN
    write(out_unit,*) 'Write Smolyak Rep 2: Basis'
    CALL Write_SmolyakRep(Srep2)
  END IF

  S(ibb,jbb)=dot_product_SmolyakRep(SRep1,SRep2,WSG)
  IF (debug) write(out_unit,*) 'dot_product: Basis',S(ibb,jbb)

  IF (Grid) THEN

    CALL BSmolyakRep_TO_GSmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_ba)
    CALL BSmolyakRep_TO_GSmolyakRep(SRep2,Smolyak_nDind%tab_ind,tab_ba)

    IF (debug) THEN
      write(out_unit,*) 'Write Smolyak Rep 1: Grid'
      CALL Write_SmolyakRep(Srep1)
      write(out_unit,*) 'Write Smolyak Rep 2: Grid'
      CALL Write_SmolyakRep(Srep2)
    END IF

    S(ibb,jbb)=dot_product_SmolyakRep(SRep1,SRep2*SRepWeight,WSG)
    IF (debug) write(out_unit,*) 'dot_product: Grid',S(ibb,jbb)
  END IF

END DO
END DO

DO ibb=1,B_nDind%MaxnD
  S(ibb,ibb) = S(ibb,ibb) - ONE
END DO
write(out_unit,*) 'max(abs(S-Id)) =',maxval(abs(S))

write(out_unit,*) 'END The overlap matrix'


END SUBROUTINE sub_main_testnD_new
SUBROUTINE sub_main_testnD_DeltaSRep()
USE EVR_system_m
USE mod_Smolyak_DInd
USE mod_Smolyak_RDP
USE mod_Smolyak_ba
USE mod_Smolyak_test
IMPLICIT NONE

integer :: D,id

integer :: LB,LG

TYPE(TypeDInd) :: B_nDind,Smolyak_nDind
TYPE(TypeBa),    allocatable :: tab_DelatBa(:,:) ! tab_DelatBa(L,D)

real(kind=Rkind), allocatable :: WSG(:) ! WSG(nG123)


real(kind=Rkind) :: R,coef

real(kind=Rkind), allocatable :: RPsi(:)
TYPE(Type_SmolyakRep)         :: SRep1,SRep2
TYPE(Type_SmolyakRep)         :: SRepWeight
TYPE(Type_SmolyakRep)         :: SBRep1

integer :: i,iG,nb_B,nb_G,ibb,jbb
integer, allocatable :: tab_i(:)
real (kind=Rkind), allocatable :: S(:,:)
logical :: Grid=.TRUE.
logical :: debug=.FALSE.
!logical :: debug=.TRUE.

D=4 !! dimension
LG=3
LB=3

read(in_unit,*) D

read(in_unit,*) LB,LG


!-------------------------------------------
!-------------------------------------------
! the 1D-basis
CALL Set_tab_DelatBa(tab_DelatBa,D,LB,LG)
!-------------------------------------------
!-------------------------------------------



write(out_unit,*)  'D,LB,LG',D,LB,LG
write(out_unit,*) 'Smolyak ind (old way): Basis'
CALL Set_Smolyak_nDInd(B_nDind,D,0,LG)

!-------------------------------------------
write(out_unit,*)
write(out_unit,*) 'Smolyak ind (new way)'
CALL Set_Smolyak_nDInd(Smolyak_nDind,D,0,LG)

IF (debug) CALL Write_TypeDInd(Smolyak_nDind)
flush(out_unit)

write(out_unit,*) 'Weight of the Smolyak grids: ONE'
allocate(WSG(Smolyak_nDind%MaxnD))
WSG = ONE
!write(out_unit,*) 'WSG',WSG
!-------------------------------------------
!-------------------------------------------

!-------------------------------------------
!-------------------------------------------
!write(out_unit,*) 'Alloc Smolyak Rep'
CALL alloc_SmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_DelatBa,delta=.TRUE.,grid=.FALSE.)
nb_B = Size_SmolyakRep(SRep1)
write(out_unit,*) 'size basis (old)',B_nDind%MaxnD
write(out_unit,*) 'size Smolyak Rep (basis)',nb_B
flush(out_unit)
!TEST B=>G
!allocate(RPsi(nb_B))
!RPsi = ZERO
!RPsi(1) = ONE
!SRep1 = RPsi
!CALL Write_SmolyakRep(Srep1)

tab_i = B_nDind%tab_ind(:,2)
DO i=1,D
  tab_i(i) = l_TO_n(tab_i(i),1)
END DO
SRep1 = ZERO
CALL R2_TO_SmolyakRep1_with_tab_i(SRep1,ONE,tab_i,tab_DelatBa,tab_ind=Smolyak_nDind%tab_ind)
SRep2 = SRep1

IF (debug) THEN
  write(out_unit,*) 'Write Smolyak Rep 1: Basis'
  CALL Write_SmolyakRep(Srep1)
END IF

CALL time_perso('sub_main_testSmat: B=>G')
CALL BSmolyakRep_TO_GSmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_DelatBa)
CALL time_perso('sub_main_testSmat: B=>G')
nb_G = Size_SmolyakRep(SRep1)

CALL time_perso('sub_main_testSmat: G=>B')
CALL GSmolyakRep_TO_BSmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_DelatBa)
CALL time_perso('sub_main_testSmat: G=>B')

write(out_unit,*) 'size basis (old)',B_nDind%MaxnD
write(out_unit,*) 'size Smolyak Rep (basis + Grid)',nb_B,nb_G

write(out_unit,*) 'Write Smolyak Rep 1: Basis (after B=>G=>B)'
IF (debug) CALL Write_SmolyakRep(Srep1)

write(out_unit,*) 'Write Smolyak Rep 1: Basis max diff',MaxVal_SmolyakRep(Srep1 - Srep2)


write(out_unit,*)
write(out_unit,*) 'nb_mult_BTOG,nb_mult_GTOB',nb_mult_BTOG,nb_mult_GTOB
STOP
IF (debug) THEN
  write(out_unit,*) 'Write Smolyak Rep 1: Grid'
  CALL Write_SmolyakRep(Srep1)
END IF




write(out_unit,*) 'Write Smolyak Rep 1 of SRepWeight: Grid'
SRepWeight = Set_weight_TO_SmolyakRep(Smolyak_nDind%tab_ind,tab_DelatBa)
!CALL Write_SmolyakRep(SRepWeight)
!STOP

allocate(S(B_nDind%MaxnD,B_nDind%MaxnD))
allocate(RPsi(nb_B))



DO ibb=1,B_nDind%MaxnD

  CALL alloc_SmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_DelatBa,delta=.TRUE.)

  RPsi      = ZERO
  RPsi(ibb) = ONE
  CALL tabR2_TO_SmolyakRep1(SRep1,RPsi)
  !SRep1     = RPsi  ! problem with ifort 16
  IF (debug) THEN
    write(out_unit,*) 'Write Smolyak Rep 1: Basis'
    CALL Write_SmolyakRep(Srep1)
  END IF

  IF (Grid) CALL BSmolyakRep_TO_GSmolyakRep(SRep1,Smolyak_nDind%tab_ind,tab_DelatBa)

  IF (debug .AND. Grid) THEN
    write(out_unit,*) 'Write Smolyak Rep 1: Grid'
    CALL Write_SmolyakRep(Srep1)
  END IF
  write(out_unit,*) 'ibb',ibb ; flush(out_unit)

  DO jbb=1,B_nDind%MaxnD



    CALL alloc_SmolyakRep(SRep2,Smolyak_nDind%tab_ind,tab_DelatBa,delta=.TRUE.)
    RPsi      = ZERO
    RPsi(jbb) = ONE
    CALL tabR2_TO_SmolyakRep1(SRep2,RPsi)
    !SRep2     = RPsi   ! problem with ifort 16
    IF (debug) THEN
      write(out_unit,*) 'Write Smolyak Rep 2: Basis'
      CALL Write_SmolyakRep(Srep2)
    END IF

    IF (.NOT. Grid) THEN
      S(ibb,jbb)=dot_product_SmolyakRep(SRep1,SRep2,WSG)
      IF (debug) write(out_unit,*) 'dot_product: Basis',S(ibb,jbb)
    ELSE

    CALL BSmolyakRep_TO_GSmolyakRep(SRep2,Smolyak_nDind%tab_ind,tab_DelatBa)

    IF (debug) THEN
      write(out_unit,*) 'Write Smolyak Rep 1: Grid'
      CALL Write_SmolyakRep(Srep1)
      write(out_unit,*) 'Write Smolyak Rep 2: Grid'
      CALL Write_SmolyakRep(Srep2)
    END IF

    S(ibb,jbb)=dot_product_SmolyakRep(SRep1,SRep2*SRepWeight,WSG)
    IF (debug) write(out_unit,*) 'dot_product: Grid',S(ibb,jbb)
  END IF

END DO
END DO

DO ibb=1,B_nDind%MaxnD
  S(ibb,ibb) = S(ibb,ibb) - ONE
END DO
write(out_unit,*) 'max(abs(S-Id)) =',maxval(abs(S))

END SUBROUTINE sub_main_testnD_DeltaSRep
