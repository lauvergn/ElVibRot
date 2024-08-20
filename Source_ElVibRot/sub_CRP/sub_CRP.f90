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
MODULE mod_CRP
  USE EVR_system_m
  IMPLICIT NONE

  integer, private :: type_LU = 3 ! LU decomposition with LAPACK
  !integer, private :: type_LU = 1 ! LU decomposition without LAPACK

  TYPE CRP_Eckart_t
    real (kind=Rkind) :: V0 = 0.0156_Rkind   ! Baloitcha values
    real (kind=Rkind) :: L  = ONE            !  //
    real (kind=Rkind) :: m  = 1060._Rkind    !  //
  END TYPE CRP_Eckart_t
  TYPE CRP_Channel_AT_TS_t
    real (kind=Rkind)              :: EneTS = 0.0105_Rkind
    real (kind=Rkind)              :: w1    = 0.015625_Rkind
    real (kind=Rkind), allocatable :: w(:)
    integer                        :: option = 1
    integer                        :: nb_channels_added = 1
  END TYPE CRP_Channel_AT_TS_t
  TYPE param_CRP
    real (kind=Rkind) :: Ene    = ZERO            ! Total energy for CRP
    real (kind=Rkind) :: DEne   = ZERO            ! Energy increment for the CRP
    integer           :: nb_Ene = 1               ! Number of CRP calculation

    integer           :: iOp_CAP_Reactif = 3      ! Operator index of reactif CAP
    integer           :: iOp_CAP_Product = 4      ! Operator index of product CAP

    integer           :: iOp_Flux_Reactif = 5      ! Operator index of reactif flux
    integer           :: iOp_Flux_Product = 6      ! Operator index of product flux

    character (len=Name_len) :: CRP_Type            = 'lanczos'
    integer                  :: KS_max_it           = 100
    real (kind=Rkind)        :: KS_accuracy         = ONETENTH**5
    character (len=Name_len) :: LinSolv_Type        = 'MatInv'
    integer                  :: LinSolv_max_it      = 100
    real (kind=Rkind)        :: LinSolv_accuracy    = ONETENTH**7
    character (len=Name_len) :: Preconditioner_Type = 'Identity'
    logical                  :: FluxOp_test         = .FALSE.

    logical                  :: With_Eckart         = .FALSE.
    logical                  :: Read_Channel_AT_TS  = .FALSE.

    logical                  :: Build_MatOp         = .FALSE.

    logical                  :: EigenVec_CAPs       = .FALSE.

    TYPE (CRP_Eckart_t)        :: Eckart
    TYPE (CRP_Channel_AT_TS_t) :: Channel_AT_TS

  END TYPE param_CRP

  INTERFACE BlockAna_Mat
    MODULE PROCEDURE BlockAna_RMat,BlockAna_CMat
  END INTERFACE
CONTAINS

SUBROUTINE read_CRP(para_CRP,ny)
USE EVR_system_m
USE mod_Constant
IMPLICIT NONE

  !----- variables pour la namelist analyse ----------------------------
  TYPE (param_CRP),     intent(inout)  :: para_CRP
  integer,              intent(in)     :: ny


  TYPE (REAL_WU) :: Ene,DEne
  integer        :: nb_Ene

  character (len=Name_len) :: CRP_Type            = 'Lanczos'
  integer                  :: KS_max_it           = 100
  real (kind=Rkind)        :: KS_accuracy         = ONETENTH**5
  character (len=Name_len) :: LinSolv_Type        = 'MatInv'
  integer                  :: LinSolv_max_it      = 100
  real (kind=Rkind)        :: LinSolv_accuracy    = ONETENTH**7
  character (len=Name_len) :: Preconditioner_Type = 'Identity'
  logical                  :: FluxOp_test         = .FALSE.
  logical                  :: EigenVec_CAPs       = .FALSE.

    TYPE (CRP_Eckart_t)    :: Eckart ! to be able to compar with Eckart CRP
    logical                :: With_Eckart

    logical                :: Read_Channel

  !----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = "read_CRP"
  !logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------

  NAMELIST /CRP/Ene,DEne,nb_Ene,CRP_Type,                               &
                KS_max_it,KS_accuracy,                                  &
                LinSolv_Type,LinSolv_max_it,LinSolv_accuracy,           &
                Preconditioner_Type,FluxOp_test,EigenVec_CAPs,          &
                Eckart,With_Eckart,Read_Channel

  Ene                 = REAL_WU(ZERO,'cm-1','E')
  DEne                = REAL_WU(ZERO,'cm-1','E')
  nb_Ene              = 1

  CRP_Type            = 'lanczos'
  KS_max_it           = 100
  KS_accuracy         = ONETENTH**5
  LinSolv_Type        = 'QMR'
  LinSolv_max_it      = 100
  LinSolv_accuracy    = ONETENTH**7
  Preconditioner_Type = 'Diag'
  FluxOp_test         = .FALSE.
  With_Eckart         = .FALSE.
  EigenVec_CAPs       = .FALSE.
  Eckart              = CRP_Eckart_t(V0=0.0156_Rkind,L=ONE,m=1060._Rkind)
  Read_Channel        = .FALSE.

  read(in_unit,CRP)
  write(out_unit,CRP)

  para_CRP%With_Eckart = With_Eckart
  IF (With_Eckart) para_CRP%Eckart = Eckart

  CALL string_uppercase_TO_lowercase(CRP_Type)
  CALL string_uppercase_TO_lowercase(LinSolv_Type)
  CALL string_uppercase_TO_lowercase(Preconditioner_Type)

  IF (print_level > 0) write(out_unit,CRP)
  write(out_unit,*)

  para_CRP%Ene                  = convRWU_TO_R_WITH_WorkingUnit(Ene)
  para_CRP%DEne                 = convRWU_TO_R_WITH_WorkingUnit(DEne)
  para_CRP%nb_Ene               = nb_Ene
  para_CRP%CRP_Type             = CRP_Type
  para_CRP%KS_max_it            = KS_max_it
  para_CRP%KS_accuracy          = KS_accuracy
  para_CRP%LinSolv_Type         = LinSolv_Type
  para_CRP%LinSolv_max_it       = LinSolv_max_it
  para_CRP%LinSolv_accuracy     = LinSolv_accuracy
  para_CRP%Preconditioner_Type  = Preconditioner_Type
  para_CRP%FluxOp_test          = FluxOp_test
  para_CRP%EigenVec_CAPs        = EigenVec_CAPs
  para_CRP%Read_Channel_AT_TS   = Read_Channel

  para_CRP%Build_MatOp          = (CRP_type == 'withmat')                  .OR. &
                                  (CRP_type == 'withmat_flux')             .OR. &
                                  (CRP_type == 'withmatspectral')          .OR. &
                                  (CRP_type == 'withmat_test')             .OR. &
         (CRP_type == 'lanczos'        .AND. LinSolv_Type == 'matinv')     .OR. &
         (CRP_type == 'lanczos'        .AND. LinSolv_Type == 'matlinsolv') .OR. &
         (CRP_type == 'lanczos_arpack' .AND. LinSolv_Type == 'matinv')     .OR. &
         (CRP_type == 'lanczos_arpack' .AND. LinSolv_Type == 'matlinsolv') .OR. &
         EigenVec_CAPs

  IF (debug) write(out_unit,*) 'E,DE,nb_E   : ',para_CRP%Ene,para_CRP%DEne,para_CRP%nb_Ene

  IF (Read_Channel) CALL Read_Channel_AT_TS(para_CRP%Channel_AT_TS,ny)

  write(out_unit,*)
  flush(out_unit)

END SUBROUTINE read_CRP
!================================================================
!     CRP
!================================================================
      SUBROUTINE sub_CRP(tab_Op,nb_Op,print_Op,para_CRP)
      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,           intent(in)    :: nb_Op
      TYPE (param_Op),   intent(inout) :: tab_Op(nb_Op)
      logical,           intent(in)    :: print_Op

      TYPE (param_CRP),  intent(in)    :: para_CRP


      integer                           :: i
      real (kind=Rkind)                 :: Ene
      complex(kind=Rkind), allocatable  :: GuessVec(:)

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP'
!-----------------------------------------------------------
      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'shape tab_op',shape(tab_Op)
        flush(out_unit)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      IF (para_CRP%FluxOp_test .AND. nb_Op < 6) Then
        write(out_unit,*) ' The number of operator is wrong'
        write(out_unit,*) ' nb_Op=',nb_Op
        write(out_unit,*) ' For testing the flux, you MUST have 6 or more operators.'
        write(out_unit,*)
        STOP ' ERROR in sub_CRP: wrong operator number'
      END IF
      IF (nb_Op < 4) THEN
        write(out_unit,*) ' The number of operator is wrong'
        write(out_unit,*) ' nb_Op=',nb_Op
        write(out_unit,*) ' You MUST have 4 or more operators.'
        write(out_unit,*) ' You HAVE to set-up: '
        write(out_unit,*) '   - nb_scalar_Op=2 in the &minimum namelist'
        write(out_unit,*) '  or '
        write(out_unit,*) '   - nb_CAP=2 in the &active namelist'
        write(out_unit,*)
        STOP ' ERROR in sub_CRP: wrong operator number'
      END IF

      IF (para_CRP%Build_MatOp) THEN
        CALL sub_MatOp(tab_Op(1),print_Op) ! H
        DO i=3,nb_Op ! for the CAP
          IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product)   &
              CALL sub_MatOp(tab_Op(i),print_Op)
        END DO
      END IF
      IF (tab_Op(1)%Partial_MatOp) STOP 'STOP the Matrices are incomplete'
      IF (para_CRP%EigenVec_CAPs) THEN
        CALL Calc_EigenVec_CAPs(tab_Op,para_CRP)
      END IF

      SELECT CASE (para_CRP%CRP_type)
      CASE ('withmat') ! old one
        CALL sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('withmat_flux')
        CALL sub_CRP_BasisRep_WithMat_flux(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('withmatspectral') ! old one
        CALL sub_CRP_BasisRep_WithMatSpectral(tab_Op,nb_Op,print_Op,para_CRP)
        !CALL sub_CRP_BasisRep_WithMatSpectral_old(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('lanczos') ! lanczos (Lucien Dupuy)

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(para_CRP,tab_Op,nb_Op) &
!$OMP   PRIVATE(i,Ene,GuessVec) &
!$OMP   NUM_THREADS(CRP_maxth)

        CALL alloc_NParray(GuessVec, [tab_Op(1)%nb_tot], 'GuessVec', name_sub)
        GuessVec(:) = ZERO

!$OMP   DO SCHEDULE(STATIC)

        DO i = 0, para_CRP%nb_Ene-1
          Ene = para_CRP%Ene+real(i,kind=Rkind)*para_CRP%DEne
          CALL calc_crp_P_lanczos(tab_Op, nb_Op,para_CRP,Ene,GuessVec)
        END DO

!$OMP   END DO

        CALL dealloc_NParray(GuessVec, 'GuessVec', name_sub)

!$OMP   END PARALLEL

      CASE ('lanczos_arpack') ! lanczos (Lucien Dupuy)

        ! OpenMP ne fonctionne pas avec Arpack
        DO i = 0, para_CRP%nb_Ene-1
          Ene = para_CRP%Ene+real(i,kind=Rkind)*para_CRP%DEne
          CALL calc_crp_IRL(tab_Op, nb_Op,para_CRP,Ene)
        END DO


      CASE ('withmat_test') ! old one
        CALL sub_CRP_BasisRep_WithMat_test(tab_Op,nb_Op,print_Op,para_CRP)
      END SELECT

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
!----------------------------------------------------------

      end subroutine sub_CRP
      SUBROUTINE sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,para_CRP)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP


!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum

!----- working variables -----------------------------
      integer       :: n,i,j,k,ie,nb_col
      real (kind=Rkind), allocatable :: EneH(:),Vec(:,:) ! for debuging


      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMat'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'shape tab_op',shape(tab_Op)
        flush(out_unit)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      write(out_unit,*) 'nb_tot of H',tab_Op(1)%nb_tot
      flush(out_unit)

      IF (debug) THEN
        write(out_unit,*) 'shape H',shape(tab_Op(1)%Rmat)
        CALL alloc_NParray(Vec,shape(tab_Op(1)%Rmat),'Vec',name_sub)
        CALL alloc_NParray(EneH,shape(tab_Op(1)%Rmat(:,1)),'EneH',name_sub)

        CALL sub_diago_H(tab_Op(1)%Rmat,EneH,Vec,tab_Op(1)%nb_tot,.TRUE.)
        write(out_unit,*) 'Ene (ua)',EneH(1:min(10,tab_Op(1)%nb_tot))

        CALL dealloc_NParray(Vec,'Vec',name_sub)
        CALL dealloc_NParray(EneH,'EneH',name_sub)
      END IF


      CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
      CALL alloc_NParray(G,shape(tab_Op(1)%Rmat),'G',name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      IF (debug) THEN
        nb_col = 5
        write(out_unit,*) 'H:'
        CALL Write_Mat(tab_Op(1)%Rmat,out_unit,nb_col)
        write(out_unit,*) 'Reactif CAP:'
        CALL Write_Mat(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,out_unit,nb_col)
        write(out_unit,*) 'Product CAP:'
        CALL Write_Mat(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,out_unit,nb_col)
      END IF

      write(out_unit,*) 'Ginv calc'
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + para_CRP%Ene-para_CRP%DEne
      END DO
      Ene = para_CRP%Ene-para_CRP%DEne

      DO ie=1,para_CRP%nb_Ene

        DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + para_CRP%DEne
        END DO
        Ene = Ene + para_CRP%DEne

        IF (debug) THEN
          nb_col = 5
          write(out_unit,*) 'Ginv:'
          CALL Write_Mat(Ginv,out_unit,nb_col)
        END IF

        G = inv_OF_Mat_TO(Ginv)

        IF (debug) THEN
          nb_col = 5
          write(out_unit,*) 'G:'
          CALL Write_Mat(G,out_unit,nb_col)
        END IF

        !Ginv = matmul(Ginv,G)
        !DO i=1,tab_Op(1)%nb_tot
        !  Ginv(i,i) = Ginv(i,i) - CONE
        !END DO
        !write(out_unit,*) 'id diff ?',maxval(abs(Ginv))

        gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
           matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

        IF (debug) THEN
          nb_col = 5
          write(out_unit,*) 'gGgG:'
          CALL Write_Mat(gGgG,out_unit,nb_col)
        END IF

        RWU_E  = REAL_WU(Ene,'au','E')

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO


        if (para_CRP%With_Eckart) then
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        flush(out_unit)

      END DO


      CALL dealloc_NParray(Ginv,'Ginv',name_sub)
      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,   'G',   name_sub)
!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMat
      SUBROUTINE sub_CRP_BasisRep_WithMat_testblock(tab_Op,nb_Op,print_Op,para_CRP)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP


!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum

!----- working variables -----------------------------
      integer       :: n,i,j,k,ie,nb_col,nbc
      real (kind=Rkind), allocatable :: EneH(:),Vec(:,:) ! for debuging
      integer       :: list_block(2)


      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMat_testblock'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum
      n = tab_Op(1)%nb_tot
      nbc = 65
      list_block=[nbc,n]
      !list_block=[n]

      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'shape tab_op',shape(tab_Op)
        flush(out_unit)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      write(out_unit,*) 'nb_tot of H',tab_Op(1)%nb_tot
      flush(out_unit)

      IF (debug) THEN
        write(out_unit,*) 'shape H',shape(tab_Op(1)%Rmat)
        CALL alloc_NParray(Vec,shape(tab_Op(1)%Rmat),'Vec',name_sub)
        CALL alloc_NParray(EneH,shape(tab_Op(1)%Rmat(:,1)),'EneH',name_sub)

        CALL sub_diago_H(tab_Op(1)%Rmat,EneH,Vec,tab_Op(1)%nb_tot,.TRUE.)
        write(out_unit,*) 'Ene (ua)',EneH(1:min(10,tab_Op(1)%nb_tot))

        CALL dealloc_NParray(Vec,'Vec',name_sub)
        CALL dealloc_NParray(EneH,'EneH',name_sub)
      END IF


      CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
      CALL alloc_NParray(G,shape(tab_Op(1)%Rmat),'G',name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      IF (debug) THEN
        nb_col = 5
        write(out_unit,*) 'H:'
        CALL Write_Mat(tab_Op(1)%Rmat,out_unit,nb_col)
        write(out_unit,*) 'Reactif CAP:'
        CALL Write_Mat(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,out_unit,nb_col)
        write(out_unit,*) 'Product CAP:'
        CALL Write_Mat(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,out_unit,nb_col)
      END IF
      CALL BlockAna_Mat(tab_Op(1)%Rmat,list_block,info='H')
      CALL BlockAna_Mat(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,list_block,info='Reactif CAP')
      CALL BlockAna_Mat(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,list_block,info='Product CAP')

      write(out_unit,*) 'Ginv calc'
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + para_CRP%Ene-para_CRP%DEne
      END DO
      Ene = para_CRP%Ene-para_CRP%DEne

      DO ie=1,para_CRP%nb_Ene

        DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + para_CRP%DEne
        END DO
        Ene = Ene + para_CRP%DEne

        IF (debug) THEN
          nb_col = 5
          write(out_unit,*) 'Ginv:'
          CALL Write_Mat(Ginv,out_unit,nb_col)
        END IF
        CALL BlockAna_Mat(Ginv,list_block,info='Ginv')

        G = inv_OF_Mat_TO(Ginv)

        IF (debug) THEN
          nb_col = 5
          write(out_unit,*) 'G:'
          CALL Write_Mat(G,out_unit,nb_col)
        END IF
        CALL BlockAna_Mat(G,list_block,info='G')

        !Ginv = matmul(Ginv,G)
        !DO i=1,tab_Op(1)%nb_tot
        !  Ginv(i,i) = Ginv(i,i) - CONE
        !END DO
        !write(out_unit,*) 'id diff ?',maxval(abs(Ginv))

        gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
           matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

        IF (debug) THEN
          nb_col = 5
          write(out_unit,*) 'gGgG:'
          CALL Write_Mat(gGgG,out_unit,nb_col)
        END IF
        CALL BlockAna_Mat(gGgG,list_block,info='gGgG')

        RWU_E  = REAL_WU(Ene,'au','E')

        CRP = ZERO
        DO i=1,nbc
          CRP = CRP + gGgG(i,i)
        END DO
        write(out_unit,*) 'CRP at (nbc)',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                  real(CRP,kind=Rkind),aimag(CRP)


        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO


        if (para_CRP%With_Eckart) then
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        flush(out_unit)

      END DO


      CALL dealloc_NParray(Ginv,'Ginv',name_sub)
      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,   'G',   name_sub)
!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMat_testblock
SUBROUTINE sub_CRP_BasisRep_WithMatSpectral(tab_Op,nb_Op,print_Op,para_CRP)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP

      !real (kind=Rkind) :: CRP_Ene,CRP_DEne
      !integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,j,k,ie

      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: VecPGinv(:,:)
      complex (kind=Rkind), allocatable :: ValPGinv(:)
      complex (kind=Rkind), allocatable :: ValPG(:)
      complex (kind=Rkind), allocatable :: Vec1(:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMatSpectral'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'shape tab_op',shape(tab_Op)
        flush(out_unit)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      write(out_unit,*) 'nb_tot of H',tab_Op(1)%nb_tot
      flush(out_unit)

      CALL alloc_NParray(Ginv,    shape(tab_Op(1)%Rmat),'Ginv',    name_sub)
      CALL alloc_NParray(VecPGinv,shape(tab_Op(1)%Rmat),'VecPGinv',name_sub)
      CALL alloc_NParray(ValPGinv,[tab_Op(1)%nb_tot],   'ValPGinv',name_sub)
      CALL alloc_NParray(ValPG,   [tab_Op(1)%nb_tot],   'ValPG',   name_sub)

      write(out_unit,*) 'Ginv calc' ; flush(out_unit)
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      write(out_unit,*) 'Ginv diago' ; flush(out_unit)
      CALL sub_diago_CH(Ginv,ValPGinv,VecPGinv,tab_Op(1)%nb_tot)
      write(out_unit,*) 'Ginv diago: done' ; flush(out_unit)

      CALL dealloc_NParray(Ginv,'Ginv',name_sub)


      CALL alloc_NParray(Vec1,[tab_Op(1)%nb_tot],'Vec1',name_sub)


      Ene = para_CRP%Ene-para_CRP%DEne

      DO ie=1,para_CRP%nb_Ene

        Ene = Ene + para_CRP%DEne

        ValPG(:) = ONE/(ValPGinv+Ene)

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot

          Vec1 = conjg(matmul(VecPGinv,ValPG * VecPGinv(i,:)))

          Vec1 = matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,Vec1)
          Vec1 = matmul(Vec1,VecPGinv) ! equvivalent to matmul(transpose(VecPGinv),Vec1)

          Vec1 = matmul(VecPGinv,ValPG * Vec1)
          CRP = CRP + sum(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat(i,:)*Vec1) ! Cannot use dot_product because of conjg of CAP_Reactif

        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        flush(out_unit)

      END DO

      CALL dealloc_NParray(VecPGinv,'VecPGinv',name_sub)
      CALL dealloc_NParray(ValPGinv,'ValPGinv',name_sub)
      CALL dealloc_NParray(ValPG,   'ValPG',   name_sub)
      CALL dealloc_NParray(Vec1,   'Vec1',   name_sub)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMatSpectral
SUBROUTINE sub_CRP_BasisRep_WithMat_test(tab_Op,nb_Op,print_Op,para_CRP)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP

      !real (kind=Rkind) :: CRP_Ene,CRP_DEne
      !integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,j,k,ie,iVecPro,iVecRea,nb_VecPro

      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E

      real (kind=Rkind), allocatable :: VecPro(:,:),ValPro(:)
      real (kind=Rkind), allocatable :: VecRea(:,:),ValRea(:)
      complex (kind=Rkind), allocatable :: Vec(:),Vec2(:)
      complex (kind=Rkind), allocatable :: CRP_Mat(:,:)
      complex (kind=Rkind), allocatable :: CRP_Mat_inv(:,:)
      complex (kind=Rkind) :: Coef_iVecRea
      real (kind=Rkind) :: Thresh = ONETENTH**6


!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMat_test'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'shape tab_op',shape(tab_Op)
        flush(out_unit)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      write(out_unit,*) 'nb_tot of H',tab_Op(1)%nb_tot
      flush(out_unit)

      CALL alloc_NParray(VecPro, shape(tab_Op(1)%Rmat),'VecPro',   name_sub)
      CALL alloc_NParray(VecRea, shape(tab_Op(1)%Rmat),'VecRea',   name_sub)
      CALL alloc_NParray(ValPro, [tab_Op(1)%nb_tot],   'ValPro',   name_sub)
      CALL alloc_NParray(ValRea, [tab_Op(1)%nb_tot],   'ValRea',   name_sub)
      CALL alloc_NParray(Vec,    [tab_Op(1)%nb_tot],   'Vec',      name_sub)
      CALL alloc_NParray(Vec2,   [tab_Op(1)%nb_tot],   'Vec2',     name_sub)

      CALL diagonalization(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
                           ValRea,VecRea,tab_Op(1)%nb_tot,3,-1,.FALSE.)


      CALL diagonalization(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,               &
                           ValPro,VecPro,tab_Op(1)%nb_tot,3,-1,.FALSE.)
      nb_VecPro = count(ValPro > Thresh)
      CALL alloc_NParray(CRP_Mat,[nb_VecPro,nb_VecPro],'CRP_Mat',name_sub)


      CALL alloc_NParray(Ginv,    shape(tab_Op(1)%Rmat),'Ginv',    name_sub)
      CALL alloc_NParray(G,       shape(tab_Op(1)%Rmat),'G',       name_sub)

      write(out_unit,*) 'Ginv calc' ; flush(out_unit)
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

        !gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
        !   matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

      Ene = para_CRP%Ene-para_CRP%DEne

      DO ie=1,para_CRP%nb_Ene

        Ene = Ene + para_CRP%DEne

        DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + para_CRP%DEne
        END DO
        G = inv_OF_Mat_TO(Ginv)

        CRP = ZERO
        DO iVecPro=1,tab_Op(1)%nb_tot
          IF (abs(ValPro(iVecPro))< Thresh) CYCLE

          Vec = matmul(conjg(G),VecPro(:,iVecPro))
          Vec = matmul(G,ValPro * matmul(transpose(VecPro),Vec))

          Vec2(:) = CZERO
          DO iVecRea=1,tab_Op(1)%nb_tot
            IF (abs(ValRea(iVecRea))< Thresh) CYCLE
            Coef_iVecRea = ValRea(iVecRea) * sum(VecRea(:,iVecRea)*Vec)
            Vec2(:) = Vec2(:) + VecRea(:,iVecRea) * Coef_iVecRea
          END DO
          CRP_Mat(:,iVecPro) = matmul(transpose(VecPro(:,1:nb_VecPro)),Vec2)

        END DO

        DO iVecPro=1,nb_VecPro
          CRP = CRP + CRP_Mat(iVecPro,iVecPro)
        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        flush(out_unit)

      END DO

      CALL dealloc_NParray(VecPro,'VecPro',   name_sub)
      CALL dealloc_NParray(VecRea,'VecRea',   name_sub)
      CALL dealloc_NParray(ValPro,'ValPro',   name_sub)
      CALL dealloc_NParray(ValRea,'ValRea',   name_sub)
      CALL dealloc_NParray(Vec,   'Vec',      name_sub)
      CALL dealloc_NParray(Vec2,  'Vec2',     name_sub)

      CALL dealloc_NParray(Ginv,  'Ginv',     name_sub)
      CALL dealloc_NParray(G,     'G',        name_sub)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMat_test

SUBROUTINE sub_CRP_BasisRep_WithMat_flux(tab_Op,nb_Op,print_Op,para_CRP)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP

      !real (kind=Rkind) :: CRP_Ene,CRP_DEne
      !integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,k,ie


      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E

      real (kind=Rkind), allocatable :: mEYE_FluxOpReactif_mat(:,:)
      real (kind=Rkind), allocatable :: mEYE_FluxOpProduct_mat(:,:)



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMat_flux'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'shape tab_op',shape(tab_Op)
        flush(out_unit)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      write(out_unit,*) 'nb_tot of H',tab_Op(1)%nb_tot


      CALL alloc_NParray(G,shape(tab_Op(1)%Rmat),'G',name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      CALL alloc_NParray(mEYE_FluxOpProduct_mat,shape(tab_Op(1)%Rmat), &
                        'mEYE_FluxOpProduct_mat',name_sub)
      CALL alloc_NParray(mEYE_FluxOpReactif_mat,shape(tab_Op(1)%Rmat), &
                        'mEYE_FluxOpReactif_mat',name_sub)

      DO i=3,nb_Op ! for the flux and the cap
        ! Here,  we don't calculate the flux operator, but -i.FluxOp = [H,HStep]
        ! Because the corresponding matrix is real.
        IF (i == para_CRP%iOp_Flux_Reactif) Then
          write(out_unit,*) 'Op name: ',tab_Op(i)%name_Op
          CALL sub_MatOp(tab_Op(i),print_Op)
          CALL FluxOp_Mat_v0(tab_Op(1),tab_Op(i),mEYE_FluxOpReactif_mat)
          STOP
        END IF
        IF (i == para_CRP%iOp_Flux_Product) Then
          write(out_unit,*) 'Op name: ',tab_Op(i)%name_Op
          CALL sub_MatOp(tab_Op(i),print_Op)
          CALL FluxOp_Mat(tab_Op(1),tab_Op(i),mEYE_FluxOpProduct_mat)
        END IF

      END DO

      Ene = para_CRP%Ene
      DO ie=1,para_CRP%nb_Ene

        CALL G_Mat(tab_Op(1),tab_Op(para_CRP%iOp_CAP_Reactif),                  &
                             tab_Op(para_CRP%iOp_CAP_Product),Ene,G)

        gGgG(:,:) = matmul(mEYE_FluxOpReactif_mat,                              &
                    matmul(G,matmul(mEYE_FluxOpProduct_mat,conjg(G))))

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)

        else
          write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        flush(out_unit)

        Ene = Ene + para_CRP%DEne


      END DO


      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,'G',name_sub)

      CALL dealloc_NParray(mEYE_FluxOpProduct_mat,'mEYE_FluxOpProduct_mat',name_sub)
      CALL dealloc_NParray(mEYE_FluxOpReactif_mat,'mEYE_FluxOpReactif_mat',name_sub)


!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMat_flux
SUBROUTINE calc_crp_p_lanczos(tab_Op,nb_Op,para_CRP,Ene,GuessVec)

      USE mod_Constant
      USE mod_Op
      implicit none

!----- Operator variables ----------------------------------------------
      integer,             intent(in)      :: nb_Op
      TYPE (param_Op),     intent(inout)   :: tab_Op(nb_Op)

      TYPE (param_CRP),    intent(in)      :: para_CRP
      real(kind=Rkind),    intent(in)      :: Ene
      complex(kind=Rkind), intent(inout)   :: GuessVec(:)

!----- working variables -----------------------------
      TYPE(REAL_WU)     :: RWU_E

      ! Calculate the inverse matrix explicitly? For debuging. It is equivalent to sub_CRP_BasisRep_WithMat
      logical, parameter :: Inv = .FALSE.

      ! Size of Hamiltonian matrix
      integer :: ncooked
      ! Loop integers
      integer :: nks, mks, i, j
      ! Vectors for Pmult
      complex(kind=Rkind),  allocatable :: Krylov_vectors(:,:),h(:,:)
      real (kind=Rkind),    allocatable :: Eigvals(:)
      complex(kind=Rkind),  allocatable :: EVec(:,:)

      complex(kind=Rkind) y, len
      real(kind=Rkind) ranr, rani

      ! Operator variables ----------------------------------------------
      logical           :: print_Op


      integer           ::    k,ie
      real (kind=Rkind) :: CRP, oldcrp,crp2,DeltaCRP
      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind), allocatable :: M1(:)

       ! temporary variables for the LU decomposition
       integer,             allocatable :: indx(:)
       complex(kind=Rkind), allocatable :: trav(:)
       complex(kind=Rkind)              :: d

       integer :: lwork,ierr
       real (kind=Rkind), ALLOCATABLE :: work(:)

       TYPE (Time_t) :: CRP_Time
       real(kind=Rkind)  :: RealTime

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc_crp_p_lanczos'
!-----------------------------------------------------------


      ncooked   = tab_Op(1)%nb_tot


    ! If need be, generate explicit matrix representation of operators
    IF ( Inv ) THEN

!$OMP SINGLE
      CALL sub_MatOp(tab_Op(1),print_Op=.FALSE.) ! H
      DO i=3,nb_Op ! for the CAP
        IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product)   &
           CALL sub_MatOp(tab_Op(i),print_Op=.FALSE.)
      END DO
!$OMP END SINGLE

      CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
      CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
      END DO

      G = inv_OF_Mat_TO(Ginv)

      gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,                 &
         matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

      crp2 = ZERO
      do i=1,tab_Op(1)%nb_tot
        crp2 = crp2 + real( gGgG(i,i), kind=Rkind)
      end do

      CALL dealloc_NParray(Ginv,'Ginv',name_sub)
      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,'G',name_sub)

      write(out_unit,*) 'CRP at E (ua)', Ene, crp,'CRP with explicit inversion =', crp2
    ELSE
      RealTime = Delta_RealTime(CRP_Time)

      CALL alloc_NParray(Krylov_vectors,[tab_Op(1)%nb_tot,para_CRP%KS_max_it], &
                        'Krylov_vectors',name_sub,tab_lb=[1,0])

      CALL alloc_NParray(h,      [para_CRP%KS_max_it,para_CRP%KS_max_it],'h',name_sub)
      CALL alloc_NParray(Eigvals,[para_CRP%KS_max_it],                   'Eigvals',name_sub)

      ! Generate first Krylov vector randomly or from a guess (previous energy iteration)
      IF (size(GuessVec) /= tab_Op(1)%nb_tot) THEN
       write(out_unit,*) ' ERROR in',name_sub
       write(out_unit,*) '  The GuessVec size is wrong: ',size(GuessVec)
       write(out_unit,*) '  H%nb_tot:                   ',tab_Op(1)%nb_tot
       write(out_unit,*) ' CHECK the fortran source !'
       STOP ' ERROR in calc_crp_p_lanczos: The GuessVec size is wrong.'
      END IF
      IF (sqrt(dot_product(GuessVec,GuessVec)) == 0) THEN
        write(out_unit,*) '  Random vector'
        CALL Random_CplxVec(GuessVec)
      END IF
      Krylov_vectors(:,0) = GuessVec


      IF (para_CRP%LinSolv_type == 'matinv') THEN

        CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
        CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
        CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

        Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                  tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

       DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + Ene
       END DO

        G = inv_OF_Mat_TO(Ginv)

        gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
           matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

       CALL dealloc_NParray(Ginv,'Ginv',name_sub)
       CALL dealloc_NParray(G,   'G',   name_sub)
     ELSE IF (para_CRP%LinSolv_type == 'matlinsolv') THEN

       CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
       CALL alloc_NParray(trav,[tab_Op(1)%nb_tot],'trav',name_sub)
       CALL alloc_NParray(indx,[tab_Op(1)%nb_tot],'indx',name_sub)

       Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                 tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

       DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + Ene
       END DO

       CALL LU_decomp(Ginv,tab_Op(1)%nb_tot,indx,d,type_LU) ! lapack
       !CALL ludcmp_cplx(Ginv,tab_Op(1)%nb_tot,trav,indx,d)
       !CALL ZGETRF(tab_Op(1)%nb_tot,tab_Op(1)%nb_tot,Ginv,tab_Op(1)%nb_tot,indx,ierr)
       !IF (ierr /= 0) STOP 'LU decomposition'
       !now in Ginv we have its LU decomposition

       CALL dealloc_NParray(trav,'trav',name_sub)
     ELSE IF (para_CRP%LinSolv_type == 'qmr' .OR. para_CRP%LinSolv_type == 'gmres') THEN
       CALL alloc_NParray(M1,[tab_Op(1)%nb_tot],'M1',name_sub)

       IF (allocated(tab_Op(1)%BasisnD%EneH0)) THEN
         M1(:) = ONE/(Ene-tab_Op(1)%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
         write(out_unit,*) 'precon /= 1. DML'
       ELSE
         M1(:)        = CONE
         write(out_unit,*) 'precon = 1. DML'
       END IF
       !M1(:)        = CONE
       !write(out_unit,*) 'precon = 1. DML'
     END IF

      ! Begin Lanczos scheme
      oldcrp = ZERO
      do nks=1,para_CRP%KS_max_it

         IF (print_level > 1) then
           write(out_unit,*) '######################'
           write(out_unit,*) '# in KS iterations, n=',nks
           write(out_unit,*) '# before p_multiply'
           flush(out_unit)
         end if

         SELECT CASE ( para_CRP%LinSolv_type )

         CASE ( 'matinv' )

            Krylov_vectors(:,nks) = matmul(gGgG,Krylov_vectors(:,nks-1))

         CASE ( 'matlinsolv')

            call p_multiplyLU(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks),    &
                              tab_Op,nb_Op,Ene,tab_Op(1)%nb_tot,Ginv,indx,      &
                              para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

         CASE ( 'qmr' )

            call p_multiplyQMR(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks),   &
                         tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

         CASE ( 'gmres' )
#if __CERFACS == 1
            call p_multiplyGMRES(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks), &
                         tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

#else
           write(out_unit,*) ' ERROR in',name_sub
           write(out_unit,*) '  CERFACS GMRES is not implemented.'
           write(out_unit,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR CERFACS GMRES is not implemented'
#endif
         CASE Default
           write(out_unit,*) ' ERROR in',name_sub
           write(out_unit,*) '  No Default for LinSolv_type:',para_CRP%LinSolv_type
           write(out_unit,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR No Default for LinSolv_type'
         END SELECT
         flush(out_unit)

         ! Calculate matrix
         IF (debug) write(out_unit,*) '# in KS iterations, buiding h'
         do mks = 0, nks-1
            h(mks+1, nks) = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
            h(nks, mks+1) = conjg(h(mks+1,nks))
         end do
         IF (debug) write(out_unit,*) '# in KS iterations, h'
         IF (debug) CALL Write_Mat(h(1:nks,1:nks),out_unit,5)

         ! Orthogonalize vectors (twice)
         IF (debug) write(out_unit,*) '# in KS iterations: Orthogonalize the vectors'
         do mks = 0, nks-1
            y = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
            Krylov_vectors(:,nks) = Krylov_vectors(:,nks) - y * Krylov_vectors(:,mks)
         end do
         ! Normalize vector
         CALL ReNorm_CplxVec(Krylov_vectors(:,nks))

         do mks = 0, nks-1
          y = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
          Krylov_vectors(:,nks) = Krylov_vectors(:,nks) - y * Krylov_vectors(:,mks)
         end do
         ! Normalize vector
         CALL ReNorm_CplxVec(Krylov_vectors(:,nks))

       if (nks > 1) then

            Eigvals(:) = ZERO
            IF (allocated(EVec)) CALL dealloc_NParray(EVec,'EVec',name_sub)
            CALL alloc_NParray(EVec,[nks,nks],'EVec',name_sub)
            CALL diagonalization(h(1:nks,1:nks),Eigvals(1:nks),EVec,3,0,.FALSE.)
            IF (debug) write(out_unit,*) '# in KS iterations, Eigvals',Eigvals(1:nks)
            !write(out_unit,*) '# in KS iterations, Eigvals',Eigvals(1:nks)

            crp = ZERO
            do i=1,nks
              crp = crp + h(i,i)
            end do

            DeltaCRP = oldcrp-crp
            if (abs(DeltaCRP) < para_CRP%KS_accuracy) EXIT
            oldcrp = crp
         endif

      !write(out_unit,*) 'Krylov_vectors(:,0)',Krylov_vectors(:,0)
      !write(out_unit,*) 'Krylov_vectors(:,1)',Krylov_vectors(:,1) ; stop

      end do

      !actual_iterations = nks
      write(out_unit,*) '# in KS iterations, n=',nks
      write(out_unit,*) 'accuracy: ',DeltaCRP
      IF (nks > para_CRP%KS_max_it .OR. DeltaCRP >= para_CRP%KS_accuracy) THEN
        write(out_unit,*) 'CRP diago, minval: ',sum(Eigvals(1:para_CRP%KS_max_it)),&
                                         minval(abs(Eigvals(1:para_CRP%KS_max_it)))
        write(out_unit,*) 'WARNING: Lanczos did not converged'
        nks = para_CRP%KS_max_it
      ELSE
        write(out_unit,*) 'CRP diago, minval: ',sum(Eigvals(1:nks)),minval(abs(Eigvals(1:nks)))
      END IF

      RWU_E  = REAL_WU(Ene,'au','E')
      if (para_CRP%With_Eckart) then
        write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                                     CRP,CRP_Eckart(Ene,para_CRP%Eckart),       &
                                     CRP-CRP_Eckart(Ene,para_CRP%Eckart)

      else
        write(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                          CRP
      end if
      RealTime = Delta_RealTime(CRP_Time)
      IF (debug .OR. print_Op .OR. print_level > 0) Then
        write(out_unit,*) 'CRP Energy iteration: Delta Real Time',RealTime
      END IF
      flush(out_unit)

      !CALL Random_CplxVec(GuessVec)
      GuessVec(:) = ZERO
      do mks = 0, nks-1
        GuessVec(:) = GuessVec(:) + Krylov_vectors(:,mks)* sum(EVec(mks+1,:))
      end do
      CALL ReNorm_CplxVec(GuessVec)
      IF (allocated(EVec)) CALL dealloc_NParray(EVec,'EVec',name_sub)


    end if

    IF (allocated(gGgG))           CALL dealloc_NParray(gGgG,'gGgG',name_sub)

    IF (allocated(M1))             CALL dealloc_NParray(M1,'M1',name_sub)
    IF (allocated(Krylov_vectors)) CALL dealloc_NParray(Krylov_vectors,'Krylov_vectors',name_sub)
    IF (allocated(h))              CALL dealloc_NParray(h,'h',name_sub)
    IF (allocated(Eigvals))        CALL dealloc_NParray(Eigvals,'Eigvals',name_sub)

    IF (allocated(indx))           CALL dealloc_NParray(indx,'indx',name_sub)
    IF (allocated(Ginv))           CALL dealloc_NParray(Ginv,'Ginv',name_sub)


END SUBROUTINE calc_crp_p_lanczos

SUBROUTINE calc_crp_IRL(tab_Op,nb_Op,para_CRP,Ene)

!---------------- Attempt at using ARPACK routine ---------------!
!                  Implicitly Restarted Lanczos                  !
!                  Scheme without shift-invert                   !
!----------------------------------------------------------------!

  USE mod_Constant
  USE mod_Op
  IMPLICIT NONE

!----- Operator variables ----------------------------------------------
  INTEGER,            INTENT(in)      :: nb_Op
  TYPE (param_Op),    INTENT(inout)   :: tab_Op(nb_Op)

  TYPE (param_CRP),   INTENT(in)      :: para_CRP
  REAL(kind=Rkind),   INTENT(in)      :: Ene

!----- working variables -----------------------------
  TYPE(REAL_WU)     :: RWU_E
  INTEGER           :: nio

      ! Calculate the inverse matrix explicitly? For debuging. It is equivalent to sub_CRP_BasisRep_WithMat
  LOGICAL, PARAMETER :: Inv = .FALSE.

      ! Size of Hamiltonian matrix
  INTEGER :: ncooked
  ! number of modes
  INTEGER :: ny
  ! number of possible excitations
  INTEGER :: nv

      ! Loop integers
  INTEGER :: i, j


  ! y harmonic frequency
  REAL (kind=Rkind) ::  wy

  ! zero point energy
  REAL (kind=Rkind) :: zpe


      !    IRL local arrays     !
  INTEGER           iparam(11), ipntr(14)
  LOGICAL, ALLOCATABLE ::  SELECT(:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: &
       &                  ax(:), d(:), &
       &                  v(:,:), workd(:), &
       &                  workev(:), resid(:), &
       &                  workl(:)
  REAL (kind=Rkind), ALLOCATABLE ::   rwork(:), rd(:,:)

     !    IRL  local scalars    !
  CHARACTER         bmat*1, which*2
  INTEGER           ido, n, nx, nev, ncv, lworkl, info, &
       &                  ierr, nconv, maxitr, ishfts, mode
  INTEGER ldv
  COMPLEX (kind=Rkind)   sigma
  REAL (kind=Rkind)   tol
  LOGICAL           rvec

     ! BLAS & LAPACK routines used by IRL !
#if __LAPACK == 1
  REAL (kind=Rkind)  :: dznrm2 , dlapy2
  EXTERNAL           :: dznrm2 , zaxpy , dlapy2
#endif

      ! Operator variables ----------------------------------------------
  LOGICAL           :: print_Op

  INTEGER           ::    k,ie
  REAL (kind=Rkind) :: CRP, oldcrp,crp2,DeltaCRP
  COMPLEX (kind=Rkind), ALLOCATABLE :: G(:,:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: Ginv(:,:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: gGgG(:,:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: M1(:)

       ! temporary variables for the LU decomposition
  INTEGER,             ALLOCATABLE :: indx(:)
  COMPLEX(kind=Rkind), ALLOCATABLE :: trav(:)
  COMPLEX(kind=Rkind)              :: dLU

  INTEGER :: lwork
  REAL (kind=Rkind), ALLOCATABLE :: work(:)

  TYPE (Time_t) :: CRP_Time
  real(kind=Rkind)  :: RealTime
  TYPE (Time_t) :: LU_Time
  real(kind=Rkind)  :: RealTime_LU

!----- for debuging --------------------------------------------------
  INTEGER   :: err
  LOGICAL, PARAMETER :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
  CHARACTER (len=*), PARAMETER :: name_sub = 'calc_crp_IRL'
!-----------------------------------------------------------

  RealTime = Delta_RealTime(CRP_Time)

!======================= END OF HEADER ==========================!

  ncooked   = tab_Op(1)%nb_tot

  IF (para_CRP%LinSolv_type == 'matinv') THEN

     CALL alloc_NParray(Ginv,SHAPE(tab_Op(1)%Rmat),'Ginv',name_sub)
     CALL alloc_NParray(G,   SHAPE(tab_Op(1)%Rmat),'G',   name_sub)
     CALL alloc_NParray(gGgG,SHAPE(tab_Op(1)%Rmat),'gGgG',name_sub)

     Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
          tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

     DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
     END DO

     G = inv_OF_Mat_TO(Ginv)

     gGgG(:,:) = MATMUL(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
          MATMUL(G,MATMUL(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,CONJG(G))))

     CALL dealloc_NParray(Ginv,'Ginv',name_sub)
     CALL dealloc_NParray(G,   'G',   name_sub)

      open(newunit=nio, file='P.dat')
      DO i=1,ncooked
         DO j=1, ncooked
            WRITE(nio,*) gGgG(i,j)
         END DO
      END DO
      CLOSE (nio)


     CRP = ZERO
     DO i=1,tab_Op(1)%nb_tot
        CRP = CRP + gGgG(i,i)
     END DO

     RWU_E  = REAL_WU(Ene,'au','E')

     WRITE(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
          CRP

     CALL dealloc_NParray(gGgG,   'gGgG',   name_sub)


  ELSE IF (para_CRP%LinSolv_type == 'matlinsolv') THEN

     CALL alloc_NParray(Ginv,SHAPE(tab_Op(1)%Rmat),'Ginv',name_sub)
     CALL alloc_NParray(trav,[tab_Op(1)%nb_tot],'trav',name_sub)
     CALL alloc_NParray(indx,[tab_Op(1)%nb_tot],'indx',name_sub)

     Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
          tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

     DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
     END DO

     RealTime_LU = Delta_RealTime(LU_Time)

     CALL LU_decomp(Ginv,tab_Op(1)%nb_tot,indx,dLU,type_LU) ! lapack
     !CALL ludcmp_cplx(Ginv,tab_Op(1)%nb_tot,trav,indx,dLU)
     !CALL ZGETRF(tab_Op(1)%nb_tot,tab_Op(1)%nb_tot,Ginv,tab_Op(1)%nb_tot,indx,ierr)
     !IF (ierr /= 0) STOP 'LU decomposition'
     !now in Ginv we have its LU decomposition

     RealTime_LU = Delta_RealTime(LU_Time)
     write(out_unit,*) 'Real Time in LU:',RealTime_LU
     write(out_unit,*) 'LU of Ginv: done' ; flush(out_unit)

     CALL dealloc_NParray(trav,'trav',name_sub)

  ELSE IF (para_CRP%LinSolv_type == 'qmr' .OR. para_CRP%LinSolv_type == 'gmres') THEN
     CALL alloc_NParray(M1,[tab_Op(1)%nb_tot],'M1',name_sub)

     IF (ALLOCATED(tab_Op(1)%BasisnD%EneH0)) THEN
        M1(:) = ONE/(Ene-tab_Op(1)%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
        WRITE(out_unit,*) 'precon /= 1. DML'
     ELSE
        M1(:)        = CONE
        WRITE(out_unit,*) 'precon = 1. DML'
     END IF
  END IF
  flush(out_unit)

!     %--------------------------------------------------%
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%

  n = ncooked
  ldv = n

      ! As we know the number of opened channels at a given energy,
      ! we know how many eigenvalues we want ???

  nev = ChannelNumber_AT_TS(Ene,para_CRP,tab_Op(1))

  ncv   = 2*nev ! recommended in manual

         ! array allocation for IRL
  ALLOCATE ( ax(n), d(ncv), &
       &     v(ldv,ncv), workd(3*n), &
       &     workev(3*ncv), resid(n), &
       &     workl(3*ncv*ncv+5*ncv), &
       &     SELECT(ncv), &
       &     rwork(ncv), rd(ncv,3) &
       & )

      ! Problem type: AX = l x
  bmat  = 'I'
      ! We want highest eigenvalues
  which = 'LM'
!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD  as         |
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         |
!     %---------------------------------------------------%
!
  lworkl  = 3*ncv**2+5*ncv
  tol    = para_CRP%KS_accuracy
  ido    = 0
  info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD .                                           |
!     %---------------------------------------------------%
!
  ishfts = 1
!  max iterations:
  maxitr = para_CRP%KS_max_it
  mode   = 1
!
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode

!----------------------------------------------------------------!
!-----------------------Begin Lanczos scheme---------------------!
!----------------------------------------------------------------!
!----------------------------------------------------------------!
  DO        ! Stop is handled by reverse communication  !
!----------------------------------------------------------------!
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD  and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
#if __ARPACK == 1
     CALL znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv,&
          &        v, ldv, iparam, ipntr, workd, workl, lworkl,&
          &        rwork,info )
#else
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) ' The ARPACK library is not present!'
        write(out_unit,*) "Use CRP_Type='lanczos' instead of CRP_Type='lanczos_Arpack'"
        write(out_unit,*) '  or recompile ElVibRot with ARPACK = 1 (makefile)'
        STOP 'ARPACK has been removed'
#endif
!
     IF (ido .EQ. -1 .OR. ido .EQ. 1) THEN
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | takes workd(ipntr(1)) as the input vector,|
!           |  and return the matrix vector product     |
!           |         to workd(ipntr(2)).               |
!           %-------------------------------------------%
!
        SELECT CASE ( para_CRP%LinSolv_type )

        CASE ( 'matinv' )

           workd(ipntr(2):ipntr(2)+n-1) = MATMUL(gGgG,workd(ipntr(1):ipntr(1)+n-1))

        CASE ( 'matlinsolv')

           CALL p_multiplyLU(workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1),    &
                &  tab_Op,nb_Op,Ene,tab_Op(1)%nb_tot,Ginv,indx,      &
                &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

        CASE ( 'qmr' )

           CALL p_multiplyQMR(workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1),   &
                &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

        CASE ( 'gmres' )
#if __CERFACS == 1
           CALL p_multiplyGMRES(workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1), &
                &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

#else
           WRITE(out_unit,*) ' ERROR in',name_sub
           WRITE(out_unit,*) '  CERFACS GMRES is not implemented.'
           WRITE(out_unit,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR CERFACS GMRES is not implemented'
#endif
        CASE Default
           WRITE(out_unit,*) ' ERROR in',name_sub
           WRITE(out_unit,*) '  No Default for LinSolv_type:',para_CRP%LinSolv_type
           WRITE(out_unit,*) '  You have to choose between: '
           WRITE(out_unit,*) '   "MatInv" or "QMR" or "MatLinSolv".'
           STOP ' ERROR No Default for LinSolv_type'
        END SELECT
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
     ELSE
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
        IF ( info .LT. 0 ) THEN
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD   |
!        %--------------------------%
!
           PRINT *, ' '
           PRINT *, ' Error with _naupd, info = ', info
           PRINT *, ' Check the documentation of _naupd'
           PRINT *, ' '
!
        ELSE
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD .                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
           rvec = .FALSE.
!
#if __ARPACK == 1
           CALL zneupd  (rvec, 'A', SELECT, d, v, ldv, sigma, &
                &        workev, bmat, n, which, nev, tol, resid, ncv, &
                &        v, ldv, iparam, ipntr, workd, workl, lworkl, &
                &        rwork, ierr)
#else
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) ' The ARPACK library is not present!'
        write(out_unit,*) "Use CRP_Type='lanczos' instead of CRP_Type='lanczos_Arpack'"
        write(out_unit,*) '  or recompile ElVibRot with ARPACK = 1 (makefile)'
        STOP 'ARPACK has been removed'
#endif

!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
           IF ( ierr .NE. 0) THEN
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD . |
!           %------------------------------------%
!
              PRINT *, ' '
              PRINT *, ' Error with _neupd, info = ', ierr
              PRINT *, ' Check the documentation of _neupd. '
              PRINT *, ' '
!
           ELSE
!
              nconv = iparam(5)
              DO j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                 SELECT CASE ( para_CRP%LinSolv_type )

                 CASE ( 'matinv' )

                    ax = MATMUL(gGgG,v(:,j))

                 CASE ( 'matlinsolv')

                    CALL p_multiplyLU(v(1,j),ax,    &
                         &  tab_Op,nb_Op,Ene,tab_Op(1)%nb_tot,Ginv,indx,      &
                         &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

                 CASE ( 'qmr' )

                    CALL p_multiplyQMR(v(1,j),ax,   &
                         &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

                 CASE ( 'gmres' )
#if __CERFACS == 1
                    CALL p_multiplyGMRES(v(1,j),ax, &
                         &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

#else
                    WRITE(out_unit,*) ' ERROR in',name_sub
                    WRITE(out_unit,*) '  CERFACS GMRES is not implemented.'
                    WRITE(out_unit,*) '  You have to choose between: "MatInv" or "QMR".'
                    STOP ' ERROR CERFACS GMRES is not implemented'
#endif
                 CASE Default
                    WRITE(out_unit,*) ' ERROR in',name_sub
                    WRITE(out_unit,*) '  No Default for LinSolv_type:',para_CRP%LinSolv_type
                    WRITE(out_unit,*) '  You have to choose between: "MatInv" or "QMR".'
                    STOP ' ERROR No Default for LinSolv_type'
                 END SELECT

#if __LAPACK == 1
                 CALL zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                 rd(j,1) = real(d(j))
                 rd(j,2) = aimag (d(j))
                 rd(j,3) = dznrm2 (n, ax, 1)
                 rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
#else
                 write(out_unit,*) ' ERROR in ',name_sub
                 write(out_unit,*) '  LAPACK is not linked (LAPACK=0 in the makfile).'
                 write(out_unit,*) "Use CRP_Type='lanczos' instead of CRP_Type='lanczos_Arpack'"
                 STOP 'ERROR in calc_crp_IRL: is not linked'
#endif
              END DO
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
#if __ARPACK == 1
              CALL dmout (6, nconv, 3, rd, ncv, -6, &
                   &            'Ritz values (Real, Imag) and relative residuals')
#else
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) ' The ARPACK library is not present!'
        write(out_unit,*) "Use CRP_Type='lanczos' instead of CRP_Type='lanczos_Arpack'"
        write(out_unit,*) '  or recompile ElVibRot with ARPACK = 1 (makefile)'
        STOP 'ARPACK has been removed'
#endif
           END IF
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
           IF ( info .EQ. 1) THEN
              PRINT *, ' '
              PRINT *, ' Maximum number of iterations reached.'
              PRINT *, ' Eigvals '
              DO i=1, nev
                 PRINT *, REAL(D(i))
              END DO
              PRINT *, ' '
              nconv = 0 ! DML 18/01/2021
           ELSE IF ( info .EQ. 3) THEN
              PRINT *, ' '
              PRINT *, ' No shifts could be applied during implicit', &
                   & ' Arnoldi update, try increasing NCV.'
              PRINT *, ' '
           END IF

           PRINT *, ' '
           PRINT *, '_NDRV1'
           PRINT *, '====== '
           PRINT *, ' '
           PRINT *, ' Size of the matrix is ', n
           PRINT *, ' The number of Ritz values requested is ', nev
           PRINT *, ' The number of Arnoldi vectors generated',&
                &   ' (NCV) is ', ncv
           PRINT *, ' What portion of the spectrum: ', which
           PRINT *, ' The number of converged Ritz values is ',&
                &   nconv
           PRINT *, ' The number of Implicit Arnoldi update',&
                &   ' iterations taken is ', iparam(3)
           PRINT *, ' The number of OP*x is ', iparam(9)
           PRINT *, ' The convergence criterion is ', tol
           PRINT *, ' '

        END IF
!
!     %---------------------------%
!     |           Done.           |
!     %---------------------------%
!
        EXIT
     END IF
!----------------------------------------------------------------!
  END DO         !        End of main loop         !
!----------------------------------------------------------------!

  CRP = ZERO
  DO i=1,nconv
     CRP = CRP + real(D(i))
  END DO

  RWU_E  = REAL_WU(Ene,'au','E')

  IF (para_CRP%With_Eckart) THEN
     WRITE(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                                     CRP,CRP_Eckart(Ene,para_CRP%Eckart),       &
                                     CRP-CRP_Eckart(Ene,para_CRP%Eckart)

  ELSE
     WRITE(out_unit,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
          CRP
  END IF
  RealTime = Delta_RealTime(CRP_Time)
  IF (debug .OR. print_Op .OR. print_level > 0) Then
    write(out_unit,*) 'CRP Energy iteration: Delta Real Time',RealTime
  END IF
  flush(out_unit)

  DEALLOCATE ( ax, d, &
       &     v, workd, &
       &     workev, resid, &
       &     workl, &
       &     SELECT, &
       &     rwork, rd &
       & )

  IF (ALLOCATED(gGgG))           CALL dealloc_NParray(gGgG,'gGgG',name_sub)
  IF (ALLOCATED(M1))             CALL dealloc_NParray(M1,'M1',name_sub)
  IF (ALLOCATED(indx))           CALL dealloc_NParray(indx,'indx',name_sub)
  IF (ALLOCATED(Ginv))           CALL dealloc_NParray(Ginv,'Ginv',name_sub)


END SUBROUTINE calc_crp_IRL

SUBROUTINE p_multiplyLU(Vin,Vut,tab_Op,nb_Op,Ene,N,Ginv_LU,indx,                &
                        iOp_CAP_Reactif,iOp_CAP_Product)
      use EVR_system_m
      USE mod_Op
      implicit none

      integer,             intent(in)    :: N
      complex(kind=Rkind), intent(in)    :: Vin(N)
      complex(kind=Rkind), intent(inout) :: Vut(N)
!----- Operator variables ----------------------------------------------
      integer,             intent(in)    :: nb_Op,iOp_CAP_Reactif,iOp_CAP_Product
      TYPE (param_Op),     intent(in)    :: tab_Op(nb_Op)
      real (kind=Rkind),   intent(in)    :: Ene
      complex(kind=Rkind), intent(in)    :: Ginv_LU(N,N)
      integer,             intent(in)    :: indx(N)



      complex(kind=Rkind) :: b(N)


!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub ='p_multiplyLU'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Vin',Vin(:)
        write(out_unit,*)
        flush(out_unit)
      END IF
!-----------------------------------------------------------

!     |b>=e_r|0>
      b(:)=Vin(:)
      call OpOnVec(b,tab_Op(iOp_CAP_Reactif),'NOC')
      IF (debug) write(out_unit,*) 'e_r |Vin>',b(:)

!     |b>=1/(H-E-ie)|b>
      IF (print_level > 1) write(out_unit,*) '# here before LU 1 '
      b(:) = conjg(b)
      !CALL lubksb_cplx(Ginv_LU,N,indx,b)
      !CALL ZGETRS('No transpose',N,1,Ginv_LU,N,indx,B,N,err)
      CALL LU_solve(Ginv_LU,N,indx,b,type_LU) ! here lapack


      b(:) = conjg(b)
      IF (debug) write(out_unit,*) '1/(H-E-ie)|b>',b(:)

!     |b>=e_p|x>
      call OpOnVec(b,tab_Op(iOp_CAP_Product),'NOC')
      IF (debug) write(out_unit,*) 'e_p |b>',b(:)

!     |b>=1/(H-E+ie)|b>
      IF (print_level > 1) write(out_unit,*) '# here before LU 2 '
      !CALL lubksb_cplx(Ginv_LU,N,indx,b)
      !CALL ZGETRS('No transpose',N,1,Ginv_LU,N,indx,B,N,err)
      CALL LU_solve(Ginv_LU,N,indx,b,type_LU) ! here lapack

      Vut(:)=b(:)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Vut',Vut(:)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!-----------------------------------------------------------

END SUBROUTINE p_multiplyLU

SUBROUTINE Gpsi(Vect,tab_Op,nb_Op,Ene,iOp_CAP_Reactif,iOp_CAP_Product,l_conjg)
!SUBROUTINE Gpsi(Vect,tab_Op,nb_Op,Ene,l_conjg)
      use EVR_system_m
      USE mod_psi,     ONLY : param_psi,alloc_psi,dealloc_psi
      USE mod_Op
      implicit none

      integer,             intent(in)           :: nb_Op
      TYPE (param_Op)                           :: tab_Op(nb_Op)
      integer,             intent(in)           :: iOp_CAP_Reactif,iOp_CAP_Product
      complex(kind=Rkind), intent(inout)        :: Vect(tab_Op(1)%nb_tot)
      real(kind=Rkind),    intent(in)           :: Ene
      character(len=3),    intent(in)           :: l_conjg

      TYPE (param_psi)   :: Psi
      TYPE (param_psi)   :: OpPsi
      integer            :: i

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      CALL init_psi(Psi,tab_Op(1),cplx=.TRUE.)
      CALL alloc_psi(Psi,BasisRep=.TRUE.,GridRep=.FALSE.)

      Psi%cvecB(:)=Vect(:)
      OpPsi = Psi

      call sub_OpPsi(Psi,OpPsi,tab_Op(1))

      Vect(:)= Vect(:)*Ene - OpPsi%cvecB(:)


      call sub_OpPsi(Psi,OpPsi,tab_Op(iOp_CAP_Reactif))

      Vect(:) = Vect(:)+EYE*HALF*OpPsi%cvecB(:)

      call sub_OpPsi(Psi,OpPsi,tab_Op(iOp_CAP_Product))

      Vect(:) = Vect(:)+EYE*HALF*OpPsi%cvecB(:)

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      call dealloc_psi(Psi, .TRUE.)
      call dealloc_psi(OpPsi, .TRUE.)

END SUBROUTINE Gpsi

SUBROUTINE G_Mat(H,CAP_Reactif,CAP_Product,Ene,G)
      use EVR_system_m
      USE mod_Op
      implicit none

      TYPE (param_Op),      intent(in)           :: H,CAP_Reactif,CAP_Product
      real(kind=Rkind),     intent(in)           :: Ene

      complex (kind=Rkind), intent(inout)        :: G(H%nb_tot,H%nb_tot)


      complex (kind=Rkind) :: Ginv(H%nb_tot,H%nb_tot)
      integer              :: i

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'G_Mat'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Ene',Ene
        write(out_unit,*)
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      Ginv(:,:) = -H%Rmat + EYE*HALF * (CAP_Reactif%Rmat+CAP_Product%Rmat)

      DO i=1,H%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
      END DO

      G = inv_OF_Mat_TO(Ginv)

      IF (debug) THEN
        Ginv = matmul(Ginv,G)
        DO i=1,H%nb_tot
          Ginv(i,i) = Ginv(i,i) - CONE
        END DO
        write(out_unit,*) 'id diff ?',maxval(abs(Ginv))
    END IF

    IF (debug) THEN
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

END SUBROUTINE G_Mat

SUBROUTINE FluxOp_Mat(H,HStep_Op,FluxOp)
      use EVR_system_m
      USE mod_Op
      implicit none

      TYPE (param_Op)                           :: H,HStep_Op
      real(kind=Rkind),    intent(inout)        :: FluxOp(H%nb_tot,H%nb_tot)



      complex(kind=Rkind)        :: Rvp(H%nb_tot,H%nb_tot)
      complex(kind=Rkind)        :: FluxOp_loc(H%nb_tot,H%nb_tot)

      real(kind=Rkind)           :: Rdiag(H%nb_tot)

      integer :: i,nb_col

      FluxOp = (matmul(H%Rmat,HStep_Op%Rmat) - matmul(HStep_Op%Rmat,H%Rmat))

      FluxOp_loc = EYE*FluxOp

      CALL diagonalization(FluxOp_loc,Rdiag,Rvp,3,2,.TRUE.)

      write(out_unit,*) 'Eigenvalues'
      DO i=1,H%nb_tot
        write(out_unit,*) i,Rdiag(i)
      END DO

      nb_col = 5
      write(out_unit,*) 'Flux eigenvectors in column'
      write(out_unit,*) nb_col,H%nb_tot,H%nb_tot
      CALL Write_Mat(Rvp,out_unit,nb_col)


      write(out_unit,*) 'Ortho ?'
      Rvp = matmul(transpose(Rvp),Rvp)
      CALL Write_Mat(Rvp,out_unit,nb_col)


      !write(out_unit,*) 'Diag',Rdiag

END SUBROUTINE FluxOp_Mat
SUBROUTINE FluxOp_Mat_old(H,HStep_Op,FluxOp)
      use EVR_system_m
      USE mod_Op
      implicit none

      TYPE (param_Op)                           :: H,HStep_Op
      real(kind=Rkind),    intent(inout)        :: FluxOp(H%nb_tot,H%nb_tot)



      real(kind=Rkind)        :: Rdiag(H%nb_tot),Rvp(H%nb_tot,H%nb_tot)
      integer :: i,nb_col


      FluxOp = matmul(H%Rmat,HStep_Op%Rmat) - matmul(HStep_Op%Rmat,H%Rmat)
      !CALL  diagonalization(FluxOp,Rdiag,Rvp,H%nb_tot,4,1,.TRUE.)

      FluxOp = matmul(FluxOp,FluxOp)
      CALL  diagonalization(FluxOp,Rdiag,Rvp,H%nb_tot,3,1,.TRUE.)

      write(out_unit,*) 'Eigenvalues'
      DO i=1,H%nb_tot
        write(out_unit,*) i,Rdiag(i)
      END DO

      nb_col = 5
      write(out_unit,*) 'Flux eigenvectors in column'
      write(out_unit,*) nb_col,H%nb_tot,H%nb_tot
      CALL Write_Mat(Rvp,out_unit,nb_col)


      !write(out_unit,*) 'Ortho ?'
      !Rvp = matmul(transpose(Rvp),Rvp)
      !CALL Write_Mat(Rvp,out_unit,nb_col)


      !write(out_unit,*) 'Diag',Rdiag

END SUBROUTINE FluxOp_Mat_old
SUBROUTINE FluxOp_Mat_v0(H,HStep_Op,FluxOp)
      use EVR_system_m
      USE mod_Op
      implicit none

      TYPE (param_Op)                           :: H,HStep_Op
      real(kind=Rkind),    intent(inout)        :: FluxOp(H%nb_tot,H%nb_tot)



      real(kind=Rkind)        :: Rdiag(H%nb_tot),Rvp(H%nb_tot,H%nb_tot)
      integer :: i,nb_col


      FluxOp = matmul(H%Rmat,HStep_Op%Rmat) - matmul(HStep_Op%Rmat,H%Rmat)
      CALL  diagonalization(FluxOp,Rdiag,Rvp,H%nb_tot,4,1,.TRUE.)

      ! WARNNING: since FluxOp is a skew matrix, they are pairs of  eigenvalues (i*wk, -iwk) ...
      ! and the eigenvectors are V(:,k) + i*V(:,k+1) and V(:,k) - i*V(:,k+1)
      ! its means the V(:,k) and V(:,k+1) are not normalized to one
      DO i=1,H%nb_tot
        Rvp(:,i) = Rvp(:,i)/sqrt(dot_product(Rvp(:,i),Rvp(:,i)))
      end do
      nb_col = 5
      write(out_unit,*) 'Flux eigenvectors in column'
      write(out_unit,*) nb_col,H%nb_tot,H%nb_tot
      CALL Write_Mat(Rvp,out_unit,nb_col)


      !write(out_unit,*) 'Ortho ?'
      !Rvp = matmul(transpose(Rvp),Rvp)
      !CALL Write_Mat(Rvp,out_unit,nb_col)


      !write(out_unit,*) 'Diag',Rdiag

END SUBROUTINE FluxOp_Mat_v0
SUBROUTINE OpOnVec(Vect,tab_Op,l_conjg)
      use EVR_system_m
      USE mod_psi,     ONLY : param_psi,alloc_psi,dealloc_psi
      USE mod_Op
      implicit none

      TYPE (param_Op)     :: tab_Op
      logical             :: print_Op
      logical, parameter  :: cplx=.TRUE.
      TYPE (param_psi)    :: Psi
      TYPE (param_psi)    :: OpPsi
      character(len=3)    :: l_conjg
      complex(kind=Rkind), dimension(tab_Op%nb_tot) :: Vect
      integer         :: i

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      CALL init_psi(Psi,tab_Op,cplx)
      CALL alloc_psi(Psi,BasisRep=.TRUE.,GridRep=.FALSE.)

      Psi%cvecB(:)=Vect(:)
      OpPsi = Psi

      call sub_OpPsi(Psi,OpPsi,tab_Op)

      Vect(:) = OpPsi%cvecB(:)

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      call dealloc_psi(Psi, .TRUE.)
      call dealloc_psi(OpPsi, .TRUE.)
END SUBROUTINE OpOnVec

SUBROUTINE Calc_EigenVec_CAPs(tab_Op,para_CRP)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      TYPE (param_Op)                     :: tab_Op(:)
      TYPE (param_CRP),   intent(in)      :: para_CRP

!----- working variables -----------------------------
      integer                        :: nb_Vec,nb_col
      real (kind=Rkind), allocatable :: Vec(:,:),Val(:),Mat(:,:)

      real (kind=Rkind)              :: Thresh = ONETENTH**10


!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMat_test'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'shape tab_op',shape(tab_Op)
        flush(out_unit)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      write(out_unit,*) 'nb_tot of H',tab_Op(1)%nb_tot
      flush(out_unit)

      CALL alloc_NParray(Mat, shape(tab_Op(1)%Rmat),'Mat',   name_sub)
      CALL alloc_NParray(Vec, shape(tab_Op(1)%Rmat),'Vec',   name_sub)
      CALL alloc_NParray(Val, [tab_Op(1)%nb_tot],   'Val',   name_sub)

      Mat(:,:) =  tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat +                       &
                  tab_Op(para_CRP%iOp_CAP_Product)%Rmat

      CALL diagonalization(Mat,Val,Vec,tab_Op(1)%nb_tot,3,-1,.FALSE.)
      nb_Vec = count(Val >= Thresh)

      IF (debug) THEN
        nb_col = 5
        Mat = matmul(transpose(Vec),matmul(mat,Vec))
        write(out_unit,*) 'Reactif+Product CAP: diago?'
        CALL Write_Mat(Mat,out_unit,nb_col)
      END IF

      write(out_unit,*) 'Val',Val
      write(out_unit,*) 'nb_Vec (Eigenvalues>E-10)',nb_Vec
      write(out_unit,*) 'nb_Vec (Eigenvalues>E-8)',count(Val >= ONETENTH**8)
      write(out_unit,*) 'nb_Vec (Eigenvalues>E-6)',count(Val >= ONETENTH**6)

      nb_col = 5
      write(out_unit,*) 'CAP eigenvectors in column'
      write(out_unit,*) nb_col,tab_Op(1)%nb_tot,tab_Op(1)%nb_tot
      CALL Write_Mat(Vec,out_unit,nb_col)

      IF (debug) THEN
        nb_col = 5
        !check if each CAP matrices are diagonal (on the grid they are)
        Mat = matmul(transpose(Vec),matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,Vec))
        write(out_unit,*) 'Reactif CAP: diago?'
        CALL Write_Mat(Mat,out_unit,nb_col)

        Mat = matmul(transpose(Vec),matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,Vec))
        write(out_unit,*) 'Product CAP: diago?'
        CALL Write_Mat(Mat,out_unit,nb_col)
      END IF

      CALL dealloc_NParray(Vec,'Vec',   name_sub)
      CALL dealloc_NParray(Mat,'Mat',   name_sub)
      CALL dealloc_NParray(Val,'Val',   name_sub)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
!----------------------------------------------------------

END SUBROUTINE Calc_EigenVec_CAPs

SUBROUTINE ReNorm_CplxVec(Vect)
      use EVR_system_m
      implicit none

      complex (kind=Rkind), intent(inout) :: Vect(:)

      Vect(:) = Vect(:)/sqrt(dot_product(Vect,Vect))

END SUBROUTINE ReNorm_CplxVec
SUBROUTINE Random_CplxVec(Vect)
      use EVR_system_m
      implicit none

      complex (kind=Rkind), intent(inout) :: Vect(:)

      integer           :: i
      real (kind=Rkind) :: ranr,rani

      ! Generate first Krylov vector randomly
      do i =1,size(Vect)
         CALL random_number(ranr)
         CALL random_number(rani)
         Vect(i) = cmplx(ranr,rani,kind=Rkind)
      end do

      Vect(:) = Vect(:)/sqrt(dot_product(Vect,Vect))

END SUBROUTINE Random_CplxVec
SUBROUTINE SchmidtProjectOut_CplxVec(Vect,tab_Vect)
      use EVR_system_m
      implicit none

      complex (kind=Rkind), intent(inout) :: Vect(:)
      complex (kind=Rkind), intent(in)    :: tab_Vect(:,:)

      integer               :: i
      complex (kind=Rkind)  :: s


    ! Orthogonalize vectors
    do i = lbound(tab_Vect,dim=2),ubound(tab_Vect,dim=2)
      s = dot_product(Vect,tab_Vect(:,i))
      Vect(:) = Vect(:) - s * tab_Vect(:,i)
    end do
    Vect(:) = Vect(:)/sqrt(dot_product(Vect,Vect))

END SUBROUTINE SchmidtProjectOut_CplxVec
FUNCTION CRP_Eckart(E,Eckart)
  USE EVR_system_m
  IMPLICIT NONE
  real (kind=Rkind)                 :: CRP_Eckart
  real (kind=Rkind),    intent(in)  :: E
  TYPE (CRP_Eckart_t),  intent(in)  :: Eckart


      real (kind=Rkind) :: b,c
      !real (kind=Rkind), parameter :: V0=0.0156_Rkind,m=1060._Rkind,L=ONE
      !real (kind=Rkind), parameter :: V0=0.015625_Rkind,m=1061._Rkind,L=ONE

       b = Eckart%L * Pi * sqrt(TWO*Eckart%m*E)
       c = (Pi/TWO) * sqrt(EIGHT * Eckart%V0*Eckart%m*Eckart%L**2 - ONE)

       CRP_Eckart = ONE / (ONE + (cosh(c)/sinh(b))**2)

END FUNCTION CRP_Eckart
FUNCTION combination(nv,ny)

  IMPLICIT NONE
  INTEGER nv, ny
  INTEGER combination
  INTEGER i


  combination=1

  DO i=nv+1,nv+ny-1
     combination = combination*i
  END DO

  DO i=2,ny-1
     combination=combination/i
  END DO

END FUNCTION combination

SUBROUTINE Read_Channel_AT_TS(Channel_AT_TS_var,ny)
USE EVR_system_m
USE mod_RealWithUnit
USE mod_dnSVM
USE mod_nDindex
USE mod_Op
IMPLICIT NONE

  TYPE (CRP_Channel_AT_TS_t),     intent(inout)      :: Channel_AT_TS_var
  integer,                        intent(in)         :: ny

  integer                         :: option,err_unit
  TYPE (REAL_WU)                  :: w1,EneTS
  character (len=Name_len)        :: w_unit = 'cm-1'
  real (kind=Rkind)               :: conv
  integer                         :: nb_channels_added

  NAMELIST / Channel_AT_TS / option,w1,EneTS,w_unit,nb_channels_added


!----- for debuging --------------------------------------------------
      integer   :: err
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'Read_Channel_AT_TS'
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'ny',ny
    flush(out_unit)
  END IF
!-----------------------------------------------------------


  w_unit              = 'cm-1'
  EneTS               = REAL_WU(0.0105_Rkind,'au','E')
  w1                  = REAL_WU(0.015625_Rkind,'au','E')
  option              = 1
  nb_channels_added   = 1

  read(in_unit,Channel_AT_TS)
  write(out_unit,Channel_AT_TS)

  Channel_AT_TS_var%EneTS             = convRWU_TO_RWU(EneTS)
  Channel_AT_TS_var%w1                = convRWU_TO_RWU(w1)
  Channel_AT_TS_var%option            = option
  Channel_AT_TS_var%nb_channels_added = nb_channels_added

  IF (option == 2) Then
    conv = get_Conv_au_TO_unit(quantity='E',Unit=w_unit,err_unit=err_unit)
    IF (err_unit /= 0) STOP 'in Read_Channel_AT_TS: Wrong w_unit !'
    write(out_unit,*) 'For w_unit= "',trim(w_unit),'", conv=',conv

    allocate(Channel_AT_TS_var%w(ny))
    read(in_unit,*) Channel_AT_TS_var%w(:)
    Channel_AT_TS_var%w(:) = Channel_AT_TS_var%w(:)/conv
  END IF


  IF (debug) THEN
    CALL Write_Channel_AT_TS(Channel_AT_TS_var)
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF


END SUBROUTINE Read_Channel_AT_TS
SUBROUTINE Write_Channel_AT_TS(Channel_AT_TS)
USE EVR_system_m
USE mod_Constant
USE mod_dnSVM
USE mod_nDindex
USE mod_Op
IMPLICIT NONE

  TYPE (CRP_Channel_AT_TS_t),     intent(in)      :: Channel_AT_TS

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'Write_Channel_AT_TS'
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  END IF
!-----------------------------------------------------------

    write(out_unit,*) 'option            ',Channel_AT_TS%option
    write(out_unit,*) 'EneTS (au)        ',Channel_AT_TS%EneTS
    write(out_unit,*) 'nb_channels_added ',Channel_AT_TS%nb_channels_added
    SELECT CASE (Channel_AT_TS%option)
    CASE(1)
      write(out_unit,*) 'w1 (au) ',Channel_AT_TS%w1
    CASE(2)
      IF (allocated(Channel_AT_TS%w)) THEN
        write(out_unit,*) 'w(:) (au) ',Channel_AT_TS%w
      ELSE
        write(out_unit,*) 'w(:) is not allocated !'
      END IF
    END SELECT

  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF


END SUBROUTINE Write_Channel_AT_TS
FUNCTION ChannelNumber_AT_TS(Ene,para_CRP,para_H) RESULT(nb_channels)
USE EVR_system_m
USE mod_dnSVM
USE mod_nDindex
USE mod_Op
IMPLICIT NONE

  INTEGER                               :: nb_channels
  TYPE (param_CRP),     intent(in)      :: para_CRP
  real (kind=Rkind),    intent(in)      :: Ene
  TYPE (param_Op),      intent(in)      :: para_H


  integer :: option = 2 ! 1: Lucien Dupuy, one degenerate frequency
                        ! 2: ny frequencies at TS
                        ! 3: general, Energy levels calculation at the TS

  ! parameters for option=1 (Lucien)
  integer             :: i,ny,nv
  real (kind=Rkind)   :: wy,zpe,EneTS

  ! more general parameters for option=2
  TYPE (Type_nDindex)             :: nDindB_Channels
  TYPE (IntVec_t), allocatable :: tab_i_TO_l(:)
  integer,            allocatable :: nbSize(:),tab_ib(:)
  integer                         :: LB,ib,nb,n
  real (kind=Rkind),  allocatable :: w(:)
  real (kind=Rkind)               :: EneChannel

  ! more general parameters for option=3
  TYPE (basis),        pointer    :: basisnD
  integer,            allocatable :: nDval(:)
  real (kind=Rkind)               :: E0_func_of_s


!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'ChannelNumber_AT_TS'
!-----------------------------------------------------------
  basisnD => para_H%para_AllBasis%BasisnD

  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'Ene',Ene
    write(out_unit,*) 'EneTS',para_CRP%Channel_AT_TS%EneTS
    write(out_unit,*)
    IF (allocated(BasisnD%EneH0)) THEN
      write(out_unit,*) 'size BasisnD%EneH0',size(BasisnD%EneH0)
      write(out_unit,*) 'BasisnD%EneH0',BasisnD%EneH0
    END IF
    flush(out_unit)
  END IF
!-----------------------------------------------------------


  ny = para_H%mole%nb_act-1

  SELECT CASE (para_CRP%Channel_AT_TS%option)
  CASE (1) ! Lucien Dupuy, one degenerate frequency
    wy    = para_CRP%Channel_AT_TS%w1
    zpe   = HALF*ny*wy
    EneTS = para_CRP%Channel_AT_TS%EneTS
    nv    = int( (Ene-zpe-EneTS)/wy )

    nb_channels = 1
    IF ( nv > 0 ) THEN
      DO i=1,nv
        nb_channels = nb_channels + combination(i,ny)
      END DO
    END IF

  CASE (2) ! ny frequencies at TS with basis set (SG4)
    allocate(nbSize(para_H%para_AllBasis%BasisnD%nb_basis-1))
    allocate(tab_ib(para_H%para_AllBasis%BasisnD%nb_basis-1))
    EneTS = para_CRP%Channel_AT_TS%EneTS

    LB = para_H%para_AllBasis%BasisnD%L_SparseBasis

    DO ib=2,para_H%para_AllBasis%BasisnD%nb_basis
      nbSize(ib-1) = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nb
    END DO


    allocate(tab_i_TO_l(para_H%para_AllBasis%BasisnD%nb_basis-1))
    DO ib=2,para_H%para_AllBasis%BasisnD%nb_basis
      nb = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nb
      CALL alloc_IntVec(tab_i_TO_l(ib-1),nb)
      IF (para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nb_basis < 1) THEN
        tab_i_TO_l(ib-1)%vec(:) = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nDindB%Tab_L(:)
      ELSE
        DO i=1,nb
          n = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nDindB%Tab_L(i)
          tab_i_TO_l(ib-1)%vec(i) = get_L_FROM_Basis_L_TO_n(para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%L_TO_nb,n)
        END DO
      END IF
    END DO


    nDindB_Channels%packed = .TRUE. ! with false the mapping is too long !!
    CALL init_nDindexPrim(nDindB_Channels,ny,nbSize,                            &
                          type_OF_nDindex=5,Lmax=LB,tab_i_TO_l=tab_i_TO_l)

    nb_channels = para_CRP%Channel_AT_TS%nb_channels_added

    CALL init_nDval_OF_nDindex(nDindB_Channels,tab_ib)
    DO ib=1,nDindB_Channels%Max_nDI
      CALL ADD_ONE_TO_nDindex(nDindB_Channels,tab_ib,iG=ib)
      tab_ib(:) = tab_ib(:)-1
      EneChannel = sum(tab_ib(:)*para_CRP%Channel_AT_TS%w(:)) +                 &
                  HALF*sum(para_CRP%Channel_AT_TS%w)
      IF (debug) write(out_unit,*) 'ib,tab_ib(:)-1',ib,tab_ib,' :',            &
                                    EneChannel,EneTS + EneChannel
      IF (Ene >= EneTS + EneChannel) nb_channels = nb_channels + 1
    END DO


    DO ib=1,size(tab_i_TO_l)
      CALL dealloc_IntVec(tab_i_TO_l(ib))
    END DO
    deallocate(tab_i_TO_l)

    CALL dealloc_nDindex(nDindB_Channels)

    deallocate(nbSize)
    deallocate(tab_ib)

  CASE (3) ! With the energy of the "inactive" basis functions

    !write(out_unit,*) 'size BasisnD%EneH0',size(BasisnD%EneH0)
    !write(out_unit,*) 'BasisnD%EneH0',BasisnD%EneH0(:)
    IF (.NOT. allocated(BasisnD%EneH0))                                         &
        STOP 'ERROR in ChannelNumber_AT_TS: EneH0 is not allocated'


    ! first the energy of BasisnD%tab
    SELECT CASE (BasisnD%SparseGrid_type)
    CASE (0) ! Direct product
      E0_func_of_s = BasisnD%tab_Pbasis(1)%Pbasis%EneH0(1)
      ! IF (debug) Then
      !   write(out_unit,*) 'BasisnD%tab_Pbasis(1)%Pbasis%EneH0',BasisnD%tab_Pbasis(1)%Pbasis%EneH0
      !   write(out_unit,*) 'BasisnD%tab_Pbasis(2)%Pbasis%EneH0',BasisnD%tab_Pbasis(2)%Pbasis%EneH0
      ! END IF
    CASE (1) ! Sparse basis
      E0_func_of_s = BasisnD%tab_basisPrimSG(1,BasisnD%L_SparseBasis)%EneH0(1)
    CASE (2,4) ! Sparse basis
      E0_func_of_s = BasisnD%tab_basisPrimSG(BasisnD%L_SparseBasis,1)%EneH0(1)
    END SELECT
    IF (debug) write(out_unit,*) 'E0_func_of_s',E0_func_of_s

    allocate(nDval(BasisnD%nb_basis))

    EneTS = para_CRP%Channel_AT_TS%EneTS

    nb_channels = para_CRP%Channel_AT_TS%nb_channels_added

    IF (allocated(BasisnD%nDindB_contracted)) THEN
      CALL init_nDval_OF_nDindex(BasisnD%nDindB_contracted,nDval)
      DO ib=1,BasisnD%nDindB_contracted%Max_nDI
        CALL ADD_ONE_TO_nDindex(BasisnD%nDindB_contracted,nDval)
        IF (nDval(1) == 1) Then
          EneChannel = BasisnD%EneH0(ib) - E0_func_of_s
          IF (debug) THEN
            write(out_unit,*) ib,'nDval',nDval
            write(out_unit,*) ib,'EneChannel',EneChannel
          END IF
          write(out_unit,*) 'Ene       ',Ene
          write(out_unit,*) 'EneChannel',EneChannel

          IF (Ene >= EneChannel) nb_channels = nb_channels + 1
        END IF
      END DO
      nb_channels = min(nb_channels,BasisnD%nDindB_contracted%Max_nDI)
    ELSE
      CALL init_nDval_OF_nDindex(BasisnD%nDindB,nDval)
      DO ib=1,BasisnD%nDindB%Max_nDI
        CALL ADD_ONE_TO_nDindex(BasisnD%nDindB,nDval)
        IF (nDval(1) == 1) Then
          EneChannel = BasisnD%EneH0(ib) - E0_func_of_s
          IF (debug) THEN
            write(out_unit,*) ib,'nDval',nDval
            write(out_unit,*) ib,'EneChannel',EneChannel
          END IF
          write(out_unit,*) 'Ene       ',Ene
          write(out_unit,*) 'EneChannel',EneChannel

          IF (Ene >= EneChannel) nb_channels = nb_channels + 1
        END IF
      END DO
      nb_channels = min(nb_channels,BasisnD%nDindB%Max_nDI)
    END IF
    deallocate(nDval)

  END SELECT

  IF (print_level > 1 .OR. debug)  write(out_unit,*) 'nb_channels',nb_channels
  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF

END FUNCTION ChannelNumber_AT_TS

!to be moved

SUBROUTINE BlockAna_RMat(f,list_block,info)
  USE mod_MPI

  character (len=*), optional :: info

  integer,          intent(in) :: list_block(:)
  real(kind=Rkind), intent(in) :: f(:,:)

  integer           :: i,j,ib1,ib2,jb1,jb2
  real(kind=Rkind)  :: valmax

  IF (present(info)) THEN
    write(out_unit,*) 'Block analysis, ',info
  ELSE
    write(out_unit,*) 'Block analysis'
  END IF

  IF (size(list_block) > 1) THEN
    DO i=1,size(list_block)
    DO j=1,size(list_block)

      ib1 = 1
      IF (i > 1) ib1 = list_block(i-1)
      ib2 = list_block(i)

      jb1 = 1
      IF (j > 1) jb1 = list_block(j-1)
      jb2 = list_block(j)

      valmax = maxval(abs(f(ib1:ib2,jb1:jb2)))
      write(out_unit,*) 'block',i,j,valmax

    END DO
    END DO
  ELSE
    valmax = maxval(abs(f))
    i=1
    j=1
    write(out_unit,*) 'block',i,j,valmax
  END IF

END SUBROUTINE BlockAna_RMat
SUBROUTINE BlockAna_CMat(f,list_block,info)
  USE mod_MPI

  character (len=*), optional :: info

  integer,          intent(in) :: list_block(:)
  complex(kind=Rkind), intent(in) :: f(:,:)

  integer           :: i,j,ib1,ib2,jb1,jb2
  real(kind=Rkind)  :: valmax

  IF (present(info)) THEN
    write(out_unit,*) 'Block analysis, ',info
  ELSE
    write(out_unit,*) 'Block analysis'
  END IF

  IF (size(list_block) > 1) THEN
    DO i=1,size(list_block)
    DO j=1,size(list_block)

      ib1 = 1
      IF (i > 1) ib1 = list_block(i-1)
      ib2 = list_block(i)

      jb1 = 1
      IF (j > 1) jb1 = list_block(j-1)
      jb2 = list_block(j)

      valmax = maxval(abs(f(ib1:ib2,jb1:jb2)))
      write(out_unit,*) 'block',i,j,valmax

    END DO
    END DO
  ELSE
    valmax = maxval(abs(f))
    i=1
    j=1
    write(out_unit,*) 'block',i,j,valmax
  END IF

END SUBROUTINE BlockAna_CMat

END MODULE mod_CRP
