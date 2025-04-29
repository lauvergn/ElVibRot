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
MODULE ReadOp_m
  USE EVR_system_m
  USE mod_OpGrid
  USE mod_PrimOp
  IMPLICIT NONE

  PRIVATE

  TYPE, EXTENDS(param_TypeOp) :: SplitOp_t
    integer                          :: nb_SplitOp        =  0
    integer                          :: nb_SplitCoord     =  0
    TYPE (param_TypeOp), allocatable :: tab_TypeOp(:,:) ! size: (nb_SplitCoord,nb_SplitOp)
  END TYPE SplitOp_t

  TYPE, EXTENDS(PrimOp_t) :: ReadOp_t ! used for transfert info from read_active to para_Op

    logical                  :: OpPsi_WithGrid      = .FALSE.

    logical                  :: pack_Op             = .FALSE.
    real (kind=Rkind)        :: tol_pack            = ONETENTH**7
    real (kind=Rkind)        :: tol_nopack          = NINE*ONETENTH
    logical                  :: read_Op             = .FALSE.
  
    logical                  :: make_Mat            = .FALSE.
    logical                  :: save_MatOp          = .FALSE.
    logical                  :: restart_MatOp       = .FALSE.
    integer                  :: Partial_MatOp_i     = 0
    integer                  :: Partial_MatOp_f     = huge(1)
    logical                  :: formatted_Mat       = .TRUE.
    character (len=Line_len) :: name_Mat            = 'MatOp'
  
    logical               :: spectral            = .FALSE.       ! IF T, spectral represention
    integer               :: spectral_Op         = 0             ! IF sepctral, we use the H (nOp=0)
    logical               :: Op_WithContracRVec  = .FALSE.
    logical               :: pot_only            = .FALSE.
    logical               :: T_only              = .FALSE.
    integer               :: nb_bRot             = 0
  
    TYPE (File_tGrid) :: para_FileGrid                       ! parameters to tranfer to OpGrid%...
    TYPE (File_t)     :: FileMat                             ! file Operator Matrix
  
    logical               :: comput_S            = .FALSE.       ! calculation of the active overlap matrix
  
    logical               :: Op_Transfo          = .FALSE.       ! true => we are using Transfo(Op) instead of Op
    real (kind=Rkind)     :: E0_Transfo          = ZERO          ! ScaledOp = Op - E0_transfo * I
    integer               :: degree_Transfo      = -1            ! degree of the transformation
    real (kind=Rkind), allocatable :: Poly_Transfo(:) ! Poly_transfo(0:degree_Transfo)
                                        ! Transfo(Op) = Sum_k Poly_transfo(k) * ScaledOp^k

    ! new operator type
    logical               :: NewOp               = .FALSE.       ! true => read the list of the operator parameters
    integer               :: nb_Op               = -1            ! number of operator
    integer               :: nb_SplitOp          = 0
    integer               :: nb_SplitCoord       = 0
    TYPE(SplitOp_t), allocatable :: tab_SplitOp(:)               ! tab_SplitOp(1:nb_Op)

  CONTAINS
    PROCEDURE, PRIVATE, PASS(ReadOp1) :: ReadOp2_TO_ReadOp1
    GENERIC,   PUBLIC  :: assignment(=) => ReadOp2_TO_ReadOp1
  END TYPE ReadOp_t

  PUBLIC :: SplitOp_t, write_SplitOp, dealloc_SplitOp
  PUBLIC :: ReadOp_t, read_ReadOp, init_ReadOp, dealloc_ReadOp, ReadOp2_TO_ReadOp1_FOR_AutoBasis

CONTAINS
  !===============================================================================
  ! SplitOp part
  !===============================================================================
  SUBROUTINE write_SplitOp(SplitOp)
    USE EVR_system_m
    IMPLICIT NONE

    TYPE (SplitOp_t), intent(in)   :: SplitOp
    integer :: iSO,iSC

    write(out_unit,*) 'BEGINNING write_SplitOp'
    write(out_unit,*) 'nb_SplitOp',SplitOp%nb_SplitOp
    write(out_unit,*) 'nb_SplitCoord',SplitOp%nb_SplitCoord

    IF (SplitOp%nb_SplitOp == 0) THEN
      CALL Write_TypeOp(SplitOp%param_TypeOp)
    ELSE
      IF (allocated(SplitOp%tab_TypeOp)) THEN
        DO iSO=1,SplitOp%nb_SplitOp
          write(out_unit,*) 'SplitOp term #',iSO
          DO iSC=1,SplitOp%nb_SplitCoord
            write(out_unit,*) 'SplitCoord part #',iSC
            CALL Write_TypeOp(SplitOp%tab_TypeOp(iSC,iSO))
          END DO
        END DO
      END IF
    END IF

    write(out_unit,*) 'END write_SplitOp'

  END SUBROUTINE write_SplitOp
  SUBROUTINE read_SplitOp(SplitOp,nb_SplitOp,nb_SplitCoord)
    USE EVR_system_m
    IMPLICIT NONE

    TYPE (SplitOp_t), intent(inout)   :: SplitOp
    integer,          intent(in)      :: nb_SplitOp,nb_SplitCoord

    integer :: iSO,iSC

    integer :: type_Op
    logical :: direct_KEO,direct_ScalOp,cplx

    namelist /TypeOp/ type_Op,direct_KEO,direct_ScalOp,cplx

    write(out_unit,*) 'BEGINNING read_SplitOp'
    write(out_unit,*) 'nb_SplitOp',nb_SplitOp
    write(out_unit,*) 'nb_SplitCoord',nb_SplitCoord

    IF ( (nb_SplitOp < 1 .AND. nb_SplitCoord > 0) .OR. &
         (nb_SplitOp > 0 .AND. nb_SplitCoord < 1) ) THEN
      write(out_unit,*) 'ERROR in read_SplitOp'
      write(out_unit,*) ' nb_SplitOp and nb_SplitCoord values are incompatible'
      write(out_unit,*) 'nb_SplitOp',nb_SplitOp
      write(out_unit,*) 'nb_SplitCoord',nb_SplitCoord
      STOP 'ERROR in read_SplitOp: nb_SplitOp and nb_SplitCoord values are incompatible'
    END IF
    SplitOp%nb_SplitOp    = nb_SplitOp
    SplitOp%nb_SplitCoord = nb_SplitCoord

    IF (SplitOp%nb_SplitOp == 0) THEN
      type_Op       = -1
      direct_KEO    = .FALSE.
      direct_ScalOp = .FALSE.
      cplx          = .FALSE.
      read(in_unit,TypeOp)
      SplitOp%param_TypeOp = param_TypeOp(type_Op=type_Op,direct_KEO=direct_KEO, &
                                          direct_ScalOp=direct_ScalOp,cplx=cplx)
    ELSE
      IF (allocated(SplitOp%tab_TypeOp)) THEN
        DO iSO=1,SplitOp%nb_SplitOp
          DO iSC=1,SplitOp%nb_SplitCoord
            type_Op       = -1
            direct_KEO    = .FALSE.
            direct_ScalOp = .FALSE.
            cplx          = .FALSE.
            read(in_unit,TypeOp)
            SplitOp%tab_TypeOp(iSC,iSO) = param_TypeOp(type_Op=type_Op,direct_KEO=direct_KEO, &
                                                direct_ScalOp=direct_ScalOp,cplx=cplx)
          END DO
        END DO
      END IF
    END IF
    write(out_unit,*) 'END read_SplitOp'

  END SUBROUTINE read_SplitOp
  SUBROUTINE dealloc_SplitOp(SplitOp)
    USE EVR_system_m
    IMPLICIT NONE

    TYPE (SplitOp_t), intent(inout)   :: SplitOp
    integer :: iSO,iSC

    CALL dealloc_TypeOp(SplitOp%param_TypeOp)

    SplitOp%nb_SplitOp    = 0
    SplitOp%nb_SplitCoord = 0
 
    IF (allocated(SplitOp%tab_TypeOp)) THEN
        DO iSO=1,SplitOp%nb_SplitOp
          DO iSC=1,SplitOp%nb_SplitCoord
            CALL dealloc_TypeOp(SplitOp%tab_TypeOp(iSC,iSO))
          END DO
        END DO
        deallocate(SplitOp%tab_TypeOp)
    END IF

  END SUBROUTINE dealloc_SplitOp
  SUBROUTINE read_ReadOp(ReadOp,mole)
    USE EVR_system_m
    USE mod_nDindex
    USE mod_Constant,  only : REAL_WU,convRWU_TO_R_WITH_WorkingUnit
    USE mod_Coord_KEO, only : CoordType
    USE mod_PrimOp
    USE mod_CAP
    USE mod_HStep
    IMPLICIT NONE
    
    !----- variables for operators ---------------------------
    TYPE (ReadOp_t), intent(inout) :: ReadOp
    
        !----- for the CoordType -----------------------------------------
    TYPE (CoordType), intent(in) :: mole

    ! local variables
    logical       :: test
    logical       :: comput_S
    
    logical       :: lect,restart
    logical       :: read_Grid,restart_Grid
    logical       :: Save_FileGrid,Keep_FileGrid,Save_MemGrid,Seq_FileGrid
  
    integer       :: num_grid_i,num_grid_f
    integer       :: JJ,Type_HamilOp
    integer       :: nb_CAP,nb_FluxOp
  
    logical       :: pot_only,T_only,pack_Op
    logical       :: read_Op,make_MatOp,save_MatOp,restart_MatOp
    integer       :: Partial_MatOp_i,Partial_MatOp_f
    logical       :: direct_KEO,direct_ScalOp
    logical       :: Op_WithContracRVec

    logical       :: NewOP
    integer       :: nb_Op
  
    character (len=Name_len) :: name0
  
    integer           :: direct,Type_FileGrid
    real (kind=Rkind) :: tol_pack,tol_nopack
  
    logical           :: Op_Transfo
    TYPE (REAL_WU)    :: E0_Transfo
    
    !     - working variables -----
    integer       :: ib
    integer       :: i,j,k,i_term
    
    !-------variables for the file names -------------------------------------
    character (len=Line_len) :: name_HADA
    logical                  :: formatted_HADA
    character (len=Line_len) :: name_Grid
    logical                  :: formatted_Grid
    character (len=Line_len) :: name_Mat
    logical                  :: formatted_Mat
    !-------------------------------------------------------------------------
  
    NAMELIST /actives/comput_S,test,                                &
                    read_Grid,lect,restart_Grid,restart,            &
                    Save_FileGrid,Keep_FileGrid,Save_MemGrid,       &
                    Seq_FileGrid,Type_FileGrid,                     &
                    num_grid_i,num_grid_f,                          &
                    pot_only,T_only,read_Op,                        &
                    name_HADA,formatted_HADA,                       &
                    name_Grid,formatted_Grid,                       &
                    name_Mat,formatted_Mat,                         &
                    JJ,Type_HamilOp,direct_KEO,direct_ScalOp,       &
                    direct,make_MatOp,save_MatOp,restart_MatOp,     &
                    Partial_MatOp_i,Partial_MatOp_f,                &
                    pack_Op,tol_pack,tol_nopack,                    &
                    Op_Transfo,E0_Transfo,nb_CAP,nb_FluxOp,         &
                    Op_WithContracRVec,                             &
                    NewOP, nb_Op
  
    !----- for debuging --------------------------------------------------
     character (len=*), parameter ::  name_sub='read_ReadOp'
    logical,parameter :: debug=.FALSE.
    !logical,parameter :: debug=.TRUE.
    !-----------------------------------------------------------

    write(out_unit,*) ' OPERATOR PARAMETERS (ACTIVE)'
  
    IF (print_level > 0) write(out_unit,*) 'BEGINNING ',name_sub

    !------- read the active namelist ----------------------------
    make_MatOp           = .FALSE.
    save_MatOp           = .FALSE.
    restart_MatOp        = .FALSE.
    Partial_MatOp_i      = 0
    Partial_MatOp_f      = huge(1)
    formatted_Mat        = .FALSE.
    name_Mat             = 'MatOp'
    tol_pack             = ONETENTH**7
    tol_nopack           = NINE/TEN
    pack_Op              = .FALSE.
    read_Op              = .FALSE.

    NewOP                = .FALSE.
    nb_Op                = 0

    test                 = .TRUE.
    comput_S             = .FALSE.
    Type_HamilOp         = 1
    direct_KEO           = .FALSE.
    direct_ScalOp        = .FALSE.
    Op_WithContracRVec  = .FALSE.
    nb_CAP              = 0
    nb_FluxOp           = 0

    direct              = 0
    Type_FileGrid       = 0
    Read_Grid           = .FALSE.
    Save_FileGrid       = .FALSE.
    Save_MemGrid        = .FALSE.
    Keep_FileGrid       = .FALSE.
    Seq_FileGrid        = .FALSE.
    lect                = .FALSE.
    Restart_Grid        = .FALSE.
    restart             = .FALSE.
    num_grid_i          = 0
    num_grid_f          = 0

    JJ                  = -1

    name_HADA           = 'SH_HADA'
    formatted_HADA      = .TRUE.
    name_Grid           = 'SH_HADA'
    formatted_Grid      = .TRUE.

    pot_only            = .FALSE.
    T_only              = .FALSE.

    Op_Transfo          = .FALSE.
    E0_Transfo          = REAL_WU(ZERO,   'cm-1','E') !ZERO with cm-1 as default unit

    read(in_unit,actives)
    IF (direct == 0 .OR. read_Op) make_MatOp = .TRUE.
    IF (print_level > 1) write(out_unit,actives)

    IF (name_HADA /= 'SH_HADA' .OR. .NOT. formatted_HADA) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) '  You should not use the parameters (name_HADA, formatted_HADA) '
      write(out_unit,*) '    in the "actives" namelist, '
      write(out_unit,*) '  name_Grid and formatted_Grid'
      STOP 'ERROR in read_ReadOp: You should not use the parameters (name_HADA, formatted_HADA)'
      name_Grid      = name_HADA
      formatted_Grid = formatted_HADA
    END IF

    IF (lect .OR. restart) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) '  You should not use the parameters (lect or restart) '
      write(out_unit,*) '    in the "actives" namelist, '
      write(out_unit,*) '  instead, you must use: Read_Grid and Restart_Grid'
      STOP 'ERROR in read_ReadOp: You should not use the parameters (lect or restart)'
    END IF

    ReadOp%Type_HamilOp       = Type_HamilOp
    ReadOp%direct_KEO         = direct_KEO
    ReadOp%direct_ScalOp      = direct_ScalOp
    ReadOp%Op_WithContracRVec = Op_WithContracRVec

    ReadOp%NewOP = NewOP
    ReadOp%nb_Op = nb_Op
    IF (NewOP .AND. nb_Op < 1) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) '  NewOp=t and  nb_Op < 1: it is not possible'
      write(out_unit,*) '  NewOp, nb_Op',NewOp, nb_Op
      write(out_unit,*) '  Either NewOp=t and  nb_Op >= 1'
      write(out_unit,*) '  Or     NewOp=f and  nb_Op <  1'
      STOP 'ERROR in read_ReadOp: NewOp=t and  nb_Op < 1: it is not possible'
    END IF
    IF (NewOP ) THEN
      allocate(ReadOp%tab_SplitOp(nb_Op))
      DO i=1,nb_Op
        CALL read_SplitOp(ReadOp%tab_SplitOp(i),nb_SplitOp=0,nb_SplitCoord=0)
        CALL write_SplitOp(ReadOp%tab_SplitOp(i))
      END DO
      STOP 'coucou'
    END IF

    IF (print_level > 0) write(out_unit,*) ' END of the "actives" namelist reading'
  
    !- Copy parameters in ReadOp --------
    SELECT CASE (direct)
    !         direct=0    => Make_Mat=T, SaveFile_Grid=T, SaveMem_Grid=F
    !         direct=1    => Make_Mat=F, SaveFile_Grid=T, SaveMem_Grid=T
    !         direct=2    => Make_Mat=F, SaveFile_Grid=F, SaveMem_Grid=T
    !         direct=3    => Make_Mat=F, SaveFile_Grid=T, SaveMem_Grid=F (for huge grid, like cHAC)
    !         direct=5    => Make_Mat=F, SaveFile_Grid=T, SaveMem_Grid=T (grid save at the end)
    CASE Default ! id case(0)
      ReadOp%Make_Mat = .TRUE.
      CALL init_FileGrid(ReadOp%para_FileGrid,                   &
                         Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                         Keep_FileGrid=.TRUE.,Save_MemGrid=.FALSE.,   &
                         Read_FileGrid=Read_Grid,                     &
                         Restart_Grid=Restart_Grid,test_Grid=test,    &
                         First_GridPoint=num_grid_i,                  &
                         Last_GridPoint=num_grid_f,                   &
                         Base_FileName_Grid=name_Grid)
    CASE (0)
      ReadOp%Make_Mat = .TRUE.
      CALL init_FileGrid(ReadOp%para_FileGrid,                   &
                         Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                         Formatted_FileGrid=formatted_Grid,           &
                         Keep_FileGrid=.TRUE.,Save_MemGrid=.FALSE.,   &
                         Read_FileGrid=Read_Grid,                     &
                         Restart_Grid=Restart_Grid,test_Grid=test,    &
                         First_GridPoint=num_grid_i,                  &
                         Last_GridPoint=num_grid_f,                   &
                         Base_FileName_Grid=name_Grid)
    CASE (1)
      ReadOp%Make_Mat = .FALSE.
      CALL init_FileGrid(ReadOp%para_FileGrid,                   &
                         Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                         Formatted_FileGrid=formatted_Grid,           &
                         Keep_FileGrid=.TRUE.,Save_MemGrid=.TRUE.,    &
                         Read_FileGrid=Read_Grid,                     &
                         Restart_Grid=Restart_Grid,test_Grid=test,    &
                         First_GridPoint=num_grid_i,                  &
                         Last_GridPoint=num_grid_f,                   &
                         Base_FileName_Grid=name_Grid)
    CASE (2)
      ReadOp%Make_Mat = .FALSE.
      CALL init_FileGrid(ReadOp%para_FileGrid,                   &
                         Type_FileGrid=0,Save_FileGrid =.FALSE.,      &
                         Formatted_FileGrid=formatted_Grid,           &
                         Keep_FileGrid=.FALSE.,Save_MemGrid=.TRUE.,   &
                         Read_FileGrid=Read_Grid,                     &
                         Restart_Grid=Restart_Grid,test_Grid=test,    &
                         First_GridPoint=num_grid_i,                  &
                         Last_GridPoint=num_grid_f,                   &
                         Base_FileName_Grid=name_Grid)
      IF (Type_FileGrid /= 0) THEN
        ReadOp%para_FileGrid%Save_MemGrid       = Save_MemGrid
        ReadOp%para_FileGrid%Type_FileGrid      = Type_FileGrid
        ReadOp%para_FileGrid%Keep_FileGrid      = Keep_FileGrid
        ReadOp%para_FileGrid%Formatted_FileGrid = .FALSE.
        ReadOp%para_FileGrid%Save_FileGrid      = .TRUE.
      END IF
    CASE (3)
      ReadOp%Make_Mat = .FALSE.
      CALL init_FileGrid(ReadOp%para_FileGrid,                   &
                         Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                         Formatted_FileGrid=formatted_Grid,           &
                         Keep_FileGrid=.TRUE.,Save_MemGrid=.FALSE.,   &
                         Read_FileGrid=Read_Grid,                     &
                         Restart_Grid=Restart_Grid,test_Grid=test,    &
                         First_GridPoint=num_grid_i,                  &
                         Last_GridPoint=num_grid_f,                   &
                         Base_FileName_Grid=name_Grid)
    CASE (4) ! for SG4, it enable to use directKEO and also direct pot
      ReadOp%Make_Mat = .FALSE.
      CALL init_FileGrid(ReadOp%para_FileGrid,                   &
                         Type_FileGrid=4,Save_FileGrid=Save_FileGrid, &
                         Formatted_FileGrid=.FALSE.,                  &
                         Keep_FileGrid=Keep_FileGrid,                 &
                         Save_MemGrid=Save_MemGrid,                   &
                         Read_FileGrid=Read_Grid,                     &
                         Restart_Grid=Restart_Grid,test_Grid=test,    &
                         First_GridPoint=num_grid_i,                  &
                         Last_GridPoint=num_grid_f,                   &
                         Base_FileName_Grid=name_Grid)
    CASE (5) ! The grid is save in memory, then it save on disk
      ReadOp%Make_Mat  = .FALSE.
      Save_FileGrid         = .NOT. Read_Grid
      CALL init_FileGrid(ReadOp%para_FileGrid,                   &
                         Type_FileGrid=5,Save_FileGrid=Save_FileGrid, &
                         Formatted_FileGrid=.FALSE.,                  &
                         Keep_FileGrid=.TRUE.,                        &
                         Save_MemGrid=.TRUE.,                         &
                         Read_FileGrid=Read_Grid,                     &
                         Restart_Grid=Restart_Grid,test_Grid=test,    &
                         First_GridPoint=num_grid_i,                  &
                         Last_GridPoint=num_grid_f,                   &
                         Base_FileName_Grid=name_Grid)
    END SELECT
  
    !CALL Write_FileGrid(ReadOp%para_FileGrid)
  
    IF (make_MatOp) THEN
      save_MatOp = save_MatOp .OR. restart_MatOp ! when restart=t, the matrix is always saved
      ReadOp%Make_Mat        = .TRUE.
      ReadOp%save_MatOp      = save_MatOp
      ReadOp%restart_MatOp   = restart_MatOp
      ReadOp%Partial_MatOp_i = Partial_MatOp_i
      ReadOp%Partial_MatOp_f = Partial_MatOp_f
      ReadOp%formatted_Mat   = formatted_Mat
      ReadOp%name_Mat        = name_Mat

      ReadOp%comput_S        = comput_S .AND. (direct < 2)

      ReadOp%pack_Op         = pack_Op .AND. ReadOp%Make_Mat
      ReadOp%read_Op         = read_Op
      ReadOp%tol_pack        = tol_pack
      ReadOp%tol_nopack      = tol_nopack
    ELSE
      ReadOp%Make_Mat        = .FALSE.
      ReadOp%save_MatOp      = .FALSE.
      ReadOp%restart_MatOp   = .FALSE.
      ReadOp%Partial_MatOp_i = 0
      ReadOp%Partial_MatOp_f = huge(1)
      ReadOp%formatted_Mat   = .TRUE.
      ReadOp%name_Mat        = name_Mat

      ReadOp%comput_S        = .FALSE.

      ReadOp%pack_Op         = .FALSE.
      ReadOp%read_Op         = .FALSE.
      ReadOp%tol_pack        = ONETENTH**7
      ReadOp%tol_nopack      = NINE/TEN
    END IF

    ReadOp%pot_only        = pot_only
    ReadOp%T_only          = T_only

    ReadOp%Op_Transfo      = Op_Transfo
    IF(MPI_id==0) write(out_unit,*) 'ReadOp%Op_Transfo',ReadOp%Op_Transfo
    IF (Op_Transfo) THEN
      !write(out_unit,*) 'E0_Transfo',E0_Transfo

      ReadOp%E0_Transfo      = convRWU_TO_R_WITH_WorkingUnit(E0_Transfo)
      ReadOp%degree_Transfo  = 2
      CALL alloc_NParray(ReadOp%Poly_Transfo,[2],            &
                        'ReadOp%Poly_Transfo','read_active',[0])
      ReadOp%Poly_Transfo = [ZERO,ZERO,ONE] ! x^2
      write(out_unit,*) 'ReadOp%E0_Transfo',ReadOp%E0_Transfo
      !write(out_unit,*) 'l u bounds',ubound(ReadOp%Poly_Transfo),lbound(ReadOp%Poly_Transfo)

    END IF

    ReadOp%nb_CAP  = nb_CAP

    IF (ReadOp%nb_CAP > 0) THEN
      allocate(ReadOp%tab_CAP(ReadOp%nb_CAP))
      DO i=1,size(ReadOp%tab_CAP)
        CALL Read_CAP(ReadOp%tab_CAP(i),mole%nb_Qtransfo)
      END DO
    END IF

    ReadOp%nb_FluxOp  = nb_FluxOp

    IF (ReadOp%nb_FluxOp > 0) THEN
      allocate(ReadOp%tab_HStep(ReadOp%nb_FluxOp))
      DO i=1,size(ReadOp%tab_HStep)
        CALL Read_HStep(ReadOp%tab_HStep(i))
      END DO
    END IF


    IF (JJ > -1) THEN
      ReadOp%nb_bRot         = 2*JJ+1
    ELSE
      ReadOp%nb_bRot         = 1
    END IF
    IF(MPI_id==0) write(out_unit,*) 'The number of rotational basis is:',ReadOp%nb_bRot

    IF(openmpi .AND. direct/=4) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'use direct=4 for namelist "actives" when running with MPI'
      STOP 'ERROR in read_ReadOp: use direct=4 for namelist "actives" when running with MPI'
    END IF
  
  END SUBROUTINE read_ReadOp


  SUBROUTINE init_ReadOp(ReadOp)
    TYPE (ReadOp_t), intent(inout) :: ReadOp

    integer :: i

    ReadOp%pack_Op             = .FALSE.
    ReadOp%tol_pack            = ONETENTH**7
    ReadOp%tol_nopack          = NINE*ONETENTH
    ReadOp%read_Op             = .FALSE.
    ReadOp%make_Mat            = .FALSE.
    ReadOp%save_MatOp          = .FALSE.
    ReadOp%restart_MatOp       = .FALSE.
    ReadOp%Partial_MatOp_i     = 0
    ReadOp%Partial_MatOp_f     = huge(1)
    ReadOp%formatted_Mat       = .TRUE.
    ReadOp%name_Mat            = 'MatOp'

    ReadOp%spectral            = .FALSE.
    ReadOp%Op_WithContracRVec  = .FALSE.

    ReadOp%spectral_Op         = 0
    ReadOp%nb_bRot             = 0

    ReadOp%pot_only            = .FALSE.
    ReadOp%T_only              = .FALSE.
    ReadOp%comput_S            = .FALSE.   ! calculation of the active overlap matrix

    CALL init_FileGrid(ReadOp%para_FileGrid)
    CALL file_dealloc(ReadOp%FileMat)

    ReadOp%Op_Transfo          = .FALSE.   ! true => we are using Transfo(Op) instead of Op
    ReadOp%E0_Transfo          = ZERO      ! ScaledOp = Op - E0_transfo * I
    ReadOp%degree_Transfo      = -1        ! degree of the transformation

    ReadOp%NewOp               = .FALSE.   ! true => read the list of operator parameters
    ReadOp%nb_Op               = -1
    ReadOp%nb_SplitOp          = 0
    ReadOp%nb_SplitCoord       = 0
    IF (allocated(ReadOp%tab_SplitOp)) THEN
      DO i=1,size(ReadOp%tab_SplitOp)
        CALL dealloc_SplitOp(ReadOp%tab_SplitOp(i))
      END DO
    END IF
  END SUBROUTINE init_ReadOp

  SUBROUTINE ReadOp2_TO_ReadOp1(ReadOp1,ReadOp2)
    CLASS (ReadOp_t), intent(inout) :: ReadOp1
    TYPE (ReadOp_t),  intent(in)    :: ReadOp2

    integer :: i

    ReadOp1%PrimOp_t           = ReadOp2%PrimOp_t

    ReadOp1%OpPsi_WithGrid     = ReadOp2%OpPsi_WithGrid

    ReadOp1%pack_Op            = ReadOp2%pack_Op
    ReadOp1%tol_pack           = ReadOp2%tol_pack
    ReadOp1%tol_nopack         = ReadOp2%tol_nopack
    ReadOp1%read_Op            = ReadOp2%read_Op
    ReadOp1%make_Mat           = ReadOp2%make_Mat
    ReadOp1%save_MatOp         = ReadOp2%save_MatOp
    ReadOp1%restart_MatOp      = ReadOp2%restart_MatOp
    ReadOp1%formatted_Mat      = ReadOp2%formatted_Mat
    ReadOp1%Partial_MatOp_i    = ReadOp2%Partial_MatOp_i
    ReadOp1%Partial_MatOp_f    = ReadOp2%Partial_MatOp_f
    ReadOp1%name_Mat           = ReadOp2%name_Mat
    ReadOp1%FileMat            = ReadOp2%FileMat

    ReadOp1%spectral           = ReadOp2%spectral
    ReadOp1%spectral_Op        = ReadOp2%spectral_Op
    ReadOp1%Op_WithContracRVec = ReadOp2%Op_WithContracRVec

    ReadOp1%nb_bRot            = ReadOp2%nb_bRot

    ReadOp1%pot_only           = ReadOp2%pot_only
    ReadOp1%T_only             = ReadOp2%T_only
    ReadOp1%comput_S           = ReadOp2%comput_S

    ReadOp1%para_FileGrid      = ReadOp2%para_FileGrid

    ReadOp1%Op_Transfo         = ReadOp2%Op_Transfo
    ReadOp1%E0_Transfo         = ReadOp2%E0_Transfo
    ReadOp1%degree_Transfo     = ReadOp2%degree_Transfo
    IF (allocated(ReadOp2%Poly_Transfo)) THEN
      ReadOp1%Poly_Transfo     = ReadOp2%Poly_Transfo
    END IF

    ReadOp1%NewOp              = ReadOp2%NewOp
    ReadOp1%nb_Op              = ReadOp2%nb_Op
    ReadOp1%nb_SplitOp         = ReadOp2%nb_SplitOp
    ReadOp1%nb_SplitCoord      = ReadOp2%nb_SplitCoord
    IF (allocated(ReadOp2%tab_SplitOp)) THEN
      ReadOp1%tab_SplitOp = ReadOp2%tab_SplitOp
    END IF

  END SUBROUTINE ReadOp2_TO_ReadOp1

  SUBROUTINE ReadOp2_TO_ReadOp1_FOR_AutoBasis(ReadOp1,ReadOp2,nq)
    TYPE (ReadOp_t),  intent(inout) :: ReadOp1
    TYPE (ReadOp_t),  intent(in)    :: ReadOp2
    integer,          intent(in)    :: nq

    ReadOp1%PrimOp_t           = ReadOp2%PrimOp_t
    ReadOp1%nb_scalar_Op       = 0
    ReadOp1%nb_CAP             = 0
    ReadOp1%nb_FluxOp          = 0
    ReadOp1%calc_scalar_Op     = .FALSE.
    ReadOp1%type_HamilOp       = 1
    ReadOp1%direct_KEO         = .FALSE.


    ReadOp1%OpPsi_WithGrid     = .FALSE.

    ReadOp1%pack_Op            = .FALSE.
    ReadOp1%read_Op            = .FALSE.
    ReadOp1%make_Mat           = .TRUE.
    ReadOp1%save_MatOp         = .FALSE.
    ReadOp1%restart_MatOp      = .FALSE.
    ReadOp1%formatted_Mat      = .TRUE.
    ReadOp1%name_Mat           = 'MatOp'
    ReadOp1%Partial_MatOp_i    = 0
    ReadOp1%Partial_MatOp_f    = huge(1)
    CALL file_dealloc(ReadOp1%FileMat)

    ReadOp1%spectral           = .FALSE.
    ReadOp1%Op_WithContracRVec = .FALSE.

    ReadOp1%nb_bRot            = 1

    ReadOp1%pot_only           = .FALSE.
    ReadOp1%T_only             = .FALSE.
    ReadOp1%comput_S           = .FALSE.

    ReadOp1%para_FileGrid      = ReadOp2%para_FileGrid
    ReadOp1%para_FileGrid%Save_FileGrid   = .FALSE.
    ReadOp1%para_FileGrid%First_GridPoint = 1
    ReadOp1%para_FileGrid%Last_GridPoint  = nq
    ReadOp1%para_FileGrid%Restart_Grid    = .FALSE.
    ReadOp1%para_FileGrid%Test_Grid       = .FALSE.
    ReadOp1%para_FileGrid%Read_FileGrid   = .FALSE.
    ReadOp1%para_FileGrid%Type_FileGrid   = 0

    ReadOp1%Op_Transfo         = .FALSE.
    ReadOp1%E0_Transfo         = ZERO
    ReadOp1%degree_Transfo     = -1

    ReadOp1%NewOp              = ReadOp2%NewOp ! ???????
    ReadOp1%nb_Op              = ReadOp2%nb_Op ! ??????? NO it should be only for H

    ReadOp1%nb_SplitOp          = ReadOp2%nb_SplitOp
    ReadOp1%nb_SplitCoord       = ReadOp2%nb_SplitCoord
    IF (allocated(ReadOp2%tab_SplitOp)) THEN
      ReadOp1%tab_SplitOp = ReadOp2%tab_SplitOp
    END IF

  END SUBROUTINE ReadOp2_TO_ReadOp1_FOR_AutoBasis

  SUBROUTINE dealloc_ReadOp(ReadOp)
    CLASS (ReadOp_t), intent(inout) :: ReadOp

    integer :: i

    CALL dealloc_PrimOp(ReadOp%PrimOp_t)

    ReadOp%OpPsi_WithGrid      = .FALSE.
    ReadOp%pack_Op             = .FALSE.
    ReadOp%tol_pack            = ONETENTH**7
    ReadOp%tol_nopack          = 0.9_Rkind
    ReadOp%read_Op             = .FALSE.
    ReadOp%make_Mat            = .FALSE.
    ReadOp%save_MatOp          = .FALSE.
    ReadOp%restart_MatOp       = .FALSE.
    ReadOp%formatted_Mat       = .TRUE.
    ReadOp%Partial_MatOp_i     = 0
    ReadOp%Partial_MatOp_f     = huge(1)
    ReadOp%name_Mat            = 'MatOp'
    CALL file_dealloc(ReadOp%FileMat)
    ReadOp%spectral            = .FALSE.
    ReadOp%spectral_Op         = 0
    ReadOp%Op_WithContracRVec  = .FALSE.

    ReadOp%nb_bRot             = 0

    ReadOp%pot_only            = .FALSE.
    ReadOp%T_only              = .FALSE.
    ReadOp%comput_S            = .FALSE.

    CALL dealloc_FileGrid(ReadOp%para_FileGrid)

    ReadOp%Op_Transfo          = .FALSE.
    ReadOp%E0_Transfo          = ZERO
    ReadOp%degree_Transfo      = -1
    IF (allocated(ReadOp%Poly_Transfo)) deallocate(ReadOp%Poly_Transfo)

    ReadOp%NewOp               = .FALSE.   ! true => read the list of operator parameters
    ReadOp%nb_Op               = -1
    ReadOp%nb_SplitOp          = 0
    ReadOp%nb_SplitCoord       = 0
    IF (allocated(ReadOp%tab_SplitOp)) THEN
      DO i=1,size(ReadOp%tab_SplitOp)
        CALL dealloc_SplitOp(ReadOp%tab_SplitOp(i))
      END DO
    END IF

  END SUBROUTINE dealloc_ReadOp

END MODULE ReadOp_m
