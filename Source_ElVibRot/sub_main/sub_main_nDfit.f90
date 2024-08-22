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
MODULE mod_nDGridFit
  USE EVR_system_m
  USE mod_nDindex,     only : Type_nDindex
  use mod_dnSVM,       only : Type_dnVec
  USE QDUtil_IntVec_m, ONLY : IntVec_t 
  IMPLICIT NONE

  TYPE param_nDGrid ! it mays change in the futur (more like "basis" type)

    TYPE (Type_nDindex)             :: nDindG  ! enable to use multidimensional index for the basis
    real (kind=Rkind), allocatable  :: Q0(:)

    ! for Grid_v1
    real (kind=Rkind),  allocatable  :: stepQ(:)
    !for Grid_v2
    TYPE (Type_dnVec),  allocatable :: Tab_stepQ(:)
    TYPE (IntVec_t),    allocatable :: tab_i_TO_l(:)


    TYPE (File_t)               :: Grid_FOR_Fit_file1
    TYPE (File_t)               :: Grid_FOR_Fit_file2
    character (len=Line_len)        :: name_Grid = ''

    integer                         :: MinCoupling = 0
    integer                         :: MaxCoupling = 4
    integer                         :: nb_DeltaQ   = 20
    integer                         :: nb_G        = 0

    real (kind=Rkind), allocatable  :: Val(:)

    real (kind=Rkind)               :: MinNorm = ZERO
    real (kind=Rkind)               :: MaxNorm = TWO

  END TYPE param_nDGrid

  PRIVATE
  PUBLIC  :: sub_nDGrid_nDfit

CONTAINS

SUBROUTINE dealloc_nDGrid(para_nDGrid)
      USE mod_nDindex, only : dealloc_nDindex
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'dealloc_nDGrid'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      CALL dealloc_nDindex(para_nDGrid%nDindG)
      IF (allocated(para_nDGrid%Q0))                                   &
                CALL dealloc_NParray(para_nDGrid%Q0,'Q0',name_sub)
      IF (allocated(para_nDGrid%stepQ))                                &
                CALL dealloc_NParray(para_nDGrid%stepQ,'stepQ',name_sub)

      para_nDGrid%MinNorm      = ZERO
      para_nDGrid%MaxNorm      = TWO
      para_nDGrid%MaxCoupling  = 4
      para_nDGrid%MinCoupling  = 0
      para_nDGrid%nb_DeltaQ    = 20
      para_nDGrid%nb_G         = 0

      CALL file_close(para_nDGrid%Grid_FOR_Fit_file1)
      CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

END SUBROUTINE dealloc_nDGrid


SUBROUTINE sub_nDGrid_nDfit()
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================
!

!----- physical and mathematical constants ---------------------------
      TYPE (constant)             :: const_phys

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)            :: mole
      TYPE (Tnum)                 :: para_Tnum
      real (kind=Rkind), pointer  :: hCC(:,:)
      TYPE(Type_dnMat)            :: dnGG

!----- variables for the construction of H ----------------------------
      TYPE (PrimOp_t)             :: PrimOp


!----- variables divers ----------------------------------------------
      logical                    :: Grid,Fit,Grid1_TO_Grid2,Transfo_fit,Analysis
      integer                    :: Grid_type,Fit_type,nb_Grid,iG,FirstGrid,LastGrid
      TYPE (param_nDGrid),allocatable  :: para_nDGrid(:)
      TYPE (param_nDFit)         :: para_nDFit
      real (kind=Rkind)          :: val_nDfit
      real (kind=Rkind), pointer :: Q(:)
      character (len=Line_len)   :: name_Grid,name_Fit
      real (kind=Rkind), allocatable :: Qact0(:)
      real (kind=Rkind)          :: auTOenergy


      namelist /Grid_Fit/ Grid,Grid_type,nb_Grid,name_Grid,Grid1_TO_Grid2,      &
                          Fit,Fit_type,FirstGrid,LastGrid,                      &
                          Transfo_fit,name_Fit,Analysis

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDGrid_nDfit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!para_mem%mem_debug = .TRUE.
!para_mem%mem_print = .FALSE.

!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================


!---------------------------------------------------------------------
!------ read or set up the physical constants ------------------------
      CALL sub_constantes(const_phys,.TRUE.)

      auTOenergy = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)
      flush(out_unit)
      write(out_unit,*) "============================================================"
      write(out_unit,*) "============================================================"
      write(out_unit,*) "=== COORDINATES (TNUM) ====================================="
      write(out_unit,*) "============================================================"
      write(out_unit,*) "============================================================"
      flush(out_unit)
!=======================================================================
!=======================================================================
      CALL versionEVRT(.TRUE.)

      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   CoordType, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_CoordType(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)
      !-----------------------------------------------------------------

      write(out_unit,*) "============================================================"
      write(out_unit,*) "============================================================"
      write(out_unit,*) "=== END COORDINATES (TNUM) ================================="
      write(out_unit,*) "============================================================"
      write(out_unit,*) "============================================================"
      flush(out_unit)


      CALL alloc_NParray(Qact0,[mole%nb_var],'Qact0',name_sub)
      CALL get_Qact0(Qact0,mole%ActiveTransfo)

!=====================================================================
        name_Grid = "Grid_FOR_Fit"
        name_Fit  = "Param_FOR_Fit-col"

        Grid           = .FALSE.
        Grid_type      = 1
        nb_Grid        = 0
        Fit            = .FALSE.
        Fit_type       = 1
        FirstGrid      = 0 ! for Fit_type=2
        LastGrid       = 0 ! for Fit_type=2
        Grid1_TO_Grid2 = .FALSE.
        Transfo_fit    = .FALSE.
        Analysis       = .FALSE.
        read(in_unit,Grid_Fit,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' End-of-file or End-of-record'
          write(out_unit,*) ' The namelist "Grid_Fit" is probably absent'
          write(out_unit,*) ' check your data!'
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Some parameter name of the namelist "Grid_Fit" are probaly wrong'
          write(out_unit,*) ' check your data!'
          write(out_unit,Grid_Fit)
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unit,Grid_Fit)

        para_nDFit%Analysis   = Analysis
        para_nDFit%name_Fit   = name_Fit

        IF (nb_Grid < 1) nb_Grid = 1
        allocate(para_nDGrid(nb_Grid))

        IF (Fit_type == 2 .AND. FirstGrid == 0) FirstGrid = min(2,nb_Grid)
        IF (Fit_type == 2 .AND. LastGrid  == 0) LastGrid  = nb_Grid

        DO iG=1,nb_Grid
          para_nDGrid(iG)%name_Grid = name_Grid

          CALL alloc_NParray(para_nDGrid(iG)%Q0,[mole%nb_act],'Q0',name_sub)
          para_nDGrid(iG)%Q0(:) = Qact0(1:mole%nb_act)
        END DO

        IF (Grid1_TO_Grid2) THEN
          CALL sub_nGrid1_TO_nGrid2(para_nDGrid(1)%Q0)
          STOP
        END IF

        IF (Grid) THEN
          SELECT CASE (Grid_type)
          CASE (1)
            CALL sub_nDGrid_v1(para_nDGrid(1),Qact0,mole,PrimOp)
          CASE (2)
            write(out_unit,*) '--------------------------------------------------'
            write(out_unit,*) '-------No coupling -------------------------------'
            CALL sub_nDGrid_coupling0_v2(para_nDGrid(1),Qact0,para_Tnum,mole,PrimOp)
            write(out_unit,*) '--------------------------------------------------'
            DO iG=2,nb_Grid
              write(out_unit,*) '-------Grid:',iG,' -----------------------------'
              CALL sub_nDGrid_coupling1_v2(para_nDGrid(iG),Qact0,para_Tnum,mole,PrimOp)
              write(out_unit,*) '--------------------------------------------------'
            END DO
          CASE DEFAULT
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' Wrong "Grid_type value"',Grid_type
            write(out_unit,*) ' check your data!'
            write(out_unit,Grid_Fit)
            write(out_unit,*) ' ERROR in ',name_sub
            STOP
          END SELECT
        ELSE
          CALL sub_nDGrid_WiTHOUT_calc(para_nDGrid(1),Qact0,mole,PrimOp)
        END IF

        IF (Fit) THEN
          SELECT CASE (Fit_type)
          CASE (1)
            CALL sub_nDFit_v1(para_nDFit,para_nDGrid(1),Qact0,mole,PrimOp)
          CASE (2)
            CALL sub_TransfoGrid_v2(para_nDGrid)
            CALL sub_nDFit_v2(para_nDFit,para_nDGrid(FirstGrid:LastGrid),Qact0,mole,PrimOp)
          CASE DEFAULT
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' Wrong "Fit_type value"',Fit_type
            write(out_unit,*) ' check your data!'
            write(out_unit,Grid_Fit)
            write(out_unit,*) ' ERROR in ',name_sub
            STOP
          END SELECT
        END IF

        IF (Transfo_fit) THEN
          CALL Read_FOR_nDFit1_TO_TnDFit2()
        END IF

        IF (Analysis) THEN
          CALL Read_Analysis(para_nDFit%nDFitAna,para_nDGrid(1)%Q0)
          CALL Analysis_nDFitW(para_nDFit,auTOenergy)
          !CALL Analysis_nDFit(para_nDFit,auTOenergy)
        END IF

        IF (PrimOp%nDfit_Op) THEN
          PrimOp%para_nDFit_V%Param_Fit_file%name = para_nDFit%Param_Fit_file%name

          nullify(Q)
          CALL alloc_array(Q,[ mole%nb_act ],'Q',name_sub)

          Q(:) = para_nDGrid(1)%Q0(:)
          Q(1) = para_nDGrid(1)%Q0(1) + ZERO

          CALL sub_nDFunc_FROM_nDFit(val_nDfit,Q,PrimOp%para_nDFit_V)
          write(out_unit,*) 'val_nDfit',val_nDfit

          CALL dealloc_array(Q,'Q',name_sub)
        END IF

!=====================================================================
      CALL dealloc_table_at(const_phys%mendeleev)

      CALL dealloc_CoordType(mole)
      IF (associated(para_Tnum%Gref)) THEN
        CALL dealloc_array(para_Tnum%Gref,"para_Tnum%Gref",name_sub)
      END IF


      write(out_unit,*) 'mem_tot',para_mem%mem_tot

      END SUBROUTINE sub_nDGrid_nDfit
      SUBROUTINE sub_nDGrid(para_nDGrid,Qact,mole,PrimOp)
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Constant,  only : get_Conv_au_TO_unit
      USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (PrimOp_t)    :: PrimOp
!=====================================================================


!----- variables divers ----------------------------------------------
      integer :: nderivE,MinCoupling,MaxCoupling,nb_DeltaQ
      real (kind=Rkind) :: MinNorm,MaxNorm,conv_ene
      real (kind=Rkind), pointer :: nDweight(:)

      integer, pointer :: nDsize(:)
      integer, pointer :: nDinit(:)
      integer :: i,ic,iGP,isign,nb_couplings,iGPtot

      integer           :: iOp,nioGrid1,nioGrid2

      integer                  :: idum
      character (len=Name_len) :: name_dum
      logical                  :: nDsize_read,nb_Gonly,add_Grid

      TYPE (param_dnMatOp) :: dnMatOp(1)


      TYPE (Type_nDindex),pointer :: Tab_PointGridSign(:) ! enable to change the sign of the DelatQ

      namelist /nDGrid/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,        &
                        nDsize_read,nb_DeltaQ,nb_Gonly,add_Grid

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDGrid'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      CALL alloc_NParray(para_nDGrid%Q0,[mole%nb_act],'para_nDGrid%Q0',name_sub)
      para_nDGrid%Q0(:) = Qact(1:mole%nb_act)
      CALL alloc_NParray(para_nDGrid%stepQ,[mole%nb_act],'para_nDGrid%stepQ',name_sub)
      para_nDGrid%stepQ(:) = ZERO

      nullify(nDsize)
      CALL alloc_array(nDsize,[mole%nb_act],'nDsize',name_sub)

      nullify(nDinit)
      CALL alloc_array(nDinit,[mole%nb_act],'nDinit',name_sub)
      nDinit(:) = 1
      nullify(nDweight)
      CALL alloc_array(nDweight,[mole%nb_act],'nDweight',name_sub)

      CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,PrimOp%nb_elec,    &
                               nderiv=0,cplx=.FALSE.,JRot=0)

!=====================================================================
      conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)


      para_nDGrid%Grid_FOR_Fit_file1%name = adjustl(trim(para_nDGrid%name_Grid)) // '1'
      para_nDGrid%Grid_FOR_Fit_file2%name = adjustl(trim(para_nDGrid%name_Grid)) // '2'

      CALL file_open(para_nDGrid%Grid_FOR_FIT_file1,nioGrid1)
      CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)

      add_Grid = .TRUE.
      DO ; IF (.NOT. add_Grid) EXIT  ! for the end of the loop

        MinNorm     = ZERO
        MaxNorm     = TWO
        MaxCoupling = 4
        MinCoupling = 0
        nb_DeltaQ   = 20
        nDsize_read = .FALSE.
        nb_Gonly    = .FALSE.
        add_Grid    = .FALSE.
        read(in_unit,nDGrid,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' End-of-file or End-of-record'
          write(out_unit,*) ' The namelist "nDGrid" is probably absent'
          write(out_unit,*) ' check your data!'
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Some parameter name of the namelist "nDGrid" are probaly wrong'
          write(out_unit,*) ' check your data!'
          write(out_unit,nDGrid)
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unit,nDGrid)
        read(in_unit,*) nDweight(:)
        read(in_unit,*) para_nDGrid%stepQ(:)
        IF (nDsize_read) read(in_unit,*) nDsize(:)

        IF (debug) write(out_unit,*) 'nDweight',nDweight(:)
        IF (debug) write(out_unit,*) 'stepQ',para_nDGrid%stepQ(:)
        flush(out_unit)


        para_nDGrid%MinNorm      = MinNorm
        para_nDGrid%MaxNorm      = MaxNorm
        para_nDGrid%MaxCoupling  = MaxCoupling
        para_nDGrid%MinCoupling  = MinCoupling
        para_nDGrid%nb_DeltaQ    = nb_DeltaQ
        IF (.NOT. nDsize_read) nDsize(:) = nb_DeltaQ
        IF (debug) write(out_unit,*) 'nDsize',nDsize(:)
        flush(out_unit)

        CALL init_nDindexPrim(para_nDGrid%nDindG,mole%nb_act,nDsize=nDsize,    &
                              nDweight=nDweight,type_OF_nDindex=0,      &
                              nDinit=nDinit,                            &
                              MaxNorm=MaxNorm,MinNorm=MinNorm,          &
                              MaxCoupling=MaxCoupling,MinCoupling=MinCoupling)
        CALL sort_nDindex(para_nDGrid%nDindG)
        para_nDGrid%nDindG%Tab_nDval(:,:) = para_nDGrid%nDindG%Tab_nDval(:,:) - 1
        CALL Write_nDindex(para_nDGrid%nDindG)

        nullify(Tab_PointGridSign)
        CALL alloc_array(Tab_PointGridSign,[MaxCoupling],             &
                        'Tab_PointGridSign',name_sub)
        DO ic=1,MaxCoupling
          CALL init_nDindexPrim(Tab_PointGridSign(ic),ic,[ (2,i=1,ic) ], &
                                type_OF_nDindex=1)
          Tab_PointGridSign(ic)%packed = .TRUE.
          CALL pack_nDindex(Tab_PointGridSign(ic))
          WHERE (Tab_PointGridSign(ic)%Tab_nDval(:,:) == 1)
            Tab_PointGridSign(ic)%Tab_nDval(:,:) = -1
          ELSEWHERE
            Tab_PointGridSign(ic)%Tab_nDval(:,:) = 1
          END WHERE
          !CALL Write_nDindex(Tab_PointGridSign(ic))
        END DO

        !number of points
        iGPtot = 0
        DO iGP=1,para_nDGrid%nDindG%Max_nDI
          nb_couplings=count(para_nDGrid%nDindG%Tab_nDval(:,iGP) > 0)
          IF (nb_couplings < MinCoupling) CYCLE

          IF (nb_couplings == 0) THEN
            iGPtot = iGPtot + 1
          ELSE
            iGPtot = iGPtot + Tab_PointGridSign(nb_couplings)%Max_nDI
          END IF
        END DO
        write(out_unit,*) "nb_G",iGPtot
        write(out_unit,*) "======================================"

        !IF (nb_Gonly) CYCLE

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        iGPtot = 0
        CALL get_Qact0(Qact,mole%ActiveTransfo)
        DO iGP=1,para_nDGrid%nDindG%Max_nDI
          nb_couplings=count(para_nDGrid%nDindG%Tab_nDval(:,iGP) > 0)
          IF (nb_couplings < MinCoupling) CYCLE

          IF (nb_couplings == 0) THEN
            Qact(1:mole%nb_act) = para_nDGrid%Q0(:)

            CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
            IF (.NOT. nb_Gonly) CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)

            iGPtot = iGPtot +1
            isign = 0
            nDinit(:) = para_nDGrid%nDindG%Tab_nDval(:,iGP)

            IF (debug) THEN
              write(out_unit,"(a,3i7,a,10i3)") 'iGPtot,iGp,isign,Qact',&
                                       iGPtot,iGp,isign,' : ',nDinit(:)
              write(out_unit,*) 'Grid',iGPtot,Qact(1:mole%nb_act),     &
                             dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            END IF

            write(nioGrid1,*) 'Grid1',iGPtot,nDinit(:),                 &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            write(nioGrid2,*) 'Grid2',iGPtot,Qact(1:mole%nb_act),       &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
          ELSE

            DO isign=1,Tab_PointGridSign(nb_couplings)%Max_nDI
              iGPtot = iGPtot +1
              ic = 0
              nDinit(:) = para_nDGrid%nDindG%Tab_nDval(:,iGP)
              DO i=1,mole%nb_act
                IF (para_nDGrid%nDindG%Tab_nDval(i,iGP) > 0) THEN
                  ic = ic +1
                  nDinit(i) = para_nDGrid%nDindG%Tab_nDval(i,iGP) *     &
                     Tab_PointGridSign(nb_couplings)%Tab_nDval(ic,isign)
                END IF
              END DO
              Qact(1:mole%nb_act) = para_nDGrid%Q0(:) +                 &
                       real(nDinit(:),kind=Rkind) * para_nDGrid%stepQ(:)

              CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
              IF (.NOT. nb_Gonly) CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)

              IF (debug) THEN
                write(out_unit,"(a,3i7,a,10i3)") 'iGPtot,iGp,isign,Qact',&
                                       iGPtot,iGp,isign,' : ',nDinit(:)
                write(out_unit,*) 'Grid',iGPtot,Qact(1:mole%nb_act),   &
                             dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
              END IF

              write(nioGrid1,*) 'Grid1',iGPtot,nDinit(:),               &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
              write(nioGrid2,*) 'Grid2',iGPtot,Qact(1:mole%nb_act),     &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            END DO
          END IF
        END DO
        para_nDGrid%nb_G = iGPtot
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"

      END DO

      CALL file_close(para_nDGrid%Grid_FOR_Fit_file1)
      CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

      IF (nb_Gonly) STOP 'nb_Gonly=t'


      CALL dealloc_array(nDsize,'nDsize',name_sub)
      CALL dealloc_array(nDinit,'nDinit',name_sub)
      DO ic=1,MaxCoupling
        CALL dealloc_nDindex(Tab_PointGridSign(ic))
      END DO
      CALL dealloc_array(Tab_PointGridSign,'Tab_PointGridSign',name_sub)

      CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
END SUBROUTINE sub_nDGrid

SUBROUTINE sub_nDGrid_coupling0_v2(para_nDGrid,Qact0,para_Tnum,mole,PrimOp)
  USE EVR_system_m
  USE mod_dnSVM
  USE mod_nDindex
  USE mod_Constant,  only : get_Conv_au_TO_unit
  USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0
  USE mod_PrimOp
  IMPLICIT NONE

  TYPE (param_nDGrid),  intent(inout) :: para_nDGrid

  !=====================================================================
  !----- for the CoordType and Tnum --------------------------------------
  TYPE (CoordType),     intent(in)    :: mole
  TYPE (Tnum),          intent(in)    :: para_Tnum
  real (kind=Rkind),    intent(in)    :: Qact0(:)

  !----- variables for the construction of H ---------------------------
  TYPE (PrimOp_t),      intent(in)    :: PrimOp
  !=====================================================================


  !----- variables divers ----------------------------------------------
  real (kind=Rkind)                 :: ValRef,DeltaVal,DeltaQ,kini,conv_ene
  integer                           :: i,j,iGP
  integer                           :: iOp,nioGrid1,nioGrid2

  TYPE (param_dnMatOp)              :: dnMatOp(1)
  real (kind=Rkind), allocatable    :: Qact(:)
  real (kind=Rkind), allocatable    :: StepQp(:)
  real (kind=Rkind), allocatable    :: StepQm(:)
  real (kind=Rkind), allocatable    :: HessRef(:,:)
  integer,           allocatable    :: nDinit(:),nDsize(:)

  integer                  :: iGdum
  character (len=Name_len) :: name_dum


  integer           :: MinCoupling,MaxCoupling,nb_DeltaQ
  real (kind=Rkind) :: MinNorm,MaxNorm
  logical           :: nb_Gonly,Read_Grid

  namelist /nDGrid/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,                    &
                        nb_DeltaQ,nb_Gonly,Read_Grid

  !----- for debuging --------------------------------------------------
  integer :: err_read
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = 'sub_nDGrid_coupling0_v2'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------


  CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,PrimOp%nb_elec,                  &
                               nderiv=2,cplx=.FALSE.,JRot=0)

  conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)


  para_nDGrid%Grid_FOR_Fit_file1%name = adjustl(trim(para_nDGrid%name_Grid)) // 'coup0_1'
  para_nDGrid%Grid_FOR_Fit_file2%name = adjustl(trim(para_nDGrid%name_Grid)) // 'coup0_2'
  CALL file_open(para_nDGrid%Grid_FOR_FIT_file1,nioGrid1)
  CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)


  MinNorm     = ZERO
  MaxNorm     = TWO
  MaxCoupling = 0
  MinCoupling = 0
  nb_DeltaQ   = 20
  nb_Gonly    = .FALSE.
  Read_Grid   = .FALSE.
  read(in_unit,nDGrid,IOSTAT=err_read)
  IF (err_read < 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' End-of-file or End-of-record'
    write(out_unit,*) ' The namelist "nDGrid" is probably absent'
    write(out_unit,*) ' check your data!'
    write(out_unit,*) ' ERROR in ',name_sub
    STOP
  ELSE IF (err_read > 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' Some parameter name of the namelist "nDGrid" are probaly wrong'
    write(out_unit,*) ' check your data!'
    write(out_unit,nDGrid)
    write(out_unit,*) ' ERROR in ',name_sub
    STOP
  END IF
  IF (debug) write(out_unit,nDGrid)
  flush(out_unit)

  para_nDGrid%MinNorm      = MinNorm
  para_nDGrid%MaxNorm      = MaxNorm
  para_nDGrid%MaxCoupling  = MaxCoupling
  para_nDGrid%MinCoupling  = MinCoupling

  write(out_unit,*) "======================================"
  write(out_unit,*) "======================================"
  CALL alloc_NParray(nDinit,[mole%nb_act],'nDinit',name_sub)
  nDinit(:) = 0

  CALL alloc_NParray(nDsize,[mole%nb_act],'nDsize',name_sub)
  nDsize(:) = 2*nb_DeltaQ+1

  CALL init_nDindexPrim(para_nDGrid%nDindG,ndim=mole%nb_act,                    &
                        nDsize=nDsize,nDinit=nDinit,                            &
                        type_OF_nDindex=5,Lmax=0)
  IF (debug) CALL Write_nDindex(para_nDGrid%nDindG)



  iGP = 1
  CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
  IF (nb_Gonly) STOP 'nb_Gonly=t'

  para_nDGrid%nb_G = 1
  allocate(para_nDGrid%Val(1))


  IF (.NOT. Read_Grid) THEN

    CALL get_dnMatOp_AT_Qact(Qact0,dnMatOp,mole,para_Tnum,PrimOp,nderiv=2)
    ValRef  = dnMatOp(1)%tab_dnMatOp(1,1,1)%d0
    HessRef = dnMatOp(1)%tab_dnMatOp(1,1,1)%d2

    IF (debug) THEN
      write(out_unit,"(a,i7,a,10i3)") 'iGP,iGp,isign,Qact',iGP,' : ',nDinit(:)
      write(out_unit,*) 'Grid',iGP,Qact0(1:mole%nb_act),                         &
                   dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene,                   &
        (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
    END IF
    write(nioGrid1,*) 'Grid1',iGP,nDinit(:),                                      &
                            dnMatOp(1)%tab_dnMatOp(:,:,1)%d0,                   &
      (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
    write(nioGrid2,*) 'Grid2',iGP,Qact0(1:mole%nb_act),                           &
                            dnMatOp(1)%tab_dnMatOp(:,:,1)%d0,                   &
      (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))


    CALL alloc_NParray(StepQp,[mole%nb_act],'StepQp',name_sub)
    CALL alloc_NParray(StepQm,[mole%nb_act],'StepQm',name_sub)

    DO i=1,mole%nb_act
      kini      = HessRef(i,i)
      DO
        DeltaQ = sqrt(MaxNorm/kini)
        write(out_unit,*) 'DeltaQ (one PES)',DeltaQ
        Qact    = Qact0
        Qact(i) = Qact0(i) + DeltaQ
        CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)
        DeltaVal = dnMatOp(1)%tab_dnMatOp(1,1,1)%d0-ValRef
        write(out_unit,*) 'i,DeltaQ,DeltaVal',i,DeltaQ,DeltaVal
        IF (DeltaVal > MaxNorm*1.1_Rkind .OR. DeltaVal < MaxNorm*0.9_Rkind) THEN
          kini = DeltaVal/DeltaQ**2
        ELSE
          StepQp(i) = DeltaQ
          EXIT
        END IF
      END DO
      kini      = HessRef(i,i)
      DO
        DeltaQ = -sqrt(MaxNorm/kini)
        write(out_unit,*) 'DeltaQ (one PES)',DeltaQ
        Qact    = Qact0
        Qact(i) = Qact0(i) + DeltaQ
        CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)
        DeltaVal = dnMatOp(1)%tab_dnMatOp(1,1,1)%d0-ValRef
        write(out_unit,*) 'i,DeltaQ,DeltaVal',i,DeltaQ,DeltaVal
        IF (DeltaVal > MaxNorm*1.1_Rkind .OR. DeltaVal < MaxNorm*0.9_Rkind) THEN
          kini = DeltaVal/DeltaQ**2
        ELSE
          StepQm(i) = DeltaQ
          EXIT
        END IF
      END DO
    END DO

    StepQm(:) = StepQm(:) / nb_DeltaQ
    StepQp(:) = StepQp(:) / nb_DeltaQ

    DO i=1,mole%nb_act
      write(out_unit,'(a,i0,100(1x,i0))')    'l(:)_',     i,0,([j,j],j=1,nb_DeltaQ)
      write(out_unit,'(a,i0,100(1x,f6.2))') 'DeltaQ(:)_',i,ZERO,([StepQp(i)*j,StepQm(i)*j],j=1,nb_DeltaQ)
    END DO
  ELSE
    read(nioGrid1,*) Name_dum,iGdum,nDinit(:),ValRef
    para_nDGrid%Val(1) = ValRef
    write(out_unit,"(a,i7,a,10i3)") 'iGP,nDval',iGP,' : ',nDinit(:)
    write(out_unit,*) 'Grid',iGP,Qact0(1:mole%nb_act),ValRef
  END IF

  CALL file_close(para_nDGrid%Grid_FOR_Fit_file1)
  CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

  CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

END SUBROUTINE sub_nDGrid_coupling0_v2

SUBROUTINE sub_nDGrid_coupling1_v2(para_nDGrid,Qact0,para_Tnum,mole,PrimOp)
  USE EVR_system_m
  use mod_dnSVM,       only : Type_dnVec, alloc_dnVec
  USE QDUtil_IntVec_m, ONLY : IntVec_t, alloc_IntVec
  USE mod_nDindex
  USE mod_Constant,    only : get_Conv_au_TO_unit
  USE mod_Coord_KEO,   only : CoordType, Tnum, get_Qact0
  USE mod_PrimOp
  IMPLICIT NONE

  TYPE (param_nDGrid),  intent(inout) :: para_nDGrid

  !=====================================================================
  !----- for the CoordType and Tnum --------------------------------------
  TYPE (CoordType),     intent(in)    :: mole
  TYPE (Tnum),          intent(in)    :: para_Tnum
  real (kind=Rkind),    intent(in)    :: Qact0(:)

  !----- variables for the construction of H ---------------------------
  TYPE (PrimOp_t),      intent(in)    :: PrimOp
  !=====================================================================


  !----- variables divers ----------------------------------------------
  !character (len=:), allocatable :: name_dum
  character (len=Line_len) :: name_dum


  real (kind=Rkind)                 :: ValRef,conv_ene
  integer                           :: i,j,iGP
  integer                           :: iOp,nioGrid1,nioGrid2

  TYPE (param_dnMatOp)              :: dnMatOp(1)
  real (kind=Rkind), allocatable    :: Qact(:)
  real (kind=Rkind), allocatable    :: StepQp(:)
  real (kind=Rkind), allocatable    :: StepQm(:)
  real (kind=Rkind), allocatable    :: DeltaQ(:)

  real (kind=Rkind), allocatable    :: HessRef(:,:)
  integer,           allocatable    :: nDinit(:),nDsize(:),nDval(:)
  integer,           allocatable    :: nDval_dum(:)
  integer                           :: iGdum,ii
  integer,           allocatable    :: tab_i_TO_l(:)
  real (kind=Rkind), allocatable    :: Tab_stepQ(:)
  real (kind=Rkind)                 :: stepQ


  integer           :: MinCoupling,MaxCoupling,nb_DeltaQ
  real (kind=Rkind) :: MinNorm,MaxNorm
  logical           :: nb_Gonly,Read_Grid

  namelist /nDGrid/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,                    &
                        nb_DeltaQ,nb_Gonly,Read_Grid

  !----- for debuging --------------------------------------------------
  integer :: err_read
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = 'sub_nDGrid_coupling1_v2'
  !logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------

  conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)

  MinNorm     = ZERO
  MaxNorm     = TWO
  MaxCoupling = 0
  MinCoupling = 0
  nb_DeltaQ   = 20
  nb_Gonly    = .FALSE.
  Read_Grid   = .FALSE.
  read(in_unit,nml=nDGrid,IOSTAT=err_read)
  IF (err_read < 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' End-of-file or End-of-record'
    write(out_unit,*) ' The namelist "nDGrid" is probably absent'
    write(out_unit,*) ' check your data!'
    write(out_unit,*) ' ERROR in ',name_sub
    STOP
  ELSE IF (err_read > 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' Some parameter name of the namelist "nDGrid" are probaly wrong'
    write(out_unit,*) ' check your data!'
    write(out_unit,nDGrid)
    write(out_unit,*) ' ERROR in ',name_sub
    STOP
  END IF
  IF (debug) write(out_unit,nml=nDGrid)
  IF (MaxCoupling == 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' This subroutine must be used with MaxCoupling > 0'
    write(out_unit,*) '    MaxCoupling',MaxCoupling
    write(out_unit,*) ' check your data!'
    write(out_unit,nml=nDGrid)
    write(out_unit,*) ' ERROR in ',name_sub
    STOP
  END IF

  flush(out_unit)

  para_nDGrid%MinNorm      = MinNorm
  para_nDGrid%MaxNorm      = MaxNorm
  para_nDGrid%MaxCoupling  = MaxCoupling
  para_nDGrid%MinCoupling  = MinCoupling
  para_nDGrid%nb_DeltaQ    = nb_DeltaQ

  allocate(para_nDGrid%Tab_stepQ(mole%nb_act))
  allocate(para_nDGrid%tab_i_TO_l(mole%nb_act))

  DO i=1,mole%nb_act
    ! for tab_i_TO_l
    IF (allocated(tab_i_TO_l)) deallocate(tab_i_TO_l)
    CALL read_name_advNo(in_unit,name_dum,err_read)
    DO
      CALL read_name_advNo(in_unit,name_dum,err_read)
      read(name_dum,*) ii
      IF (allocated(tab_i_TO_l)) THEN
        tab_i_TO_l = [tab_i_TO_l,ii]
      ELSE
        tab_i_TO_l = [ii]
      END IF
      IF (err_read /= 0) EXIT
    END DO
    CALL alloc_IntVec(para_nDGrid%tab_i_TO_l(i),SizeVec=size(tab_i_TO_l))
    para_nDGrid%tab_i_TO_l(i)%vec = tab_i_TO_l

    ! for Tab_stepQ
    IF (allocated(Tab_stepQ))  deallocate(Tab_stepQ)
    CALL read_name_advNo(in_unit,name_dum,err_read)
    DO
      CALL read_name_advNo(in_unit,name_dum,err_read)
      read(name_dum,*) stepQ
      IF (allocated(Tab_stepQ)) THEN
        Tab_stepQ = [Tab_stepQ,stepQ]
      ELSE
        Tab_stepQ = [stepQ]
      END IF
      IF (err_read /= 0) EXIT
    END DO
    CALL alloc_dnVec(para_nDGrid%Tab_stepQ(i),nb_var_vec=size(Tab_stepQ))
    para_nDGrid%Tab_stepQ(i)%d0 = Tab_stepQ

    write(out_unit,*) i,'tab_i_TO_l',para_nDGrid%tab_i_TO_l(i)%vec
    write(out_unit,*) i,'Tab_stepQ',para_nDGrid%Tab_stepQ(i)%d0

    !read(in_unit,*) name_dum,para_nDGrid%tab_i_TO_l(i)%vec(:)
    !read(in_unit,*) name_dum,para_nDGrid%Tab_stepQ(i)%d0(:)
  END DO



  write(out_unit,*) "======================================"
  write(out_unit,*) "======================================"
  CALL alloc_NParray(nDinit,[mole%nb_act],'nDinit',name_sub)
  nDinit(:) = 0
  CALL alloc_NParray(nDsize,[mole%nb_act],'nDsize',name_sub)
  !nDsize(:) = 2*nb_DeltaQ+1
  DO i=1,mole%nb_act
    nDsize(i) = size(para_nDGrid%tab_i_TO_l(i)%vec)
  END DO
  CALL alloc_NParray(nDVal,[mole%nb_act],'nDVal',name_sub)
  CALL alloc_NParray(DeltaQ,[mole%nb_act],'DeltaQ',name_sub)
  CALL alloc_NParray(nDVal_dum,[mole%nb_act],'nDVal_dum',name_sub)

  CALL init_nDindexPrim(para_nDGrid%nDindG,ndim=mole%nb_act,                    &
                        nDsize=nDsize,nDinit=nDinit,                            &
                        type_OF_nDindex=5,Lmax=para_nDGrid%nb_DeltaQ,           &
                        MinCoupling=para_nDGrid%MinCoupling,                    &
                        MaxCoupling=para_nDGrid%MaxCoupling,                    &
                        tab_i_TO_l=para_nDGrid%tab_i_TO_l)
  IF (debug) CALL Write_nDindex(para_nDGrid%nDindG)

  CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,PrimOp%nb_elec,                  &
                           nderiv=0,cplx=.FALSE.,JRot=0)

  para_nDGrid%nb_G = para_nDGrid%nDindG%Max_nDI
  allocate(para_nDGrid%Val(para_nDGrid%nDindG%Max_nDI))


  para_nDGrid%Grid_FOR_Fit_file1%name = adjustl(trim(para_nDGrid%name_Grid)) // &
        'coup' // TO_string(MaxCoupling) // '_1'
  para_nDGrid%Grid_FOR_Fit_file2%name = adjustl(trim(para_nDGrid%name_Grid)) // &
        'coup' // TO_string(MaxCoupling) // '_2'
  CALL file_open(para_nDGrid%Grid_FOR_FIT_file1,nioGrid1)
  CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)


  CALL init_nDval_OF_nDindex(para_nDGrid%nDindG,nDval)
  DO iGP=1,para_nDGrid%nDindG%Max_nDI
    CALL ADD_ONE_TO_nDindex(para_nDGrid%nDindG,nDval)

    DO i=1,mole%nb_act
      DeltaQ(i) = para_nDGrid%Tab_stepQ(i)%d0(nDval(i)+1)
    END DO
    Qact = Qact0
    Qact(1:mole%nb_act) = Qact(1:mole%nb_act) + DeltaQ

    IF (.NOT. Read_Grid) THEN
      IF (nb_Gonly) THEN
        ValRef = ZERO
      ELSE
        CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)
        ValRef  = dnMatOp(1)%tab_dnMatOp(1,1,1)%d0
      END IF

      IF (debug) THEN
        write(out_unit,"(a,i7,a,10i3)") 'iGP,nDval',iGP,' : ',nDval(:)
        write(out_unit,*) 'Grid',iGP,Qact(1:mole%nb_act),                        &
                       dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene,                 &
          (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
      END IF
      write(nioGrid1,*) 'Grid1',iGP,nDval(:),                                     &
                                dnMatOp(1)%tab_dnMatOp(:,:,1)%d0,                 &
          (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
      write(nioGrid2,*) 'Grid2',iGP,Qact(1:mole%nb_act),                          &
                                dnMatOp(1)%tab_dnMatOp(:,:,1)%d0,                 &
          (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))

      IF (debug) write(out_unit,*) 'iGP,nDval',iGP,nDval,DeltaQ,ValRef
    ELSE
      read(nioGrid1,*) Name_dum,iGdum,nDval_dum(:),ValRef
      para_nDGrid%Val(iGP) = ValRef
      write(out_unit,"(a,i7,a,10i3)") 'iGP,nDval',iGP,' : ',nDval(:)
      write(out_unit,*) 'Grid',iGP,Qact(1:mole%nb_act),ValRef
    END IF
  END DO


  CALL file_close(para_nDGrid%Grid_FOR_Fit_file1)
  CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

  CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

END SUBROUTINE sub_nDGrid_coupling1_v2

SUBROUTINE sub_TransfoGrid_v2(para_nDGrid)
  USE EVR_system_m
  USE mod_dnSVM
  USE mod_nDindex
  USE mod_Constant
  USE mod_Coord_KEO
  USE mod_PrimOp
  IMPLICIT NONE

  TYPE (param_nDGrid), intent(inout)   :: para_nDGrid(:)

  !----- variables divers ----------------------------------------------
  integer, allocatable            :: nDvalG(:)
  real (kind=Rkind)               :: valGP
  integer                         :: iGP,iG


  !----- for debuging --------------------------------------------------
  integer :: err_read
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = 'sub_TransfoGrid_v2'
  !      logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------

  CALL alloc_NParray(nDvalG,[para_nDGrid(1)%nDindG%ndim],'nDvalG',name_sub)

  write(out_unit,*) "======================================"
  write(out_unit,*) "=== TRANSFORM the Grids =============="
  write(out_unit,*) "======================================"
  DO iG=2,size(para_nDGrid)
  DO iGP=1,para_nDGrid(iG)%nb_G
      CALL calc_nDindex(para_nDGrid(iG)%nDindG,iGP,nDvalG)
      valGP = ValGridPoint(para_nDGrid(1:iG-1),nDvalG)
      write(out_unit,*) 'val',iG,iGP,':',nDvalG,para_nDGrid(iG)%Val(iGP)
      para_nDGrid(iG)%Val(iGP) = para_nDGrid(iG)%Val(iGP) - valGP
      write(out_unit,*) 'new val',iG,iGP,':',nDvalG,para_nDGrid(iG)%Val(iGP)
      write(out_unit,*)
      write(out_unit,*)
  END DO
  END DO

  CALL dealloc_NParray(nDvalG,'nDvalG',name_sub)

END SUBROUTINE sub_TransfoGrid_v2

FUNCTION ValGridPoint(para_nDGrid,nDvalG)
  USE EVR_system_m
  USE mod_nDindex
  IMPLICIT NONE

  real (kind=Rkind)    :: ValGridPoint

  TYPE (param_nDGrid),  intent(in) :: para_nDGrid(:)
  integer,              intent(in) :: nDvalG(:)



  !----- variables divers ----------------------------------------------
  real (kind=Rkind)                 :: Val
  integer                           :: i,j,iGP,iG
  integer                           :: nb_coupling

  integer,           allocatable    :: nDval(:)
  integer                           :: err_sub


  !----- for debuging --------------------------------------------------
  integer :: err_read
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = 'ValGridPoint'
  !logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------

  CALL alloc_NParray(nDVal,shape(nDvalG),'nDVal',name_sub)

  nb_coupling = count(nDvalG > 0)

  !zero coupling (allways)
  iG  = 1
  Val = para_nDGrid(iG)%Val(1)

  ! 1 mode coupling
  IF (nb_coupling > 1) THEN ! because when nb_coupling=1, the 1D-coupling cannot be removed
    iG    = iG + 1
    DO i=1,size(nDvalG)
      nDval    = 0
      IF (nDvalG(i) /= 0) THEN
        nDval(i) = nDvalG(i)
        CALL calc_nDI(iGP,nDval,para_nDGrid(iG)%nDindG,err_sub)
        Val = Val + para_nDGrid(iG)%Val(iGP)
      END IF
    END DO
  END IF

  ! 2 mode couplings
  IF (nb_coupling > 2) THEN
    iG    = iG + 1
    DO i=1,size(nDvalG)
    DO j=i+1,size(nDvalG)

      nDval    = 0

      IF (nDvalG(i) /= 0 .AND. nDvalG(j) /= 0) THEN
        nDval(i) = nDvalG(i)
        nDval(j) = nDvalG(j)
        CALL calc_nDI(iGP,nDval,para_nDGrid(iG)%nDindG,err_sub)
        Val = Val + para_nDGrid(iG)%Val(iGP)
      END IF
    END DO
    END DO

  END IF
  ValGridPoint = Val

  CALL dealloc_NParray(nDVal,'nDVal',name_sub)

END FUNCTION ValGridPoint

      SUBROUTINE sub_nDGrid_v1(para_nDGrid,Qact,mole,PrimOp)
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Constant,  only : get_Conv_au_TO_unit
      USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)      :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (PrimOp_t)    :: PrimOp
!=====================================================================


!----- variables divers ----------------------------------------------
      integer :: nderivE,MinCoupling,MaxCoupling,nb_DeltaQ
      real (kind=Rkind) :: MinNorm,MaxNorm,conv_ene
      real (kind=Rkind), pointer :: nDweight(:)

      integer, pointer :: nDsize(:)
      integer, pointer :: nDinit(:)
      integer :: i,ic,iGP,isign,nb_couplings,iGPtot

      integer           :: iOp,nioGrid1,nioGrid2

      integer                  :: idum
      character (len=Name_len) :: name_dum
      logical                  :: nDsize_read,nb_Gonly,add_Grid

      TYPE (param_dnMatOp) :: dnMatOp(1)


      TYPE (Type_nDindex),pointer :: Tab_PointGridSign(:) ! enable to change the sign of the DelatQ

      namelist /nDGrid/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,        &
                        nDsize_read,nb_DeltaQ,nb_Gonly,add_Grid

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDGrid'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      CALL alloc_NParray(para_nDGrid%Q0,[mole%nb_act],'para_nDGrid%Q0',name_sub)
      para_nDGrid%Q0(:) = Qact(1:mole%nb_act)
      CALL alloc_NParray(para_nDGrid%stepQ,[mole%nb_act],'para_nDGrid%stepQ',name_sub)
      para_nDGrid%stepQ(:) = ZERO

      nullify(nDsize)
      CALL alloc_array(nDsize,[mole%nb_act],'nDsize',name_sub)

      nullify(nDinit)
      CALL alloc_array(nDinit,[mole%nb_act],'nDinit',name_sub)
      nDinit(:) = 1
      nullify(nDweight)
      CALL alloc_array(nDweight,[mole%nb_act],'nDweight',name_sub)

      CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,PrimOp%nb_elec,    &
                               nderiv=0,cplx=.FALSE.,JRot=0)

!=====================================================================
      conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)


      para_nDGrid%Grid_FOR_Fit_file1%name = adjustl(trim(para_nDGrid%name_Grid)) // '1'
      para_nDGrid%Grid_FOR_Fit_file2%name = adjustl(trim(para_nDGrid%name_Grid)) // '2'

      CALL file_open(para_nDGrid%Grid_FOR_FIT_file1,nioGrid1)
      CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)

      add_Grid = .TRUE.
      DO ; IF (.NOT. add_Grid) EXIT  ! for the end of the loop

        MinNorm     = ZERO
        MaxNorm     = TWO
        MaxCoupling = 4
        MinCoupling = 0
        nb_DeltaQ   = 20
        nDsize_read = .FALSE.
        nb_Gonly    = .FALSE.
        add_Grid    = .FALSE.
        read(in_unit,nDGrid,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' End-of-file or End-of-record'
          write(out_unit,*) ' The namelist "nDGrid" is probably absent'
          write(out_unit,*) ' check your data!'
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Some parameter name of the namelist "nDGrid" are probaly wrong'
          write(out_unit,*) ' check your data!'
          write(out_unit,nDGrid)
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unit,nDGrid)
        read(in_unit,*) nDweight(:)
        read(in_unit,*) para_nDGrid%stepQ(:)
        IF (nDsize_read) read(in_unit,*) nDsize(:)

        IF (debug) write(out_unit,*) 'nDweight',nDweight(:)
        IF (debug) write(out_unit,*) 'stepQ',para_nDGrid%stepQ(:)
        flush(out_unit)


        para_nDGrid%MinNorm      = MinNorm
        para_nDGrid%MaxNorm      = MaxNorm
        para_nDGrid%MaxCoupling  = MaxCoupling
        para_nDGrid%MinCoupling  = MinCoupling
        para_nDGrid%nb_DeltaQ    = nb_DeltaQ
        IF (.NOT. nDsize_read) nDsize(:) = nb_DeltaQ
        IF (debug) write(out_unit,*) 'nDsize',nDsize(:)
        flush(out_unit)

        CALL init_nDindexPrim(para_nDGrid%nDindG,mole%nb_act,nDsize,    &
                              nDweight=nDweight,type_OF_nDindex=0,      &
                              nDinit=nDinit,                            &
                              MaxNorm=MaxNorm,MinNorm=MinNorm,          &
                              MaxCoupling=MaxCoupling,MinCoupling=MinCoupling)
        CALL sort_nDindex(para_nDGrid%nDindG)
        para_nDGrid%nDindG%Tab_nDval(:,:) = para_nDGrid%nDindG%Tab_nDval(:,:) - 1
        CALL Write_nDindex(para_nDGrid%nDindG)

        nullify(Tab_PointGridSign)
        CALL alloc_array(Tab_PointGridSign,[MaxCoupling],             &
                        'Tab_PointGridSign',name_sub)
        DO ic=1,MaxCoupling
          CALL init_nDindexPrim(Tab_PointGridSign(ic),ic,[ (2,i=1,ic) ], &
                                type_OF_nDindex=1)
          Tab_PointGridSign(ic)%packed = .TRUE.
          CALL pack_nDindex(Tab_PointGridSign(ic))
          WHERE (Tab_PointGridSign(ic)%Tab_nDval(:,:) == 1)
            Tab_PointGridSign(ic)%Tab_nDval(:,:) = -1
          ELSEWHERE
            Tab_PointGridSign(ic)%Tab_nDval(:,:) = 1
          END WHERE
          !CALL Write_nDindex(Tab_PointGridSign(ic))
        END DO

        !number of points
        iGPtot = 0
        DO iGP=1,para_nDGrid%nDindG%Max_nDI
          nb_couplings=count(para_nDGrid%nDindG%Tab_nDval(:,iGP) > 0)
          IF (nb_couplings < MinCoupling) CYCLE

          IF (nb_couplings == 0) THEN
            iGPtot = iGPtot + 1
          ELSE
            iGPtot = iGPtot + Tab_PointGridSign(nb_couplings)%Max_nDI
          END IF
        END DO
        write(out_unit,*) "nb_G",iGPtot
        write(out_unit,*) "======================================"

        !IF (nb_Gonly) CYCLE

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        iGPtot = 0
        CALL get_Qact0(Qact,mole%ActiveTransfo)
        DO iGP=1,para_nDGrid%nDindG%Max_nDI
          nb_couplings=count(para_nDGrid%nDindG%Tab_nDval(:,iGP) > 0)
          IF (nb_couplings < MinCoupling) CYCLE

          IF (nb_couplings == 0) THEN
            Qact(1:mole%nb_act) = para_nDGrid%Q0(:)

            CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
            IF (.NOT. nb_Gonly) CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)

            iGPtot = iGPtot +1
            isign = 0
            nDinit(:) = para_nDGrid%nDindG%Tab_nDval(:,iGP)

            IF (debug) THEN
              write(out_unit,"(a,3i7,a,10i3)") 'iGPtot,iGp,isign,Qact',&
                                       iGPtot,iGp,isign,' : ',nDinit(:)
              write(out_unit,*) 'Grid',iGPtot,Qact(1:mole%nb_act),     &
                             dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            END IF

            write(nioGrid1,*) 'Grid1',iGPtot,nDinit(:),                 &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            write(nioGrid2,*) 'Grid2',iGPtot,Qact(1:mole%nb_act),       &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
          ELSE

            DO isign=1,Tab_PointGridSign(nb_couplings)%Max_nDI
              iGPtot = iGPtot +1
              ic = 0
              nDinit(:) = para_nDGrid%nDindG%Tab_nDval(:,iGP)
              DO i=1,mole%nb_act
                IF (para_nDGrid%nDindG%Tab_nDval(i,iGP) > 0) THEN
                  ic = ic +1
                  nDinit(i) = para_nDGrid%nDindG%Tab_nDval(i,iGP) *     &
                     Tab_PointGridSign(nb_couplings)%Tab_nDval(ic,isign)
                END IF
              END DO
              Qact(1:mole%nb_act) = para_nDGrid%Q0(:) +                 &
                       real(nDinit(:),kind=Rkind) * para_nDGrid%stepQ(:)

              CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
              IF (.NOT. nb_Gonly) CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp,nderiv=0)

              IF (debug) THEN
                write(out_unit,"(a,3i7,a,10i3)") 'iGPtot,iGp,isign,Qact',&
                                       iGPtot,iGp,isign,' : ',nDinit(:)
                write(out_unit,*) 'Grid',iGPtot,Qact(1:mole%nb_act),   &
                             dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
              END IF

              write(nioGrid1,*) 'Grid1',iGPtot,nDinit(:),               &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
              write(nioGrid2,*) 'Grid2',iGPtot,Qact(1:mole%nb_act),     &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            END DO
          END IF
        END DO
        para_nDGrid%nb_G = iGPtot
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"

      END DO

      CALL file_close(para_nDGrid%Grid_FOR_Fit_file1)
      CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

      IF (nb_Gonly) STOP 'nb_Gonly=t'


      CALL dealloc_array(nDsize,'nDsize',name_sub)
      CALL dealloc_array(nDinit,'nDinit',name_sub)
      DO ic=1,MaxCoupling
        CALL dealloc_nDindex(Tab_PointGridSign(ic))
      END DO
      CALL dealloc_array(Tab_PointGridSign,'Tab_PointGridSign',name_sub)

      CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

    END SUBROUTINE sub_nDGrid_v1
      SUBROUTINE sub_nDGrid_WiTHOUT_calc(para_nDGrid,Qact,mole,PrimOp)
      USE EVR_system_m
      USE mod_Constant,  only : get_Conv_au_TO_unit
      USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)      :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (PrimOp_t)    :: PrimOp
!=====================================================================


!----- variables divers ----------------------------------------------
      real(kind=Rkind), pointer :: Q(:),val(:)
      integer :: iGPtot

      integer           :: nioGrid1,nioGrid2

      integer                  :: idum
      character (len=Name_len) :: name_dum
      real(kind=Rkind)         :: auTOenergy


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDGrid_WiTHOUT_calc'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       auTOenergy = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)

       CALL get_Qact0(Qact,mole%ActiveTransfo)


       CALL alloc_NParray(para_nDGrid%Q0,[mole%nb_act],                 &
                         'para_nDGrid%Q0',name_sub)
        para_nDGrid%Q0(:) = Qact(1:mole%nb_act)

        nullify(Q)
        CALL alloc_array(Q,[mole%nb_act],'Q',name_sub)
        nullify(val)
        CALL alloc_array(val,[4],'val',name_sub)

        para_nDGrid%Grid_FOR_Fit_file1%name = adjustl(trim(para_nDGrid%name_Grid)) // '1'
        para_nDGrid%Grid_FOR_Fit_file2%name = adjustl(trim(para_nDGrid%name_Grid)) // '2'

        CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        iGPtot = 0
        DO
          !read(nioGrid2,*,IOSTAT=err_read) name_dum,idum,Q(:),val(:)
          read(nioGrid2,*,IOSTAT=err_read) name_dum,idum,Q(:)

          IF (err_read < 0) EXIT
          !write(out_unit,*) name_dum,idum,Q(:),val(1)/auTOenergy,val(2:4)
          iGPtot = iGPtot + 1
        END DO
        para_nDGrid%nb_G = iGPtot

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

        CALL dealloc_array(Q,'Q',name_sub)
        CALL dealloc_array(val,'val',name_sub)


      END SUBROUTINE sub_nDGrid_WiTHOUT_calc
      SUBROUTINE sub_nDFit_v2(para_nDFit,para_nDGrid,Qact,mole,PrimOp)
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(in)      :: para_nDGrid(:)
      TYPE (param_nDFit),  intent(inout)   :: para_nDFit

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (PrimOp_t)    :: PrimOp
!=====================================================================


!----- variables divers ----------------------------------------------

      integer, allocatable            :: nDinit(:)
      integer, allocatable            :: nDvalG(:),nDvalB(:)
      real (kind=Rkind)               :: valGP
      real (kind=Rkind), allocatable  :: Q(:)
      integer                  :: i,ic,iGP,iB,jB,idum,nb_B,nb_WB,iG,nb_coupling
      character (len=Name_len) :: name_dum
      character (len=Name_len) :: name_file
      character (len=Name_len) :: name_int

      integer           :: nioGrid1,nioGrid2
      integer           :: nioFit

      real (kind=Rkind), allocatable :: F(:)   ! size nb_B ! intermediate
      real (kind=Rkind), allocatable :: A(:,:) ! size nb_B,nb_B
      real (kind=Rkind), allocatable :: B(:)   ! size nb_B ! result
      real (kind=Rkind), allocatable :: w(:)
      real (kind=Rkind), allocatable :: vv(:,:)
      real (kind=Rkind)              :: wmin,wmax,conv_ene,Weight_iGP


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDFit_v2'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------


  write(out_unit,*) "======================================"
  write(out_unit,*) "=== FIT the Grids ===================="
  write(out_unit,*) "======================================"
  write(out_unit,*) "lbound,ubound",lbound(para_nDGrid),ubound(para_nDGrid)

  conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)

  CALL Read_nDFit(para_nDFit,Qact(1:mole%nb_act))

  CALL alloc_NParray(Q,[mole%nb_act],'Q',name_sub)
  CALL alloc_NParray(nDvalG,[mole%nb_act],'nDvalG',name_sub)
  CALL alloc_NParray(nDvalB,[mole%nb_act],'nDvalB',name_sub)
  CALL alloc_NParray(nDinit,[mole%nb_act],'nDinit',name_sub)
  nDinit(:) = 0

  CALL init_nDindexPrim(para_nDFit%nDindB,mole%nb_act,para_nDFit%nDsize,        &
                        nDweight=para_nDFit%nDweight,type_OF_nDindex=0,         &
                        nDinit=nDinit,                                          &
                        MinNorm=para_nDFit%MinNorm,MaxNorm=para_nDFit%MaxNorm,  &
                        MinCoupling=para_nDFit%MinCoupling,                     &
                        MaxCoupling=para_nDFit%MaxCoupling)
  CALL sort_nDindex(para_nDFit%nDindB)
  CALL Write_nDindex(para_nDFit%nDindB)
  nb_B = para_nDFit%nDindB%Max_nDI

  CALL alloc_NParray(F,[nb_B],'F',name_sub)
  CALL alloc_NParray(B,[nb_B],'B',name_sub)
  B(:) = ZERO
  CALL alloc_NParray(A,[nb_B,nb_B],'A',name_sub)
  A(:,:) = ZERO
  CALL alloc_NParray(w,[nb_B],'w',name_sub)
  CALL alloc_NParray(vv,[nb_B,nb_B],'vv',name_sub)


  DO iG=lbound(para_nDGrid,dim=1),ubound(para_nDGrid,dim=1)
  DO iGP=1,para_nDGrid(iG)%nb_G

    CALL calc_nDindex(para_nDGrid(iG)%nDindG,iGP,nDvalG)

    ValGP = para_nDGrid(iG)%Val(iGP)

    Q(:) = para_nDGrid(iG)%Q0
    IF (allocated(para_nDGrid(iG)%Tab_stepQ)) THEN
      DO i=1,mole%nb_act
        Q(i) = Q(i) + para_nDGrid(iG)%Tab_stepQ(i)%d0(nDvalG(i)+1)
      END DO
    END IF

    IF (para_nDFit%Col_FOR_WeightOFFit > 0) THEN
      Weight_iGP = exp(-para_nDFit%Scal_FOR_WeightOFFit * ValGP)
    ELSE
      Weight_iGP = ONE
    END IF

    DO iB=1,nb_B
      F(iB) = nDFunct_WITH_Q(Q,iB,para_nDFit)
    END DO
    !CALL Write_Vec(F, out_unit, 5, name_info='F')

    DO iB=1,nb_B
       DO jB=1,nb_B
         A(jB,iB) = A(jB,iB) + Weight_iGP*F(iB)*F(jB)
       END DO
       B(iB)      = B(iB)    + Weight_iGP*F(iB)*ValGP
    END DO

  END DO
  END DO

  !CALL Write_Mat(A, out_unit, 5, name_info='A')
  !CALL Write_Vec(B, out_unit, 5, name_info='B')

  IF (para_nDFit%svd) THEN
    !une facon.... SVD
    !write(out_unit,*) 'a',a
    !write(out_unit,*) 'b',b
    !CALL SVDCMP(A,nb_B,nb_B,w,vv,nb_B)
    !Find maximum singular value
    !wmax = maxval(w)
    !Define "small"
    !wmin = wmax * para_nDFit%epsi
    !Zero the "small" singular values
    !WHERE (w<wmin) w = ZERO
    !F(:) = B(:)
    !CALL SVBKSB(a,w,vv,nb_B,nb_B,F,b,nb_B)

    F(:) = B(:)
    B = LinearSys_Solve(A,F,LS_type=1,epsi=para_nDFit%epsi)

  ELSE
    STOP 'not SVD, not yet'
  END IF
  write(out_unit,*) "======================================"
  write(out_unit,*) "======================================"

  para_nDFit%ndim = mole%nb_act
  CALL ReadWrite_nDFitW(para_nDFit,.FALSE.,B=B,conv_ene=conv_ene)

  ! check the fit
  CALL sub_ChecknDFit2_tabGrid(para_nDFit,para_nDGrid)

    STOP 'coucou'


        ! check the fit
        CALL Analysis_nDFit(para_nDFit,conv_ene=conv_ene)
        ! check the fit
        !CALL Analysis_nDFitW(para_nDFit,conv_ene=conv_ene)


        CALL dealloc_NParray(Q,'Q',name_sub)
        CALL dealloc_NParray(nDinit,'nDinit',name_sub)
        CALL dealloc_NParray(nDvalB,'nDvalB',name_sub)
        CALL dealloc_NParray(nDvalG,'nDvalG',name_sub)

        CALL dealloc_NParray(F,'F',name_sub)
        CALL dealloc_NParray(B,'B',name_sub)
        CALL dealloc_NParray(A,'A',name_sub)
        CALL dealloc_NParray(w,'w',name_sub)
        CALL dealloc_NParray(vv,'vv',name_sub)


      END SUBROUTINE sub_nDFit_v2

      SUBROUTINE sub_nDFit_v1(para_nDFit,para_nDGrid,Qact,mole,PrimOp)
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid
      TYPE (param_nDFit),  intent(inout) :: para_nDFit

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (PrimOp_t)    :: PrimOp
!=====================================================================


!----- variables divers ----------------------------------------------

      integer, pointer         :: nDinit(:)
      integer, pointer         :: nDvalG(:),nDvalB(:)
      real (kind=Rkind), pointer :: Q(:)
      integer                  :: i,ic,iGP,iB,jB,idum,nb_B,nb_WB
      character (len=Name_len) :: name_dum
      character (len=Name_len) :: name_file
      character (len=Name_len) :: name_int

      integer           :: nioGrid1,nioGrid2
      integer           :: nioFit

      real (kind=Rkind), pointer :: F(:)   ! size nb_B ! intermediate
      real (kind=Rkind), pointer :: A(:,:) ! size nb_B,nb_B
      real (kind=Rkind), pointer :: B(:)   ! size nb_B ! result
      real (kind=Rkind), pointer :: val(:) ! number of data which can be fit
      real (kind=Rkind), pointer :: w(:)
      real (kind=Rkind), pointer :: vv(:,:)
      real (kind=Rkind)          :: wmin,wmax,conv_ene,Weight_iGP


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDFit_v1'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!
!=====================================================================
        conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)

        CALL Read_nDFit(para_nDFit,Qact(1:mole%nb_act))

        nullify(val)
        CALL alloc_array(val,[para_nDFit%nb_val],'val',name_sub)
        val(:) = ZERO

        nullify(Q)
        CALL alloc_array(Q,[mole%nb_act],'Q',name_sub)

        nullify(nDinit)
        CALL alloc_array(nDinit,[mole%nb_act],'nDinit',name_sub)
        nDinit(:) = 1

        nullify(nDvalB)
        CALL alloc_array(nDvalB,[mole%nb_act],'nDvalB',name_sub)
        nullify(nDvalG)
        CALL alloc_array(nDvalG,[mole%nb_act],'nDvalG',name_sub)


        CALL init_nDindexPrim(para_nDFit%nDindB,mole%nb_act,para_nDFit%nDsize,     &
                              nDweight=para_nDFit%nDweight,type_OF_nDindex=0,      &
                              nDinit=nDinit,                            &
                              MinNorm=para_nDFit%MinNorm,MaxNorm=para_nDFit%MaxNorm,&
                              MinCoupling=para_nDFit%MinCoupling,MaxCoupling=para_nDFit%MaxCoupling)
        CALL sort_nDindex(para_nDFit%nDindB)
        para_nDFit%nDindB%Tab_nDval(:,:) = para_nDFit%nDindB%Tab_nDval(:,:) - 1
        CALL Write_nDindex(para_nDFit%nDindB)
        nb_B = para_nDFit%nDindB%Max_nDI

        nullify(F)
        CALL alloc_array(F,[nb_B],'F',name_sub)
        nullify(B)
        CALL alloc_array(B,[nb_B],'B',name_sub)
        B(:) = ZERO
        nullify(A)
        CALL alloc_array(A,[nb_B,nb_B],'A',name_sub)
        A(:,:) = ZERO
        nullify(w)
        CALL alloc_array(w,[nb_B],'w',name_sub)
        nullify(vv)
        CALL alloc_array(vv,[nb_B,nb_B],'vv',name_sub)

        write(out_unit,*) "======================================"
        write(out_unit,*) "=== FIT =============================="
        write(out_unit,*) "======================================"
        CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)

        DO iGP=1,para_nDGrid%nb_G

          read(nioGrid2,*) name_dum,idum,Q(:),val(:)

          IF (para_nDFit%Col_FOR_WeightOFFit > 0) THEN
            Weight_iGP = exp(-para_nDFit%Scal_FOR_WeightOFFit *         &
                                    val(para_nDFit%Col_FOR_WeightOFFit))
          ELSE
            Weight_iGP = ONE
          END IF

          DO iB=1,nb_B
            F(iB) = nDFunct_WITH_Q(Q,iB,para_nDFit)
          END DO

           DO iB=1,nb_B
             DO jB=1,nb_B
               A(jB,iB) = A(jB,iB) + Weight_iGP*F(iB)*F(jB)
             END DO
             B(iB)      = B(iB)    + Weight_iGP*F(iB)*val(para_nDFit%ind_val)
           END DO

        END DO

        IF (para_nDFit%svd) THEN
          !une facon.... SVD
          !write(out_unit,*) 'a',a
          !write(out_unit,*) 'b',b
          !CALL SVDCMP(A,nb_B,nb_B,w,vv,nb_B)
          !Find maximum singular value
          !wmax = maxval(w)
          !Define "small"
          !wmin = wmax * para_nDFit%epsi
          !Zero the "small" singular values
          !WHERE (w<wmin) w = ZERO

          !F(:) = B(:)
          !CALL SVBKSB(a,w,vv,nb_B,nb_B,F,b,nb_B)

          F(:) = B(:)
          B = LinearSys_Solve(A,F,LS_type=1,epsi=para_nDFit%epsi)
        ELSE
          STOP 'not SVD, not yet'
        END IF
        CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"


        para_nDFit%ndim = mole%nb_act
        CALL ReadWrite_nDFitW(para_nDFit,.FALSE.,B=B,conv_ene=conv_ene)

        ! check the fit
        CALL sub_ChecknDFit2(para_nDFit,para_nDGrid)


        ! check the fit
        CALL Analysis_nDFit(para_nDFit,conv_ene=conv_ene)
        ! check the fit
        !CALL Analysis_nDFitW(para_nDFit,conv_ene=conv_ene)


        CALL dealloc_array(Q,'Q',name_sub)
        CALL dealloc_array(nDinit,'nDinit',name_sub)
        CALL dealloc_array(nDvalB,'nDvalB',name_sub)
        CALL dealloc_array(nDvalG,'nDvalG',name_sub)

        CALL dealloc_array(F,'F',name_sub)
        CALL dealloc_array(B,'B',name_sub)
        CALL dealloc_array(A,'A',name_sub)
        CALL dealloc_array(w,'w',name_sub)
        CALL dealloc_array(vv,'vv',name_sub)


      END SUBROUTINE sub_nDFit_v1

      SUBROUTINE sub_ChecknDFit2(para_nDFit,para_nDGrid)
      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid
      TYPE (param_nDFit),  intent(inout) :: para_nDFit


      integer                    :: iGP,iB
      real (kind=Rkind), pointer :: val(:),Q(:)

      integer                    :: idum
      character (len=Name_len)   :: name_dum

      integer           :: nioGrid1,nioGrid2
      real (kind=Rkind) :: val_fit,max_err,RMS_err


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_ChecknDFit2'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

!
!=====================================================================
        nullify(val)
        CALL alloc_array(val,[para_nDFit%ind_val],'val',name_sub)
        val(:) = ZERO

        nullify(Q)
        CALL alloc_array(Q,[para_nDFit%nDindB%ndim],'Q',name_sub)


        write(out_unit,*) "======================================"
        write(out_unit,*) "=== CHECKING FIT ====================="
        write(out_unit,*) "======================================"

        CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)

        max_err = ZERO
        RMS_err = ZERO

        DO iGP=1,para_nDGrid%nb_G

          read(nioGrid2,*) name_dum,idum,Q(:),val(:)

          CALL sub_nDFunc_FROM_nDFit(val_fit,Q,para_nDFit)
          IF (iGP == 1) write(out_unit,*) 'Q (1st point),val_fit,val',Q,val_fit,val


          IF (abs(val(para_nDFit%ind_val)-val_fit) > FIVE*ONETENTH**4) THEN
            IF (para_nDFit%ind_val == 1) THEN
              write(out_unit,*) 'iG,val,val_fit,diff (cm-1)',iGP,        &
                                      val(para_nDFit%ind_val),val_fit,    &
                                    abs(val(para_nDFit%ind_val)-val_fit)* &
                                      get_Conv_au_TO_unit('E','cm-1')
            ELSE
              write(out_unit,*) 'iG,val,val_fit,diff',iGP,               &
                                      val(para_nDFit%ind_val),val_fit,    &
                                      abs(val(para_nDFit%ind_val)-val_fit)
            END IF
          END IF


          val_fit = abs(val(para_nDFit%ind_val)-val_fit)
          IF (val_fit > max_err) max_err = val_fit
          RMS_err = RMS_err + val_fit**2
        END DO
        CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)


        RMS_err = sqrt(RMS_err/para_nDGrid%nb_G)

        IF (para_nDFit%ind_val == 1) THEN
          write(out_unit,*) 'Max_error (cm-1): ',max_err*              &
                                         get_Conv_au_TO_unit('E','cm-1')

          write(out_unit,*) 'RMS_error (cm-1): ',RMS_err*              &
                                         get_Conv_au_TO_unit('E','cm-1')

        ELSE
          write(out_unit,*) 'Max_error: ',max_err
          write(out_unit,*) 'RMS_error: ',RMS_err
        END IF

        write(out_unit,*) "======================================"
        write(out_unit,*) "=== END CHECKING FIT ================="
        write(out_unit,*) "======================================"


        CALL dealloc_array(val,'val',name_sub)
        CALL dealloc_array(Q,'Q',name_sub)


      END SUBROUTINE sub_ChecknDFit2

      SUBROUTINE sub_ChecknDFit2_tabGrid(para_nDFit,para_nDGrid)
      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      USE mod_nDindex
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(in)   :: para_nDGrid(:)
      TYPE (param_nDFit), intent(inout) :: para_nDFit


      integer                         :: i,iG,iGP,iB
      real (kind=Rkind), allocatable  :: Q(:)
      integer,           allocatable  :: nDvalG(:)

      real (kind=Rkind) :: val_fit,max_err,RMS_err,valGP


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_ChecknDFit2_tabGrid'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------


        CALL alloc_NParray(Q,[para_nDFit%nDindB%ndim],'Q',name_sub)
        CALL alloc_NParray(nDvalG,[para_nDFit%nDindB%ndim],'nDvalG',name_sub)


        write(out_unit,*) "======================================"
        write(out_unit,*) "=== CHECKING FIT ====================="
        write(out_unit,*) "======================================"
        write(out_unit,*) "lbound,ubound",lbound(para_nDGrid),ubound(para_nDGrid)
        max_err = ZERO
        RMS_err = ZERO

        DO iG=lbound(para_nDGrid,dim=1),ubound(para_nDGrid,dim=1)
        DO iGP=1,para_nDGrid(iG)%nb_G

          CALL calc_nDindex(para_nDGrid(iG)%nDindG,iGP,nDvalG)

          ValGP = para_nDGrid(iG)%Val(iGP)

          Q(:) = para_nDGrid(iG)%Q0
          IF (allocated(para_nDGrid(iG)%Tab_stepQ)) THEN
            DO i=1,size(para_nDGrid(iG)%Q0)
              Q(i) = Q(i) + para_nDGrid(iG)%Tab_stepQ(i)%d0(nDvalG(i)+1)
            END DO
          END IF

          CALL sub_nDFunc_FROM_nDFit(val_fit,Q,para_nDFit)

          IF (abs(ValGP-val_fit) > FIVE*ONETENTH**4) THEN
            IF (para_nDFit%ind_val == 1) THEN
              write(out_unit,*) 'iG,val,val_fit,diff (cm-1)',iGP,ValGP,val_fit,&
                              abs(ValGP-val_fit)*get_Conv_au_TO_unit('E','cm-1')
            ELSE
              write(out_unit,*) 'iG,val,val_fit,diff',iGP, ValGP,val_fit,abs(ValGP-val_fit)
            END IF
          END IF


          val_fit = abs(ValGP-val_fit)
          IF (val_fit > max_err) max_err = val_fit
          RMS_err = RMS_err + val_fit**2
        END DO
        END DO


        RMS_err = sqrt(RMS_err/sum(para_nDGrid(:)%nb_G))

        IF (para_nDFit%ind_val == 1) THEN
          write(out_unit,*) 'Max_error (cm-1): ',max_err*              &
                                         get_Conv_au_TO_unit('E','cm-1')

          write(out_unit,*) 'RMS_error (cm-1): ',RMS_err*              &
                                         get_Conv_au_TO_unit('E','cm-1')

        ELSE
          write(out_unit,*) 'Max_error: ',max_err
          write(out_unit,*) 'RMS_error: ',RMS_err
        END IF

        write(out_unit,*) "======================================"
        write(out_unit,*) "=== END CHECKING FIT ================="
        write(out_unit,*) "======================================"


        CALL dealloc_NParray(Q,'Q',name_sub)
        CALL dealloc_NParray(nDvalG,'nDvalG',name_sub)


      END SUBROUTINE sub_ChecknDFit2_tabGrid
      SUBROUTINE sub_nGrid1_TO_nGrid2(Q0)
      USE EVR_system_m
      IMPLICIT NONE

      real (kind=Rkind) :: Q0(:)

      integer, pointer           :: nDvalG1(:),nDvalG2(:)
      integer, pointer           :: list1(:)
      real (kind=Rkind), pointer :: val(:),Q1(:),Q2(:)



      integer           :: nioGrid1,nioGrid2,i
      integer                    :: idum
      character (len=Name_len)   :: name_dum

      ! for the namelist
      integer                    :: nb_act1,nb_act2
      character (len=Line_len)   :: name_Grid1,name_Grid2
      character (len=Line_len)   :: name_file

      namelist /nGrid1_TO_nGrid2/ nb_act1,name_Grid1,name_Grid2

      ! for the namelist


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nGrid1_TO_nGrid2'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

!
!=====================================================================
        nb_act2 = size(Q0)

        write(out_unit,*) 'Q0',Q0(:)

        nb_act1    = 0
        name_Grid1 = ''
        name_Grid2 = ''
        read(in_unit,nGrid1_TO_nGrid2,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' End-of-file or End-of-record'
          write(out_unit,*) ' The namelist "nGrid1_TO_nGrid2" is probably absent'
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Some parameter name of the namelist "nGrid1_TO_nGrid2" are probaly wrong'
          write(out_unit,*) ' It should never append !!'
          write(out_unit,*) ' Check the fortran'
          write(out_unit,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unit,nGrid1_TO_nGrid2)

        nullify(val)
        CALL alloc_array(val,[4],'val',name_sub)
        nullify(Q1)
        CALL alloc_array(Q1,[nb_act1],'Q1',name_sub)
        nullify(Q2)
        CALL alloc_array(Q2,[nb_act2],'Q2',name_sub)

        nullify(nDvalG1)
        CALL alloc_array(nDvalG1,[nb_act1],'nDvalG1',name_sub)
        nullify(list1)
        CALL alloc_array(list1,[nb_act1],'list1',name_sub)
        nullify(nDvalG2)
        CALL alloc_array(nDvalG2,[nb_act2],'nDvalG2',name_sub)


        read(in_unit,*) list1(:) ! for list of active variables of grid 1
        IF (debug) write(out_unit,*) list1(:)
        flush(out_unit)


        write(out_unit,*) "======================================"
        write(out_unit,*) "=== READING THE GRIDS ================"
        write(out_unit,*) "======================================"
        flush(out_unit)

        ! for the file with the grid points
        name_file = adjustl(trim(name_Grid1)) // '2'
        CALL file_open2(name_file,nioGrid1)
        name_file = adjustl(trim(name_Grid2)) // '2'
        CALL file_open2(name_file,nioGrid2)

        DO

          read(nioGrid1,*,IOSTAT=err_read) name_dum,idum,Q1(:),val(:)
          IF (err_read /= 0) EXIT

          Q2(:) = Q0(:)
          DO i=1,nb_act1
            Q2(list1(i)) = Q1(i)
          END DO
          write(nioGrid2,*) name_dum,idum,Q2(:),val(:)


        END DO
        close(nioGrid1) ! CALL file_close cannot be used
        close(nioGrid2) ! CALL file_close cannot be used

        ! for the file with the basis indices
        name_file = adjustl(trim(name_Grid1)) // '1'
        CALL file_open2(name_file,nioGrid1)
        name_file = adjustl(trim(name_Grid2)) // '1'
        CALL file_open2(name_file,nioGrid2)

        DO

          read(nioGrid1,*,IOSTAT=err_read) name_dum,idum,nDvalG1(:),val(:)
          IF (err_read /= 0) EXIT

          nDvalG2(:) = 1
          DO i=1,nb_act1
            nDvalG2(list1(i)) = nDvalG1(i)
          END DO
          write(nioGrid2,*) name_dum,idum,nDvalG2(:),val(:)


        END DO
        close(nioGrid1) ! CALL file_close cannot be used
        close(nioGrid2) ! CALL file_close cannot be used


        CALL dealloc_array(val,    'val',name_sub)
        CALL dealloc_array(Q1,     'Q1',name_sub)
        CALL dealloc_array(Q2,     'Q2',name_sub)
        CALL dealloc_array(nDvalG1,'nDvalG1',name_sub)
        CALL dealloc_array(list1,  'list1',name_sub)
        CALL dealloc_array(nDvalG2,'nDvalG2',name_sub)

      END SUBROUTINE sub_nGrid1_TO_nGrid2

END MODULE mod_nDGridFit
