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
      SUBROUTINE sub_Optimization_OF_VibParam(max_mem)
      USE mod_system
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0
      USE mod_PrimOp
      USE mod_basis
      USE BasisMakeGrid
      USE mod_psi
      USE mod_propa
      USE mod_Op
      USE mod_analysis
      USE mod_Auto_Basis
      USE mod_Optimization
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================
!

!----- variables for the dynamical memory allocation -----------------
      logical  :: intensity_only
      integer  :: max_mem



!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- variables for the construction of H ----------------------------
      TYPE (param_AllOp), target  :: para_AllOp
      TYPE (param_Op),    pointer :: para_H      => null()


!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)           :: para_ana
      TYPE (param_intensity)     :: para_intensity

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa


!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis) :: para_AllBasis
      TYPE (basis)          :: BasisnD_Save

!----- variables divers ----------------------------------------------
      TYPE (param_Optimization)      :: para_Optimization
      real (kind=Rkind)              :: Energ,P1,Norm_min
      integer                        :: err_io,i,i0,i1,iOpt,nb_Opt,n,ib
      real (kind=Rkind), allocatable :: xOpt_min(:),SQ(:),SQini(:),QA(:),QB(:)

      real (kind=Rkind), allocatable :: freq(:),grad(:)
      real (kind=Rkind), allocatable :: Qact(:),Qdyn(:)
      real (kind=Rkind), allocatable :: Qopt(:)
      TYPE (param_dnMatOp)           :: dnMatOp(1)
      character (len=Name_len)       :: name_dum

      TYPE (Type_dnVec) :: dnx

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_Optimization_OF_VibParam'
!para_mem%mem_debug = .TRUE.
!para_mem%mem_print = .FALSE.

!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================


      write(out_unitp,*)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: BEGINNING ini_data'
      CALL time_perso('ini_data')
      write(out_unitp,*)
      write(out_unitp,*)
      CALL     ini_data(const_phys,para_Tnum,mole,                      &
                        para_AllBasis,para_AllOp,                       &
                        para_ana,para_intensity,intensity_only,         &
                        para_propa)

      para_H => para_AllOp%tab_Op(1)
      CALL basis2TObasis1(BasisnD_Save,para_AllBasis%BasisnD) ! basis saved here

      write(out_unitp,*)
      write(out_unitp,*)
      CALL time_perso('ini_data')
      write(out_unitp,*) ' VIB: END ini_data'
      write(out_unitp,*) '================================================='
      write(out_unitp,*)
!=====================================================================

      CALL Read_param_Optimization(para_Optimization,mole,para_AllBasis%BasisnD,read_nml=.TRUE.)
      CALL Write_param_Optimization(para_Optimization)

      ! get the values from mole or Basis or Qopt ...
      allocate(Qopt(para_Optimization%nb_Opt))
      CALL Set_ALL_para_FOR_optimization(mole,para_AllBasis%BasisnD,Qopt,1)

      CALL Sub_Optimization(para_AllBasis%BasisnD,     &
                            para_Tnum,mole,para_H%para_ReadOp%PrimOp_t,&
                            Qopt,para_Optimization)

      write(out_unitp,*) '============ FINAL ANLYSIS =================='

      allocate(Qact(mole%nb_var))
      CALL get_Qact0(Qact,mole%ActiveTransfo)
      Qact(1:para_Optimization%nb_Opt) = Qopt
      IF (para_Optimization%FinalEnergy) THEN
        write(out_unitp,*) '============ ENERGY:'
        CALL set_print_level(0,force=.TRUE.)  ! print_level = 0
        CALL Sub_Energ_OF_ParamBasis(Energ,para_Optimization%xOpt_min,para_Optimization%nb_Opt,BasisnD_Save,&
                                     para_Tnum,mole,                    &
                                     para_H%para_ReadOp%PrimOp_t,Qact)
        write(out_unitp,*) 'Optimal param',xOpt_min,' Energy',Energ

      END IF
      IF (para_Optimization%Freq) THEN
        write(out_unitp,*) '============ Freq:'
        CALL alloc_NParray(freq,[para_Optimization%nb_Opt],'freq',name_sub)
        CALL sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,                 &
                           para_H%para_ReadOp%PrimOp_t,print_freq=.TRUE.)
        CALL dealloc_NParray(freq,'freq',name_sub)
      END IF

      IF (para_Optimization%Grad) THEN
        write(out_unitp,*) '============ GRAD:'
        CALL alloc_NParray(Grad,[nb_Opt],'Grad',name_sub)

        CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,                   &
                           para_H%para_ReadOp%PrimOp_t%nb_elec,nderiv=1)

        CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,           &
                                 para_H%para_ReadOp%PrimOp_t)

        CALL Get_Grad_FROM_Tab_OF_dnMatOp(Grad,dnMatOp) ! for the first electronic state
        write(out_unitp,*) 'Grad: size, RMS',size(grad),sqrt(sum(grad**2)/size(grad))
        write(out_unitp,*) 'Grad',Grad

        CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
        CALL dealloc_NParray(Grad,'Grad',name_sub)

      END IF

      IF (para_Optimization%Optimization_param == 'geometry') THEN
        write(out_unitp,*) '============ GEOM:'
        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv=0)

        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv=0,Gcenter=.FALSE.,WriteCC=.TRUE.)

        CALL dealloc_dnSVM(dnx)
        write(out_unitp,*) 'Qact geometry:'
        DO i=1,size(Qact)
          write(out_unitp,*) Qact(i)
        END DO
        write(out_unitp,*) 'END Qact geometry:'


        CALL alloc_NParray(Qdyn,[mole%nb_var],'Qdyn',name_sub)
        CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
        write(out_unitp,*) 'Qdyn geometry:'
        DO i=1,size(Qdyn)
          write(out_unitp,*) Qdyn(i)
        END DO
        write(out_unitp,*) 'END Qdyn geometry:'

        CALL dealloc_NParray(Qdyn,'Qdyn',name_sub)
        CALL dealloc_NParray(Qact,'Qact',name_sub)

      END IF
      write(out_unitp,*) '========= END FINAL ANLYSIS =================='


      CALL dealloc_NParray(para_FOR_optimization%Val_RVec,'para_FOR_optimization%Val_RVec',name_sub)
      CALL dealloc_NParray(para_FOR_optimization%opt_RVec,'para_FOR_optimization%opt_RVec',name_sub)

!=====================================================================
      CALL dealloc_table_at(const_phys%mendeleev)

      CALL dealloc_CoordType(mole)
      IF (associated(para_Tnum%Gref)) THEN
        CALL dealloc_array(para_Tnum%Gref,"para_Tnum%Gref",name_sub)
      END IF
      CALL dealloc_Tnum(para_Tnum)

      CALL dealloc_para_AllOp(para_AllOp)


      CALL dealloc_param_propa(para_propa)

      CALL dealloc_AllBasis(para_AllBasis)

      write(out_unitp,*) 'mem_tot',para_mem%mem_tot

  END SUBROUTINE sub_Optimization_OF_VibParam

  SUBROUTINE Sub_Energ_OF_ParamBasis(Energ,xOpt,nb_Opt,BasisnD,     &
                                         para_Tnum,mole,PrimOp,Qact)
      USE mod_system
      USE mod_nDindex
      USE mod_Constant
      use mod_Coord_KEO, only: CoordType, tnum
      USE mod_PrimOp

      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

      real(kind=Rkind)  :: Energ ! parameter to be minimized
      integer           :: nb_Opt
      real (kind=Rkind) :: xOpt(nb_Opt)

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
      logical        :: Cart_Transfo_save


      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout) :: BasisnD

!----- variables pour la namelist minimum ----------------------------
      TYPE (PrimOp_t)  :: PrimOp
      integer          :: nb_scalar_Op
      logical          :: calc_scalar_Op

!----- variables for the construction of H ---------------------------
      TYPE (param_ReadOp)         :: para_ReadOp
      logical                     :: Save_FileGrid,Save_MemGrid


!----- local variables -----------------------------------------------
!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp)          :: para_AllOp_loc

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis)       :: para_AllBasis_loc
      TYPE (basis)                :: basis_temp

      integer            :: i,i0,i1,iOpt
      real(kind=Rkind)   :: ene0
      TYPE(REAL_WU)      :: RWU_E,RWU_DE

      TYPE (param_dnMatOp) :: dnMatOp(1)
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Sub_Energ_OF_ParamBasis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        !CALL Write_CoordType(mole)
        write(out_unitp,*) 'basis to be optimized'
        write(out_unitp,*)
        !CALL RecWrite_basis(basis_temp)
      END IF
!---------------------------------------------------------------------


    ! For cubature, it works directly with xOpt
    IF (para_FOR_optimization%Optimization_param /= 'cubature') THEN
      para_FOR_optimization%Val_RVec(:) = xOpt(:)
      CALL Set_ALL_para_FOR_optimization(mole,BasisnD,Qact,1)
    END IF



    SELECT CASE (para_FOR_optimization%Optimization_param)

    CASE ('coordbasis','coordbasis_avene')

      CALL basis2TObasis1(basis_temp,BasisnD)

      ! save some parameters of PrimOp, Tnum, mole
      nb_scalar_Op              = PrimOp%nb_scalar_Op
      PrimOp%nb_scalar_Op       = 0

      calc_scalar_Op            = PrimOp%calc_scalar_Op
      PrimOp%calc_scalar_Op     = .FALSE.

      Cart_Transfo_save         = mole%Cart_transfo
      mole%Cart_transfo         = .FALSE.

      ! init0 of para_H
      para_AllOp_loc%nb_Op = 2 ! just H and S
      CALL alloc_array(para_AllOp_loc%tab_Op,[para_AllOp_loc%nb_Op],  &
                      "para_AllOp_loc%tab_Op",name_sub)


      !CALL RecWrite_basis(basis_temp,write_all=.TRUE.)

      CALL basis_TO_Allbasis(basis_temp,para_AllBasis_loc,mole)
      CALL dealloc_basis(basis_temp)

      CALL init_ReadOp(para_ReadOp)
      para_ReadOp%nb_bRot                       = 1
      para_ReadOp%make_Mat                      = .TRUE.
      para_ReadOp%para_FileGrid%Save_MemGrid    = .TRUE.
      para_ReadOp%para_FileGrid%Save_FileGrid   = .FALSE.
      para_ReadOp%comput_S                      = .FALSE.
      para_ReadOp%para_FileGrid%First_GridPoint = 1
      para_ReadOp%para_FileGrid%Last_GridPoint  = get_nq_FROM_basis(para_AllBasis_loc%BasisnD)
      para_ReadOp%para_FileGrid%Restart_Grid    = .FALSE.
      para_ReadOp%para_FileGrid%Test_Grid       = .FALSE.


      CALL Auto_basis(para_Tnum,mole,para_AllBasis_loc,para_ReadOp)


      !---------------------------------------------------------------
      ! make Operators: H and S
      !i=1 => for H
      CALL All_param_TO_para_H(para_Tnum,mole,para_AllBasis_loc,        &
                               para_ReadOp,para_AllOp_loc%tab_Op(1))

      !i=2 => for S
      CALL param_Op1TOparam_Op2(para_AllOp_loc%tab_Op(1),               &
                                para_AllOp_loc%tab_Op(2))
      para_AllOp_loc%tab_Op(2)%name_Op = 'S'
      para_AllOp_loc%tab_Op(2)%n_Op    = -1
      para_AllOp_loc%tab_Op(2)%cplx    = .FALSE.
      para_AllOp_loc%tab_Op(2)%nb_Term = 1

      CALL alloc_NParray(para_AllOp_loc%tab_Op(2)%derive_termQact,        &
              [2,1],"para_AllOp_loc%tab_Op(2)%derive_termQact",name_sub)
      para_AllOp_loc%tab_Op(2)%derive_termQact(:,:) = 0

      CALL alloc_NParray(para_AllOp_loc%tab_Op(2)%derive_termQdyn,        &
              [2,1],"para_AllOp_loc%tab_Op(2)%derive_termQdyn",name_sub)
      para_AllOp_loc%tab_Op(2)%derive_termQdyn(:,:) = 0

      IF (debug) THEN
        DO i=1,para_AllOp_loc%nb_Op
          CALL write_param_Op(para_AllOp_loc%tab_Op(i))
        END DO
      END IF

      !---------------------------------------------------------------
      ! make the Grid
      CALL sub_qa_bhe(para_AllOp_loc)

      !---------------------------------------------------------------
      ! make the matrix of H
      CALL sub_MatOp(para_AllOp_loc%tab_Op(1),debug)

      !---------------------------------------------------------------
      ! digonalization of H
      para_AllOp_loc%tab_Op(1)%diago = .TRUE.
      CALL alloc_para_Op(para_AllOp_loc%tab_Op(1),Grid=.FALSE.,Mat=.TRUE.)
      CALL sub_diago_H(para_AllOp_loc%tab_Op(1)%Rmat,                   &
                       para_AllOp_loc%tab_Op(1)%Rdiag,                  &
                       para_AllOp_loc%tab_Op(1)%Rvp,                    &
                       para_AllOp_loc%tab_Op(1)%nb_tot,                 &
                       para_AllOp_loc%tab_Op(1)%sym_Hamil)

      IF (print_level > -1) THEN
        ! Energy levels
        ene0 = para_AllOp_loc%tab_Op(1)%Rdiag(1)
        write(out_unitp,*) 'levels: '
        DO i=1,para_AllOp_loc%tab_Op(1)%nb_tot

          RWU_E  = REAL_WU(para_AllOp_loc%tab_Op(1)%Rdiag(i),     'au','E')
          RWU_DE = REAL_WU(para_AllOp_loc%tab_Op(1)%Rdiag(i)-ene0,'au','E')

          write(out_unitp,*) i,                                         &
               RWU_Write(RWU_E, WithUnit=.FALSE.,WorkingUnit=.FALSE.),  &
               RWU_Write(RWU_DE,WithUnit=.TRUE., WorkingUnit=.FALSE.)
        END DO
        flush(out_unitp)
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'Eigenvectors of ',name_sub
        CALL Write_Mat(para_AllOp_loc%tab_Op(1)%Rvp,out_unitp,5)
      END IF
      !---------------------------------------------------------------

      Energ = para_AllOp_loc%tab_Op(1)%Rdiag(1)
      Energ = sum(para_AllOp_loc%tab_Op(1)%Rdiag(:)) /                  &
                       real(para_AllOp_loc%tab_Op(1)%nb_tot,kind=Rkind)

      !-----------------------------------------------------------------
      ! deallocation ....
      CALL dealloc_AllBasis(para_AllBasis_loc)
      CALL dealloc_para_AllOp(para_AllOp_loc)


      ! restor some parameters of PrimOp, Tnum, mole
      PrimOp%nb_scalar_Op                   = nb_scalar_Op
      PrimOp%calc_scalar_Op                 = calc_scalar_Op
      mole%Cart_transfo                     = Cart_Transfo_save

    CASE ('geometry') ! potential (find a minimum)

      IF (debug) write(out_unitp,*) 'xOpt',xOpt
      IF (debug) write(out_unitp,*) 'Qact',Qact


      CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,PrimOp%nb_elec,nderiv=0)

      CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp)

      Energ = Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp) ! for the first electronic state

      CALL dealloc_Tab_OF_dnMatOp(dnMatOp)


    CASE ('cubature') ! find grid points

      CALL Sub_Energ_FOR_cubature(Energ,xOpt,nb_Opt,BasisnD)
    END SELECT

    IF (debug) THEN
      write(out_unitp,*) 'xOpt,Energ',xOpt,':',Energ
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF


    END SUBROUTINE Sub_Energ_OF_ParamBasis
      SUBROUTINE Set_ALL_para_FOR_optimization(mole,BasisnD,Qact,Set_Val)

      USE mod_system
      use mod_Coord_KEO, only: CoordType
      USE mod_basis
      IMPLICIT NONE


!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(inout)    :: mole
      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)
      TYPE (basis), intent(in)         :: BasisnD
      integer, intent(in)              :: Set_Val

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Set_ALL_para_FOR_optimization'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------
      para_FOR_optimization%nb_OptParam = 0
      para_FOR_optimization%i_OptParam  = 0

      SELECT CASE (para_FOR_optimization%Optimization_param)
      CASE ('cubature')
        para_FOR_optimization%nb_OptParam = BasisnD%nqc*BasisnD%ndim
      CASE ('geometry')
        CALL Set_paramQ_FOR_optimization(Qact,mole,Set_Val)
      CASE ('coordbasis','coordbasis_avene')
        CALL Set_OptimizationPara_FROM_CoordType(mole,Set_Val)
        CALL Set_basis_para_FOR_optimization(BasisnD,Set_Val)
      CASE DEFAULT
        CALL Set_OptimizationPara_FROM_CoordType(mole,Set_Val)
        CALL Set_basis_para_FOR_optimization(BasisnD,Set_Val)
      END SELECT



      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nb_OptParam ',para_FOR_optimization%nb_OptParam

        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

      END SUBROUTINE Set_ALL_para_FOR_optimization
      SUBROUTINE Sub_Energ_FOR_cubature(Energ,xOpt,nb_Opt,BasisnD)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      USE BasisMakeGrid
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout) :: BasisnD

      real(kind=Rkind)  :: Energ ! parameter to be minimized
      integer           :: nb_Opt
      real (kind=Rkind) :: xOpt(nb_Opt)
      real (kind=Rkind) :: x0(BasisnD%ndim,BasisnD%nqc)
      integer           :: deg ! degree


      integer            :: ib,iq,id
      logical            :: err_cuba

      real (kind=Rkind) :: x
      real (kind=Rkind) :: poly_Hermite ! function

      TYPE (Type_nDindex) :: nDindB
      integer            :: nDsize(BasisnD%ndim),nDval(BasisnD%ndim)
      integer            :: nDinit(BasisnD%ndim)

      real (kind=Rkind), allocatable :: M(:,:),W(:),R(:),MM(:,:),MR(:),ER(:)
      logical :: FirstTime = .TRUE.
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Sub_Energ_FOR_cubature'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'basis to be optimized'
        write(out_unitp,*)
        CALL RecWrite_basis(BasisnD)
      END IF
!---------------------------------------------------------------------
      x0(:,:) = reshape(xOpt,shape(x0))
      CALL ReCentered_grid(x0,BasisnD%ndim,BasisnD%nqc)
      CALL ReOriented_grid(x0,BasisnD%ndim,BasisnD%nqc)
      Xopt = reshape(x0,shape(xOpt))

      deg = int(BasisnD%Norm_OF_nDindB)
      !write(out_unitp,*) 'ndim,deg,nqc',BasisnD%ndim,deg,BasisnD%nqc

      nDsize(:) = deg+1
      nDinit(:) = 1
      CALL init_nDindexPrim(nDindB,BasisnD%ndim,nDsize,nDinit=nDinit,   &
                            type_OF_nDindex=0,                          &
                            MaxNorm=BasisnD%Norm_OF_nDindB)
      nDindB%Tab_nDval(:,:) = nDindB%Tab_nDval(:,:) - 1
      IF (FirstTime) CALL Write_nDindex(nDindB)
      FirstTime = .FALSE.

      CALL alloc_NParray(M ,[nDindB%Max_nDI,BasisnD%nqc],'M', name_sub)
      CALL alloc_NParray(W ,[BasisnD%nqc],               'W', name_sub)
      CALL alloc_NParray(R ,[nDindB%Max_nDI],            'R', name_sub)
      CALL alloc_NParray(ER,[nDindB%Max_nDI],            'ER',name_sub)
      CALL alloc_NParray(MM,[BasisnD%nqc,BasisnD%nqc],   'MM',name_sub)
      CALL alloc_NParray(MR,[BasisnD%nqc],               'MR',name_sub)


      R(:)   = ZERO
      M(:,:) = ONE

      DO ib=1,nDindB%Max_nDI
        CALL calc_nDindex(nDindB,ib,nDval)
        IF (count(nDval == 0) == BasisnD%ndim) R(ib) = sqrt(sqrt(PI**BasisnD%ndim))
        DO iq=1,BasisnD%nqc
        DO id=1,BasisnD%ndim
          x = x0(id,iq)
          M(ib,iq) = M(ib,iq)*poly_Hermite(x,nDval(id))*exp(-x*x)
        END DO
        END DO

      END DO
      !CALL Write_VecMat(R,out_unitp,5)
      !CALL Write_VecMat(M,out_unitp,5)


      ! matrices of the linear (square nq*nq) system (Mt.M).W=(Mt.R)
      MM(:,:) = matmul(transpose(M),M)
      MR(:)   = matmul(transpose(M),R)


      ! Solve the linear system AtA.W=AtB
      W = LinearSys_Solve(MM,MR)

      WHERE (W < ZERO) W = W * 0.9_Rkind

      !CALL Write_VecMat(W,out_unitp,5)


      ER(:) = matmul(M,W)-R(:)
      !CALL Write_VecMat(ER,out_unitp,5)

      Energ = sqrt(dot_product(ER,ER)/real(nDindB%Max_nDI,kind=Rkind))

      !write(out_unitp,*) 'Energ',Energ

      CALL dealloc_NParray(M,'M',name_sub)
      CALL dealloc_NParray(W,'W',name_sub)
      CALL dealloc_NParray(R,'R',name_sub)
      CALL dealloc_NParray(ER,'ER',name_sub)
      CALL dealloc_NParray(MM,'MM',name_sub)
      CALL dealloc_NParray(MR,'MR',name_sub)

      CALL dealloc_nDindex(nDindB)
!---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
!---------------------------------------------------------------------

    END SUBROUTINE Sub_Energ_FOR_cubature
      SUBROUTINE Sub_Energ_FOR_cubatureWeight(Energ,xOpt,nb_Opt,BasisnD)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      USE BasisMakeGrid
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout) :: BasisnD

      real(kind=Rkind)  :: Energ ! parameter to be minimized
      integer           :: nb_Opt
      real (kind=Rkind) :: xOpt(nb_Opt)
      real (kind=Rkind) :: x0(BasisnD%ndim,BasisnD%nqc),w(BasisnD%nqc)
      integer           :: deg ! degree


      integer            :: ib,iq,id
      logical            :: err_cuba

      real (kind=Rkind) :: x
      real (kind=Rkind) :: poly_Hermite ! function

      TYPE (Type_nDindex) :: nDindB
      integer            :: nDsize(BasisnD%ndim),nDval(BasisnD%ndim)
      integer            :: nDinit(BasisnD%ndim)

      real (kind=Rkind), allocatable :: M(:,:),R(:),MM(:,:),MR(:),ER(:)
      logical :: FirstTime = .TRUE.
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Sub_Energ_FOR_cubatureWeight'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'basis to be optimized'
        write(out_unitp,*)
        CALL RecWrite_basis(BasisnD)
      END IF
!---------------------------------------------------------------------

      x0(:,:) = reshape(xOpt(1:size(x0)),shape(x0))
      w(:)    = xOpt(size(x0)+1:BasisnD%nqc)
      CALL ReCentered_grid(x0,BasisnD%ndim,BasisnD%nqc)
      CALL ReOriented_grid(x0,BasisnD%ndim,BasisnD%nqc)
      Xopt(1:size(x0)) = reshape(x0,[size(x0)] )

      deg = int(BasisnD%Norm_OF_nDindB)
      !write(out_unitp,*) 'ndim,deg,nqc',BasisnD%ndim,deg,BasisnD%nqc

      nDsize(:) = deg+1
      nDinit(:) = 1
      CALL init_nDindexPrim(nDindB,BasisnD%ndim,nDsize,nDinit=nDinit,   &
                            type_OF_nDindex=0,                          &
                            MaxNorm=BasisnD%Norm_OF_nDindB)
      nDindB%Tab_nDval(:,:) = nDindB%Tab_nDval(:,:) - 1
      IF (FirstTime) CALL Write_nDindex(nDindB)
      FirstTime = .FALSE.

      CALL alloc_NParray(M ,[nDindB%Max_nDI,BasisnD%nqc],'M', name_sub)
      CALL alloc_NParray(R ,[nDindB%Max_nDI],            'R', name_sub)
      CALL alloc_NParray(ER,[nDindB%Max_nDI],            'ER',name_sub)

      R(:)   = ZERO
      M(:,:) = ONE

      DO ib=1,nDindB%Max_nDI
        CALL calc_nDindex(nDindB,ib,nDval)
        IF (count(nDval == 0) == BasisnD%ndim) R(ib) = sqrt(sqrt(PI**BasisnD%ndim))
        DO iq=1,BasisnD%nqc
        DO id=1,BasisnD%ndim
          x = x0(id,iq)
          M(ib,iq) = M(ib,iq)*poly_Hermite(x,nDval(id))*exp(-x*x)
        END DO
        END DO

      END DO
      !CALL Write_VecMat(R,out_unitp,5)
      !CALL Write_VecMat(M,out_unitp,5)


      ER(:) = matmul(M,W)-R(:)
      !CALL Write_VecMat(ER,out_unitp,5)

      Energ = sqrt(dot_product(ER,ER)/real(nDindB%Max_nDI,kind=Rkind))

      !write(out_unitp,*) 'Energ',Energ

      CALL dealloc_NParray(M,'M',name_sub)
      CALL dealloc_NParray(R,'R',name_sub)
      CALL dealloc_NParray(ER,'ER',name_sub)

      CALL dealloc_nDindex(nDindB)
!---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
!---------------------------------------------------------------------

    END SUBROUTINE Sub_Energ_FOR_cubatureWeight
