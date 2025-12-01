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
MODULE mod_Optimization
  USE EVR_system_m
  USE mod_SimulatedAnnealing
  USE mod_BFGS
  IMPLICIT NONE


  TYPE param_Optimization

        character (len=Name_len) :: Optimization_method = 'SimulatedAnnealing'
        character (len=Name_len) :: Optimization_param  = 'coordbasis'
        logical                  :: FinalEnergy         = .TRUE.
        logical                  :: Grad                = .FALSE.
        logical                  :: Freq                = .FALSE.
        logical                  :: ReadRange           = .FALSE.

        TYPE (param_SimulatedAnnealing) :: para_SimulatedAnnealing
        TYPE (param_BFGS)               :: para_BFGS

        integer                         :: nb_Opt = 0
        real (kind=Rkind), allocatable  :: xOpt_min(:)
        real (kind=Rkind), allocatable  :: SQ(:)
        real (kind=Rkind), allocatable  :: QA(:),QB(:)

  END TYPE param_Optimization

CONTAINS

  SUBROUTINE Read_param_Optimization(para_Optimization,mole,BasisnD,read_nml)
    USE mod_Coord_KEO, only : CoordType, get_Qact0
    USE mod_basis
    USE BasisMakeGrid
    IMPLICIT NONE

    TYPE (param_Optimization), intent(inout) :: para_Optimization
    TYPE (CoordType),          intent(in)    :: mole
    TYPE (basis),              intent(in)    :: BasisnD
    logical,                   intent(in)    :: read_nml


    real (kind=Rkind), allocatable :: Qact(:)
    integer :: i
    character (len=Name_len)       :: name_dum

    character (len=Name_len) :: Optimization_method = 'SimulatedAnnealing'
    character (len=Name_len) :: Optimization_param  = 'coordbasis'
    logical                  :: FinalEnergy         = .TRUE.
    logical                  :: Freq                = .FALSE.
    logical                  :: Grad                = .FALSE.
    logical                  :: ReadRange           = .FALSE.

    integer :: err_io
    NAMELIST /Optimization/ Optimization_method,Optimization_param, &
                            FinalEnergy,Freq,Grad,ReadRange

    character (len=*), parameter :: name_sub = 'Read_param_Optimization'

    CALL alloc_NParray(Qact,[mole%nb_var],'Qact',name_sub)
    CALL get_Qact0(Qact,mole%tab_Qtransfo(mole%itActive)%ActiveTransfo)
    CALL Set_ALL_para_FOR_optimization(mole,BasisnD,Qact,0)

    para_Optimization%nb_Opt = para_FOR_optimization%nb_OptParam
    IF (para_Optimization%nb_Opt < 1) THEN
      write(out_unit,*) 'ERROR in sub_Optimization_OF_VibParam'
      write(out_unit,*) 'Wrong number  nb_OptParam',para_Optimization%nb_Opt
      write(out_unit,*) ' Check your data!!'
      STOP
    END IF
    IF (print_level > 1) write(out_unit,*) 'para_Optimization%nb_Opt',para_Optimization%nb_Opt

    IF (read_nml) THEN
      Optimization_method = 'SimulatedAnnealing'
      Optimization_param  = 'coordbasis'
      FinalEnergy         = .TRUE.  ! redo the optimal energy
      Freq                = .FALSE. ! Calculation of the frequencies (geometry optimization)
      Grad                = .FALSE. ! Calculation of the gradiant
      ReadRange           = .FALSE.

      read(in_unit,Optimization,IOSTAT=err_io)
      IF (err_io /= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the "Optimization" namelist'
        write(out_unit,*) ' end of file or end of record'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Read_param_Optimization while reading the "Optimization" namelist'
      END IF
      IF (print_level > 1) write(out_unit,Optimization)

      para_Optimization%Optimization_method    = TO_lowercase(Optimization_method)
      para_Optimization%Optimization_param     = TO_lowercase(Optimization_param)

      para_FOR_optimization%Optimization_param = Optimization_param
      IF (Optimization_param /= 'geometry') Freq = .FALSE.

      para_Optimization%FinalEnergy            = FinalEnergy
      para_Optimization%Freq                   = Freq
      para_Optimization%Grad                   = Grad
      para_Optimization%ReadRange              = ReadRange

      SELECT CASE (para_Optimization%Optimization_method)
      CASE ('simulatedannealing','sa')
        !CALL Read_param_SimulatedAnnealing(para_Optimization%para_SimulatedAnnealing)
STOP   'Compilation pb in Read_param_Optimization: !!! ?????'
      CASE ('bfgs')
        CALL Read_param_BFGS(para_Optimization%para_BFGS,para_Optimization%nb_Opt)

      CASE DEFAULT
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'Wrong Optimization_method: ',para_Optimization%Optimization_method
        write(out_unit,*) ' Check your data!!'
        STOP 'ERROR in Read_param_Optimization: Wrong Optimization_method'
      END SELECT
    END IF

    IF (allocated(para_FOR_optimization%Val_RVec)) &
       CALL dealloc_NParray(para_FOR_optimization%Val_RVec,'para_FOR_optimization%Val_RVec',name_sub)

    IF (allocated(para_FOR_optimization%opt_RVec)) &
       CALL dealloc_NParray(para_FOR_optimization%opt_RVec,'para_FOR_optimization%opt_RVec',name_sub)

    IF (para_Optimization%nb_Opt > 0) THEN
      !write(out_unit,*) 'para_FOR_optimization%nb_Opt',para_Optimization%nb_Opt
      CALL alloc_NParray(para_FOR_optimization%Val_RVec,[para_Optimization%nb_Opt],     &
                        'para_FOR_optimization%Val_RVec',name_sub)
      para_FOR_optimization%Val_RVec(:) = ZERO
      CALL alloc_NParray(para_FOR_optimization%opt_RVec,[para_Optimization%nb_Opt],     &
                        'para_FOR_optimization%opt_RVec',name_sub)
      para_FOR_optimization%opt_RVec(:) = 1
    END IF
    CALL Set_ALL_para_FOR_optimization(mole,BasisnD,Qact,-1)

    ! transfert xOpt to para_FOR_optimization%tab_OF_RVec(1)%P_RVec
    CALL alloc_NParray(para_Optimization%xOpt_min,[para_Optimization%nb_Opt],'xOpt_min',name_sub)
    CALL alloc_NParray(para_Optimization%SQ,      [para_Optimization%nb_Opt],'SQ',      name_sub)
    CALL alloc_NParray(para_Optimization%QA,      [para_Optimization%nb_Opt],'QA',      name_sub)
    CALL alloc_NParray(para_Optimization%QB,      [para_Optimization%nb_Opt],'QB',      name_sub)

    para_Optimization%SQ(:) = ZERO
    IF (para_FOR_optimization%Optimization_param /= 'cubature') THEN
      para_Optimization%xOpt_min(:) = para_FOR_optimization%Val_RVec(:)
      IF (para_Optimization%para_SimulatedAnnealing%With_RangeInit) THEN
        para_Optimization%SQ(:) = para_Optimization%para_SimulatedAnnealing%RangeInit
      ELSE
        para_Optimization%SQ(:)       = ONETENTH * real(para_FOR_optimization%opt_RVec(:),kind=Rkind)
      END IF
    ELSE
      IF (para_Optimization%para_SimulatedAnnealing%With_RangeInit) THEN
        para_Optimization%SQ(:) = para_Optimization%para_SimulatedAnnealing%RangeInit
        para_Optimization%xOpt_min(:) = ZERO
      ELSE
        para_Optimization%SQ(:) = (maxval(BasisnD%x)-minval(BasisnD%x))
        para_Optimization%xOpt_min(:) = ZERO
      END IF
      CALL ReCentered_grid(reshape(para_Optimization%xOpt_min,                          &
          [BasisnD%ndim,BasisnD%nqc]),&
                  BasisnD%ndim,BasisnD%nqc)
      CALL ReOriented_grid(reshape(para_Optimization%xOpt_min,                          &
          [BasisnD%ndim,BasisnD%nqc]),&
                  BasisnD%ndim,BasisnD%nqc)
      IF (para_Optimization%ReadRange) THEN
         DO i=1,size(para_Optimization%QA)
           read(in_unit,*,IOSTAT=err_io) name_dum,para_Optimization%QA(i),para_Optimization%QB(i)
           IF (err_io /= 0) THEN
             write(out_unit,*) ' WARNING in name_sub'
             write(out_unit,*) '  while reading the variable range'
             write(out_unit,*) ' Check your data !!'
             STOP
           END IF
         END DO
         para_Optimization%SQ(:) = para_Optimization%QA-para_Optimization%QB
      ELSE
        para_Optimization%QA(:) = para_Optimization%xOpt_min(:) - para_Optimization%SQ(:)
        para_Optimization%QB(:) = para_Optimization%xOpt_min(:) + para_Optimization%SQ(:)
      END IF

      !BasisnD_Save%nb  = BasisnD%nb
      !BasisnD_Save%nbc = BasisnD%nbc
      !CALL Set_nq_OF_basis(BasisnD,BasisnD%nqc)
      !BasisnD_Save%nqc = BasisnD%nqc

    END IF

  END SUBROUTINE Read_param_Optimization
  SUBROUTINE Write_param_Optimization(para_Optimization)
    IMPLICIT NONE

    TYPE (param_Optimization), intent(in) :: para_Optimization

    write(out_unit,*) 'WRITE param_Optimization'
    write(out_unit,*)
    write(out_unit,*) 'Optimization_method        ',para_Optimization%Optimization_method
    write(out_unit,*) 'Optimization_param         ',para_Optimization%Optimization_param
    write(out_unit,*) 'ReadRange                  ',para_Optimization%ReadRange
    write(out_unit,*) 'FinalEnergy                ',para_Optimization%FinalEnergy
    write(out_unit,*) 'Grad                       ',para_Optimization%Grad
    write(out_unit,*) 'Freq                       ',para_Optimization%Freq
    write(out_unit,*) 'nb_Opt                     ',para_Optimization%nb_Opt

    IF (allocated(para_Optimization%xOpt_min)) &
    write(out_unit,*) 'xOpt_min(:)                 ',para_Optimization%xOpt_min

    IF (allocated(para_Optimization%SQ)) &
    write(out_unit,*) 'SQ(:)                       ',para_Optimization%SQ

    IF (allocated(para_Optimization%QA)) &
    write(out_unit,*) 'QA(:)                       ',para_Optimization%QA

    IF (allocated(para_Optimization%QB)) &
    write(out_unit,*) 'QB(:)                       ',para_Optimization%QB

    SELECT CASE (para_Optimization%Optimization_method)
    CASE ('simulatedannealing','sa')
      CALL Write_param_SimulatedAnnealing(para_Optimization%para_SimulatedAnnealing)
    CASE ('bfgs')
      CALL Write_param_BFGS(para_Optimization%para_BFGS)
    CASE DEFAULT
      write(out_unit,*) 'WARNING in Write_param_param_Optimization'
      write(out_unit,*) 'No parameter for this Optimization_method',para_Optimization%Optimization_method
    END SELECT
    write(out_unit,*) 'END WRITE param_Optimization'
    flush(out_unit)
  END SUBROUTINE Write_param_Optimization

  SUBROUTINE dealloc_param_Optimization(para_Optimization)
    IMPLICIT NONE

    TYPE (param_Optimization), intent(inout) :: para_Optimization

    character (len=*), parameter :: name_sub='dealloc_param_Optimization'

    para_Optimization%Optimization_method = 'SimulatedAnnealing'
    para_Optimization%Optimization_param  = 'coordbasis'
    para_Optimization%FinalEnergy         = .TRUE.
    para_Optimization%Grad                = .FALSE.
    para_Optimization%Freq                = .FALSE.
    para_Optimization%ReadRange           = .FALSE.

    CALL dealloc_param_SimulatedAnnealing(para_Optimization%para_SimulatedAnnealing)
    CALL dealloc_param_BFGS(para_Optimization%para_BFGS)

    para_Optimization%nb_Opt = 0
    CALL dealloc_NParray(para_Optimization%xOpt_min,'xOpt_min',name_sub)
    CALL dealloc_NParray(para_Optimization%SQ,      'SQ',      name_sub)
    CALL dealloc_NParray(para_Optimization%QA,      'QA',      name_sub)
    CALL dealloc_NParray(para_Optimization%QB,      'QB',      name_sub)

  END SUBROUTINE dealloc_param_Optimization
  SUBROUTINE Sub_Optimization(BasisnD,para_Tnum,mole,PrimOp,Qopt,para_Optimization)
    USE EVR_system_m
    USE mod_dnSVM
    use mod_Coord_KEO, only: CoordType, Tnum, alloc_array, dealloc_array
    USE mod_PrimOp
    USE mod_basis
    USE mod_Op
    USE mod_Auto_Basis
    IMPLICIT NONE

    !----- for the optimization -------------------------------------------
    TYPE (param_Optimization),  intent(inout)    :: para_Optimization
    real (kind=Rkind),          intent(inout)    :: Qopt(:)

    !----- for the CoordType and Tnum --------------------------------------
    TYPE (CoordType) :: mole
    TYPE (Tnum)      :: para_Tnum
    
    !----- for the basis set ----------------------------------------------
    TYPE (basis)          :: BasisnD
    
    !----- variables pour la namelist minimum ----------------------------
    TYPE (PrimOp_t)  :: PrimOp

    !----- local variables ----------------------------
    real (kind=Rkind)           :: SQini(para_Optimization%nb_Opt),Norm_min
    integer :: i


    SQini(:) = para_Optimization%SQ(:)
    !write(out_unit,*) 'SQini(1:10)',SQini(1:min(10,size(SQini)))

    SELECT CASE (para_Optimization%Optimization_method)
    CASE ('simulatedannealing','sa') ! simulated annealing
      DO i=0,para_Optimization%para_SimulatedAnnealing%Restart_Opt
        para_Optimization%SQ(:) = SQini(:) * para_Optimization%para_SimulatedAnnealing%RangeScalInit
        !write(out_unit,*) 'SQini(1:10)',SQini(1:min(10,size(SQ)))
        !write(out_unit,*) 'SQ(1:10)',SQ(1:min(10,size(SQ)))

        IF (para_FOR_optimization%Optimization_param == 'cubature') THEN
          CALL Sub_SimulatedAnnealing_cuba(BasisnD,para_Optimization%xOpt_min,           &
                                           Norm_min,para_Optimization%SQ,para_Optimization%QA,&
                                           para_Optimization%QB,para_Optimization%nb_Opt,   &
                                           para_Tnum,mole,PrimOp,Qopt, &
                            para_Optimization%para_SimulatedAnnealing)
          SQini(:) = Norm_min*TEN**1
        ELSE
          CALL Sub_SimulatedAnnealing(BasisnD,para_Optimization%xOpt_min,Norm_min, &
          para_Optimization%SQ,para_Optimization%nb_Opt,para_Tnum,mole,PrimOp,Qopt,&
                            para_Optimization%para_SimulatedAnnealing)
          SQini(:) = SQini(:) * para_Optimization%para_SimulatedAnnealing%RangeScal
        END IF

      END DO
      write(out_unit,*) 'Norm_min',i,Norm_min

    CASE ('bfgs')
      IF (para_Optimization%para_BFGS%calc_hessian_always) THEN
        CALL Sub_Newton(BasisnD,para_Optimization%xOpt_min,para_Optimization%SQ, &
                        para_Optimization%nb_Opt,para_Tnum,mole,   &
                        PrimOp,Qopt,                 &
                        para_Optimization%para_BFGS)
      ELSE
        CALL Sub_BFGS(BasisnD,para_Optimization%xOpt_min,para_Optimization%SQ, &
                      para_Optimization%nb_Opt,para_Tnum,mole,   &
                      PrimOp,Qopt,                 &
                      para_Optimization%para_BFGS)
      END IF

    CASE DEFAULT
      write(out_unit,*) 'ERROR in sub_Optimization_OF_VibParam'
      write(out_unit,*) 'Wrong Optimization_method: ',                 &
                                 para_Optimization%Optimization_method
      write(out_unit,*) ' Check your data!!'
      STOP
    END SELECT

  END SUBROUTINE Sub_Optimization
END MODULE mod_Optimization

