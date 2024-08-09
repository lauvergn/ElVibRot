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
      MODULE mod_SimulatedAnnealing
      USE EVR_system_m

      IMPLICIT NONE

        TYPE param_SimulatedAnnealing

        integer           :: nb_mc_tot           =  100000
        integer           :: nb_mc_partial       =  100

        integer           :: TempInit_type       =  1
        real (kind=Rkind) :: Tmax                = -ONE
        real (kind=Rkind) :: Tmin                =  ONETENTH**7
        real (kind=Rkind) :: DeltaT              =  ZERO

        real (kind=Rkind) :: RangeScal           =  0.8_Rkind
        real (kind=Rkind) :: RangeScalInit       =  1._Rkind

        logical           :: With_RangeInit      = .FALSE.
        real (kind=Rkind) :: RangeInit           = 1._Rkind

        integer           :: TempScheduling_type =  2 ! 1: linear, 2: geometrical ...
        real (kind=Rkind) :: ExpCoolParam        =  0.95_Rkind

        logical           :: ResetTemp           = .TRUE.
        real (kind=Rkind) :: ResetTempScal       =  ONE/THREE

        integer           :: Restart_Opt         =  0


        END TYPE param_SimulatedAnnealing

      CONTAINS

      SUBROUTINE Read_param_SimulatedAnnealing(para_SimulatedAnnealing)
      TYPE (param_SimulatedAnnealing), intent(inout) :: para_SimulatedAnnealing

        integer :: nb_mc_tot     = 1000
        integer :: nb_mc_partial = 100

        real (kind=Rkind) :: Tmax          = -ONE
        real (kind=Rkind) :: Tmin          =  ONETENTH**7
        real (kind=Rkind) :: DeltaT        =  ZERO
        real (kind=Rkind) :: ResetTempScal = ONE/THREE


        real (kind=Rkind) :: ExpCoolParam = 0.95_Rkind

        real (kind=Rkind) :: RangeScal     = 0.8_Rkind
        real (kind=Rkind) :: RangeScalInit = 1._Rkind
        logical           :: With_RangeInit = .FALSE.
        real (kind=Rkind) :: RangeInit     = 1._Rkind

        logical :: ResetTemp               = .TRUE.
        integer :: TempScheduling_type     = 2 ! 1: linear, 2: geometrical ...
        integer :: TempInit_type           = 1

        integer :: Restart_Opt     = 0

        integer :: err_io
        NAMELIST /SimulatedAnnealing/nb_mc_tot,nb_mc_partial,           &
                                        Tmax,Tmin,DeltaT,TempInit_type, &
                                               RangeScal,RangeScalInit, &
                                              With_RangeInit,RangeInit, &
                                      TempScheduling_type,ExpCoolParam, &
                                    ResetTemp,ResetTempScal,Restart_Opt

        nb_mc_tot           =  1000
        nb_mc_partial       =  100

        TempInit_type       =  1
        Tmax                = -ONE
        Tmin                =  ONETENTH**7
        DeltaT              =  ZERO

        RangeScal           =  0.8_Rkind
        RangeScalInit       =  1._Rkind

        With_RangeInit      = .FALSE.
        RangeInit           = 1._Rkind

        TempScheduling_type =  2 ! 1: linear, 2: geometrical ...
        ExpCoolParam        =  0.95_Rkind

        ResetTemp           = .TRUE.
        ResetTempScal       =  ONE/THREE

        Restart_Opt         =  0

        read(in_unit,SimulatedAnnealing,IOSTAT=err_io)
        IF (err_io /= 0) THEN
           write(out_unit,*) ' WARNING in Read_param_SimulatedAnnealing'
           write(out_unit,*) '  while reading the "SimulatedAnnealing" namelist'
           write(out_unit,*) ' end of file or end of record'
           write(out_unit,*) ' Check your data !!'
           STOP
        END IF
        IF (print_level > 1) write(out_unit,SimulatedAnnealing)

        para_SimulatedAnnealing%nb_mc_tot           =  nb_mc_tot
        para_SimulatedAnnealing%nb_mc_partial       =  nb_mc_partial

        para_SimulatedAnnealing%TempInit_type       =  TempInit_type
        para_SimulatedAnnealing%Tmax                =  Tmax
        para_SimulatedAnnealing%Tmin                =  Tmin
        para_SimulatedAnnealing%DeltaT              =  DeltaT

        para_SimulatedAnnealing%With_RangeInit      =  With_RangeInit
        para_SimulatedAnnealing%RangeInit           =  RangeInit
        IF (With_RangeInit) RangeScalInit           =  ONE

        para_SimulatedAnnealing%RangeScal           =  RangeScal
        para_SimulatedAnnealing%RangeScalInit       =  RangeScalInit



        para_SimulatedAnnealing%TempScheduling_type =  TempScheduling_type
        para_SimulatedAnnealing%ExpCoolParam        =  ExpCoolParam

        para_SimulatedAnnealing%ResetTemp           =  ResetTemp
        para_SimulatedAnnealing%ResetTempScal       =  ResetTempScal

        para_SimulatedAnnealing%Restart_Opt         =  Restart_Opt

      END SUBROUTINE Read_param_SimulatedAnnealing
  SUBROUTINE dealloc_param_SimulatedAnnealing(para_SimulatedAnnealing)
    TYPE (param_SimulatedAnnealing), intent(inout)   :: para_SimulatedAnnealing
  

    para_SimulatedAnnealing%nb_mc_tot           =  100000
    para_SimulatedAnnealing%nb_mc_partial       =  100

    para_SimulatedAnnealing%TempInit_type       =  1
    para_SimulatedAnnealing%Tmax                = -ONE
    para_SimulatedAnnealing%Tmin                =  ONETENTH**7
    para_SimulatedAnnealing%DeltaT              =  ZERO

    para_SimulatedAnnealing%RangeScal           =  0.8_Rkind
    para_SimulatedAnnealing%RangeScalInit       =  1._Rkind

    para_SimulatedAnnealing%With_RangeInit      = .FALSE.
    para_SimulatedAnnealing%RangeInit           = 1._Rkind

    para_SimulatedAnnealing%TempScheduling_type =  2 ! 1: linear, 2: geometrical ...
    para_SimulatedAnnealing%ExpCoolParam        =  0.95_Rkind

    para_SimulatedAnnealing%ResetTemp           = .TRUE.
    para_SimulatedAnnealing%ResetTempScal       =  ONE/THREE

    para_SimulatedAnnealing%Restart_Opt         =  0

  END SUBROUTINE dealloc_param_SimulatedAnnealing
      SUBROUTINE Write_param_SimulatedAnnealing(para_SimulatedAnnealing)
      TYPE (param_SimulatedAnnealing), intent(in)   :: para_SimulatedAnnealing

      write(out_unit,*) '  WRITE param_SimulatedAnnealing'
      write(out_unit,*)
      write(out_unit,*) '  nb_mc_tot          ',para_SimulatedAnnealing%nb_mc_tot
      write(out_unit,*) '  nb_mc_partial      ',para_SimulatedAnnealing%nb_mc_partial
      write(out_unit,*)
      write(out_unit,*) '  TempInit_type      ',para_SimulatedAnnealing%TempInit_type
      write(out_unit,*) '  Tmax               ',para_SimulatedAnnealing%Tmax
      write(out_unit,*) '  Tmin               ',para_SimulatedAnnealing%Tmin
      write(out_unit,*) '  DeltaT             ',para_SimulatedAnnealing%DeltaT
      write(out_unit,*)
      write(out_unit,*) '  RangeScal          ',para_SimulatedAnnealing%RangeScal
      write(out_unit,*) '  RangeScalInit      ',para_SimulatedAnnealing%RangeScalInit

      write(out_unit,*) '  With_RangeInit     ',para_SimulatedAnnealing%With_RangeInit
      write(out_unit,*) '  RangeInit          ',para_SimulatedAnnealing%RangeInit

      write(out_unit,*)
      write(out_unit,*) '  TempScheduling_type',para_SimulatedAnnealing%TempScheduling_type
      write(out_unit,*) '  ExpCoolParam       ',para_SimulatedAnnealing%ExpCoolParam
      write(out_unit,*)
      write(out_unit,*) '  ResetTemp          ',para_SimulatedAnnealing%ResetTemp
      write(out_unit,*) '  ResetTempScal      ',para_SimulatedAnnealing%ResetTempScal
      write(out_unit,*)
      write(out_unit,*) '  Restart_Opt        ',para_SimulatedAnnealing%Restart_Opt
      write(out_unit,*)
      write(out_unit,*) '  END WRITE param_SimulatedAnnealing'

      END SUBROUTINE write_param_SimulatedAnnealing

      SUBROUTINE Sub_SimulatedAnnealing(BasisnD,xOpt_min,Norm_min,      &
                                        SQ,nb_Opt,                      &
                                        para_Tnum,mole,PrimOp,Qact,     &
                                        para_SimulatedAnnealing)

      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
      logical        :: Cart_Transfo_save
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: BasisnD

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

!----- for the optimization -------------------------------------------
      TYPE (param_SimulatedAnnealing) :: para_SimulatedAnnealing
      integer, intent(in) :: nb_Opt
      real (kind=Rkind), intent(inout) :: xOpt_min(nb_Opt),SQ(nb_Opt)

!---------- working variables ----------------------------------------
      real (kind=Rkind) :: Norm_min,Norm_max

      real (kind=Rkind) :: DNorm,NormA,NormB
      real (kind=Rkind) :: Temp,Temp_max,PTemp,DTemp


      real (kind=Rkind) :: xi,xav
      real (kind=Rkind), pointer :: x(:)
      real (kind=Rkind), pointer :: x0(:)
      real (kind=Rkind), pointer :: QA(:)
      real (kind=Rkind), pointer :: QB(:)
      real (kind=Rkind), pointer :: QAv(:)

      real (kind=Rkind), pointer :: xOpt(:)

      integer           :: i,imc,nb_Norm_min,nb_block_WithoutMin
      integer           :: nb_Accepted_DNorm,nb_Accepted_proba,nb_NotAccepted_proba


!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Sub_SimulatedAnnealing'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------

      ! allocation x ...
      nullify(x)
      nullify(x0)
      nullify(QA)
      nullify(QB)
      nullify(QAv)
      nullify(xOpt)

      CALL alloc_array(x,   [nb_Opt],'x',   name_sub)
      CALL alloc_array(x0,  [nb_Opt],'x0',  name_sub)
      CALL alloc_array(xOpt,[nb_Opt],'xOpt',name_sub)
      CALL alloc_array(QA,  [nb_Opt],'QA',  name_sub)
      CALL alloc_array(QB,  [nb_Opt],'QB',  name_sub)
      CALL alloc_array(QAv, [nb_Opt],'Qav', name_sub)
      ! end allocation

      CALL Write_param_SimulatedAnnealing(para_SimulatedAnnealing)

      CALL set_print_level(0,force=.TRUE.)
      CALL Sub_Energ_OF_ParamBasis(Norm_min,xOpt_min,nb_Opt,BasisnD,    &
                                   para_Tnum,mole,PrimOp,Qact)
      IF (size(xOpt_min) > 20) THEN
        write(out_unit,*) ' Initial param',xOpt_min(1:10),' ... Energy',Norm_min
      ELSE
        write(out_unit,*) ' Initial param',xOpt_min,' Energy',Norm_min
      END IF


      QA(:) = xOpt_min(:) - SQ(:)*HALF
      QB(:) = xOpt_min(:) + SQ(:)*HALF
      QAv(:)= xOpt_min(:)
      IF (size(QA) > 20) THEN
        write(out_unit,*) 'SQ(1:10)',SQ(1:10),' ...'
        write(out_unit,*) 'QA(1:10)',QA(1:10),' ...'
        write(out_unit,*) 'QB(1:10)',QB(1:10),' ...'
      ELSE
        write(out_unit,*) 'SQ',SQ
        write(out_unit,*) 'QA',QA
        write(out_unit,*) 'QB',QB
      END IF


      Norm_max    = Norm_min
      NormA       = ZERO

      CALL set_print_level(-1,force=.TRUE.)
      ! first find the average Energy (Norm), then Temp
      DO imc=1,para_SimulatedAnnealing%nb_mc_tot/10
          CALL random_number(x)
          xOpt(:) = QA(:) + x(:)*SQ(:)

          CALL Sub_Energ_OF_ParamBasis(NormB,xOpt,nb_Opt,BasisnD,       &
                                       para_Tnum,mole,PrimOp,Qact)

          !write(out_unit,*) 'NormB',imc,NormB
          NormA = NormA + NormB
          IF (NormB > Norm_max) Norm_max = NormB

          IF (NormB < Norm_min) THEN
            xOpt_min(:) = xOpt(:)
            Norm_min    = NormB
            write(out_unit,*) 'NormB',imc,NormB
          END IF
      END DO
      NormA = NormA / real(para_SimulatedAnnealing%nb_mc_tot/10,kind=Rkind)
      write(out_unit,*) 'Min, Average, Max Energy',Norm_min,NormA,Norm_max

      x0(:)  = xOpt_min(:)
      IF (size(xOpt_min) > 20) THEN
        write(out_unit,*) ' Energy_min',xOpt_min(1:10),' ...',Norm_min
      ELSE
        write(out_unit,*) ' Energy_min',xOpt_min,Norm_min
      END IF

      Temp_max = Norm_max-Norm_min
      Temp     = Temp_max
      write(out_unit,*) 'Average Energy, Temp',NormA,Temp
      flush(out_unit)

      DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot),kind=Rkind)

      imc                 = 1
      nb_Norm_min         = 0
      nb_block_WithoutMin = 0
      DO

        DO i=1,nb_Opt
          DO
            CALL random_number(xi)
            xi = TWO*xi-ONE
            xOpt(i) = x0(i) + xi*SQ(i)
            IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
          END DO
        END DO

        CALL Sub_Energ_OF_ParamBasis(NormA,xOpt,nb_Opt,BasisnD,         &
                                     para_Tnum,mole,PrimOp,Qact)

        DNorm = NormA - NormB
        !write(out_unit,*) 'Norm',imc,xOpt,NormA
        !write(out_unit,*) 'Norm',imc,sqrt(sum((xOpt-xOpt_min)**2)),NormA

        flush(out_unit)

        IF ( NormA < Norm_min) THEN
          nb_Norm_min   = nb_Norm_min + 1
          xOpt_min(:)   = xOpt(:)
          Norm_min      = NormA
          !write(out_unit,*) ' imc, Temp, nb_Energy_min, Energy_min',imc,Temp,nb_Norm_min,Norm_min
        END IF


        IF (DNorm < ZERO) THEN
          x0(:)       = xOpt(:)
          NormB       = NormA
        ELSE
          CALL random_number(PTemp)
          IF ( PTemp < exp(-DNorm/Temp) ) THEN
            x0(:)       = xOpt(:)
            NormB       = NormA
            !write(out_unit,*) ' imc, Temp, Norm (propa)',imc,Temp,NormA
            !write(out_unit,*) ' accepted, Proba',PTemp,exp(-DNorm/Temp)
          END IF
        END IF
        imc = imc + 1
        SELECT CASE (para_SimulatedAnnealing%TempScheduling_type)
        CASE (1)
          !write(out_unit,*) 'LinCoolParam'
          Temp = Temp - DTemp
        CASE (2)
          !write(out_unit,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        CASE DEFAULT
          !write(out_unit,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        END SELECT

        IF (para_SimulatedAnnealing%ResetTemp .AND.                     &
             Temp < para_SimulatedAnnealing%ResetTempScal*Temp_max) THEN

          write(out_unit,*) 'Temp_max,nb_Energy_min,Energy_min',           &
                                           Temp_max,nb_Norm_min,Norm_min
          flush(out_unit)
          SQ(:)       = SQ(:)*para_SimulatedAnnealing%RangeScal
          Temp_max    = (ONE-para_SimulatedAnnealing%ResetTempScal)*Temp_max
          Temp        = Temp_max
          DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot-imc),kind=Rkind)
          x0(:)       = xOpt_min(:)

          IF (nb_Norm_min == 0) THEN
            nb_block_WithoutMin = nb_block_WithoutMin + 1
          ELSE
            nb_block_WithoutMin = 0
            nb_Norm_min         = 0
          END IF
        END IF

        IF (nb_block_WithoutMin > 20) THEN
          write(out_unit,*) 'restart because, nb_block_WithoutMin > 20'
          flush(out_unit)
          EXIT
        END IF

        IF (Temp < para_SimulatedAnnealing%Tmin .OR. imc > para_SimulatedAnnealing%nb_mc_tot) EXIT
      END DO
      IF (size(xOpt_min) > 20) THEN
        write(out_unit,*) ' Optimal param',xOpt_min(1:10),' ... Energy',Norm_min
      ELSE
        write(out_unit,*) ' Optimal param',xOpt_min,' Energy',Norm_min
      END IF

      !deallocation
      CALL dealloc_array(x,   'x',   name_sub)
      CALL dealloc_array(x0,  'x0',  name_sub)
      CALL dealloc_array(xOpt,'xOpt',name_sub)
      CALL dealloc_array(QA,  'QA',  name_sub)
      CALL dealloc_array(QB,  'QB',  name_sub)
      CALL dealloc_array(QAv, 'Qav', name_sub)
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Sub_SimulatedAnnealing

      SUBROUTINE Sub_SimulatedAnnealing_cuba(BasisnD,xOpt_min,Norm_min, &
                                             SQ,QA,QB,nb_Opt,           &
                                            para_Tnum,mole,PrimOp,Qact, &
                                             para_SimulatedAnnealing)

      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
      logical        :: Cart_Transfo_save
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: BasisnD

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

!----- for the optimization -------------------------------------------
      TYPE (param_SimulatedAnnealing)  :: para_SimulatedAnnealing
      integer, intent(in)              :: nb_Opt
      real (kind=Rkind), intent(inout) :: xOpt_min(nb_Opt),SQ(nb_Opt),  &
                                          QA(nb_Opt),QB(nb_Opt)
!---------- working variables ----------------------------------------
      real (kind=Rkind) :: Norm_min,Norm_max

      real (kind=Rkind) :: DNorm,NormA,NormB
      real (kind=Rkind) :: Temp,Temp_max,PTemp,DTemp


      real (kind=Rkind) :: xi,xav
      real (kind=Rkind), pointer :: x(:)
      real (kind=Rkind), pointer :: x0(:)
      real (kind=Rkind), pointer :: QAv(:)

      real (kind=Rkind), pointer :: xOpt(:)

      integer           :: i,k,imc,imc_block,nb_Norm_min,nb_block_WithoutMin
      integer           :: nb_Accepted_DNorm,nb_Accepted_proba,nb_NotAccepted_proba

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Sub_SimulatedAnnealing_cuba'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
      END IF
!---------------------------------------------------------------------

      ! allocation x ...
      nullify(x)
      nullify(x0)
      nullify(QAv)
      nullify(xOpt)

      CALL alloc_array(x,   [nb_Opt],'x',   name_sub)
      CALL alloc_array(x0,  [nb_Opt],'x0',  name_sub)
      CALL alloc_array(xOpt,[nb_Opt],'xOpt',name_sub)
      CALL alloc_array(QAv, [nb_Opt],'Qav', name_sub)
      ! end allocation

      CALL Write_param_SimulatedAnnealing(para_SimulatedAnnealing)

      CALL set_print_level(0,force=.TRUE.)
      QAv(:)= xOpt_min(:)
      x0(:) = xOpt_min(:)
      DO i=1,nb_Opt
        DO k=1,1000
          CALL random_number(xi)
          xi = TWO*xi-ONE
          xOpt(i) = x0(i) + xi*SQ(i)/TWO
          IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
        END DO
      END DO
      CALL Sub_Energ_OF_ParamBasis(Norm_min,xOpt,nb_Opt,BasisnD,        &
                                   para_Tnum,mole,PrimOp,Qact)

      write(out_unit,*) 'Initial param',xOpt,' Energy',Norm_min


      write(out_unit,*) 'SQ',SQ
      write(out_unit,*) 'QA',QA
      write(out_unit,*) 'QB',QB


      Norm_max = Norm_min
      NormA    = ZERO
      CALL set_print_level(-1,force=.TRUE.)
      ! first find the average Energy (Norm), then Temp
      DO imc=1,para_SimulatedAnnealing%nb_mc_tot/10
        DO i=1,nb_Opt
          DO k=1,1000
            CALL random_number(xi)
            xi = TWO*xi-ONE
            xOpt(i) = x0(i) + xi*SQ(i)
            IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
          END DO
        END DO

        CALL Sub_Energ_OF_ParamBasis(NormB,xOpt,nb_Opt,BasisnD,         &
                                     para_Tnum,mole,PrimOp,Qact)

        !write(out_unit,*) 'NormB',imc,NormB
        NormA = NormA + NormB
        IF (NormB > Norm_max) Norm_max = NormB

        IF (NormB < Norm_min) THEN
          xOpt_min(:) = xOpt(:)
          Norm_min = NormB
        END IF
      END DO
      NormA = NormA / real(para_SimulatedAnnealing%nb_mc_tot/10,kind=Rkind)
      write(out_unit,*) 'Min, Average, Max Energy',Norm_min,NormA,Norm_max

      x0(:)   = xOpt_min(:)
      Xopt(:) = xOpt_min(:)
      CALL Sub_Energ_OF_ParamBasis(NormB,xOpt,nb_Opt,BasisnD,         &
                                   para_Tnum,mole,PrimOp,Qact)
      Norm_min = NormB
      write(out_unit,*) ' Energy_min',xOpt_min,Norm_min

      Temp_max = Norm_max-Norm_min
      Temp_max = NormA
      Temp     = Temp_max
      write(out_unit,*) 'Average Energy, Temp',NormA,Temp
      flush(out_unit)

      DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot),kind=Rkind)

      nb_Accepted_DNorm    = 0
      nb_Accepted_proba    = 0
      nb_NotAccepted_proba = 0
      imc                  = 0
      nb_Norm_min          = 0
      nb_block_WithoutMin  = 0
      imc_block            = 0
      !SQ(:) = SQ(:)/TEN
      write(out_unit,*) 'NormB,SQ(1)',NormB,SQ(1)
      DO
        imc       = imc       + 1
        imc_block = imc_block + 1

        DO i=1,nb_Opt
          DO k=1,1000
            CALL random_number(xi)
            xi = TWO*xi-ONE
            xOpt(i) = x0(i) + xi*SQ(i)
            IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
          END DO
        END DO


        CALL Sub_Energ_OF_ParamBasis(NormA,xOpt,nb_Opt,BasisnD,         &
                                     para_Tnum,mole,PrimOp,Qact)

        !write(out_unit,*) 'imc,Xopt,NormA,SQ',imc,xOpt,NormA,SQ(1)

        DNorm = NormA - NormB
        !write(out_unit,*) 'Norm',imc,xOpt,NormA
        flush(out_unit)

        IF ( NormA < Norm_min) THEN
          nb_Norm_min   = nb_Norm_min + 1
          xOpt_min(:)   = xOpt(:)
          Norm_min      = NormA
          !write(out_unit,*) ' imc, Temp, nb_Energy_min, Energy_min',imc,Temp,nb_Norm_min,Norm_min
        END IF

        IF (DNorm < ZERO) THEN
          x0(:)       = xOpt(:)
          NormB       = NormA
          !write(out_unit,*) ' accepted, DNorm',imc,DNorm
          nb_Accepted_DNorm = nb_Accepted_DNorm + 1
        ELSE
          CALL random_number(PTemp)
          IF ( PTemp < exp(-DNorm/Temp) ) THEN
            x0(:)       = xOpt(:)
            NormB       = NormA
            !write(out_unit,*) ' imc, Temp, Norm (propa)',imc,Temp,NormA
            !write(out_unit,*) ' accepted, Proba',imc,PTemp,exp(-DNorm/Temp)
            nb_Accepted_proba = nb_Accepted_proba + 1
          ELSE
            nb_NotAccepted_proba = nb_NotAccepted_proba + 1
            !write(out_unit,*) ' Not accepted, Proba',imc,PTemp,exp(-DNorm/Temp)
          END IF
        END IF


        IF (imc_block > para_SimulatedAnnealing%nb_mc_partial) THEN
        SELECT CASE (para_SimulatedAnnealing%TempScheduling_type)
        CASE (1)
          !write(out_unit,*) 'LinCoolParam'
          Temp = Temp - DTemp
        CASE (2)
          !write(out_unit,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        CASE DEFAULT
          !write(out_unit,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        END SELECT
        END IF

        !write(out_unit,*) 'imc,Temp',imc,Temp,Temp_max,para_SimulatedAnnealing%ResetTempScal
        !write(out_unit,*) 'imc,Temp',imc,Temp,Temp_max*para_SimulatedAnnealing%ResetTempScal

        !write(out_unit,*) 'test ?',(Temp < para_SimulatedAnnealing%ResetTempScal*Temp_max)


        IF (para_SimulatedAnnealing%ResetTemp .AND.                     &
             Temp < para_SimulatedAnnealing%ResetTempScal*Temp_max) THEN

          write(out_unit,*) 'imc,Temp_max,nb_Energy_min,Energy_min',   &
                                       imc,Temp_max,nb_Norm_min,Norm_min
          write(out_unit,*) 'imc,Temp_max,nb_Accepted_...',imc,Temp_max,&
                nb_Accepted_DNorm,nb_Accepted_proba,nb_NotAccepted_proba
          flush(out_unit)
          SQ(:)       = SQ(:)*para_SimulatedAnnealing%RangeScal
          Temp_max    = (ONE-para_SimulatedAnnealing%ResetTempScal)*Temp_max
          DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot-imc),kind=Rkind)
          Temp        = Temp_max
          x0(:)       = xOpt_min(:)
          IF (nb_Norm_min == 0) THEN
            nb_block_WithoutMin = nb_block_WithoutMin + 1
          ELSE
            nb_block_WithoutMin = 0
            nb_Norm_min         = 0
          END IF

          nb_Accepted_DNorm    = 0
          nb_Accepted_proba    = 0
          nb_NotAccepted_proba = 0
          imc_block            = 0

        END IF

        IF (Temp < para_SimulatedAnnealing%Tmin .OR. imc > para_SimulatedAnnealing%nb_mc_tot) EXIT
      END DO
      write(out_unit,*) 'Optimal param',xOpt_min,' Energy',Norm_min

      !deallocation
      CALL dealloc_array(x,   'x',   name_sub)
      CALL dealloc_array(x0,  'x0',  name_sub)
      CALL dealloc_array(xOpt,'xOpt',name_sub)
      CALL dealloc_array(QAv, 'Qav', name_sub)
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Sub_SimulatedAnnealing_cuba


      END MODULE mod_SimulatedAnnealing

