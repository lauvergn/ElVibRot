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
    PROGRAM ElVibRot
      USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT
      USE EVR_system_m
!$    USE omp_lib, only : omp_get_max_threads
      USE mod_nDGridFit
      IMPLICIT NONE

      INTERFACE
        SUBROUTINE read_arg(input_filename)
          character(len=:), allocatable, intent(inout) :: input_filename
        END SUBROUTINE read_arg
      END INTERFACE

      logical  :: intensity_only,analysis_only,Grid_only
      logical  :: Popenmp,Popenmpi

      integer  :: PMatOp_omp,POpPsi_omp,PBasisTOGrid_omp,PGrid_omp,optimization
      integer  :: maxth,PMatOp_maxth,POpPsi_maxth,PBasisTOGrid_maxth,PGrid_maxth
      integer  :: PCRP_omp,PCRP_maxth,PAna_omp,PAna_maxth
      logical  :: PTune_SG4_omp,PTune_Grid_omp
      integer  :: PSG4_omp,PSG4_maxth
      integer (kind=ILkind)  :: max_mem
      integer  :: printlevel,err
      logical  :: test,EVR,cart,nDfit,nDGrid,Opt_CAP_Basis,mem_debug
      logical  :: GridTOBasis_test,OpPsi_test,main_test
      character (len=Name_longlen)   :: EneFormat
      character (len=Name_longlen)   :: RMatFormat
      character (len=Name_longlen)   :: CMatFormat
      character (len=Line_len)       :: base_FileName = ''
      character (len=Line_len)       :: File_path = ''

      character(len=:), allocatable  :: input_filename

      ! parameters for system setup
      namelist /system/ max_mem,mem_debug,test,printlevel,              &

                          Popenmp,Popenmpi,                             &
                          PSG4_omp,PSG4_maxth,                          &
                          PMatOp_omp,PMatOp_maxth,                      &
                          POpPsi_omp,POpPsi_maxth,                      &
                          PBasisTOGrid_omp,PBasisTOGrid_maxth,          &
                          PGrid_omp,PGrid_maxth,                        &
                          PCRP_omp,PCRP_maxth,                          &
                          PAna_omp,PAna_maxth,                          &
                          PTune_SG4_omp,PTune_Grid_omp,                 &

                          RMatFormat,CMatFormat,EneFormat,              &

                          intensity_only,analysis_only,Grid_only,EVR,   &
                          cart,                                         &
                          GridTOBasis_test,OpPsi_test,                  &
                          optimization,nDfit,nDGrid,                    &
                          main_test,Opt_CAP_Basis,                      &

                          EVRT_path,File_path,base_FileName,            &
                          Srep_MPI,MPI_scheme,MPI_mc,MPI_iGs_auto,      &
                          MPI_fake_nodes,MPI_mem_node


        !-------------------------------------------------------------------------------
        ! set parallelization
        maxth              = 1
        !$ maxth           = omp_get_max_threads()

        Popenmp            = .FALSE.
        Popenmpi           = .FALSE.
#if(run_MPI)
        Popenmpi           = .TRUE.  !< True to run with MPI
        Popenmp            = .FALSE.
#else
        Popenmpi           = .FALSE.
        Popenmp            = (maxth > 1)
#endif


        intensity_only     = .FALSE.
        analysis_only      = .FALSE.
        Grid_only          = .FALSE.
        test               = .FALSE.
        cart               = .FALSE.
        GridTOBasis_test   = .FALSE.
        OpPsi_test         = .FALSE.   !< True for test of action
        EVR                = .FALSE.   ! ElVibRot (default)
        nDfit              = .FALSE.
        nDGrid             = .FALSE.
        main_test          = .FALSE.
        optimization       = 0
        Opt_CAP_Basis      = .FALSE.

        PMatOp_omp         = 0
        PMatOp_maxth       = maxth
        POpPsi_omp         = 0
        POpPsi_maxth       = maxth
        PBasisTOGrid_omp   = 0
        PBasisTOGrid_maxth = maxth
        PGrid_omp          = 1
        PGrid_maxth        = maxth
        PCRP_omp           = 0
        PCRP_maxth         = maxth
        PAna_omp           = 0
        PAna_maxth         = maxth

        PTune_SG4_omp      = .FALSE.
        PTune_Grid_omp     = .FALSE.

        PSG4_omp           = 1
        PSG4_maxth         = maxth

        max_mem          = 4000000000_ILkind/Rkind ! 4GO
        mem_debug        = .FALSE.
        printlevel       = -1

        EneFormat        = "f18.10"
        RMatFormat       = "f18.10"
        CMatFormat       = "f15.7"

        ! version and copyright statement
        CALL versionEVRT(.TRUE.)
        write(out_unit,*)
        IF(Popenmpi) THEN
          CALL ini_MPI()
          CALL time_perso('MPI start, initial time')
        ENDIF

        !read the file name for the command arguments
        CALL read_arg(input_filename)
        !> automatically decide the reading of namelist, from file or shell
        !> NOTE: remember to use vib to ensure "rm namelist" to prevent the
        !> reading of old namelist
        CALL file_open2(input_filename,in_unit,old=.TRUE.,err_file=err)
        IF(err/=0) THEN
          write(out_unit,*) input_filename,' file does not exist or error.'
          write(out_unit,*) '   => reading input data from shell.'
          in_unit=INPUT_UNIT
        ELSE
          write(out_unit,*) input_filename,' file does exist.'
          write(out_unit,*) '   => reading input data from the file.'
        ENDIF
        read(in_unit,system,IOSTAT=err)

        IF (err < 0) THEN
          write(out_unit,*) ' ERROR in ElVibRot (main program)'
          write(out_unit,*) ' End-of-file or End-of-record'
          write(out_unit,*) ' The namelist "system" is probably absent'
          write(out_unit,*) ' check your data!'
          write(out_unit,*) ' ERROR in ElVibRot (main program)'
          STOP
        ELSE IF (err > 0) THEN
          write(out_unit,*) ' ERROR in ElVibRot (main program)'
          write(out_unit,*) ' Some parameter name of the namelist "system" are probaly wrong'
          write(out_unit,*) ' check your data!'
          write(out_unit,system)
          write(out_unit,*) ' ERROR in ElVibRot (main program)'
          STOP
        END IF

        IF (base_FileName /= "" .AND. File_path /= "") THEN
          write(out_unit,*) ' ERROR in ElVibRot (main program)'
          write(out_unit,*) ' base_FileName and File_path are both set!!'
          write(out_unit,*) ' You MUST define only File_path.'
          write(out_unit,*) ' check your data!'
          write(out_unit,system)
          write(out_unit,*) ' ERROR in ElVibRot (main program)'
          STOP
        ELSE IF (base_FileName /= "") THEN
          File_path = base_FileName
        END IF
        Current_Path = trim(File_path) ! Current_Path is in FOR_EVRT library

        para_mem%mem_debug = mem_debug

        EVR = .NOT. (analysis_only .OR. GridTOBasis_test .OR.            &
                     OpPsi_test .OR. cart .OR. main_test .OR. nDfit .OR. &
                     nDGrid .OR. Opt_CAP_Basis .OR. optimization /= 0 .OR. analysis_only)

        IF (printlevel > 1) write(out_unit,system)

        para_EVRT_calc%optimization     = optimization
        para_EVRT_calc%EVR              = EVR
        para_EVRT_calc%analysis_only    = analysis_only
        para_EVRT_calc%intensity_only   = intensity_only
        para_EVRT_calc%Grid_only        = Grid_only

        para_EVRT_calc%cart             = cart
        para_EVRT_calc%GridTOBasis_test = GridTOBasis_test
        para_EVRT_calc%OpPsi_test       = OpPsi_test

        para_EVRT_calc%nDfit            = nDfit
        para_EVRT_calc%nDGrid           = nDGrid
        para_EVRT_calc%main_test        = main_test
        para_EVRT_calc%Opt_CAP_Basis    = Opt_CAP_Basis

        CALL set_print_level(printlevel,force=.TRUE.) ! print_level = printlevel ! print_level is in EVR_system_m.mod

        EneIO_format  = EneFormat
        RMatIO_format = RMatFormat
        CMatIO_format = "'('," // trim(adjustl(CMatFormat)) //      &
                    ",' +i'," // trim(adjustl(CMatFormat)) // ",')'"


        openmp                = (Popenmp .AND. maxth > 1)
        openmpi               = Popenmpi

        IF (.NOT. openmp) THEN
           MatOp_omp          = PMatOp_omp
           OpPsi_omp          = 0
           BasisTOGrid_omp    = 0
           Grid_omp           = 0
           CRP_omp            = 0
           Ana_omp            = 0
           SG4_omp            = 0

           MatOp_maxth        = 1
           OpPsi_maxth        = 1
           BasisTOGrid_maxth  = 1
           Grid_maxth         = 1
           CRP_maxth          = 1
           Ana_maxth          = 1
           SG4_maxth          = 1
           maxth              = 1

           Tune_SG4_omp       = .FALSE.
           Tune_Grid_omp      = .FALSE.
        ELSE
           MatOp_omp          = PMatOp_omp
           OpPsi_omp          = POpPsi_omp
           BasisTOGrid_omp    = PBasisTOGrid_omp
           Grid_omp           = PGrid_omp
           CRP_omp            = PCRP_omp
           SG4_omp            = PSG4_omp

           IF (MatOp_omp > 0) THEN
             MatOp_maxth        = min(PMatOp_maxth,maxth)
           ELSE
             MatOp_maxth        = 1
           END IF

           IF (OpPsi_omp > 0) THEN
             OpPsi_maxth        = min(POpPsi_maxth,maxth)
           ELSE
             OpPsi_maxth        = 1
           END IF

           IF (BasisTOGrid_omp > 0) THEN
             BasisTOGrid_maxth  = min(PBasisTOGrid_maxth,maxth)
           ELSE
             BasisTOGrid_maxth  = 1
           END IF

           IF (Grid_omp > 0) THEN
             Grid_maxth         = min(PGrid_maxth,maxth)
             Tune_Grid_omp      = PTune_Grid_omp
           ELSE
             Grid_maxth         = 1
             Tune_Grid_omp      = .FALSE.
           END IF

           IF (CRP_omp > 0) THEN
             CRP_maxth          = min(PCRP_maxth,maxth)
           ELSE
             CRP_maxth         = 1
           END IF

           IF (Ana_omp > 0) THEN
            Ana_maxth          = min(PAna_maxth,maxth)
          ELSE
            Ana_maxth         = 1
          END IF

          IF (SG4_omp > 0) THEN
             SG4_maxth         = min(PSG4_maxth,maxth)
             Tune_SG4_omp      = PTune_SG4_omp
           ELSE
             SG4_maxth         = 1
             Tune_SG4_omp       = .FALSE.
           END IF

        END IF
        MatOp_maxth_init        = MatOp_maxth
        OpPsi_maxth_init        = OpPsi_maxth
        BasisTOGrid_maxth_init  = BasisTOGrid_maxth
        Grid_maxth_init         = Grid_maxth
        CRP_maxth_init          = CRP_maxth
        Ana_maxth_init          = Ana_maxth
        SG4_maxth_init          = SG4_maxth

        write(out_unit,*) '========================================='
        write(out_unit,*) 'OpenMP parameters:',openmp
        write(out_unit,*) 'Max number of threads:           ',maxth
        write(out_unit,*) 'MatOp_omp,      MatOp_maxth      ',MatOp_omp,MatOp_maxth
        write(out_unit,*) 'OpPsi_omp,      OpPsi_maxth      ',OpPsi_omp,OpPsi_maxth
        write(out_unit,*) 'BasisTOGrid_omp,BasisTOGrid_maxth',BasisTOGrid_omp,BasisTOGrid_maxth
        write(out_unit,*) 'Grid_omp,       Grid_maxth       ',Grid_omp,Grid_maxth
        write(out_unit,*) 'CRP_omp,        CRP_maxth        ',CRP_omp,CRP_maxth
        write(out_unit,*) 'Ana_omp,        Ana_maxth        ',Ana_omp,Ana_maxth
        write(out_unit,*) 'SG4_omp,        SG4_maxth        ',SG4_omp,SG4_maxth
        write(out_unit,*) '========================================='

        write(out_unit,*) '========================================='
        write(out_unit,*) 'File_path: ',trim(adjustl(File_path))
        write(out_unit,*) '========================================='

        para_mem%max_mem    = max_mem/Rkind
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='

        IF (para_EVRT_calc%optimization /= 0) THEN
          IF(MPI_id==0) write(out_unit,*) ' Optimization calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL sub_Optimization_OF_VibParam(max_mem)

        ELSE IF (para_EVRT_calc%nDfit .OR. para_EVRT_calc%nDGrid) THEN
          IF(MPI_id==0) write(out_unit,*) ' nDfit or nDGrid calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL sub_nDGrid_nDfit()

        ELSE IF (para_EVRT_calc%EVR) THEN
          IF(MPI_id==0) write(out_unit,*) ' ElVibRot calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL vib(max_mem,test,intensity_only)

        ELSE IF (para_EVRT_calc%cart) THEN
          IF(MPI_id==0) write(out_unit,*) ' cart calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL sub_cart(max_mem)

        ELSE IF (para_EVRT_calc%GridTOBasis_test) THEN
          IF(MPI_id==0) write(out_unit,*) ' sub_GridTOBasis calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL sub_GridTOBasis_test(max_mem)

        ELSE IF (para_EVRT_calc%OpPsi_test) THEN
          IF(MPI_id==0) write(out_unit,*) ' OpPsi calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL Sub_OpPsi_test(max_mem)

        ELSE IF (para_EVRT_calc%analysis_only) THEN
          IF(MPI_id==0) write(out_unit,*) ' WP analysis calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL sub_analysis_only(max_mem)

        ELSE IF (para_EVRT_calc%main_test) THEN
          IF(MPI_id==0) write(out_unit,*) ' Smolyat test calculation'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL sub_main_Smolyak_test()

        ELSE IF (para_EVRT_calc%Opt_CAP_Basis) THEN
          IF(MPI_id==0) write(out_unit,*) ' Optimization of the CAP and basis parameters (CRP)'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL sub_Opt_CAP_basis()

        ELSE
          IF(MPI_id==0) write(out_unit,*) ' ElVibRot calculation (default)'
          IF(MPI_id==0) write(out_unit,*) '========================================='
          CALL vib(max_mem,test,intensity_only)
        END IF

        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='

        close(in_unit)   ! CALL file_close cannot be used
        IF(openmpi) THEN
          CALL time_perso('MPI closed')
          CALL end_MPI()
        ENDIF

END PROGRAM ElVibRot
SUBROUTINE read_arg(input_filename)
  USE EVR_system_m
  IMPLICIT NONE

  character(len=:), allocatable, intent(inout) :: input_filename


  character(len=:), allocatable :: arg,arg2
  integer :: i,arg_len

  IF (COMMAND_ARGUMENT_COUNT() /= 0 .AND. COMMAND_ARGUMENT_COUNT() /= 2) THEN
    write(out_unit,*) ' ERROR in read_arg'
    write(out_unit,*) ' Wrong ElVibRot argument number!'
    write(out_unit,*) 'argument number',COMMAND_ARGUMENT_COUNT()
    write(out_unit,*) ' You can have 0 or 2 arguments.'
    STOP 'Wrong ElVibRot argument number'
  END IF


  DO i=1, COMMAND_ARGUMENT_COUNT(),2

    CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=arg_len )
    allocate( character(len=arg_len) :: arg )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=arg_len )
    allocate( character(len=arg_len) :: arg2 )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

    SELECT CASE(arg)
    CASE("-i","--input")
      input_filename = arg2
    CASE Default
      write(out_unit,*) ' ERROR in read_arg'
      write(out_unit,*) ' Wrong ElVibRot argument!'
      write(out_unit,*) '   arg: "',arg,'"'
      write(out_unit,*) ' The possibilities are:'
      write(out_unit,*) '    -i or --input'
      STOP 'Wrong ElVibRot argument'
    END SELECT

    write(out_unit,*) 'Argument number: ',i,' ==> arg: "',arg,'", arg2: "',arg2,'"'

    deallocate(arg)
    deallocate(arg2)
  END DO
  IF (.NOT. allocated(input_filename)) THEN
    write(out_unit,*) ' WARNING in read_arg'
    write(out_unit,*) ' No input file name argument'
    write(out_unit,*) '    => the file name is "namelist".'
    input_filename = 'namelist'
  END IF

  write(out_unit,*) '=================================='
  write(out_unit,*) '=================================='

END SUBROUTINE read_arg
