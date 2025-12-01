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
MODULE EVR_system_m
  USE QDUtil_m
  USE mod_MPI
  USE FOR_EVRT_system_m, ONLY : param_FOR_optimization, para_FOR_optimization, &
                                Write_Mat_MPI, Write_Vec_MPI
  IMPLICIT NONE

#if defined(__EVR_VER)
  character (len=Name_len) :: EVR_version = __EVR_VER
#else
  character (len=Name_len) :: EVR_version = "unknown: -D__EVR_VER=?"
#endif

#if defined(__EVRTPATH)
  character (len=Line_len) :: EVRT_path   =                         &
       __EVRTPATH
#else
  character (len=Line_len) :: EVRT_path   = '~/ElVibRot'
#endif

#if defined(__COMPILE_DATE)
  character (len=Line_len) :: compile_date = __COMPILE_DATE
#else
  character (len=Line_len) :: compile_date = "unknown: -D__COMPILE_DATE=?"
#endif

#if defined(__COMPILE_HOST)
  character (len=Line_len) :: compile_host = __COMPILE_HOST
#else
  character (len=Line_len) :: compile_host = "unknown: -D__COMPILE_HOST=?"
#endif

  logical :: openmp = .FALSE.
  logical :: openmpi= .FALSE.
 
  integer :: MatOp_omp,MatOp_maxth,MatOp_maxth_init
  integer :: OpPsi_omp,OpPsi_maxth,OpPsi_maxth_init
  integer :: BasisTOGrid_omp,BasisTOGrid_maxth,BasisTOGrid_maxth_init
  integer :: Grid_omp,Grid_maxth,Grid_maxth_init
  integer :: SG4_omp,SG4_maxth,SG4_maxth_init
  integer :: CRP_omp,CRP_maxth,CRP_maxth_init
  integer :: Ana_omp,Ana_maxth,Ana_maxth_init

  logical :: Tune_SG4_omp  = .FALSE.
  logical :: Tune_Grid_omp = .FALSE.

  integer (kind=ILkind) :: nb_mult_BTOG  = 0
  integer (kind=ILkind) :: nb_mult_GTOB  = 0
  integer (kind=ILkind) :: nb_mult_OpPsi = 0

  integer, parameter :: max_HADA = 5000
  integer, parameter :: max_nb_G_FOR_print = 2000
  !integer, parameter :: max_nb_G_FOR_print = 20000

  integer :: SGtype               = -1
  integer :: FilePsiVersion       = 0
  logical :: NewBasisEl           = .FALSE.

  character (len=:), allocatable :: Current_Path


  character (len=Name_longlen) :: EneIO_format = "f20.5"

  TYPE param_EVRT_calc
    integer :: optimization     = 0
    logical :: EVR              = .TRUE.   ! ElVibRot (default)
    logical :: analysis_only    = .FALSE.
    logical :: intensity_only   = .FALSE.
    logical :: Grid_only        = .FALSE.
    logical :: cart             = .FALSE.
    logical :: GridTOBasis_test = .FALSE.
    logical :: OpPsi_test       = .FALSE.
    logical :: nDfit            = .FALSE.
    logical :: nDGrid           = .FALSE.
    logical :: main_test        = .FALSE.
    logical :: Opt_CAP_Basis    = .FALSE.
  END TYPE param_EVRT_calc

  TYPE (param_EVRT_calc),        save :: para_EVRT_calc

CONTAINS

!=======================================================================================
!> @brief subroutine for recording time
!> @param time_sum should be initialized before calling this function
!=======================================================================================
      SUBROUTINE time_record(time_sum,time1,time2,point)
        USE mod_MPI
        IMPLICIT NONE

        Integer,                        intent(inout) :: time_sum
        Integer,                        intent(inout) :: time1
        Integer,                        intent(inout) :: time2
        Integer,                           intent(in) :: point

        IF(point==1) THEN
          CALL system_clock(time1,time_rate,time_max)
        ELSEIF(point==2) THEN
          CALL system_clock(time2,time_rate,time_max)
          time_sum=time_sum+merge(time2-time1,time2-time1+time_max,time2>=time1)
        ELSE
          STOP 'error when calling time_record'
        ENDIF
  END SUBROUTINE

  FUNCTION nom_i(nom1,i1)
    IMPLICIT NONE
  
    character (len=14) :: nom_i
    character (len=10) :: nom1
    character (len=14) :: nom2
    integer            :: j,i1
  
    write(out_unit,*) nom1,i1
    IF (i1 .GT. 100 ) STOP ' in nom_i: i1 too big'
  
    write(nom2,'(a10,i2)') nom1,i1
    DO j=1,12  ! it has to be 12 and not 14
      IF (nom2(j:j) .EQ. ' ') nom2(j:j)='_'
    END DO
    nom_i=nom2
  
  END FUNCTION nom_i
  
  FUNCTION nom_ii(nom1,i1,i2)
    IMPLICIT NONE
  
    character (len=14) :: nom_ii
  
    character (len=10) :: nom1
    character (len=14) :: nom2
    integer            :: j,i1,i2
  
    !write(out_unit,*) nom1,i1,i2
    IF (i1 .GT. 100 .OR. i2 .GT. 100) STOP ' in nom_ii: i1 or i2 too big'
  
    write(nom2,'(a10,2i2)') nom1,i1,i2
    DO j=1,14
      IF (nom2(j:j) .EQ. ' ') nom2(j:j)='_'
    END DO
    nom_ii=nom2
  
  END FUNCTION nom_ii
  !================================================================
  !       analysis of a string
  !       output : nb_word,word (the i th word of name)
  !================================================================
  SUBROUTINE analysis_name(name,word,i,nb_word)
    IMPLICIT NONE

    !- analysis of the basis name --------------------------
    character (len=*)  :: word,name
    character          :: ch
    integer            :: i,nb_word
    integer            :: iw,ic,icw
    logical            :: blank


    iw = 0
    icw = 0
    blank = .TRUE.
    !write(out_unit,*) 'analysis_name: ',name,len(name)
    DO ic=1,len(name)
      ch = name(ic:ic)
      IF (ch .EQ. " ") THEN
        IF (.NOT. blank) THEN
          iw = iw + 1
          blank = .TRUE.
        END IF
      ELSE
        IF (iw .EQ. i-1) THEN
          icw = icw + 1
          word(icw:icw) = ch
        END IF
        blank = .FALSE.
      END IF
     !write(out_unit,*) 'analysis_name: ',ic,ch,blank,iw
    END DO

    nb_word = iw
    !write(out_unit,*) 'analysis_name: ',name,':',nb_word
    !write(out_unit,*) 'analysis_name: ',i,word


  END SUBROUTINE analysis_name
  FUNCTION make_EVRTInternalFileName(FileName,FPath) RESULT(make_FileName)
    USE QDUtil_m, ONLY : err_FileName
    IMPLICIT NONE

    character (len=:), allocatable          :: make_FileName

    character(len=*), intent(in)            :: FileName
    character(len=*), intent(in), optional  :: FPath

    character (len=:), allocatable          :: FPath_loc


    integer :: ilast_char,err

    IF (present(FPath)) THEN
      FPath_loc = FPath
    ELSE
      FPath_loc = EVRT_path
    END IF

    err = err_FileName(FileName,name_sub='make_EVRTInternalFileName')
    IF (err /= 0) STOP 'ERROR in make_EVRTInternalFileName: problem with the FileName'

    ilast_char = len_trim(FPath_loc)

    IF (FileName(1:1) == "/" .OR. FileName(1:1) == "~" .OR. ilast_char == 0) THEN
      make_FileName = trim(adjustl(FileName))
    ELSE
      IF (FPath_loc(ilast_char:ilast_char) == "/") THEN
        make_FileName = trim(adjustl(FPath_loc)) // trim(adjustl(FileName))
      ELSE
        make_FileName = trim(adjustl(FPath_loc)) // '/' // trim(adjustl(FileName))
      END IF
    END IF

    IF (allocated(FPath_loc)) deallocate(FPath_loc)

  END FUNCTION make_EVRTInternalFileName
  FUNCTION make_EVRTFileName(FileName,FPath) RESULT(make_FileName)
    USE QDUtil_m, ONLY : err_FileName
    IMPLICIT NONE

    character (len=:), allocatable          :: make_FileName

    character(len=*), intent(in)            :: FileName
    character(len=*), intent(in), optional  :: FPath

    character (len=:), allocatable          :: FPath_loc


    integer :: ilast_char,err

    IF (present(FPath)) THEN
      FPath_loc = FPath
    ELSE IF (allocated(Current_Path)) THEN
      FPath_loc = Current_Path
    ELSE
      FPath_loc = ''
    END IF

    err = err_FileName(FileName,name_sub='make_EVRTFileName')
    IF (err /= 0) THEN 
      !write(out_unit,*) 'ERROR in make_EVRTFileName: problem with the FileName'
      !write(out_unit,*) 'FileName: "',trim(FileName),'"'
      STOP 'ERROR in make_EVRTFileName: problem with the FileName'
    END IF

    ilast_char = len_trim(FPath_loc)

    IF (FileName(1:1) == "/" .OR. FileName(1:1) == "~" .OR. ilast_char == 0) THEN
      make_FileName = trim(adjustl(FileName))
    ELSE
      IF (FPath_loc(ilast_char:ilast_char) == "/") THEN
        make_FileName = trim(adjustl(FPath_loc)) // trim(adjustl(FileName))
      ELSE
        make_FileName = trim(adjustl(FPath_loc)) // '/' // trim(adjustl(FileName))
      END IF
    END IF

    IF (allocated(FPath_loc)) deallocate(FPath_loc)

  END FUNCTION make_EVRTFileName
  SUBROUTINE versionEVRT(write_version)
    USE iso_fortran_env,   ONLY : compiler_version, compiler_options
    USE TnumTana_system_m, ONLY : TnumTana_version
    IMPLICIT NONE
  
    logical, intent(in) :: write_version
  
    character (len=*), parameter :: EVR_name='ElVibRot'  
  
    IF (write_version .AND. MPI_id==0) THEN
      write(out_unit,*) '==============================================='
      write(out_unit,*) '==============================================='
      write(out_unit,*) 'Working with ',EVR_name,trim(adjustl(EVR_version))

      write(out_unit,*) 'Compiled on "',trim(compile_host), '" the ',trim(compile_date)
      write(out_unit,*) 'Compiler:         ',compiler_version()
      write(out_unit,*) 'Compiler options: ',compiler_options()

      write(out_unit,*) 'EVRT_path: ',trim(EVRT_path)

      write(out_unit,*) '-----------------------------------------------'

      write(out_unit,*) EVR_name,' is written by David Lauvergnat [1] '
      write(out_unit,*) '  with contributions of'
      write(out_unit,*) '     Josep Maria Luis (optimization) [2]'
      write(out_unit,*) '     Ahai Chen (MPI) [1,4]'
      write(out_unit,*) '     Lucien Dupuy (CRP) [5]'

      write(out_unit,*) EVR_name,' is under MIT license.'

      write(out_unit,*)
      write(out_unit,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
      write(out_unit,*) '[2]: Institut de Química Computacional and Departament de Química',&
                                 ' Universitat de Girona, Catalonia, Spain'
      write(out_unit,*) '[3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark'
      write(out_unit,*) '[4]: Maison de la Simulation USR 3441, CEA Saclay, France'
      write(out_unit,*) '[5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,', &
                                 ' Université de Montpellier, France'
      write(out_unit,*) '==============================================='
      write(out_unit,*) '==============================================='
      CALL TnumTana_version(write_version)
    END IF
  END SUBROUTINE versionEVRT

END MODULE EVR_system_m