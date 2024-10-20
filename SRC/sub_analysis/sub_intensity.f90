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

!================================================================
!
!     Intensity of transitions
!
!================================================================
      SUBROUTINE sub_intensity(para_Dip,print_Op,para_H,nb_ana,         &
                               para_intensity,const_phys,               &
                               intensity_only,nio_res_int)

      USE EVR_system_m
      use mod_Constant, only: constant, rwu_write, real_wu, get_conv_au_to_unit
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      USE mod_analysis
      IMPLICIT NONE


!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Dip(3),para_H
      logical          :: print_Op

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_intensity) :: para_intensity

!----- physical and mathematical constants ---------------------------
      TYPE (constant)            :: const_phys

      logical :: intensity_only
      integer :: nio_res_int

      integer :: nb_aie,nb_ana,nio

      real (kind=Rkind) ::    conv
      integer           ::    DE_pow


!---- variable for the Z-matrix ------------------------------------------
      TYPE (CoordType), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum


      real (kind=Rkind), allocatable ::    Mat_Aif(:,:)

      integer       ::    i,j,k

      real (kind=Rkind) ::    val,DE,Aif,e0
      real (kind=Rkind) ::    emin,emax
      real (kind=Rkind) ::    Ai,Bi,Ci,Af,Bf,Cf,Imax,Iif
      integer           ::    JJ,KK,Jmax,ind_f
      real (kind=Rkind) ::    zpe,Q
      real (kind=Rkind) ::    pop(para_H%nb_tot)
      real (kind=Rkind) ::    Ewidth,auTOenergy

      integer   :: err
!----- function ------------------------------------------------------
      real (kind=Rkind) ::    tcmc_ij
      real (kind=Rkind) ::    part_func,pop2_i
!----- function ------------------------------------------------------



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_intensity'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum

      nb_aie = para_H%nb_tot

      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'nb_aie',nb_aie
        write(out_unit,*) 'nb_ana',nb_ana
        write(out_unit,*) 'intensity_only,nio_res_int',                        &
                    intensity_only,nio_res_int
        write(out_unit,*)
        write(out_unit,*) ' l_Int,l_IntVR:',                                   &
                     para_intensity%l_Int,para_intensity%l_IntVR
        write(out_unit,*) ' l_Aif,l_Tau:',                                     &
                     para_intensity%l_Aif,para_intensity%l_Tau
        write(out_unit,*)
        write(out_unit,*) 'convIDif',const_phys%convIDif
        write(out_unit,*) 'convAif',const_phys%convAif
        write(out_unit,*)
        write(out_unit,*) 'Rvp',shape(para_H%Rvp)
!       CALL Write_Mat_MPI(para_H%Rvp,out_unit,5)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------


      IF (nb_ana > 0) THEN
        CALL alloc_NParray(Mat_Aif,[nb_ana,nb_ana],'Mat_Aif',name_sub)
        Mat_Aif(:,:) = ZERO
        IF (para_intensity%l_IntVR .AND. .NOT. allocated(para_intensity%ABC)) THEN
          CALL alloc_NParray(para_intensity%ABC,[3,nb_ana],           &
                            'para_intensity%ABC',name_sub)
          para_intensity%ABC(:,:) = ZERO
        END IF
      ELSE
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_ana <=0',nb_ana
        STOP
      END IF
      flush(out_unit)

      IF (intensity_only) THEN

        read(nio_res_int,*)
        CALL Read_Mat(Mat_Aif,nio_res_int,5,err)
        IF (err /= 0) THEN
           write(out_unit,*) 'ERROR in ',name_sub
           write(out_unit,*) ' reading the matrix "Mat_Aif"'
           STOP
        END IF
        IF (allocated(para_intensity%ABC)) THEN
          read(nio_res_int,*)
          CALL Read_Mat(para_intensity%ABC,nio_res_int,5,err)
          IF (err /= 0) THEN
            write(out_unit,*) 'ERROR in ',name_sub
            write(out_unit,*) ' reading "para_intensity%ABC"'
            STOP
          END IF
        END IF

      ELSE

!       -------------------------------------------
!       write(out_unit,*) ' ene (ua): ',nb_ana
!       DO i=1,nb_ana
!         write(out_unit,*) i,para_H%Rdiag(i)
!       END DO
!       write(out_unit,*) ' END ene',nb_ana
        write(out_unit,*)
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        write(out_unit,*) ' Calculation and saving of "Mat_Aif(:,:)": '
!       -- initialization -------------------------
        Mat_Aif(:,:) = ZERO
        DO k=1,3
          IF (para_intensity%pola_xyz(k)) THEN
            write(out_unit,*) ' sub_intensity: ',trim(para_Dip(k)%name_Op),'+ <psi|dip|psi>'
            write(out_unit,*) ' alloc Rmat: ',trim(para_Dip(k)%name_Op),' ',allocated(para_Dip(k)%Rmat)

            IF (.NOT. para_Dip(k)%spectral .OR. .NOT. para_Dip(k)%mat_done) THEN
              write(out_unit,*) ' ERROR in sub_intensity'
              write(out_unit,*) ' The Dipolar operator MUST have a spectral representation'
              write(out_unit,*) ' or mat_done is NOT .TRUE.'
              STOP
            END IF
            Mat_Aif(:,:) = Mat_Aif(:,:) + para_Dip(k)%Rmat(1:nb_ana,1:nb_ana)**2
            !write(out_unit,*) 'para_Dip(k)%Rmat',k
            !CALL Write_Mat_MPI(para_Dip(k)%Rmat,out_unit,5)
            CALL dealloc_para_Op(para_Dip(k))
          END IF

          flush(out_unit)
        END DO
        write(nio_res_int,*) 'Mat_Aif'
        CALL Write_Mat_MPI(Mat_Aif,nio_res_int,5,Rformat='e30.23')
        IF (allocated(para_intensity%ABC)) THEN
          write(nio_res_int,*) 'ABC'
          CALL Write_Mat_MPI(para_intensity%ABC,nio_res_int,5,Rformat='e30.23')
        END IF
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        flush(out_unit)

      END IF


      IF (para_intensity%l_Aif .OR. para_intensity%l_Tau) THEN
!       coef Einstein
        DE_pow = 3
        conv = const_phys%convAif
      ELSE IF (para_intensity%l_int) THEN
!       intensity (dipole)
        DE_pow = 1
        conv = const_phys%convIDif
      ELSE IF (para_intensity%l_CrossSec) THEN
!       Cross Section (m^2)
        DE_pow = 1
        conv = ONE
        STOP
      ELSE IF (para_intensity%l_intVR) THEN
!       intensity (dipole)
        DE_pow = 0
        conv = ONE
      END IF


      DO i=1,nb_ana
      DO j=i+1,nb_ana

        DE = (para_H%Rdiag(j)-para_H%Rdiag(i))
        Mat_Aif(i,j) = DE**DE_pow * Mat_Aif(i,j) * conv
        Mat_Aif(j,i) = Mat_Aif(i,j)

      END DO
      END DO

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        write(out_unit,*) ' DE(i-f)*S(i-f):'
        CALL Write_Mat_MPI(Mat_Aif,out_unit,5,Rformat='e15.8')
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='

      END IF

      write(out_unit,*)
      write(out_unit,*) '==================================================='
      write(out_unit,*) '==================================================='
      IF (para_intensity%l_Tau) THEN
        DO i=1,nb_ana
        DO j=i+1,nb_ana
          Mat_Aif(i,j) = ONE/ Mat_Aif(i,j)
          Mat_Aif(j,i) = Mat_Aif(i,j)
        END DO
        END DO
        write(out_unit,*) 'Life time (1/Aif) in s'
        CALL Write_Mat_MPI(Mat_Aif,out_unit,5,Rformat='e15.8')
        write(out_unit,*) '==================================================='

      ELSE IF (para_intensity%l_Aif) THEN
        write(out_unit,*) 'Einstein coef (Aif) in s-1'
        CALL Write_Mat_MPI(Mat_Aif,out_unit,5,Rformat='e15.8')
        write(out_unit,*) '==================================================='

      ELSE IF (para_intensity%l_Int) THEN
        write(out_unit,*) 'intensity (with dipolar moments)'
        ! add temperature
        zpe = minval(para_H%Rdiag)
        Q=part_func(para_H%Rdiag,nb_aie,para_intensity%Temp)
        write(out_unit,*)
        write(out_unit,*) 'population at Temp,Q',para_intensity%Temp,Q
        write(out_unit,*)
        DO i=1,nb_ana
          pop(i) = pop2_i(para_H%Rdiag(i),zpe,para_intensity%Temp)/Q
        END DO
        write(out_unit,*)
        DO i=1,nb_ana
        DO j=i+1,nb_ana
          Mat_Aif(i,j) = Mat_Aif(i,j)*(pop(i)-pop(j))
          Mat_Aif(j,i) = Mat_Aif(i,j)
        END DO
        END DO
      ELSE IF (para_intensity%l_IntVR) THEN
        write(out_unit,*) 'intensity Vib+Rot (dip)'
        IF (debug) CALL Write_Mat_MPI(Mat_Aif,out_unit,5,Rformat='e15.8')
      END IF
      flush(out_unit)

      IF (para_intensity%l_Int .OR. para_intensity%l_CrossSec) THEN

        emin = para_intensity%Emin

        IF (para_intensity%Emax == ZERO) THEN
          emax = (para_H%Rdiag(nb_ana)-para_H%Rdiag(1))
        ELSE
          emax = para_intensity%Emax
        END IF

        Ewidth = para_intensity%Ewidth

        write(out_unit,*) 'spectrum Temp:',para_intensity%Temp
        write(out_unit,*) 'spectrum Ewidth: ',                         &
          RWU_Write(REAL_WU(Ewidth,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emin: ',                           &
          RWU_Write(REAL_WU(emin,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emax: ',                      &
          RWU_Write(REAL_WU(emax,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        CALL sub_spectre(Mat_Aif,para_H%Rdiag,nb_ana,                   &
                         Ewidth,emin,emax,para_intensity%file_spectrum,para_intensity%Min_relativeI0)

      ELSE IF (para_intensity%l_IntVR) THEN
        emin = para_intensity%Emin


!       - calculation of the partition function Q
        Q = ZERO
        Jmax = para_intensity%JJmax
        write(out_unit,*) ' Jmax : ',Jmax
        flush(out_unit)

        IF (.NOT. allocated(para_intensity%ABC)) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) ' The table of rotational constants (ABC) is not allocated'
          write(out_unit,*) ' CHECK the fortran!!'
          STOP
        END IF
        flush(out_unit)

!       - calculation of the partition function Q
        DO i=1,nb_ana
          Ai = para_intensity%ABC(1,i)
          Bi = para_intensity%ABC(2,i)
          Ci = para_intensity%ABC(3,i)
          CALL calc_Q_VR(Q,pop(i),para_H%Rdiag(1),para_H%Rdiag(i),      &
                         Ai,Bi,Ci,para_intensity%Temp,Jmax,const_phys)

        END DO

!       - population sum over rotational levels
        write(out_unit,*) ' Population :'
        DO i=1,nb_ana
          pop(i) = pop(i) /Q
          write(out_unit,*) i-1,(para_H%Rdiag(i)-para_H%Rdiag(1)) *            &
                                    const_phys%auTOenergy,pop(i)
        END DO


!       - intensity maximal
        Imax = ZERO
        Iif= ZERO
        ind_f=0
        DO i=1,nb_ana
          Iif = (para_H%Rdiag(i)-para_H%Rdiag(1))*Mat_Aif(1,i)
          write(out_unit,*) 'calc max values i,Aij,Iif',i,Mat_Aif(1,i),Iif
          IF (Iif > Imax) THEN
            Imax = Iif
            ind_f = i
          END IF
        END DO
        write(out_unit,*) 'Vibrational max values : Imax,ind_f',Imax,ind_f

        Ai= para_intensity%ABC(1,1)
        Bi= para_intensity%ABC(2,1)
        Ci= para_intensity%ABC(3,1)
        Af= para_intensity%ABC(1,ind_f)
        Bf= para_intensity%ABC(2,ind_f)
        Cf= para_intensity%ABC(3,ind_f)
        CALL calc_Imax(Imax,para_H%Rdiag(1),Ai,Bi,Ci,                   &
                       para_H%Rdiag(ind_f),Af,Bf,Cf,                    &
                       Mat_Aif(1,ind_f),para_intensity%Temp,            &
                       Q,Jmax,const_phys)

        write(out_unit,*) 'Imax',Imax
        flush(out_unit)
!       - spectrum : RoVib
        emin = para_intensity%Emin

        IF (para_intensity%Emax == ZERO) THEN
          emax = (para_H%Rdiag(nb_ana)-para_H%Rdiag(1))
        ELSE
          emax = para_intensity%Emax
        END IF

        Ewidth = para_intensity%Ewidth

        write(out_unit,*) 'spectrum Temp,Q:',para_intensity%Temp,Q
        write(out_unit,*) 'spectrum Ewidth: ',                         &
          RWU_Write(REAL_WU(Ewidth,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emin: ',                           &
          RWU_Write(REAL_WU(emin,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emax: ',                      &
          RWU_Write(REAL_WU(emax,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)




        CALL sub_spectreRV(Mat_Aif,para_H%Rdiag,nb_ana,para_intensity%ABC,pop,&
                           para_intensity%Temp,Q,Jmax,Imax,const_phys,  &
                           Ewidth,emin,emax,                            &
                           para_intensity%file_spectrum,                &
                           para_intensity%file_intensity)

      END IF
      write(out_unit,*) '==================================================='
      write(out_unit,*) '==================================================='

      CALL dealloc_NParray(Mat_Aif,'Mat_Aif','sub_intensity')

      IF (para_intensity%l_IntVR) THEN
        CALL dealloc_NParray(para_intensity%ABC,"para_intensity%ABC",'sub_intensity')
      END IF
!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END sub_intensity'
!----------------------------------------------------------


      END SUBROUTINE sub_intensity


!================================================================
!
!     sub_AnalysePsy_ScalOp
!
!================================================================
      SUBROUTINE sub_AnalysePsy_ScalOp(para_ScalOp,nb_ScalOp,para_H,nb_ana)

      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      USE mod_analysis
      IMPLICIT NONE


!----- Operator variables --------------------------------------------
      integer          :: nb_ScalOp
      TYPE (param_Op)  :: para_ScalOp(nb_ScalOp),para_H


      integer :: nb_aie,nb_ana,nio

      real (kind=Rkind) ::    conv
      integer           ::    DE_pow


!---- variable for the Z-matrix ------------------------------------------
      TYPE (CoordType), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum


      real (kind=Rkind), allocatable ::    Mat_Aif(:,:)

      integer       ::    i,j,k

      real (kind=Rkind) ::    val,DE,Aif,e0
      real (kind=Rkind) ::    emin,emax
      real (kind=Rkind) ::    Ai,Bi,Ci,Af,Bf,Cf,Imax,Iif
      integer           ::    JJ,KK,Jmax,ind_f
      real (kind=Rkind) ::    zpe,Q
      real (kind=Rkind) ::    pop(para_H%nb_tot)
      real (kind=Rkind) ::    Ewidth,auTOenergy

      integer   :: err
!----- function ------------------------------------------------------
      real (kind=Rkind) ::    tcmc_ij
      real (kind=Rkind) ::    part_func,pop2_i
!----- function ------------------------------------------------------



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_AnalysePsy_ScalOp'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum

      nb_aie = para_H%nb_tot

      write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*) 'nb_ana',nb_ana
        write(out_unit,*) 'nb_ScalOp',nb_ScalOp
        write(out_unit,*) 'size(para_ScalOp)',size(para_ScalOp)

        write(out_unit,*) 'Rvp',shape(para_H%Rvp)
        !CALL Write_Mat_MPI(para_H%Rvp,out_unit,5)
        write(out_unit,*)
        flush(out_unit)
      END IF
!-----------------------------------------------------------


      IF (nb_ana > 0) THEN
        CALL alloc_NParray(Mat_Aif,[nb_ana,nb_ana],'Mat_Aif',name_sub)
        Mat_Aif(:,:) = ZERO
      ELSE
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_ana <=0',nb_ana
        STOP
      END IF
      flush(out_unit)

      DO k=1,nb_ScalOp
        write(out_unit,*)
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        write(out_unit,*) ' Calculation and saving of "Mat_Aif(:,:)": ',k

        !para_ScalOp(k)%spectral = .TRUE.
        !CALL sub_MatOp(para_ScalOp(k),.TRUE.)

        !Mat_Aif(:,:) = Mat_Aif(:,:) + para_ScalOp(k)%Rmat(1:nb_ana,1:nb_ana)**2
        !write(out_unit,*) 'para_ScalOp(k)%Rmat',k
        !CALL Write_Mat_MPI(para_ScalOp(k)%Rmat,out_unit,5)
        !write(out_unit,*) 'Mat_Aif'
        !CALL Write_Mat_MPI(Mat_Aif,out_unit,5,Rformat='e30.23')
        !write(out_unit,*) 'Mat_Aif(:,1)',Mat_Aif(:,1)

        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        flush(out_unit)
      END DO

      write(out_unit,*)
      write(out_unit,*) '==================================================='
      write(out_unit,*) '==================================================='
      DO i=1,nb_ana
        write(out_unit,'(i0,20(1x,f0.4))') i,(abs(para_ScalOp(k)%Rmat(i,1)),k=1,nb_ScalOp)
      END DO
      write(out_unit,*) '==================================================='
      write(out_unit,*) '==================================================='
      write(out_unit,*)

      CALL dealloc_NParray(Mat_Aif,'Mat_Aif','sub_intensity')

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
!----------------------------------------------------------


      END SUBROUTINE sub_AnalysePsy_ScalOp

!================================================================
!
!     spectrum RV
!
!================================================================
      SUBROUTINE sub_spectreRV(Mat_Aif,ene,nb_ana,ABC,pop,              &
                               Temp,Q,Jmax,Imax,const_phys,             &
                               width,emin,emax,                         &
                               file_spectrum,file_intensity)
      USE EVR_system_m
      use mod_Constant, only: constant, RWU_Write,REAL_WU,get_Conv_au_TO_unit
      IMPLICIT NONE

!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys


      integer           :: nb_ana
      real (kind=Rkind) :: Mat_Aif(nb_ana,nb_ana)
      real (kind=Rkind) :: ene(nb_ana),ABC(3,nb_ana),pop(nb_ana)

      real (kind=Rkind) :: Ai,Bi,Ci,Aj,Bj,Cj,Temp,Q,Imax
      real (kind=Rkind) :: ImaxV,IntV
      integer           :: Jmax

      real (kind=Rkind) :: width,emin,emax,pas
      TYPE (File_t) :: file_spectrum,file_intensity
      integer           :: nio,nio_int
      integer           :: n
      real (kind=Rkind),allocatable :: spectre(:)


      real (kind=Rkind) :: e,e0,auTOenergy
      integer           :: i,j

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_spectreRV'
      integer :: err_mem,memory
!     logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_ana',nb_ana
        write(out_unit,*) 'spectrum Ewidth: ',                         &
          RWU_Write(REAL_WU(width,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emin: ',                           &
          RWU_Write(REAL_WU(emin,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emax: ',                      &
          RWU_Write(REAL_WU(emax,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)
        DO i=1,nb_ana
          write(out_unit,*) 'ABC (cm-1)',i,ABC(:,i)*get_Conv_au_TO_unit('E','cm-1')
        END DO
      END IF
!-----------------------------------------------------------

!----- calculation of ImaxV --------------------------------
      ImaxV = ZERO
      DO i=1,nb_ana
      DO j=i+1,nb_ana

        IntV = (ene(j)-ene(i)) * (pop(i)-pop(j)) * Mat_Aif(i,j)

        IF (IntV > ImaxV) ImaxV = IntV
      END DO
      END DO
      write(out_unit,*) 'ImaxV',ImaxV
!-----------------------------------------------------------


!----- intialisation de spectre ----------------------------
      auTOenergy = get_Conv_au_TO_unit('E',WorkingUnit=.TRUE.)
      pas = width / TEN
      n = (emax-emin)/pas
      CALL alloc_NParray(spectre,[n],"spectre",name_sub)
      spectre(:) = ZERO
!-----------------------------------------------------------

      CALL file_open(file_intensity,nio_int)

      DO i=1,nb_ana
      DO j=i+1,nb_ana

        Ai= ABC(1,i)
        Bi= ABC(2,i)
        Ci= ABC(3,i)
        Aj= ABC(1,j)
        Bj= ABC(2,j)
        Cj= ABC(3,j)

        e0 = (ene(j)-ene(i))
        IntV = (ene(j)-ene(i)) * (pop(i)-pop(j)) * Mat_Aif(i,j)
!       write(out_unit,*) 'i,j,e0,Intv',i,j,e0,Intv,ImaxV

        IF (e0 < emax .AND. IntV > ImaxV/TEN**4) THEN
          write(out_unit,*) 'sub_spectreRV',i,j,e0
          CALL calc_Int(Imax,ene(i),Ai,Bi,Ci,ene(j),Aj,Bj,Cj,           &
                        ene(1),Mat_Aif(i,j),Temp,                       &
                        Q,Jmax,const_phys,                              &
                        width,emin,pas,spectre,n,nio_int)
        END IF
      END DO
      END DO

      CALL file_close(file_intensity)

      CALL file_open(file_spectrum,nio)

      e =   emin * auTOenergy
      pas = pas  * auTOenergy
      DO i=1,n
        write(nio,*) e,spectre(i)
        e = e + pas
      END DO
      CALL file_close(file_spectrum)

      CALL dealloc_NParray(spectre,"spectre",name_sub)


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      end subroutine sub_spectreRV
!================================================================
!
!     calculation of the intensity
!
!================================================================
      SUBROUTINE calc_Int(Imax,Evi,Ai,Bi,Ci,Evf,Af,Bf,Cf,               &
                          E0,Svif,Temp,Q,Jmax,const_phys,               &
                          width,emin,pas,spectre,n,nio_int)
      USE EVR_system_m
      USE mod_Constant, only: constant
      IMPLICIT NONE

!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

      integer :: nio

      integer :: J,Jmax,K,Jf,Jop,Kop

      real (kind=Rkind) :: Ei,Evi,Pop_i,Ai,Bi,Ci
      real (kind=Rkind) :: Ef,Evf,Pop_f,Af,Bf,Cf
      real (kind=Rkind) :: Q,Temp,Svif,E0,Iif,Imax

      real (kind=Rkind) :: width,emin,pas
      integer       :: n
      real (kind=Rkind) :: spectre(n)
      integer       :: nio_int

!----- functions -------------------
      real (kind=Rkind) :: pop2_i,ene_VR
!----- functions -------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING calc_Int'
        write(out_unit,*) 'Evi,Ai,Bi,Ci',Evi,Ai,Bi,Ci
        write(out_unit,*) 'Evf,Af,Bf,Cf',Evf,Af,Bf,Cf
        write(out_unit,*) 'Imax,E0,Svif,Temp,Q',Imax,E0,Svif,Temp,Q
        write(out_unit,*) 'Jmax',Jmax
        write(out_unit,*) 'width,emin,pas (cm-1)i,n',width,emin,pas,n
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      DO J=0,Jmax
      DO Jf=J-1,J+1
        IF (Jf <0) CYCLE
        DO K=-Jf,Jf


          IF (K <-J .OR. K > J) CYCLE
!         write(out_unit,*) 'J,K,Jf',J,K,Jf

          Ei = ene_VR(Evi,Ai,Bi,Ci,J ,K)
          Ef = ene_VR(Evf,Af,Bf,Cf,Jf,K)

          Pop_i = pop2_i(Ei,E0,Temp)
          Pop_f = pop2_i(Ef,E0,Temp)


          Iif = Svif * const_phys%convIDif * (Ef-Ei)*(Pop_i-Pop_f)/Q
          Iif = Iif * real(J,kind=Rkind)*(real(J,kind=Rkind)+ONE)

          IF (Iif > Imax/TEN**4)                                        &
            CALL add_transition((Ef-Ei),Iif,width,emin,pas,spectre,n)

          IF (Iif > Imax/TEN**4)                                        &
          write(nio_int,*) (Ef-Ei)*const_phys%auTOenergy,Iif


        END DO

      END DO
      END DO

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'END calc_Int'
      END IF
!-----------------------------------------------------------

      end subroutine calc_Int
!================================================================
!
!     Imax
!
!================================================================
      SUBROUTINE calc_Imax(Imax,Evi,Ai,Bi,Ci,Evf,Af,Bf,Cf,              &
                           Svif,Temp,Q,Jmax,const_phys)
      USE EVR_system_m
      USE mod_Constant, only: constant
      IMPLICIT NONE

!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

      integer :: J,Jmax,K,Jf,Jop,Kop

      real (kind=Rkind) :: Ei,Evi,Pop_i,Ai,Bi,Ci
      real (kind=Rkind) :: Ef,Evf,Pop_f,Af,Bf,Cf
      real (kind=Rkind) :: Q,Temp,Svif,E0,Iif,Imax

!----- functions -------------------
      real (kind=Rkind) :: pop2_i,ene_VR
!----- functions -------------------

      E0 = ene_VR(Evi,Ai,Bi,Ci,0,0)
      Jop = 0
      kop = 0
      DO J=0,Jmax
        Jf = J
        K = 0

        Ei = ene_VR(Evi,Ai,Bi,Ci,J ,K)
        Ef = ene_VR(Evf,Af,Bf,Cf,Jf,K)

        Pop_i = pop2_i(Ei,E0,Temp)
        Pop_f = pop2_i(Ef,E0,Temp)


        Iif = Svif * const_phys%convIDif * (Ef-Ei)*(Pop_i-Pop_f)/Q
        Iif = Iif * real(J,kind=Rkind)*(real(J,kind=Rkind)+ONE)

        IF (Iif > Imax) Imax = Iif
        IF (Iif > Imax) Jop = Jf

      END DO


      write(out_unit,*) 'Imax,J,K',Imax,Jop,Kop

      end subroutine calc_Imax
!================================================================
!
!     Rotational energy + Evi
!
!================================================================
      FUNCTION ene_VR(Ev,A,B,C,J,K)
      USE EVR_system_m
      IMPLICIT NONE

      real (kind=Rkind) :: ene_VR
      real (kind=Rkind) :: Ev,A,B,C
      integer       :: J,K
      real (kind=Rkind) :: RJ,RK,BC,EVR


      BC = sqrt(B*C)
      RJ = real(J,kind=Rkind)
      RK = real(K,kind=Rkind)


!     for a prolate symetric top
      EVR = Ev + BC * RJ*(RJ+ONE) + (A-BC)*RK*RK
!     for a oblate symetric top
!     EVR = Ev + BC * RJ*(RJ+ONE) + (A-BC)*RK*RK


!     write(out_unit,*) 'EVR,Ev,A,BC,RJ,RK',EVR,Ev,A,BC,RJ,RK
      ene_VR = EVR

      end function ene_VR
!================================================================
!
!     Partition function : part_func(ene,nb_aie,Temp,k)
!
!================================================================

!     Q : partial Partition function
!     P : population sums over rotational levels
      SUBROUTINE calc_Q_VR(Q,P,E0,Ev,A,B,C,Temp,Jmax,const_phys)
      USE EVR_system_m
      USE mod_Constant, only: constant
      IMPLICIT NONE


!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys


      integer :: J,Jmax,K

      real (kind=Rkind) :: E,Ev,P,A,B,C
      real (kind=Rkind) :: Q,Temp,E0

!----- functions -------------------
      real (kind=Rkind) :: pop2_i,ene_VR
!----- functions -------------------




      P = 0
      DO J=0,Jmax
      DO K=-J,J

        E = ene_VR(Ev,A,B,C,J ,K)

        P = P + (real(J,kind=Rkind)+ONE)*real(J,kind=Rkind)*pop2_i(E,E0,Temp)


      END DO
      END DO
      Q = Q + P

      end subroutine calc_Q_VR
      FUNCTION part_func(ene,nb_aie,Temp)
      USE EVR_system_m
      use mod_Constant, only:real_wu,convRWU_TO_R_WITH_WorkingUnit
      implicit none

      real (kind=Rkind) :: part_func


!----- variables pour la namelist analyse ----------------------------
      real (kind=Rkind) :: Temp

!----- eigenvalues ---------------------------------------------------
      integer           :: nb_aie
      real (kind=Rkind) :: zpe,ene(nb_aie)

!----- working variables ---------------------------------------------
      integer           :: i
      real (kind=Rkind) :: Q,Etemp
      TYPE(REAL_WU)     :: RWU_Temp,RWU_E,RWU_DE

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING part_func'
        write(out_unit,*) 'nb_aie',nb_aie
        write(out_unit,*) 'ene(1)',ene(1)
        write(out_unit,*) 'Temperature (K) :',Temp
        write(out_unit,*)
      END IF
!-----------------------------------------------------------
      IF (Temp == ZERO) THEN
        Q = ONE
      ELSE
        RWU_Temp = REAL_WU(Temp,'°K','E')
        Etemp    = convRWU_TO_R_WITH_WorkingUnit(RWU_Temp)

        zpe = minval(ene)

        Q = ZERO
        DO i=1,nb_aie
          Q = Q + exp(-(ene(i)-zpe)/Etemp)
        END DO
      END IF


      part_func = Q

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Q,Temp',Q,Temp
        write(out_unit,*) 'END part_func'
      END IF
!----------------------------------------------------------


      end function part_func
!================================================================
!
!     population in the level i
!
!================================================================
      FUNCTION pop2_i(enei,ene0,Temp)
      USE EVR_system_m
      use mod_Constant, only: real_wu,convRWU_TO_R_WITH_WorkingUnit
      IMPLICIT NONE

      real (kind=Rkind) :: pop2_i


!----- variables pour la namelist analyse ----------------------------
      real (kind=Rkind) :: Temp

!----- eigenvalues ---------------------------------------------------
      real (kind=Rkind) :: enei,ene0

!----- working variables ---------------------------------------------
      real (kind=Rkind) :: Q,Etemp
      TYPE(REAL_WU)     :: RWU_Temp,RWU_E,RWU_DE

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING pop2_i'
        write(out_unit,*) 'Temperature (K) :',Temp
      END IF
!-----------------------------------------------------------
      IF (Temp == ZERO) THEN
        IF (enei == ene0) THEN
            pop2_i = ONE
        ELSE
            pop2_i = ZERO
        END IF
      ELSE
        RWU_Temp = REAL_WU(Temp,'°K','E')
        Etemp    = convRWU_TO_R_WITH_WorkingUnit(RWU_Temp)

        pop2_i = exp(-(enei-ene0)/Etemp)
      END IF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'pop2_i,Temp',Q,Temp
        write(out_unit,*) 'END pop2_i'
      END IF
!----------------------------------------------------------
      end function pop2_i

!================================================================
!
!     spectrum
!
!================================================================
      SUBROUTINE sub_spectre(Mat_Aif,ene,nb_ana,                        &
                             Ewidth,emin,emax,file_spectrum,Min_relativeI0)
      USE EVR_system_m
      use mod_Constant, only: REAL_WU,get_Conv_au_TO_unit,RWU_Write
      IMPLICIT NONE

      integer,           intent(in) :: nb_ana
      real (kind=Rkind), intent(in) :: Mat_Aif(nb_ana,nb_ana)
      real (kind=Rkind), intent(in) :: ene(nb_ana)
      real (kind=Rkind), intent(in) :: Ewidth,emin,emax
      real (kind=Rkind), intent(in) :: Min_relativeI0
      TYPE (File_t) :: file_spectrum

      integer           :: nio
      integer           :: n
      real (kind=Rkind), allocatable :: spectre(:)
      real (kind=Rkind) :: pas



      real (kind=Rkind) :: e,e0,I0,limit_I0 = ZERO
      real (kind=Rkind) :: auTOenergy

      integer           :: i,j

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      auTOenergy = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)
      limit_I0 = maxval(Mat_Aif)*Min_relativeI0

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING sub_spectre'
        write(out_unit,*) 'nb_ana',nb_ana
        write(out_unit,*) 'spectrum Ewidth: ',                         &
          RWU_Write(REAL_WU(Ewidth,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emin: ',                           &
          RWU_Write(REAL_WU(emin,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'spectrum emax: ',                      &
          RWU_Write(REAL_WU(emax,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.)

        write(out_unit,*) 'limit_I0',limit_I0
      END IF
!-----------------------------------------------------------

      pas      = Ewidth / TEN
      n        = (emax-emin)/pas
 


      CALL alloc_NParray(spectre,[n],"spectre","sub_spectre")

      spectre(:) = ZERO
      DO i=1,nb_ana
      DO j=i+1,nb_ana
        I0 = Mat_Aif(i,j)
        IF (I0 > limit_I0) THEN
          e0 = (ene(j)-ene(i))
          IF (e0 < emax) THEN
            CALL add_transition(e0,I0,Ewidth,emin,pas,spectre,n)
            write(out_unit,11) 'i j levels: ',i,j,' Ej-Ei: ',                     &
              RWU_Write(REAL_WU(e0,'au','E'),WithUnit=.TRUE.,WorkingUnit=.FALSE.), &
              ' Intensity:',I0/TEN**3,' (km.mol-1)'
 11         format(a,i0,1x,i0,a,a,a,f20.10,a)
          END IF
        END IF
      END DO
      END DO


      CALL file_open(file_spectrum,nio)

      e = emin*auTOenergy
      DO i=1,n
        write(nio,*) e,spectre(i)
        e = e + pas*auTOenergy
      END DO

      CALL file_close(file_spectrum)


      CALL dealloc_NParray(spectre,"spectre","sub_spectre")

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END sub_spectre'
      END IF
!-----------------------------------------------------------

      end subroutine sub_spectre
!================================================================
!
!     add a transition to the spectrum
!
!================================================================
      SUBROUTINE add_transition(e0,I0,a,emin,pas,spectre,n)
      USE EVR_system_m
      IMPLICIT NONE


      real (kind=Rkind) :: e0,I0,a
      real (kind=Rkind) :: emin,pas
      integer           :: n
      real (kind=Rkind) :: spectre(n)


      real (kind=Rkind) :: e,b
      integer           :: i,imin,imax

      real (kind=Rkind) :: func_lorentzian

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING add_transition'
        write(out_unit,*) 'e0 I0',e0,I0
        write(out_unit,*)
      END IF
!-----------------------------------------------------------


      imin = (e0-30._Rkind*a-emin)/pas+1
      imin = max(1,imin)
      imax = (e0+30._Rkind*a-emin)/pas+1
      imax = min(n,imax)


      DO i=imin,imax
        e = emin + real(i-1,kind=Rkind)*pas
        b = I0 * func_lorentzian(e,e0,a*HALF)
!       IF ( b > spectre(i) ) spectre(i) = spectre(i) + b
        spectre(i) = spectre(i) + b
      END DO

!-----------------------------------------------------------
      IF (debug) THEN
!       write(out_unit,*)
!       e = emin
!       DO i=1,n
!         write(out_unit,*) e,spectre(i)
!         e = e + pas
!       END DO
!       write(out_unit,*)
        write(out_unit,*) 'END add_transition'
      END IF
!-----------------------------------------------------------


      end subroutine add_transition
!================================================================
!
!     Function lorentzian
!     Inte[-inf,+inf] func_lorentzian(x,x0,a) = 1.
!
!================================================================
      FUNCTION func_lorentzian(x,x0,a)
      USE EVR_system_m
      IMPLICIT NONE


      real (kind=Rkind) :: func_lorentzian
      real (kind=Rkind) :: x,x0,a

      func_lorentzian = ONE/(a*pi) * ONE/(ONE+((x-x0)/a)**2)

      end function func_lorentzian
