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
MODULE mod_ana_psi
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sub_analyze_tab_psi,sub_analyze_psi
  PUBLIC :: Channel_weight,Channel_weight_contracHADA
  PUBLIC :: norm2_psi,renorm_psi_With_norm2,renorm_psi
  PUBLIC :: string_t

  TYPE string_t
    character (len=:), allocatable :: str
  END TYPE string_t
CONTAINS

!================================================================
!
!     sub_analyze_psi :
!      ....
!     Rho1D, Rho2D
!
!================================================================
SUBROUTINE sub_analyze_tab_Psi(tab_psi,ana_psi,adia,Write_psi)
  USE EVR_system_m
  USE mod_psi_set_alloc
  USE mod_type_ana_psi
  IMPLICIT NONE

  TYPE (param_psi),     intent(inout)        :: tab_psi(:)
  TYPE (param_ana_psi), intent(inout)        :: ana_psi
  logical,              intent(in)           :: adia
  logical,              intent(in), optional :: Write_psi


!-- working parameters --------------------------------
  integer           :: i

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_analyze_tab_Psi'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*)
    flush(out_unit)
   END IF
!-------------------------------------------------------

  IF (present(Write_psi)) THEN
    DO i=1,size(tab_psi)
      ana_psi%num_psi = i
      ana_psi%Ene     = real(tab_psi(i)%CAvOp,kind=Rkind)
      CALL sub_analyze_psi(tab_psi(i),ana_psi,adia=adia,Write_psi=Write_psi)
    END DO
  ELSE
    DO i=1,size(tab_psi)
      ana_psi%num_psi = i
      ana_psi%Ene     = real(tab_psi(i)%CAvOp,kind=Rkind)
      CALL sub_analyze_psi(tab_psi(i),ana_psi,adia=adia)
    END DO
  END IF


!----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
!----------------------------------------------------------
END SUBROUTINE sub_analyze_tab_Psi

SUBROUTINE sub_analyze_psi(psi,ana_psi,adia,Write_psi,PsiAna)
  USE EVR_system_m
  USE mod_Constant
  USE mod_psi_set_alloc
  USE mod_type_ana_psi
  USE mod_psi_io
  USE mod_psi_B_TO_G
  IMPLICIT NONE

  !----- variables for the WP -------------------------------
  TYPE (param_psi),               intent(inout)           :: psi
  TYPE (param_ana_psi),           intent(inout)           :: ana_psi
  logical,                        intent(in)              :: adia
  logical,                        intent(in),    optional :: Write_psi
  character (len=:), allocatable, intent(inout), optional :: PsiAna


  integer                           :: iE,i,j,i_be,i_bi
  real (kind=Rkind)                 :: pop,Etemp
  character(len=:), allocatable     :: lformat
  TYPE(REAL_WU)                     :: RWU_Temp,RWU_E,RWU_DE,RWU_T
  real (kind=Rkind)                 :: E,DE

  real (kind=Rkind), allocatable    :: tab_WeightChannels(:,:)
  integer                           :: Dominant_Channel(2)
  real (kind=Rkind)                 :: w,Psi_norm2

  character(len=:), allocatable     :: info
  character(len=:), allocatable     :: psi_line
  character(len=:), allocatable     :: lines

  integer                           :: nioPsi
  character (len=Line_len)          :: name_filePsi      = " "     ! name of the file

  logical :: Grid,Basis,Write_psi_loc,GridDone

  !----- dynamic allocation memory ------------------------------------
  real (kind=Rkind), allocatable :: moy_Qba(:)
  real (kind=Rkind), allocatable :: Mij(:,:,:)
!----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_analyze_psi'
  integer :: err_mem,memory
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'ana_psi%GridDone',ana_psi%GridDone
    CALL ecri_psi(psi=psi)
    flush(out_unit)
  END IF
  IF (adia .AND. .NOT. ana_psi%GridDone) STOP 'adia=t and GridDone=f'

  !write(6,*) ' in sub_analyze_psi: 0' ; flush(6)

  ! save the GridRep and BasisRep values to be able to deallocate the unused representation
  Grid  = psi%GridRep
  Basis = psi%BasisRep

  IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
    ana_psi%AvQ = .FALSE.
  END IF

  IF (present(Write_psi)) THEN
    Write_psi_loc = Write_psi .AND. ana_psi%Write_psi
  ELSE
    Write_psi_loc = ana_psi%Write_psi
  END IF

  !----------------------------------------------------------------------
  ! Boltzmann population
  IF (ana_psi%Boltzmann_pop) THEN
    IF (ana_psi%Temp > ZERO .AND. ana_psi%Part_func > ZERO) THEN
      RWU_Temp = REAL_WU(ana_psi%Temp,'°K','E')
      Etemp    = convRWU_TO_R_WITH_WorkingUnit(RWU_Temp)
      pop = exp(-(ana_psi%Ene-ana_psi%ZPE)/Etemp) / ana_psi%Part_func
    ELSE
      IF (abs(ana_psi%Ene-ana_psi%ZPE) < ONETENTH**10) THEN
        pop = ONE
      ELSE
        pop = ZERO
      END IF
    END IF
  END IF

  !----------------------------------------------------------------------
  ! tab_WeightChannels
  IF (ana_psi%propa) THEN
    !write(6,*) ' in sub_analyze_psi0: propa' ; flush(6)

    RWU_E  = REAL_WU(ana_psi%Ene,'au','E')
    E      = convRWU_TO_R_WITH_WritingUnit(RWU_E)

    RWU_T  = REAL_WU(ana_psi%T,'au','t')

    IF (adia) THEN
      CALL Channel_weight(tab_WeightChannels,psi,GridRep=.TRUE.,BasisRep=.FALSE.)
      CALL SET_string(info,'#WPadia ',TO_string(ana_psi%num_psi))
    ELSE
      CALL Channel_weight(tab_WeightChannels,psi,GridRep=.FALSE.,BasisRep=.TRUE.)
      CALL SET_string(info,'#WP ',TO_string(ana_psi%num_psi))
    END IF
    Psi_norm2 = sum(tab_WeightChannels)

    CALL SET_string(psi_line,'norm^2-WP ',info,' ',  &
            RWU_Write(RWU_T,WithUnit=.TRUE.,WorkingUnit=.FALSE.))
     !write(6,*) ' in sub_analyze_psi: 1' ; flush(6)

    ! add the energy
    iE = int(log10(abs(E)+ONETENTH**8)) ! add 1e-8 to avoid zero
    GridDone = ana_psi%GridDone
    IF (15-iE < 5 .OR. 7-iE < 0) THEN
      CALL modif_ana_psi(ana_psi,EFormat='e15.5')
    ELSE
      CALL modif_ana_psi(ana_psi,EFormat='f' // TO_string(15-iE) // '.' // TO_string(7-iE) )
    END IF
    ana_psi%GridDone = GridDone

    !write(6,*) E,iE,'ana_psi%Eformat: ',ana_psi%Eformat
    !write(6,*) ' in sub_analyze_psi: 2' ; flush(6)

    CALL ADD_TO_string(psi_line,' ',TO_string(E,Rformat=ana_psi%Eformat))
    !write(6,*) ' in sub_analyze_psi: 2.5' ; flush(6)

    ! add the field (if necessary)
    IF (ana_psi%With_field) THEN
      CALL ADD_TO_string(psi_line,' ',TO_string(ana_psi%field,Rformat='f8.5'))
    END IF

    CALL ADD_TO_string(psi_line,' ',TO_string(Psi_norm2,Rformat='f10.7'))
    DO i_be=1,psi%nb_be
      CALL ADD_TO_string(psi_line,' ',TO_string(tab_WeightChannels(:,i_be),Rformat='f10.7'))
    END DO

    write(out_unit,*) psi_line
    IF (present(PsiAna)) CALL ADD_TO_string(PsiAna,NEW_LINE('nl')) ! first line

    IF (ana_psi%Coherence > 0) THEN
      CALL alloc_NParray(Mij,[2, psi%nb_bi*psi%nb_be, psi%nb_bi*psi%nb_be],'Mij',name_sub)

      CALL sub_Rhoi_Rhoj_Over_Rho(psi,Mij,ana_psi)

      CALL SET_string(psi_line,'Mij ',info,' ',TO_string(ana_psi%T,Rformat='f12.2'))

      DO i=1,size(Mij(1,:,1))
      DO j=i+1,size(Mij(1,:,1))

        IF (present(PsiAna)) THEN
          CALL ADD_TO_string(PsiAna,'M-',TO_string(i),'-',TO_string(j),' ',info,' ', &
                             TO_string(ana_psi%T,Rformat='f12.2'),': ')
          CALL ADD_TO_string(PsiAna,TO_string(Mij(1:2,i,j),Rformat='f12.8'), &
                             NEW_LINE('nl'))
        ELSE
          write(out_unit,*) 'M-' // TO_string(i) // '-' //                  &
                            TO_string(j) // ' ' // info // ' ' //            &
                            TO_string(ana_psi%T,Rformat='f12.2') // ': ' //  &
                            TO_string(Mij(1,i,j),Rformat='f12.8'),' ',       &
                            TO_string(Mij(2,i,j),Rformat='f12.8')
        END IF
      END DO
      END DO
      CALL dealloc_NParray(Mij,'Mij',name_sub)

    END IF
    !write(6,*) ' in sub_analyze_psif: propa' ; flush(6)

  ELSE ! not propa
    CALL Channel_weight(tab_WeightChannels,psi,                         &
                        GridRep=.FALSE.,BasisRep=.TRUE.,                &
                        Dominant_Channel=Dominant_Channel)

    RWU_E  = REAL_WU(ana_psi%Ene,'au','E')
    RWU_DE = REAL_WU(ana_psi%Ene-ana_psi%ZPE,'au','E')
    E      = convRWU_TO_R_WITH_WritingUnit(RWU_E)
    DE     = convRWU_TO_R_WITH_WritingUnit(RWU_DE)

    IF (present(PsiAna)) THEN
      !$OMP  CRITICAL (sub_analyze_psi_CRIT)
      CALL SET_string(PsiAna,'lev: ',TO_string(ana_psi%num_psi),' ', & 
                         TO_string(Dominant_Channel(2)),' ',TO_string(psi%convAvOp))
      CALL ADD_TO_string(PsiAna,' ',TO_string([E,DE,pop],Rformat=EneIO_format))
      !$OMP  END CRITICAL (sub_analyze_psi_CRIT)
    END IF                        
    IF (ana_psi%AvQ) THEN
      CALL alloc_NParray(moy_Qba,[2*Psi%BasisnD%ndim],"moy_Qba",name_sub)
      CALL sub_Qmoy(psi,moy_Qba,ana_psi)

      IF (present(PsiAna)) THEN
        !$OMP  CRITICAL (sub_analyze_psi_CRIT)
        CALL ADD_TO_string(PsiAna,' ',TO_string(moy_Qba,Rformat=EneIO_format))
        !$OMP  END CRITICAL (sub_analyze_psi_CRIT)
      ELSE
        IF (ana_psi%num_psi < 10000 .AND. Dominant_Channel(2) < 10000) THEN
          CALL SET_string(lformat,'("lev: ",i4,i4,l3,',TO_string(3+size(moy_Qba)),      &
                                  '(1x,',trim(adjustl(EneIO_format)),'))')
        ELSE
          CALL SET_string(lformat,'("lev: ",i0,i0,l3,',TO_string(3+size(moy_Qba)),      &
                                  '(1x,',trim(adjustl(EneIO_format)),'))')
        END IF
        write(out_unit,lformat) ana_psi%num_psi,Dominant_Channel(2),psi%convAvOp,E,DE,pop,moy_Qba(:)
      END IF

      CALL dealloc_NParray(moy_Qba,"moy_Qba",name_sub)
    ELSE

      IF (ana_psi%num_psi < 10000 .AND. Dominant_Channel(2) < 10000) THEN
        CALL SET_string(lformat,'("lev: ",i4,i4,l3,3(1x,',trim(adjustl(EneIO_format)),'))' )
      ELSE
        CALL SET_string(lformat,'("lev: ",i0,i0,l3,3(1x,',trim(adjustl(EneIO_format)),'))' )
      END IF

      IF (.NOT. present(PsiAna)) write(out_unit,lformat) ana_psi%num_psi,Dominant_Channel(2),psi%convAvOp,E,DE,pop
    END IF

    CALL SET_string(info," ",TO_string(E,"f12.6" ),' : ')
  END IF
  IF (present(PsiAna)) THEN
    !$OMP  CRITICAL (sub_analyze_psi_CRIT)
    CALL ADD_TO_string(PsiAna,NEW_LINE('nl'))
    !$OMP  END CRITICAL (sub_analyze_psi_CRIT)
  END IF

  !----------------------------------------------------------------------

  !----------------------------------------------------------------------

  IF (psi%nb_bi > 1 .AND. .NOT. ana_psi%propa) THEN
    IF (present(PsiAna)) THEN
      !$OMP  CRITICAL (sub_analyze_psi_CRIT)
      CALL ADD_TO_string(PsiAna,'% HACbis: ',&
        TO_string(tab_WeightChannels(:,1)*TEN**2,Rformat='f4.0',max_col=10),NEW_LINE('nl'))
      !$OMP  END CRITICAL (sub_analyze_psi_CRIT)
    ELSE
      CALL SET_string(lformat,'("% HAC: ",',TO_string(psi%nb_bi),"(1x,f4.0) )" )
      write(out_unit,lformat) (tab_WeightChannels(i_bi,1)*TEN**2,i_bi=1,psi%nb_bi)
      flush(out_unit)
    END IF
  END IF

  IF (.NOT. adia) THEN 
    IF (present(PsiAna)) THEN 
      CALL calc_1Dweight(psi,ana_psi,tab_WeightChannels,20,real(ana_psi%num_psi,kind=Rkind),info,.TRUE.,PsiAna)
    ELSE
      CALL calc_1Dweight(psi,ana_psi,tab_WeightChannels,20,real(ana_psi%num_psi,kind=Rkind),info,.TRUE.)
    END IF
  END IF

  IF (.NOT. adia .AND. .NOT. ana_psi%propa) THEN 
    IF (present(PsiAna)) THEN 
      CALL calc_MaxCoef_psi(psi,ana_psi%T,info,PsiAna)
    ELSE
      CALL calc_MaxCoef_psi(psi,ana_psi%T,info,lines)
      write(out_unit,*) lines
      flush(out_unit)
      deallocate(lines)
    END IF
  END IF

  IF (ana_psi%propa) THEN
    IF (present(PsiAna)) THEN
      CALL psi_Qba_ie_psi(ana_psi%T,psi,ana_psi,info,PsiAna)
    ELSE
      CALL psi_Qba_ie_psi(ana_psi%T,psi,ana_psi,info,lines)
      IF (allocated(lines)) THEN
        write(out_unit,*) lines
        flush(out_unit)
        deallocate(lines)
      END IF
    END IF
  END IF

  IF (allocated(psi%BasisnD%nDindG%nDsize)) THEN
    IF (present(PsiAna)) THEN
      CALL Rho1D_Rho2D_psi(psi,ana_psi,adia,PsiAna)
    ELSE
      CALL Rho1D_Rho2D_psi(psi,ana_psi,adia)
    END IF
    CALL write1D2D_psi(psi,ana_psi,adia)
  ELSE
    IF (present(PsiAna)) THEN
      !$OMP CRITICAL (sub_analyze_psi_CRIT)
      CALL ADD_TO_string(PsiAna, &
        ' WARNING Rho1D or Rho2D or 1Dcut or 2Dcut are not possible!',NEW_LINE('nl'))
      !$OMP  END CRITICAL (sub_analyze_psi_CRIT)
    ELSE
      write(out_unit,*) ' WARNING Rho1D or Rho2D or 1Dcut or 2Dcut are not possible!'
    END IF
  END IF

  !---------------------------------------------------------------------------
  IF (Write_psi_loc) THEN

    IF (string_IS_empty(ana_psi%file_Psi%name)) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' The file name in "file_Psi%name" is empty !'
      STOP 'ERROR in sub_analyze_psi: The file name in "file_Psi%name" is empty'
    END IF

    ! write one the grid
    IF (ana_psi%Write_psi2_Grid) THEN
      name_filePsi = ana_psi%file_Psi%name
      IF (adia) THEN
        name_filePsi = trim(name_filePsi) // '_GridAdiaPsi2'
      ELSE
        name_filePsi = trim(name_filePsi) // '_GridPsi2'
      END IF

      IF (ana_psi%propa) name_filePsi = trim(name_filePsi) // '-' // TO_string(ana_psi%num_psi)

      IF (ana_psi%propa .AND. ana_psi%T == ZERO) THEN
        IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi)
      ELSE
        IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi,append=.TRUE.)
      END IF

      IF(MPI_id==0) CALL ecri_psi(T=ana_psi%T,psi=psi,nioWP=nioPsi,     &
                         ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.,     &
                         ecri_psi2=.TRUE.)
      IF(MPI_id==0) close(nioPsi) ! CALL file_close cannot be used
    END IF

    IF (.NOT. adia .AND. ana_psi%Write_psi_Grid) THEN
      name_filePsi = trim(ana_psi%file_Psi%name) // '_GridPsi'
      IF (ana_psi%propa) name_filePsi = trim(name_filePsi) // '-' // TO_string(ana_psi%num_psi)

      IF (ana_psi%propa .AND. ana_psi%T == ZERO) THEN
        IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi)
      ELSE
        IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi,append=.TRUE.)
      END IF

      IF(MPI_id==0) CALL ecri_psi(T=ana_psi%T,psi=psi,nioWP=nioPsi,     &
                         ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.,     &
                         ecri_psi2=.FALSE.)
      IF(MPI_id==0) close(nioPsi) ! CALL file_close cannot be used
    END IF
    !---------------------------------------------------------------------------

    IF (.NOT. adia .AND. ana_psi%Write_psi2_Basis) THEN
      name_filePsi = trim(ana_psi%file_Psi%name) // '_BasisPsi2'
      IF (ana_psi%propa) name_filePsi = trim(name_filePsi) // '-' // TO_string(ana_psi%num_psi)

      IF (ana_psi%propa) THEN
        IF (ana_psi%T == ZERO) THEN
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi)
        ELSE
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi,append=.TRUE.)
        END IF
        IF(MPI_id==0) CALL ecri_psi(T=ana_psi%T,psi=psi,nioWP=nioPsi,   &
                           ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.,   &
                           ecri_psi2=.FALSE.)
      ELSE
        IF (ana_psi%num_psi == 1) THEN
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi)
        ELSE
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi,append=.TRUE.)
        END IF
        IF(MPI_id==0) CALL ecri_psi(psi=psi,nioWP=nioPsi,               &
                           ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.,   &
                           ecri_psi2=.FALSE.)
      END IF
      IF(MPI_id==0) close(nioPsi) ! CALL file_close cannot be used
    END IF

    IF (.NOT. adia .AND. ana_psi%Write_psi_Basis) THEN
      name_filePsi = trim(ana_psi%file_Psi%name) // '_BasisPsi'
      IF (ana_psi%propa) name_filePsi = trim(name_filePsi) // '-' // TO_string(ana_psi%num_psi)

      IF (ana_psi%propa) THEN
        IF (ana_psi%T == ZERO) THEN
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi)
        ELSE
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi,append=.TRUE.)
        END IF
        IF(MPI_id==0) THEN 
          CALL Write_Psi_nDBasis(Psi,nioPsi,iPsi=ana_psi%num_psi, &
                                 epsi=ZERO,lformated=.TRUE.,version=2,T=ana_psi%T)
          !CALL ecri_psi(T=ana_psi%T,psi=psi,nioWP=nioPsi,   &
          !                 ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.,   &
          !                 ecri_psi2=.FALSE.)
        END IF
      ELSE
        IF (ana_psi%num_psi == 1) THEN
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi)
        ELSE
          IF(MPI_id==0) CALL file_open2(name_filePsi,nioPsi,append=.TRUE.)
        END IF
        IF(MPI_id==0) THEN 
          CALL Write_Psi_nDBasis(Psi,nioPsi,iPsi=ana_psi%num_psi, &
                epsi=ZERO,lformated=.TRUE.,version=0)
          !CALL ecri_psi(psi=psi,nioWP=nioPsi,   &
          !                 ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.,   &
          !                 ecri_psi2=.FALSE.)
        END IF
      END IF
      IF(MPI_id==0) close(nioPsi) ! CALL file_close cannot be used
    END IF
  END IF

!         CALL calc_DM(WP(i),max_ecri,T,info,.TRUE.)
!         C12 = WP(i)%CvecB(1)*conjg(WP(i)%CvecB(2))
!         write(out_unit,31) T,i,C12,abs(C12)
!31       format('C12',f12.1,1x,i2,3(1x,f12.6))
          !CALL calc_nDTk(WP(i),T)

  IF (allocated(tab_WeightChannels)) THEN
    CALL dealloc_NParray(tab_WeightChannels,"tab_WeightChannels",name_sub)
  END IF
  IF (allocated(info))     deallocate(info)
  IF (allocated(psi_line)) deallocate(psi_line)
  IF (allocated(lformat))  deallocate(lformat)
  IF (allocated(moy_Qba))  deallocate(moy_Qba)

  ! enable to deallocate the unsed representation.
  CALL alloc_psi(psi,BasisRep=Basis,GridRep=Grid)

  IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
  END IF

END SUBROUTINE sub_analyze_psi
!================================================================
!
!     calculation of <psi | Qact(j) | psi>
!
!================================================================
      SUBROUTINE sub_Qmoy(Psi,moy_Qba,ana_psi)
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_type_ana_psi
      USE mod_psi_B_TO_G
      IMPLICIT NONE

      TYPE (param_psi),     intent(inout)   :: Psi
      TYPE (param_ana_psi), intent(inout)   :: ana_psi

!------ means value ---------------------------------------
      real (kind=Rkind),    intent(inout)   :: moy_Qba(2*Psi%nb_act1)


!------ working variables ---------------------------------
      integer           :: i,i_q,i_ie,iqie,dnErr
      real (kind=Rkind) :: psi2,WrhonD
      real (kind=Rkind) :: x(Psi%BasisnD%ndim)
      TYPE (Type_dnS)   :: dnt
      real (kind=Rkind) :: cte(20)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING sub_Qmoy'
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      cte(:) = ZERO
      IF (allocated(ana_psi%Qtransfo_type)) CALL alloc_dnSVM(dnt,1,3)

      IF (.NOT. ana_psi%GridDone) CALL sub_PsiBasisRep_TO_GridRep(Psi)

      moy_Qba(:)  = ZERO

      DO i_q=1,Psi%nb_qa
        psi2 = ZERO
        DO i_ie=1,Psi%nb_bi*Psi%nb_be

          iqie = (i_ie-1)*Psi%nb_qa + i_q

          IF (Psi%cplx) THEN
            psi2 = psi2 + abs(Psi%CvecG(iqie))**2
          ELSE
            psi2 = psi2 + abs(Psi%RvecG(iqie))**2
          END IF
        END DO

        !- calculation of x -------------------------------
        CALL Rec_x(x,psi%BasisnD,i_q)

        IF (allocated(ana_psi%Qtransfo_type)) THEN
          DO i=1,size(x)
            CALL sub_dntf(ana_psi%Qtransfo_type(i),dnt,x(i),cte,dnErr)
            IF (dnErr /= 0) THEN
              write(out_unit,*) ' ERROR in sub_Qmoy'
              write(out_unit,*) '   ERROR in the sub_dntf call for the coordinates, iQdyn:',psi%BasisnD%iQdyn(i)
              STOP 'ERROR in sub_dntf called from sub_Qmoy'
            END IF
            x(i) = dnt%d0
          END DO
        END IF
        !write(out_unit,*) 'x',x

        !- calculation of WrhonD ------------------------------
        WrhonD = Rec_WrhonD(Psi%BasisnD,i_q)
        ! write(out_unit,*) i,'WrhonD',WrhonD

        psi2 = psi2 * WrhonD

        moy_Qba(1:Psi%nb_act1) = moy_Qba(1:Psi%nb_act1) + psi2 * x(:)
        moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) =                         &
                  moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) + psi2 * x(:)**2
      END DO

      moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) =                           &
                                moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) - &
                                              moy_Qba(1:Psi%nb_act1)**2

      IF (allocated(ana_psi%Qtransfo_type)) CALL dealloc_dnSVM(dnt)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END sub_Qmoy'
      END IF
!----------------------------------------------------------

      END SUBROUTINE sub_Qmoy

      SUBROUTINE sub_Rhoi_Rhoj_Over_Rho(Psi,Mij,ana_psi)
      USE EVR_system_m
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      USE mod_type_ana_psi
      IMPLICIT NONE

      TYPE (param_psi)    , intent(inout) :: Psi
      TYPE (param_ana_psi), intent(inout) :: ana_psi

!------ value ---------------------------------------
      real (kind=Rkind),    intent(inout) :: Mij(2,Psi%nb_bi*Psi%nb_be,Psi%nb_bi*Psi%nb_be)


!------ working variables ---------------------------------
      TYPE(OldParam)       :: OldPara
      integer              :: i,j,i_q,i_ie,j_ie,iqie
      real (kind=Rkind)    :: WrhonD,Rho,Rho_i(Psi%nb_bi*Psi%nb_be)
      complex (kind=Rkind) :: CVec_bie(psi%nb_bi*psi%nb_be)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'sub_Rhoi_Rhoj_Over_Rho'
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      IF (.NOT. ana_psi%GridDone) CALL sub_PsiBasisRep_TO_GridRep(Psi)

      Mij(:,:,:)  = ZERO


      DO i_q=1,Psi%nb_qa

        WrhonD =  Rec_WrhonD(Psi%BasisnD,i_q,OldPara=OldPara)

        IF (psi%cplx) THEN
          CALL get_CVec_OF_psi_AT_ind_a(CVec_bie,psi,i_q,OldPara=OldPara)
          Rho_i(:) = conjg(CVec_bie(:))*CVec_bie(:)
        ELSE
          CALL get_RVec_OF_psi_AT_ind_a(Rho_i,psi,i_q,OldPara=OldPara)
          Rho_i(:) = Rho_i(:)**2
        END IF

        !- calculation of total density, Rho ----------------------------
        Rho = sum(Rho_i)

        DO i_ie=1,Psi%nb_bi*Psi%nb_be
        DO j_ie=1,Psi%nb_bi*Psi%nb_be

          IF (Rho > ana_psi%coherence_epsi) THEN
            Mij(1,i_ie,j_ie) = Mij(1,i_ie,j_ie) + Rho_i(i_ie)*Rho_i(j_ie) / Rho * WrhonD
          END IF
          Mij(2,i_ie,j_ie) = Mij(2,i_ie,j_ie) + Rho_i(i_ie)*Rho_i(j_ie)*WrhonD

        END DO
        END DO

      END DO


!----------------------------------------------------------
      IF (debug) THEN

        DO i_ie=1,Psi%nb_bi*Psi%nb_be
        DO j_ie=i_ie+1,Psi%nb_bi*Psi%nb_be

         write(out_unit,*) 'Mij ',i_ie,j_ie,Mij(:,i_ie,j_ie)

        END DO
        END DO

        write(out_unit,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE sub_Rhoi_Rhoj_Over_Rho

!================================================================
!
!     means of Qact1(i) psi  : <psi|Qact(i)|psi>
!
!================================================================
  SUBROUTINE psi_Qba_ie_psi(T,psi,ana_psi,info,PsiAna)
      USE EVR_system_m
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_type_ana_psi
      USE mod_psi_B_TO_G
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi),               intent(inout) :: psi
      TYPE (param_ana_psi),           intent(in)    :: ana_psi
      character (len=*),              intent(in)    :: info
      real (kind=Rkind),              intent(in)    :: T
      character (len=:), allocatable, intent(inout) :: PsiAna

!------ working variables ---------------------------------
      TYPE(OldParam)        :: OldPara

      real (kind=Rkind)     :: Qmean(psi%nb_act1)
      real (kind=Rkind)     :: Qmean_ie(psi%nb_act1,psi%nb_bi,psi%nb_be)
      real (kind=Rkind)     :: Qmean2(psi%nb_act1,psi%nb_act1)
      real (kind=Rkind)     :: Qmean2_ie(psi%nb_act1,psi%nb_act1,psi%nb_bi,psi%nb_be)
      real (kind=Rkind)     :: pop_ie(psi%nb_bi,psi%nb_be)

      real (kind=Rkind)     :: x(Psi%BasisnD%ndim)
      real (kind=Rkind)     :: xy(Psi%BasisnD%ndim,Psi%BasisnD%ndim)

      integer               :: i_qa
      integer               :: i_be,i_bi,i_ba,i_bie
      integer               :: i,j
      real (kind=Rkind)     :: WrhonD
      logical               :: psiN,norm2GridRep,norm2BasisRep
      real (kind=Rkind)     :: RVec_bie(psi%nb_bi*psi%nb_be)
      complex (kind=Rkind)  :: CVec_bie(psi%nb_bi*psi%nb_be)
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'psi_Qba_ie_psi'
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'psi'
        CALL ecri_psi(Psi=psi)
      END IF
!-----------------------------------------------------------

      IF (psi%nb_baie > psi%nb_tot) RETURN

      pop_ie(:,:)        = ZERO

      Qmean(:)           = ZERO
      Qmean_ie(:,:,:)    = ZERO

      Qmean2(:,:)        = ZERO
      Qmean2_ie(:,:,:,:) = ZERO

      IF (.NOT. ana_psi%GridDone) CALL sub_PsiBasisRep_TO_GridRep(psi)


      DO i_qa=1,psi%nb_qa

        !- calculation of WrhonD ------------------------------
        WrhonD = Rec_WrhonD(psi%BasisnD,i_qa,OldPara)

        !- calculation of x -------------------------------
        CALL Rec_x(x,psi%BasisnD,i_qa,OldPara)
        xy(:,:) = ZERO
        DO i=1,size(x)
        DO j=i,size(x)
          xy(j,i) = x(i) * x(j)
        END DO
        END DO

        IF (psi%cplx) THEN
          CALL get_CVec_OF_psi_AT_ind_a(CVec_bie,psi,i_qa,OldPara=OldPara)
          RVec_bie(:) = conjg(CVec_bie(:))*CVec_bie(:) * WrhonD
        ELSE
          CALL get_RVec_OF_psi_AT_ind_a(RVec_bie,psi,i_qa,OldPara=OldPara)
          RVec_bie(:) = RVec_bie(:)**2 * WrhonD
        END IF


        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_bie = (i_bi-1)+(i_be-1)*psi%nb_bi + 1

          Qmean(:)                 = Qmean(:)                 + x(:) * RVec_bie(i_bie)
          Qmean_ie(:,i_bi,i_be)    = Qmean_ie(:,i_bi,i_be)    + x(:) * RVec_bie(i_bie)

          Qmean2(:,:)              = Qmean2(:,:)              + xy(:,:) * RVec_bie(i_bie)
          Qmean2_ie(:,:,i_bi,i_be) = Qmean2_ie(:,:,i_bi,i_be) + xy(:,:) * RVec_bie(i_bie)

          pop_ie(i_bi,i_be)        = pop_ie(i_bi,i_be)        + RVec_bie(i_bie)

        END DO
        END DO

      END DO

      Qmean(:)    = Qmean(:)    / sum(pop_ie)
      Qmean2(:,:) = Qmean2(:,:) / sum(pop_ie)

      DO i_be=1,psi%nb_be
      DO i_bi=1,psi%nb_bi
        IF (pop_ie(i_bi,i_be) > ONETENTH**7) THEN
          Qmean_ie(:,i_bi,i_be)    = Qmean_ie(:,i_bi,i_be)    / pop_ie(i_bi,i_be)
          Qmean2_ie(:,:,i_bi,i_be) = Qmean2_ie(:,:,i_bi,i_be) / pop_ie(i_bi,i_be)
        ELSE
          Qmean_ie(:,i_bi,i_be)    = ZERO
          Qmean2_ie(:,:,i_bi,i_be) = ZERO
        END IF
      END DO
      END DO

      DO i=1,psi%nb_act1
        !write(out_unit,11) 'T <Q',TO_string(i),'>_ie ',info,T,Qmean_ie(i,:,:)
        !write(out_unit,11) 'T <Q',TO_string(i),'>    ',info,T,Qmean(i)
  11    format(4a,' ',f0.5,100(' ',f0.4))

        CALL ADD_TO_string(PsiAna,'T <Q',TO_string(i),'>_ie ',info,' ',TO_string(T,'f0.5'))
        DO i_be=1,psi%nb_be
          CALL ADD_TO_string(PsiAna,' ',TO_string(Qmean_ie(i,:,i_be),'f0.4'))
        END DO
        CALL ADD_TO_string(PsiAna,new_line('nl'))

        CALL ADD_TO_string(PsiAna,'T <Q',TO_string(i),'>_ie ',info,' ',TO_string(T,'f0.5'), &
                          ' ',TO_string(Qmean(i),'f0.4'),new_line('nl'))
      END DO
      DO i=1,psi%nb_act1
      DO j=i,psi%nb_act1
        !write(out_unit,21) 'T <Q',TO_string(j),'*Q',TO_string(i),'>_ie ',info,T,Qmean2_ie(j,i,:,:)
        !write(out_unit,21) 'T <Q',TO_string(j),'*Q',TO_string(i),'>    ',info,T,Qmean2(j,i)

        CALL ADD_TO_string(PsiAna,'T <Q',TO_string(j),'*Q',TO_string(i),'>_ie ',info,' ',TO_string(T,'f0.5'))

        DO i_be=1,psi%nb_be
          CALL ADD_TO_string(PsiAna,' ',TO_string(Qmean2_ie(j,i,:,i_be),'f0.4'))
        END DO
        CALL ADD_TO_string(PsiAna,new_line('nl'))

        CALL ADD_TO_string(PsiAna,'T <Q',TO_string(i),'>_ie ',info,' ',TO_string(T,'f0.5'), &
                           ' ',TO_string(Qmean2(j,i),'f0.4'),new_line('nl'))
      END DO
      END DO
 21   format(6a,' ',f0.5,100(' ',f0.4))
      flush(out_unit)

      CALL dealloc_OldParam(OldPara)

    !----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'Qmean',Qmean
      DO i=1,psi%nb_act1
        write(out_unit,*) 'Qmean_ie',i,Qmean_ie(i,:,:)
      END DO
      write(out_unit,*) '<Qi*Qj>',Qmean2(:,:)
      DO i=1,psi%nb_act1
      DO j=i,psi%nb_act1
        write(out_unit,21) 'T iQbasis <Qi*Qj>_ie ',info,T,j,i,Qmean2_ie(j,i,:,:)
        write(out_unit,21) 'T iQbasis <Qi*Qj>    ',info,T,j,i,Qmean2(j,i)
      END DO
      END DO
      write(out_unit,*) 'END ',name_sub
    END IF
  END SUBROUTINE psi_Qba_ie_psi
  SUBROUTINE write1D2D_psi(psi,ana_psi,adia)
      USE EVR_system_m
      USE mod_basis
      USE mod_dnSVM
      USE mod_psi_set_alloc
      USE mod_type_ana_psi
      USE mod_psi_B_TO_G
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi),     intent(inout)  :: psi
      TYPE (param_ana_psi), intent(inout)  :: ana_psi
      logical,              intent(in)     :: adia



!----- working variables ------------------------------------------
      real (kind=Rkind) :: val,val_mini,val_psi2,val_Rpsi,val_Cpsi

      integer           :: i_qa,i_qaie,iq,jq
      integer           :: ib,jb,i_bie


      TYPE (Type_dnVec) :: Tab1D_Qact(psi%BasisnD%nb_basis)
      TYPE (Type_dnVec) :: Tab1D_psi2(psi%BasisnD%nb_basis)
      TYPE (Type_dnVec) :: Tab1D_Rpsi(psi%BasisnD%nb_basis)
      TYPE (Type_dnVec) :: Tab1D_Cpsi(psi%BasisnD%nb_basis)

      TYPE (Type_dnMat) :: Tab2D_psi2(psi%BasisnD%nb_basis,psi%BasisnD%nb_basis)
      TYPE (Type_dnMat) :: Tab2D_Rpsi(psi%BasisnD%nb_basis,psi%BasisnD%nb_basis)
      TYPE (Type_dnMat) :: Tab2D_Cpsi(psi%BasisnD%nb_basis,psi%BasisnD%nb_basis)

      integer :: nDval0(psi%BasisnD%nb_basis)
      integer :: nDval(psi%BasisnD%nb_basis)
      integer :: nDvalib(psi%BasisnD%nb_basis)
      integer :: nDval0ib(psi%BasisnD%nb_basis)

!------ working variables ---------------------------------
      TYPE (File_t)              :: file_psi
      integer                        :: nio,nqi,nqj
      character (len=:), allocatable  :: state_name

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='write1D2D_psi'
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------

      IF (.NOT.ana_psi%psi1D_Q0 .AND. .NOT.ana_psi%psi2D_Q0) RETURN
      IF (adia) RETURN

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      IF (psi%nb_act1 /= psi%BasisnD%nb_basis) THEN
        write(out_unit,*) 'WARNNING in',name_sub
        write(out_unit,*) 'nb_act1 /= BasisnD%nb_basis',psi%nb_act1,psi%BasisnD%nb_basis
        RETURN
      END IF

      IF (.NOT. ana_psi%GridDone) CALL sub_PsiBasisRep_TO_GridRep(psi)

      nDval0(:) = 0
      DO ib=1,psi%BasisnD%nb_basis
        nqi = get_nq_FROM_basis(psi%BasisnD%tab_Pbasis(ib)%Pbasis)

        CALL alloc_dnVec(Tab1D_Qact(ib),nqi)
        CALL alloc_dnVec(Tab1D_psi2(ib),nqi)
        CALL alloc_dnVec(Tab1D_Rpsi(ib),nqi)
        CALL alloc_dnVec(Tab1D_Cpsi(ib),nqi)

        Tab1D_Qact(ib)%d0(:) = psi%BasisnD%tab_Pbasis(ib)%Pbasis%x(1,:)

        ! set the table nDval0 from Qana and the Tab1D_Qact(ib)
        val_mini  = huge(ONE)
        DO iq=1,Tab1D_Qact(ib)%nb_var_vec
          val = abs(ana_psi%Qana(ib)-Tab1D_Qact(ib)%d0(iq))
          IF (val < val_mini) THEN
            val_mini = val
            nDval0(ib) = iq
          END IF
        END DO

        DO jb=1,psi%BasisnD%nb_basis
          nqj = get_nq_FROM_basis(psi%BasisnD%tab_Pbasis(jb)%Pbasis)
          CALL alloc_dnMat(Tab2D_psi2(jb,ib),nb_var_Matl=nqj,nb_var_Matc=nqi)
          CALL alloc_dnMat(Tab2D_Rpsi(jb,ib),nb_var_Matl=nqj,nb_var_Matc=nqi)
          CALL alloc_dnMat(Tab2D_Cpsi(jb,ib),nb_var_Matl=nqj,nb_var_Matc=nqi)
        END DO

      END DO
      !write(out_unit,*) 'nDval0,Qana',nDval0,ana_psi%Qana

      i_bie = 1

      IF (allocated(psi%CvecG)) THEN

        DO i_qa=1,psi%nb_qa
          i_qaie = i_qa + (i_bie-1) * psi%nb_qa
          val_psi2 = abs(psi%CvecG(i_qaie))**2
          val_Rpsi = Real(psi%CvecG(i_qaie),kind=Rkind)
          val_Cpsi = aimag(psi%CvecG(i_qaie))

          CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

          DO ib=1,psi%BasisnD%nb_basis
            nDvalib(:)   = nDval(:)
            nDvalib(ib)  = 0

            nDval0ib(:)  = nDval0(:)
            nDval0ib(ib) = 0

            IF (compare_tab(nDval0ib,nDvalib)) THEN
              Tab1D_psi2(ib)%d0(nDval(ib)) = val_psi2
              Tab1D_Rpsi(ib)%d0(nDval(ib)) = val_Rpsi
              Tab1D_Cpsi(ib)%d0(nDval(ib)) = val_Cpsi
            END IF

            DO jb=ib+1,psi%BasisnD%nb_basis

              nDvalib(jb)  = 0
              nDval0ib(jb) = 0

              IF (compare_tab(nDval0ib,nDvalib)) THEN
                Tab2D_psi2(ib,jb)%d0(nDval(ib),nDval(jb)) = val_psi2
                Tab2D_Rpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Rpsi
                Tab2D_Cpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Cpsi
              END IF

            END DO

          END DO

        END DO
      ELSE IF (allocated(psi%RvecG)) THEN

        DO i_qa=1,psi%nb_qa
          i_qaie = i_qa + (i_bie-1) * psi%nb_qa
          val_psi2 = psi%RvecG(i_qaie)**2
          val_Rpsi = psi%RvecG(i_qaie)
          val_Cpsi = ZERO

          CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

          DO ib=1,psi%BasisnD%nb_basis
            nDvalib(:)   = nDval(:)
            nDvalib(ib)  = 0

            nDval0ib(:)  = nDval0(:)
            nDval0ib(ib) = 0

            IF (compare_tab(nDval0ib,nDvalib)) THEN
              Tab1D_psi2(ib)%d0(nDval(ib)) = val_psi2
              Tab1D_Rpsi(ib)%d0(nDval(ib)) = val_Rpsi
              Tab1D_Cpsi(ib)%d0(nDval(ib)) = val_Cpsi
            END IF

            DO jb=ib+1,psi%BasisnD%nb_basis

              nDvalib(jb)  = 0
              nDval0ib(jb) = 0

              IF (compare_tab(nDval0ib,nDvalib)) THEN
                Tab2D_psi2(ib,jb)%d0(nDval(ib),nDval(jb)) = val_psi2
                Tab2D_Rpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Rpsi
                Tab2D_Cpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Cpsi
              END IF

            END DO
          END DO


        END DO

      ELSE
        write(out_unit,*) 'WARNING in ',name_sub
        write(out_unit,*) ' Impossible to write CvecG or RvecG'
        RETURN
      END IF


      IF (ana_psi%psi1D_Q0) THEN

        DO ib=1,psi%BasisnD%nb_basis

          IF (adia) THEN
            state_name = make_EVRTFileName('psiAdia1D_')
          ELSE
            state_name = make_EVRTFileName('psi1D_')
          END IF
          CALL SET_string(file_psi%name,state_name,TO_string(ana_psi%num_psi),'-',TO_string(ib))

          IF (ana_psi%propa .AND. ana_psi%T > ZERO) THEN
            CALL file_open(file_psi,nio,append=.TRUE.)
          ELSE
            CALL file_open(file_psi,nio)
          END IF

          IF (ana_psi%propa) THEN

            DO iq=1,Tab1D_Qact(ib)%nb_var_vec
              write(nio,*) ana_psi%T,iq,Tab1D_Qact(ib)%d0(iq),Tab1D_Rpsi(ib)%d0(iq),  &
                             Tab1D_Cpsi(ib)%d0(iq),Tab1D_psi2(ib)%d0(iq)
            END DO
            write(nio,*)

          ELSE

            DO iq=1,Tab1D_Qact(ib)%nb_var_vec
              write(nio,*) iq,Tab1D_Qact(ib)%d0(iq),Tab1D_Rpsi(ib)%d0(iq),  &
                             Tab1D_Cpsi(ib)%d0(iq),Tab1D_psi2(ib)%d0(iq)
            END DO
          END IF

          CALL file_close(file_psi)

        END DO
      END IF

      IF (ana_psi%psi2D_Q0) THEN

        DO ib=1,psi%BasisnD%nb_basis
        DO jb=ib+1,psi%BasisnD%nb_basis

          IF (ana_psi%adia) THEN
            state_name = make_EVRTFileName('psiAdia2D_')
          ELSE
            state_name = make_EVRTFileName('psi2D_')
          END IF
          CALL SET_string(file_psi%name,state_name,TO_string(ana_psi%num_psi), &
                          '-',TO_string(ib),'-',TO_string(jb))

          IF (ana_psi%propa .AND. ana_psi%T > ZERO) THEN
            CALL file_open(file_psi,nio,append=.TRUE.)
          ELSE
            CALL file_open(file_psi,nio)
          END IF

          IF (ana_psi%propa) THEN

            DO iq=1,Tab1D_Qact(ib)%nb_var_vec
            DO jq=1,Tab1D_Qact(jb)%nb_var_vec

              write(nio,*) ana_psi%T,iq,jq,Tab1D_Qact(ib)%d0(iq),Tab1D_Qact(jb)%d0(jq),&
                 Tab2D_Rpsi(ib,jb)%d0(iq,jq),Tab2D_Cpsi(ib,jb)%d0(iq,jq),  &
                                             Tab2D_psi2(ib,jb)%d0(iq,jq)
            END DO
            write(nio,*)
            END DO
            write(nio,*)

          ELSE

            DO iq=1,Tab1D_Qact(ib)%nb_var_vec
            DO jq=1,Tab1D_Qact(jb)%nb_var_vec

              write(nio,*) iq,jq,Tab1D_Qact(ib)%d0(iq),Tab1D_Qact(jb)%d0(jq),&
              Tab2D_Rpsi(ib,jb)%d0(iq,jq),Tab2D_Cpsi(ib,jb)%d0(iq,jq),  &
                                             Tab2D_psi2(ib,jb)%d0(iq,jq)
            END DO
            write(nio,*)
            END DO
          END IF

          CALL file_close(file_psi)

        END DO
        END DO
      END IF


      DO ib=1,psi%BasisnD%nb_basis
        CALL dealloc_dnVec(Tab1D_Qact(ib))
        CALL dealloc_dnVec(Tab1D_psi2(ib))
        CALL dealloc_dnVec(Tab1D_Rpsi(ib))
        CALL dealloc_dnVec(Tab1D_Cpsi(ib))
        DO jb=1,psi%BasisnD%nb_basis
          CALL dealloc_dnMat(Tab2D_psi2(jb,ib))
          CALL dealloc_dnMat(Tab2D_Rpsi(jb,ib))
          CALL dealloc_dnMat(Tab2D_Cpsi(jb,ib))
        END DO
      END DO

      IF (allocated(state_name)) deallocate(state_name)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------


      END SUBROUTINE write1D2D_psi
      SUBROUTINE Rho1D_Rho2D_psi(psi,ana_psi,adia,PsiAna)
      USE EVR_system_m
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      USE mod_type_ana_psi
      IMPLICIT NONE
!----- variables for the WP ----------------------------------------
      TYPE (param_psi),               intent(inout)           :: psi
      TYPE (param_ana_psi),           intent(inout)           :: ana_psi
      logical,                        intent(in)              :: adia
      character (len=:), allocatable, intent(inout), optional :: PsiAna


!------ working variables ---------------------------------
      TYPE (File_t)                  :: file_Rho
      integer                        :: i_basis_act1,j_basis_act1,nioRho
      real (kind=Rkind), allocatable :: rho2D(:,:,:,:)
      real (kind=Rkind), allocatable :: rho1D(:,:,:)


      TYPE(OldParam)        :: OldPara
      integer               :: i_qa,i_qaie
      integer               :: i_be,i_bi,i_ba,i_bie
      integer               :: i_baie,f_baie
      real (kind=Rkind)     :: WrhonD,Weight
      complex (kind=Rkind)  :: temp
      real (kind=Rkind)     :: Rtemp,Qix
      real (kind=Rkind)     :: RVec_bie(psi%nb_bi*psi%nb_be)
      complex (kind=Rkind)  :: CVec_bie(psi%nb_bi*psi%nb_be)
      real (kind=Rkind), allocatable     :: x(:)

      integer :: iq_act1,jq_act1,i,j,iq,ix,ix_TO_iQdyn
      integer :: nDval(psi%BasisnD%nDindG%ndim)
      character (len=:), allocatable  :: state_name

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Rho1D_Rho2D_psi'
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'psi%basis'
        CALL RecWrite_basis(psi%BasisnD)
        write(out_unit,*) 'psi'
        CALL ecri_psi(psi=psi)
        write(out_unit,*) 'ana_psi%ana_level',ana_psi%ana_level
        write(out_unit,*) 'ana_psi%Rho1D,ana_psi%Rho2D',ana_psi%Rho1D,ana_psi%Rho2D
      END IF
!-----------------------------------------------------------
      IF (ana_psi%ana_level > 0 .AND. (ana_psi%Rho1D .OR. ana_psi%Rho2D)) THEN
        IF (.NOT. ana_psi%GridDone) CALL sub_PsiBasisRep_TO_GridRep(psi)

        IF (ana_psi%Rho1D) THEN
          !- loop on the basis ----------------------------------
          DO i_basis_act1=1,psi%BasisnD%nb_basis
            IF (psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis%nb_basis > 1) CYCLE

            IF (present(PsiAna)) THEN
              CALL ADD_TO_string(PsiAna,'rho1D of ',   &
                   TO_string(ana_psi%num_psi),' ',TO_string(i_basis_act1))
            ELSE
              write(out_unit,*) 'rho1D of ',ana_psi%num_psi,i_basis_act1
            END IF

            IF (adia) THEN
              state_name = make_EVRTFileName('RhoAdia1D_')
            ELSE
              state_name = make_EVRTFileName('Rho1D_')
            END IF
            CALL SET_string(file_Rho%name,state_name,TO_string(ana_psi%num_psi),'-',TO_string(i_basis_act1))
            IF (ana_psi%propa .AND. ana_psi%T > ZERO) THEN
              CALL file_open(file_Rho,nioRho,append=.TRUE.)
            ELSE
              CALL file_open(file_Rho,nioRho)
            END IF

            CALL alloc_NParray(rho1D,                                   &
                 [psi%BasisnD%nDindG%nDsize(i_basis_act1),psi%nb_bi,psi%nb_be], &
                              'rho1D',name_sub)
            rho1D(:,:,:) = ZERO


            DO i_qa=1,psi%nb_qa

              CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

              !- calculation of WrhonD ------------------------------
              WrhonD   = ONE
              DO i=1,psi%BasisnD%nb_basis
                iq  = nDval(i)

                IF (allocated(ana_psi%Weight_Rho) .AND.            &
                              allocated(ana_psi%Qana_weight) ) THEN

                  DO ix=1,psi%BasisnD%tab_Pbasis(i)%Pbasis%ndim
                    ix_TO_iQdyn = psi%BasisnD%tab_Pbasis(i)%Pbasis%iQdyn(ix)
                    IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == 1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix > ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    ELSE IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == -1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix < ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    END IF
                  END DO

                END IF

                IF (WrhonD == ZERO) CYCLE

                IF (i == i_basis_act1) THEN
                  IF (ana_psi%Rho_type == 1) THEN
                    WrhonD = WrhonD * Rec_WrhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                  ELSE IF (ana_psi%Rho_type == 2) THEN
                    WrhonD = WrhonD * Rec_rhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                  END IF
                ELSE
                  WrhonD = WrhonD * Rec_WrhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                END IF
              END DO

              IF (psi%cplx) THEN
                CALL get_CVec_OF_psi_AT_ind_a(CVec_bie,psi,i_qa,OldPara=OldPara)
                RVec_bie(:) = conjg(CVec_bie(:))*CVec_bie(:) * WrhonD
              ELSE
                CALL get_RVec_OF_psi_AT_ind_a(RVec_bie,psi,i_qa,OldPara=OldPara)
                RVec_bie(:) = RVec_bie(:)**2 * WrhonD
              END IF

              DO i_be=1,psi%nb_be
              DO i_bi=1,psi%nb_bi
                i_bie = (i_bi-1)+(i_be-1)*psi%nb_bi + 1
                iq_act1 = nDval(i_basis_act1)

                rho1D(iq_act1,i_bi,i_be) = rho1D(iq_act1,i_bi,i_be) + RVec_bie(i_bie)
              END DO
              END DO

            END DO
            allocate(x(psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis%ndim))
            IF (ana_psi%propa) THEN
              DO i=1,psi%BasisnD%nDindG%nDsize(i_basis_act1)
                CALL Rec_x(x,psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis,i)
                write(nioRho,*) ana_psi%T,i,x,rho1D(i,:,:)
              END DO
              write(nioRho,*)
            ELSE
            
              DO i=1,psi%BasisnD%nDindG%nDsize(i_basis_act1)
                CALL Rec_x(x,psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis,i)
                write(nioRho,*) i,x,rho1D(i,:,:)
              END DO
            END IF
            deallocate(x)

            CALL dealloc_NParray(rho1D,'rho1D',name_sub)
            CALL file_close(file_Rho)
          END DO

        END IF

        IF (ana_psi%Rho2D) THEN

          !- loop on coordinates ----------------------------------
          DO i_basis_act1=1,psi%BasisnD%nb_basis
          DO j_basis_act1=i_basis_act1+1,psi%BasisnD%nb_basis

            IF (adia) THEN
              state_name = make_EVRTFileName('RhoAdia2D_')
            ELSE
              state_name = make_EVRTFileName('Rho2D_')
            END IF
            CALL SET_string(file_Rho%name,state_name ,TO_string(ana_psi%num_psi), &
                            '-',TO_string(i_basis_act1),'-',TO_string(j_basis_act1))

            IF (ana_psi%propa .AND. ana_psi%T > ZERO) THEN
              CALL file_open(file_Rho,nioRho,append=.TRUE.)
            ELSE
              CALL file_open(file_Rho,nioRho)
            END IF

            CALL alloc_NParray(rho2D,                                   &
                            [psi%BasisnD%nDindG%nDsize(i_basis_act1),  &
                              psi%BasisnD%nDindG%nDsize(j_basis_act1),  &
                              psi%nb_bi,psi%nb_be],'rho2D',name_sub)
            rho2D(:,:,:,:) = ZERO

            DO i_qa=1,psi%nb_qa

              CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

              !- calculation of WrhonD ------------------------------
              WrhonD   = ONE
              DO i=1,psi%BasisnD%nb_basis
                iq  = nDval(i)

                IF (allocated(ana_psi%Weight_Rho) .AND.            &
                              allocated(ana_psi%Qana_weight) ) THEN
                  DO ix=1,psi%BasisnD%tab_Pbasis(i)%Pbasis%ndim
                    ix_TO_iQdyn = psi%BasisnD%tab_Pbasis(i)%Pbasis%iQdyn(ix)

                    IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == 1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix > ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    ELSE IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == -1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix < ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    END IF
                  END DO

                  IF (WrhonD == ZERO) CYCLE
                END IF


                IF (i == i_basis_act1 .OR. i == j_basis_act1) THEN
                  IF (ana_psi%Rho_type == 1) THEN
                    WrhonD = WrhonD * Rec_WrhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                  ELSE IF (ana_psi%Rho_type == 2) THEN
                    WrhonD = WrhonD * Rec_rhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                  END IF
                ELSE
                  WrhonD = WrhonD * Rec_WrhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                END IF
              END DO

              IF (psi%cplx) THEN
                CALL get_CVec_OF_psi_AT_ind_a(CVec_bie,psi,i_qa,OldPara=OldPara)
                RVec_bie(:) = conjg(CVec_bie(:))*CVec_bie(:) * WrhonD
              ELSE
                CALL get_RVec_OF_psi_AT_ind_a(RVec_bie,psi,i_qa,OldPara=OldPara)
                RVec_bie(:) = RVec_bie(:)**2 * WrhonD
              END IF

              DO i_be=1,psi%nb_be
              DO i_bi=1,psi%nb_bi
                i_bie = (i_bi-1)+(i_be-1)*psi%nb_bi + 1

                iq_act1 = nDval(i_basis_act1)
                jq_act1 = nDval(j_basis_act1)
                rho2D(iq_act1,jq_act1,i_bi,i_be) = rho2D(iq_act1,jq_act1,i_bi,i_be) + &
                                                                 RVec_bie(i_bie)
              END DO
              END DO
            END DO

            IF (present(PsiAna)) THEN
              CALL ADD_TO_string(PsiAna,'rho2D of ',TO_string(ana_psi%num_psi),' ', &
                                 TO_string(i_basis_act1),' ',TO_string(j_basis_act1))
            ELSE
              write(out_unit,*) 'rho2D of ',ana_psi%num_psi,i_basis_act1,j_basis_act1
            END IF
            IF (ana_psi%propa) THEN
              DO i=1,psi%BasisnD%nDindG%nDsize(i_basis_act1)
              DO j=1,psi%BasisnD%nDindG%nDsize(j_basis_act1)
                write(nioRho,*) ana_psi%T,i,j,                          &
                 psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis%x(:,i),    &
                 psi%BasisnD%tab_Pbasis(j_basis_act1)%Pbasis%x(:,j),    &
                 rho2D(i,j,:,:)
              END DO
              write(nioRho,*)
              END DO
              write(nioRho,*)
            ELSE
              DO i=1,psi%BasisnD%nDindG%nDsize(i_basis_act1)
              DO j=1,psi%BasisnD%nDindG%nDsize(j_basis_act1)
                write(nioRho,*) i,j,                                    &
                 psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis%x(:,i),    &
                 psi%BasisnD%tab_Pbasis(j_basis_act1)%Pbasis%x(:,j),    &
                 rho2D(i,j,:,:)
              END DO
              write(nioRho,*)
              END DO
            END IF

            CALL dealloc_NParray(rho2D,'rho2D',name_sub)
            CALL file_close(file_Rho)
          END DO
          END DO
        END IF

      END IF

      IF (allocated(state_name)) deallocate(state_name)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE Rho1D_Rho2D_psi

!==============================================================
!
!      calc_1Dweight(psi)
!
!==============================================================
  SUBROUTINE calc_1Dweight(psi,ana_psi,tab_WeightChannels,max_1D,T,info,print_w,PsiAna)
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_type_ana_psi
      IMPLICIT NONE

      TYPE (param_psi),               intent(in)              :: psi
      TYPE (param_ana_psi),           intent(inout)           :: ana_psi
      real (kind=Rkind), allocatable, intent(in)              :: tab_WeightChannels(:,:)
      integer,                        intent(in)              :: max_1D
      real (kind=Rkind),              intent(in)              :: T ! time
      character (len=*),              intent(in)              :: info
      logical,                        intent(in)              :: print_w
      character (len=:), allocatable, intent(inout), optional :: PsiAna

      character (len=:), allocatable :: PsiAna_loc

!----- for debuging --------------------------------------------------
      integer :: err_mem
      character (len=*), parameter :: name_sub='calc_1Dweight'
      logical, parameter :: debug =.FALSE.
!      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(psi%BasisnD,write_all=.TRUE.)
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      IF (present(PsiAna)) THEN
        IF (psi%nb_baie*psi%nb_bRot == psi%nb_tot) THEN
          CALL calc_1Dweight_act1(psi,ana_psi,max_1D,T,info,print_w,PsiAna)
        END IF

        CALL calc_1Dweight_inact2n_elec(psi,ana_psi,tab_WeightChannels,max_1D,T,info,print_w,PsiAna)
      ELSE
        IF (psi%nb_baie*psi%nb_bRot == psi%nb_tot) THEN
          CALL calc_1Dweight_act1(psi,ana_psi,max_1D,T,info,print_w,PsiAna_loc)
        END IF

        CALL calc_1Dweight_inact2n_elec(psi,ana_psi,tab_WeightChannels,max_1D,T,info,print_w,PsiAna_loc)
        write(out_unit,*) PsiAna_loc
        deallocate(PsiAna_loc)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE calc_1Dweight

      SUBROUTINE calc_1Dweight_inact2n_elec(psi,ana_psi,tab_WeightChannels,&
                                            max_1D,T,info,print_w,PsiAna)
      USE EVR_system_m
      USE mod_nDindex
      USE mod_psi_set_alloc
      USE mod_type_ana_psi
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),               intent(in)     :: psi
      TYPE (param_ana_psi),           intent(in)     :: ana_psi
      real (kind=Rkind), allocatable, intent(in)     :: tab_WeightChannels(:,:)
      integer,                        intent(in)     :: max_1D
      real (kind=Rkind),              intent(in)     :: T ! time
      character (len=*),              intent(in)     :: info
      logical,                        intent(in)     :: print_w
      character (len=:), allocatable, intent(inout)  :: PsiAna

      integer          :: max_herm
      real (kind=Rkind),allocatable :: weight1D(:,:)
      integer          :: tab(psi%Basis2n%nb_basis+1)

      real (kind=Rkind),allocatable :: weight1Dact(:,:)
      real (kind=Rkind)    :: a



      integer              :: i,j,j_herm,beg_e
      integer              :: ie,ii,ib,ibie,iq,ibiq
      integer              :: max_dim
      integer, allocatable :: nDval(:)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_1Dweight_inact2n_elec'
      integer :: err_mem,memory
      logical, parameter :: debug =.FALSE.
!      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_inact2n',psi%Basis2n%nb_basis
        write(out_unit,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
        write(out_unit,*) 'nb_bie',psi%nb_bi*psi%nb_be
        !write(out_unit,*) 'tab_WeightChannels',tab_WeightChannels
        flush(out_unit)
      END IF
!-----------------------------------------------------------
      !IF (ana_psi%adia) RETURN

      IF (.NOT. allocated(tab_WeightChannels)) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'tab_WeightChannels is not allocated!!'
        write(out_unit,*) ' It should be done in sub_analyze or sub_analyze_WP_forPropa'
        STOP
      END IF

      IF (psi%Basis2n%nb_basis == 0) THEN
        max_herm = psi%nb_be-1
      ELSE
        CALL alloc_NParray(nDval,[psi%Basis2n%nb_basis],'nDval',name_sub)
        max_herm = 0
        DO i=1,psi%Basis2n%nb_basis
          IF (psi%Basis2n%tab_Pbasis(i)%Pbasis%nb > max_herm)           &
                      max_herm = psi%Basis2n%tab_Pbasis(i)%Pbasis%nb
        END DO
        max_herm = max_herm - 1
      END IF

      CALL alloc_NParray(weight1D,[psi%Basis2n%nb_basis+1,max_herm],    &
                        "weight1D",name_sub,[1,0])

      IF (psi%nb_bi > 1 .OR. psi%nb_be > 1) THEN
        weight1D(:,:) = ZERO
        i = 0
        DO ie=1,psi%nb_be
        DO ii=1,psi%nb_bi
          i = i + 1
          IF (psi%Basis2n%nb_basis > 0) THEN
            CALL calc_nDindex(psi%Basis2n%nDindB,i,nDval)
            tab(1:psi%Basis2n%nb_basis) = nDval(:)-1
            tab(psi%Basis2n%nb_basis+1) = ie
          ELSE
            tab(psi%Basis2n%nb_basis+1) = ie-1
          END IF


          DO j=1,psi%Basis2n%nb_basis+1
            j_herm = tab(j)
            weight1D(j,j_herm) = weight1D(j,j_herm) + tab_WeightChannels(ii,ie)
          END DO
        END DO
        END DO

        IF (print_w .OR. debug) THEN
          DO i=1,psi%Basis2n%nb_basis
            !write(out_unit,11) 'harm T W ',trim(info),i,T,             &
            !        weight1D(i,0:psi%Basis2n%tab_Pbasis(i)%Pbasis%nb-1)
            CALL ADD_TO_string(PsiAna,'harm T W ',trim(info),' ',TO_string(i), &
                               ' ',TO_string(T,'f15.2'),' ',                   &
                               TO_string(weight1D(i,0:psi%Basis2n%tab_Pbasis(i)%Pbasis%nb-1),'f6.3'), &
                               new_line('nl'))
          END DO
          beg_e = 1
          IF (psi%Basis2n%nb_basis == 0) beg_e = 0
          i = psi%Basis2n%nb_basis+1
          !write(out_unit,11) 'elec T W ',trim(info),i,T,               &
          !      weight1D(i,beg_e:beg_e+psi%nb_be-1)

 11       format(a,a,1x,i3,1x,f15.2,20(1x,f6.3))
          CALL ADD_TO_string(PsiAna,'elec T W ',trim(info),' ',TO_string(i),  &
                             ' ',TO_string(T,'f15.2'),' ',                    &
                             ' ',TO_string(weight1D(i,beg_e:beg_e+psi%nb_be-1),'f6.3'), &
                             new_line('nl'))
        END IF
      END IF

      CALL dealloc_NParray(weight1D,"weight1D",name_sub)
      IF (allocated(nDval)) CALL dealloc_NParray(nDval,'nDval',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


    END SUBROUTINE calc_1Dweight_inact2n_elec

    SUBROUTINE calc_MaxCoef_psi(psi,T,info,PsiAna)
      USE EVR_system_m
      USE mod_basis
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),               intent(in)    :: psi
      real (kind=Rkind),              intent(in)    :: T ! time
      character (len=*),              intent(in)    :: info
      character (len=:), allocatable, intent(inout) :: PsiAna

      real (kind=Rkind)  :: C,maxC1,maxC2
      integer            :: i_bhe,i_b,i_e,i_h,i_R
      integer            :: i_bhe1,i_bhe2,i_b_maxC1,i_b_maxC2
      integer            :: i_R_maxC1,i_R_maxC2
      integer            :: i_e_maxC1,i_e_maxC2
      integer            :: i_h_maxC1,i_h_maxC2
      integer            :: ndim_index(psi%BasisnD%ndim)
      character (len=:), allocatable :: First_info,Second_info


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_MaxCoef_psi'
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'asso psi%BasisnD',associated(psi%BasisnD)
        flush(out_unit)
      END IF
!-----------------------------------------------------------
      IF (psi%nb_baie*psi%nb_bRot /= psi%nb_tot .AND.                   &
          .NOT. psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) RETURN   ! should be spectral WP

      CALL SET_string(First_info, 'First largest coefficient  ',TO_string(T),' ',trim(info),' ')
      CALL SET_string(Second_info,'Second largest coefficient ',TO_string(T),' ',trim(info),' ')

      maxC1 = ZERO
      maxC2 = ZERO
      i_b_maxC1 = 0
      i_b_maxC2 = 0
      i_h_maxC1 = 0
      i_h_maxC2 = 0
      i_e_maxC1 = 0
      i_e_maxC2 = 0
      i_R_maxC1 = 0
      i_R_maxC2 = 0
      i_bhe1    = 0
      i_bhe2    = 0

      i_bhe = 0
      DO i_R=1,psi%nb_bRot
      DO i_e=1,psi%nb_be
      DO i_h=1,psi%nb_bi
      DO i_b=1,psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
      !DO i_b=1,psi%nb_ba

        i_bhe = i_bhe + 1

        IF (psi%cplx) THEN
          C = abs(psi%CvecB(i_bhe))
        ELSE
          C = abs(psi%RvecB(i_bhe))
        END IF

        IF ( C > maxC1) THEN
          ! test if maxC1 > maxC2 (the old maxC1)
          IF ( maxC1 > maxC2) THEN
            maxC2     = maxC1
            i_bhe2    = i_bhe1
            i_b_maxC2 = i_b_maxC1
            i_h_maxC2 = i_h_maxC1
            i_e_maxC2 = i_e_maxC1
            i_R_maxC2 = i_R_maxC1
          END IF
          maxC1     = C
          i_bhe1    = i_bhe
          i_b_maxC1 = i_b
          i_h_maxC1 = i_h
          i_e_maxC1 = i_e
          i_R_maxC1 = i_R
        ELSE IF ( C > maxC2) THEN
          maxC2     = C
          i_bhe2    = i_bhe
          i_b_maxC2 = i_b
          i_h_maxC2 = i_h
          i_e_maxC2 = i_e
          i_R_maxC2 = i_R
        END IF
        !write(out_unit,*) i_bhe,C,'i_b_maxC1,i_b_maxC2',i_b_maxC1,i_b_maxC2
      END DO
      END DO
      END DO
      END DO

      IF (.NOT. allocated(PsiAna)) PsiAna = ''

      IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN

        IF (psi%cplx) THEN
          CALL ADD_TO_string(PsiAna,First_info, &
              TO_string([i_b_maxC1,i_h_maxC1,i_e_maxC1]),' ',TO_string(psi%CvecB(i_bhe1)),new_line('nl'))
          IF (i_b_maxC2 > 0) THEN
            CALL ADD_TO_string(PsiAna,Second_info, &
              TO_string([i_b_maxC2,i_h_maxC2,i_e_maxC2]),' ',TO_string(psi%CvecB(i_bhe2)),new_line('nl'))
          END IF
        ELSE
          CALL ADD_TO_string(PsiAna,First_info, &
              TO_string([i_b_maxC1,i_h_maxC1,i_e_maxC1]),' ',TO_string(psi%RvecB(i_bhe1)),new_line('nl'))
          
          IF (i_b_maxC2 > 0) THEN
            CALL ADD_TO_string(PsiAna,Second_info, &
              TO_string([i_b_maxC2,i_h_maxC2,i_e_maxC2]),' ',TO_string(psi%RvecB(i_bhe2)),new_line('nl'))
          END IF
        END IF
      ELSE
        IF (i_b_maxC1 < 1 .OR. i_b_maxC1 > psi%nb_ba) THEN
          CALL ADD_TO_string(PsiAna,'   WARNING, i_b_maxC1 is out-of-range !! ',TO_string(i_b_maxC1),new_line('nl'))
          CALL ADD_TO_string(PsiAna,'    it should > 0 and <= nb_ba',new_line('nl'))
        ELSE
          CALL Rec_ndim_index(psi%BasisnD,ndim_index,i_b_maxC1)
          IF (psi%cplx) THEN
            CALL ADD_TO_string(PsiAna,First_info,TO_string(ndim_index),' ', &
              TO_string([i_h_maxC1,i_e_maxC1]),' ',TO_string(psi%CvecB(i_bhe1)),new_line('nl'))

            IF (i_b_maxC2 > 0) THEN 
              CALL Rec_ndim_index(psi%BasisnD,ndim_index,i_b_maxC2)
              CALL ADD_TO_string(PsiAna,Second_info,TO_string(ndim_index),' ', &
                 TO_string([i_h_maxC2,i_e_maxC2]),' ',TO_string(psi%CvecB(i_bhe2)),new_line('nl'))
            END IF
          ELSE
            CALL ADD_TO_string(PsiAna,First_info,TO_string(ndim_index),' ', &
              TO_string([i_h_maxC1,i_e_maxC1]),' ',TO_string(psi%RvecB(i_bhe1)),new_line('nl'))

            IF (i_b_maxC2 > 0) THEN 
              CALL Rec_ndim_index(psi%BasisnD,ndim_index,i_b_maxC2)
              CALL ADD_TO_string(PsiAna,Second_info,TO_string(ndim_index),' ', &
                 TO_string([i_h_maxC2,i_e_maxC2]),' ',TO_string(psi%RvecB(i_bhe2)),new_line('nl'))
            END IF
          END IF
        END IF


      END IF
      CALL ADD_TO_string(PsiAna,'Abelian symmetry (symab): ',TO_string(psi%symab),new_line('nl'))
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      end subroutine calc_MaxCoef_psi

  SUBROUTINE calc_1Dweight_act1(psi,ana_psi,max_1D,T,info,print_w,PsiAna)
    USE EVR_system_m
    USE mod_nDindex
    USE mod_psi_set_alloc
    USE mod_type_ana_psi
    USE mod_param_RD
    IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
    TYPE (param_psi),               intent(in)    :: psi
    TYPE (param_ana_psi),           intent(inout) :: ana_psi
    character (len=:), allocatable, intent(inout) :: PsiAna

    character (len=*),              intent(in)    :: info
    logical,                        intent(in)    :: print_w
    real (kind=Rkind),              intent(in)    :: T ! time
    integer,                        intent(in)    :: max_1D


    real (kind=Rkind), allocatable :: weight1Dact(:,:)
    real (kind=Rkind)    :: a

    integer          :: i,ie,ii,ib,ibie,iq,ibiq,n,ii_baie,if_baie
    integer          :: max_dim
    integer          :: max_indGr(psi%BasisnD%nDindB%ndim)
    integer          :: ndim_AT_ib(psi%BasisnD%nDindB%ndim)
    integer          :: nDval(psi%BasisnD%nDindB%ndim)

    character (len=:), allocatable :: state_name

    real (kind=Rkind)                 :: r2
    real (kind=Rkind),    allocatable :: RDcontrac(:,:)   ! reduced density matrix with the contracted basis set (nbc,nbc)
    real (kind=Rkind),    allocatable :: RD(:,:)          ! reduced density matrix (nb,nb)
    complex (kind=Rkind), allocatable :: CRDcontrac(:,:)  ! reduced density matrix with the contracted basis set (nbc,nbc)
    complex (kind=Rkind), allocatable :: CRD(:,:)         ! reduced density matrix (nb,nb)

    TYPE (param_RD),   allocatable :: para_RD(:) ! it is allocated only for BasisnD. The size is nb_basis

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='calc_1Dweight_act1'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
!---------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nb_inact2n',psi%Basis2n%nb_basis
      write(out_unit,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
      write(out_unit,*) 'nb_baie,nb_tot',psi%nb_baie,psi%nb_tot
      flush(out_unit)
    END IF
!---------------------------------------------------------
    IF (psi%nb_baie /= psi%nb_tot) RETURN

    IF (allocated(Psi%BasisnD%nDindB%Tab_nDval)) THEN
      max_dim = maxval(Psi%BasisnD%nDindB%Tab_nDval)
    ELSE
      max_dim = maxval(psi%BasisnD%nDindB%nDsize(1:psi%BasisnD%nDindB%ndim))
    END IF
    CALL alloc_NParray(weight1Dact,[psi%BasisnD%nDindB%ndim,max_dim],"weight1Dact",name_sub)

    !---- RD analysis (for contrac_analysis=t) -------------------------
    ! initialization for RD analysis
    IF (psi%BasisnD%nb_basis > 1) THEN
      allocate(para_RD(psi%BasisnD%nb_basis))

      DO ib=1,size(para_RD)
        para_RD(ib)%RD_analysis =                                       &
                      psi%BasisnD%tab_Pbasis(ib)%Pbasis%contrac_analysis

        !para_RD(ib)%RD_analysis = .TRUE.

        para_RD(ib)%basis_index = ib
        IF (allocated(psi%BasisnD%tab_Pbasis(ib)%Pbasis%Rvec)) THEN
          CALL init_RD(para_RD(ib),psi%BasisnD%nDindB,            &
                       psi%BasisnD%tab_Pbasis(ib)%Pbasis%Rvec)
        ELSE
          CALL init_RD(para_RD(ib),psi%BasisnD%nDindB)
        END IF

      END DO
    END IF

    ibie = 0
    DO ie=1,psi%nb_be
    DO ii=1,psi%nb_bi
      weight1Dact(:,:) = ZERO
      ndim_AT_ib(:)    = 0
      IF (ie == 1 .AND. ii == 1) THEN
        CALL SET_string(state_name,'Grd Channel')
      ELSE
        CALL SET_string(state_name,'State_Se', TO_string(ie),'_Cha',TO_string(ii))
      END IF

      IF (psi%cplx) THEN

        DO ib=1,psi%nb_ba
          ibie = ibie + 1
          a = abs(psi%CvecB(ibie))**2

          CALL calc_nDindex(psi%BasisnD%nDindB,ib,nDval)

          DO iq=1,psi%BasisnD%nDindB%ndim
            ibiq = nDval(iq)
            ndim_AT_ib(iq) = max(ibiq,ndim_AT_ib(iq))
            !write(out_unit,*) 'calc_1Dweight',iq,ibiq,ib
            weight1Dact(iq,ibiq) = weight1Dact(iq,ibiq) + a
          END DO
        END DO
      ELSE

        DO ib=1,psi%nb_ba
          ibie = ibie + 1
          a = psi%RvecB(ibie)**2

          CALL calc_nDindex(psi%BasisnD%nDindB,ib,nDval)

          DO iq=1,psi%BasisnD%nDindB%ndim
            ibiq = nDval(iq)
            ndim_AT_ib(iq) = max(ibiq,ndim_AT_ib(iq))
            !write(out_unit,*) 'calc_1Dweight',iq,ibiq,ib
            weight1Dact(iq,ibiq) = weight1Dact(iq,ibiq) + a
          END DO
        END DO
      END IF

      !write(out_unit,*) 'ndim_AT_ib',ndim_AT_ib(:)
      IF (print_w .OR. debug) THEN
        DO iq=1,psi%BasisnD%nDindB%ndim
          IF (sum(weight1Dact(iq,1:ndim_AT_ib(iq)))-ONE > ONETENTH**7) THEN
            CALL ADD_TO_string(PsiAna,state_name,' Sum(RD)/=1',trim(info),' ', &
                 TO_string(iq),' ',TO_string(T,'f17.4'),' ',TO_string(T,'e10.3'),new_line('nl'))
          END IF
          CALL ADD_TO_string(PsiAna,state_name,' ',trim(info),' ', &
               TO_string(iq),' ',TO_string(T,'f17.4'),' ',TO_string(weight1Dact(iq,:),'e10.3'),new_line('nl'))

          !---- RD analysis --------------------------------------------
          ! RD analysis
          IF (allocated(para_RD)) THEN
          IF (para_RD(iq)%RD_analysis) THEN

            ii_baie = 1 + ( (ii-1)+ (ie-1)*psi%nb_bi ) * psi%nb_ba
            if_baie = ii_baie -1 + psi%nb_ba

            IF (psi%cplx) THEN
              CALL calc_CRD(para_RD(iq),psi%CvecB(ii_baie:if_baie),     &
                                           CRD=CRD,CRDcontrac=CRDcontrac)

              IF (allocated(CRD)) THEN
                n = min(ndim_AT_ib(iq),size(CRD,dim=1))
                weight1Dact(iq,1:n) = real([(CRD(i,i),i=1,n)],kind=Rkind)
                CALL ADD_TO_string(PsiAna,state_name,' new',trim(info),' ', &
                    TO_string(iq),' ',TO_string(T,'f17.4'),TO_string(weight1Dact(iq,:),'e10.3'),new_line('nl'))
                
                CALL dealloc_NParray(CRD,'CRD',name_sub)
              END IF

              IF (allocated(CRDcontrac)) THEN
                n = min(ndim_AT_ib(iq),size(CRDcontrac,dim=1))
                weight1Dact(iq,1:n) = real([(CRDcontrac(i,i),i=1,n)],kind=Rkind)

                CALL ADD_TO_string(PsiAna,state_name,' c',trim(info),' ', &
                  TO_string(iq),' ',TO_string(T,'f17.4'),TO_string(weight1Dact(iq,:),'e10.3'),new_line('nl'))

                CALL dealloc_NParray(CRDcontrac,'CRDcontrac',name_sub)
              END IF
            ELSE
              CALL calc_RD(para_RD(iq),psi%RvecB(ii_baie:if_baie),RD=RD,RDcontrac=RDcontrac)

              IF (allocated(RD)) THEN
                !n = min(ndim_AT_ib(iq),size(RD,dim=1))
                !write(out_unit,21) state_name // 'new',trim(info),iq,T,  &
                !                      [(RD(i,i),i=1,min(max_1D,n))]
                !flush(out_unit)
                CALL dealloc_NParray(RD,'RD',name_sub)
              END IF

              IF (allocated(RDcontrac)) THEN
                n = min(ndim_AT_ib(iq),size(RDcontrac,dim=1))
                weight1Dact(iq,1:n) = [(RDcontrac(i,i),i=1,n)]

                CALL ADD_TO_string(PsiAna,state_name,' c',trim(info),' ', &
                TO_string(iq),' ',TO_string(T,'f17.4'),TO_string(weight1Dact(iq,:),'e10.3'),new_line('nl'))

                CALL dealloc_NParray(RDcontrac,'RDcontrac',name_sub)
              END IF
            END IF

          END IF
          END IF
          !---- END RD analysis ----------------------------------------

          max_indGr(iq) = sum(maxloc(weight1Dact(iq,1:ndim_AT_ib(iq))))
        END DO
        CALL ADD_TO_string(PsiAna,'max_ind ',state_name ,' at ',trim(info), &
                           ' ',TO_string(max_indGr(1:size(max_indGr))),new_line('nl'))
      END IF

      !write(out_unit,*) 'max_RedDensity
      IF (.NOT. allocated(ana_psi%max_RedDensity)) THEN
        CALL alloc_NParray(ana_psi%max_RedDensity,[psi%BasisnD%nDindB%ndim], &
                          "ana_psi%max_RedDensity",name_sub)
      END IF
      DO iq=1,psi%BasisnD%nDindB%ndim
        ana_psi%max_RedDensity(iq) = ZERO
        DO ib=1,ndim_AT_ib(iq)
          r2 = real(ndim_AT_ib(iq)-ib,kind=Rkind)**2
          ana_psi%max_RedDensity(iq) = ana_psi%max_RedDensity(iq) + weight1Dact(iq,ib)*exp(-r2)
        END DO
      END DO

      CALL ADD_TO_string(PsiAna,'max_RedDensity of ',state_name,' ', &
          TO_string(ana_psi%max_RedDensity(:),Rformat='e9.2'),new_line('nl'))
    END DO
    END DO

    CALL dealloc_NParray(weight1Dact,"weight1Dact",name_sub)
    CALL dealloc_tab_RD(para_RD)

!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) PsiAna
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
!-----------------------------------------------------------

    END SUBROUTINE calc_1Dweight_act1


!=======================================================================================
!     norm^2 of psi (BasisRep or GridRep)
!=======================================================================================
      SUBROUTINE norm2_psi(psi,GridRep,BasisRep,ReNorm)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi), intent(inout) :: psi
      logical, intent(in), optional   :: ReNorm,GridRep,BasisRep

!------ working variables ---------------------------------
      integer       :: i_qa,i_qaie
      integer       :: i_be,i_bi,i_ba,i_baie
      integer       :: ii_baie,if_baie
      real (kind=Rkind) :: WrhonD,temp
      integer       :: i_max_w

      logical       :: psiN,norm2GridRep,norm2BasisRep

      real (kind=Rkind), allocatable :: tab_WeightChannels(:,:)


!----- for debuging --------------------------------------------------
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING norm2_psi'
        IF (present(ReNorm)) write(out_unit,*) 'Renormalization of psi',ReNorm
        IF (present(GridRep)) write(out_unit,*) 'norm2GridRep',GridRep
        IF (present(BasisRep)) write(out_unit,*) 'norm2BasisRep',BasisRep
        write(out_unit,*) 'psi'
        CALL ecri_psi(psi=psi)
      END IF
!-----------------------------------------------------------

      IF (present(ReNorm)) THEN
        psiN = ReNorm
      ELSE
        psiN = .FALSE.
      END IF

      IF (present(GridRep)) THEN
        IF (present(BasisRep)) THEN
          norm2GridRep  = GridRep
          norm2BasisRep = BasisRep
        ELSE
          norm2GridRep  = GridRep
          norm2BasisRep = .FALSE.
        END IF
      ELSE
        IF (present(BasisRep)) THEN
          norm2BasisRep = BasisRep
          norm2GridRep  = .FALSE.
        ELSE
          IF (psi%BasisRep .AND. psi%GridRep) THEN
            norm2BasisRep = .TRUE.
            norm2GridRep  = .FALSE.
          ELSE
            norm2BasisRep = psi%BasisRep
            norm2GridRep  = psi%GridRep
          END IF
        END IF
      END IF

      IF (debug) write(out_unit,*) 'nGridRep,nBasisRep,psiN',norm2GridRep,norm2BasisRep,psiN
      IF (norm2GridRep .AND. norm2BasisRep) THEN
        write(out_unit,*) ' ERROR in norm2_psi'
        write(out_unit,*) ' norm2GridRep=t and norm2BasisRep=t !'
        write(out_unit,*) ' BasisRep,GridRep',psi%BasisRep,psi%GridRep
        STOP
      END IF

      CALL Channel_weight(tab_WeightChannels,psi,norm2GridRep,norm2BasisRep)

      IF (debug)  write(out_unit,*) 'tab_WeightChannels : ',tab_WeightChannels

      psi%norm2 = sum(tab_WeightChannels)

      IF (debug) THEN
        write(out_unit,*) 'norm2 : ',psi%norm2,tab_WeightChannels
      END IF


      IF (psiN) THEN

        IF (psi%norm2 .EQ. ZERO ) THEN
          write(out_unit,*) ' ERROR in norm2_psi'
          write(out_unit,*) ' the norm2 is zero !',psi%norm2
          STOP
        END IF
        temp = ONE/psi%norm2
        tab_WeightChannels = tab_WeightChannels * temp
        temp = sqrt(temp)

        IF (norm2GridRep) THEN
!         - normalization of psiGridRep -------------------------
          IF (psi%cplx) THEN
            psi%CvecG(:) = psi%CvecG(:) *cmplx(temp,ZERO,kind=Rkind)
          ELSE
            psi%RvecG(:) = psi%RvecG(:) * temp
          END IF
        ELSE
!         - normalization of psiBasisRep -------------------------
          IF (psi%cplx) THEN
            psi%CvecB(:) = psi%CvecB(:) *cmplx(temp,ZERO,kind=Rkind)
          ELSE
            psi%RvecB(:) = psi%RvecB(:) * temp
          END IF
        END IF
        psi%norm2 = ONE
      END IF

      IF (allocated(tab_WeightChannels)) THEN
        CALL dealloc_NParray(tab_WeightChannels,"tab_WeightChannels","norm2_psi (alloc from Channel_weight)")
      END IF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'norm2 : ',psi%norm2,tab_WeightChannels
        write(out_unit,*) 'END norm2_psi'
      END IF
!----------------------------------------------------------

      END SUBROUTINE norm2_psi
!=======================================================================================

!=======================================================================================
      SUBROUTINE renorm_psi(psi,GridRep,BasisRep)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi), intent(inout)          :: psi
      logical,          intent(in),   optional :: GridRep,BasisRep

!------ working variables ---------------------------------
      logical                        :: norm2GridRep,norm2BasisRep
      real (kind=Rkind)              :: temp
      real (kind=Rkind), allocatable :: tab_WeightChannels(:,:)


!----- for debuging --------------------------------------------------
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING renorm_psi'
        IF (present(GridRep)) write(out_unit,*) 'norm2GridRep',GridRep
        IF (present(BasisRep)) write(out_unit,*) 'norm2BasisRep',BasisRep
        write(out_unit,*) 'psi'
        CALL ecri_psi(psi=psi)
      END IF
!-----------------------------------------------------------

      IF (present(GridRep)) THEN
        IF (present(BasisRep)) THEN
          norm2GridRep  = GridRep
          norm2BasisRep = BasisRep
        ELSE
          norm2GridRep  = GridRep
          norm2BasisRep = .FALSE.
        END IF
      ELSE
        IF (present(BasisRep)) THEN
          norm2BasisRep = BasisRep
          norm2GridRep  = .FALSE.
        ELSE
          IF (psi%BasisRep .AND. psi%GridRep) THEN
            norm2BasisRep = .TRUE.
            norm2GridRep  = .FALSE.
          ELSE
            norm2BasisRep = psi%BasisRep
            norm2GridRep  = psi%GridRep
          END IF
       END IF
     END IF

     IF (debug) write(out_unit,*) 'nGridRep,nBasisRep',norm2GridRep,norm2BasisRep
      IF (norm2GridRep .AND. norm2BasisRep) THEN
        write(out_unit,*) ' ERROR in renorm_psi'
        write(out_unit,*) ' norm2GridRep=t and norm2BasisRep=t !'
        write(out_unit,*) ' BasisRep,GridRep',psi%BasisRep,psi%GridRep
        STOP
      END IF

      CALL Channel_weight(tab_WeightChannels,psi,norm2GridRep,norm2BasisRep)

      IF (debug)  write(out_unit,*) 'tab_WeightChannels : ',tab_WeightChannels

      psi%norm2 = sum(tab_WeightChannels)

      IF (debug) THEN
        write(out_unit,*) 'norm2 : ',psi%norm2,tab_WeightChannels
      END IF

      IF (psi%norm2 .EQ. ZERO ) THEN
        write(out_unit,*) ' ERROR in renorm_psi'
        write(out_unit,*) ' the norm2 is zero !',psi%norm2
        STOP
      END IF
      temp = sqrt(ONE/psi%norm2)

      IF (norm2GridRep) THEN
        !- normalization of psiGridRep -------------------------
        IF (psi%cplx) THEN
          psi%CvecG(:) = psi%CvecG(:) *cmplx(temp,ZERO,kind=Rkind)
        ELSE
          psi%RvecG(:) = psi%RvecG(:) * temp
        END IF
      ELSE
        !- normalization of psiBasisRep -------------------------
        IF (psi%cplx) THEN
          psi%CvecB(:) = psi%CvecB(:) *cmplx(temp,ZERO,kind=Rkind)
        ELSE
          psi%RvecB(:) = psi%RvecB(:) * temp
        END IF
      END IF
      psi%norm2 = ONE

      IF (allocated(tab_WeightChannels)) THEN
        CALL dealloc_NParray(tab_WeightChannels,"tab_WeightChannels","renorm_psi (alloc from Channel_weight)")
      END IF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'norm2 : ',psi%norm2
        write(out_unit,*) 'END renorm_psi'
      END IF
!----------------------------------------------------------

      END SUBROUTINE renorm_psi
!=======================================================================================

!=======================================================================================
      SUBROUTINE renorm_psi_With_norm2(psi,GridRep,BasisRep)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi), intent(inout)          :: psi
      logical,          intent(in),   optional :: GridRep,BasisRep

!------ working variables ---------------------------------
      logical                        :: norm2GridRep,norm2BasisRep
      real (kind=Rkind)              :: temp

!----- for debuging --------------------------------------------------
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING renorm_psi_With_norm2'
        IF (present(GridRep))  write(out_unit,*) 'norm2GridRep',GridRep
        IF (present(BasisRep)) write(out_unit,*) 'norm2BasisRep',BasisRep
        write(out_unit,*) 'psi'
        CALL ecri_psi(psi=psi)
      END IF
!-----------------------------------------------------------

      IF (present(GridRep)) THEN
        IF (present(BasisRep)) THEN
          norm2GridRep  = GridRep
          norm2BasisRep = BasisRep
        ELSE
          norm2GridRep  = GridRep
          norm2BasisRep = .FALSE.
        END IF
      ELSE
        IF (present(BasisRep)) THEN
          norm2BasisRep = BasisRep
          norm2GridRep  = .FALSE.
        ELSE
          IF (psi%BasisRep .AND. psi%GridRep) THEN
            norm2BasisRep = .TRUE.
            norm2GridRep  = .FALSE.
          ELSE
            norm2BasisRep = psi%BasisRep
            norm2GridRep  = psi%GridRep
          END IF
       END IF
     END IF

     IF (debug) write(out_unit,*) 'nGridRep,nBasisRep',norm2GridRep,norm2BasisRep
      IF (norm2GridRep .AND. norm2BasisRep) THEN
        write(out_unit,*) ' ERROR in renorm_psi_With_norm2'
        write(out_unit,*) ' norm2GridRep=t and norm2BasisRep=t !'
        write(out_unit,*) ' BasisRep,GridRep',psi%BasisRep,psi%GridRep
        STOP
      END IF

      IF (debug) THEN
        write(out_unit,*) 'norm2 : ',psi%norm2
      END IF

      IF (psi%norm2 .EQ. ZERO ) THEN
        write(out_unit,*) ' ERROR in renorm_psi_With_norm2'
        write(out_unit,*) ' the norm2 is zero !',psi%norm2
        STOP
      END IF
      temp = sqrt(ONE/psi%norm2)

      IF (norm2GridRep) THEN
        !- normalization of psiGridRep -------------------------
        IF (psi%cplx) THEN
          IF(keep_MPI) psi%CvecG(:) = psi%CvecG(:) *cmplx(temp,ZERO,kind=Rkind)
        ELSE
          IF(keep_MPI) psi%RvecG(:) = psi%RvecG(:) * temp
        END IF
      ELSE
        !- normalization of psiBasisRep -------------------------
        IF (psi%cplx) THEN
          IF(keep_MPI) psi%CvecB(:) = psi%CvecB(:) *cmplx(temp,ZERO,kind=Rkind)
        ELSE
          IF(keep_MPI) psi%RvecB(:) = psi%RvecB(:) * temp
        END IF
      END IF
      psi%norm2 = ONE

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'norm2 : ',psi%norm2
        write(out_unit,*) 'END renorm_psi_With_norm2'
      END IF
!----------------------------------------------------------

      END SUBROUTINE renorm_psi_With_norm2

  SUBROUTINE Channel_weight(tab_WeightChannels,psi,                 &
                            GridRep,BasisRep,Dominant_Channel)
  USE EVR_system_m
  USE mod_basis
  USE mod_psi_set_alloc
  IMPLICIT NONE

!- variables for the WP ----------------------------------------
  real (kind=Rkind), intent(inout), allocatable :: tab_WeightChannels(:,:)
  TYPE (param_psi),  intent(inout)              :: psi
  logical,           intent(in)                 :: GridRep,BasisRep
  integer,           intent(inout), optional    :: Dominant_Channel(2)

!-- working variables ---------------------------------
  integer           :: i_qa,i_qaie
  integer           :: i_be,i_bi,i_ba,i_baie
  integer           :: ii_baie,if_baie
  real (kind=Rkind) :: WrhonD,temp,max_w
  integer           :: nb_be,nb_bi
  integer           :: iSG,iqSG,err_sub ! for SG4
  TYPE (OldParam)   :: OldPara
  real (kind=Rkind) :: WeightSG
  real (kind=Rkind) :: DSG



!- for debuging --------------------------------------------------
  character(len=*), parameter :: name_sub='Channel_weight'
  logical,parameter :: debug = .FALSE.
  !logical,parameter :: debug = .TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'GridRep',GridRep
    write(out_unit,*) 'BasisRep',BasisRep
    write(out_unit,*) 'alloc tab_WeightChannels',allocated(tab_WeightChannels)
    !write(out_unit,*) 'psi'
    !CALL ecri_psi(psi=psi)
    flush(out_unit)
  END IF
!-------------------------------------------------------

  nb_be = get_nb_be_FROM_psi(psi)
  nb_bi = get_nb_bi_FROM_psi(psi)

  IF (GridRep .AND. BasisRep) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' GridRep=t and BasisRep=t !'
    STOP
  END IF
  IF (.NOT. GridRep .AND. .NOT. BasisRep) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' GridRep=f and BasisRep=f !'
    STOP
  END IF
  IF (.NOT. allocated(tab_WeightChannels) .AND. nb_bi > 0 .AND. nb_be > 0) THEN
    CALL alloc_NParray(tab_WeightChannels,[nb_bi,nb_be],                        &
                      "tab_WeightChannels",name_sub)
  END IF
  tab_WeightChannels(:,:) = ZERO

  IF (debug) THEN
    write(out_unit,*) 'nb_bi,nb_be',nb_bi,nb_be
    write(out_unit,*) 'shape tab_WeightChannels',shape(tab_WeightChannels)
    write(out_unit,*) 'SG4',(psi%BasisnD%SparseGrid_type == 4)
  END IF

IF (psi%BasisnD%SparseGrid_type == 4) THEN
  !IF (GridRep)  CALL Channel_weight_SG4_grid(tab_WeightChannels,psi)
  IF (GridRep)  CALL Channel_weight_SG4_grid_old(tab_WeightChannels,psi)
  IF (BasisRep) CALL Channel_weight_SG4_basis(tab_WeightChannels,psi)
ELSE
  IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN

    CALL Channel_weight_contracHADA(tab_WeightChannels(:,1),psi)

  ELSE IF (psi%nb_baie == psi%nb_tot) THEN

    IF (BasisRep .AND. (allocated(psi%CvecB) .OR. allocated(psi%RvecB)) ) THEN

      DO i_be=1,psi%nb_be
      DO i_bi=1,psi%nb_bi
        ii_baie = 1 + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_ba
        if_baie = ii_baie -1 + psi%nb_ba

        IF (psi%cplx) THEN
          tab_WeightChannels(i_bi,i_be) =                               &
             real(dot_product(psi%CvecB(ii_baie:if_baie),               &
                              psi%CvecB(ii_baie:if_baie)) ,kind=Rkind)
        ELSE
          tab_WeightChannels(i_bi,i_be) =                               &
                              dot_product(psi%RvecB(ii_baie:if_baie),   &
                                          psi%RvecB(ii_baie:if_baie))
        END IF
      END DO
      END DO

    ELSE IF (GridRep .AND. (allocated(psi%CvecG) .OR. allocated(psi%RvecG)) ) THEN

      !- initialization ----------------------------------
      tab_WeightChannels(:,:) = ZERO

      CALL dealloc_OldParam(OldPara)

      DO i_qa=1,psi%nb_qa

        !- calculation of WrhonD ------------------------------
        WrhonD = Rec_WrhonD(psi%BasisnD,i_qa,OldPara)

        !write(6,*) 'i_qa,wrho*D',i_qa,WrhonD

        DO i_be=1,nb_be
        DO i_bi=1,nb_bi
          i_qaie = i_qa + ( (i_bi-1)+(i_be-1)*nb_bi ) * psi%nb_qa

          IF (psi%cplx) THEN
            temp = abs( psi%CvecG(i_qaie))
          ELSE
            temp = psi%RvecG(i_qaie)
          END IF
          tab_WeightChannels(i_bi,i_be) = tab_WeightChannels(i_bi,i_be) + &
                                          WrhonD *temp**2
        END DO
        END DO
      END DO
      CALL dealloc_OldParam(OldPara)

    ELSE

      write(out_unit,*) ' ERROR in ',name_sub,' from ',MPI_id
      IF (GridRep)  write(out_unit,*) ' impossible to calculate the weights with the Grid'
      IF (BasisRep) write(out_unit,*) ' impossible to calculate the weights with the Basis'
      write(out_unit,*) 'GridRep',GridRep
      write(out_unit,*) 'BasisRep',BasisRep
      write(out_unit,*)  'allo CvecG',allocated(psi%CvecG)
      write(out_unit,*)  'allo RvecG',allocated(psi%RvecG)
      write(out_unit,*)  'allo CvecB',allocated(psi%CvecB)
      write(out_unit,*)  'allo RvecB',allocated(psi%RvecB)
      STOP
    END IF

  ELSE
    !- To deal with a spectral representation ------

    tab_WeightChannels(:,:) = ZERO

    IF (psi%cplx) THEN
      tab_WeightChannels(1,1) = real(dot_product(psi%CvecB,psi%CvecB),kind=Rkind)
    ELSE
      tab_WeightChannels(1,1) = dot_product(psi%RvecB,psi%RvecB)
    END IF
  END IF
END IF

  IF (present(Dominant_Channel)) THEN
    Dominant_Channel(:) = 1
    max_w               = ZERO
    DO i_be=1,nb_be
    DO i_bi=1,nb_bi
      IF (tab_WeightChannels(i_bi,i_be) > max_w) THEN
        max_w = tab_WeightChannels(i_bi,i_be)
        Dominant_Channel(:) = [ i_be,i_bi ]
      END IF
    END DO
    END DO
  END IF

!------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'tab_WeightChannels : ',tab_WeightChannels
    write(out_unit,*) 'END ',name_sub
  END IF
!------------------------------------------------------

  END SUBROUTINE Channel_weight

  SUBROUTINE Channel_weight_SG4_grid(tab_WeightChannels,psi)
  USE EVR_system_m
  USE mod_psi_set_alloc
  USE mod_basis_BtoG_GtoB_SGType4
  USE mod_MPI_aux
  IMPLICIT NONE

!- variables for the WP ----------------------------------------
  real (kind=Rkind), intent(inout), allocatable :: tab_WeightChannels(:,:)
  TYPE (param_psi),  intent(inout)              :: psi

!-- working variables ---------------------------------
  real (kind=Rkind),    allocatable :: wrho(:)

  integer                           :: i,ib0,nb0,i_be,nb_be,i_bi,nb_bi
  integer                           :: iG,iq,nq_AT_iG
  complex (kind=Rkind), allocatable :: CVecG(:,:)
  real (kind=Rkind),    allocatable :: RVecG(:,:)
  Integer                           :: d1,d2

  !- for debuging --------------------------------------------------
  character(len=*), parameter :: name_sub='Channel_weight_SG4_grid'
  logical,parameter :: debug = .FALSE.
  !logical,parameter :: debug = .TRUE.
  !-------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'psi'
    CALL ecri_psi(psi=psi)
  END IF
!-------------------------------------------------------
  IF (.NOT. psi%GridRep) THEN
    write(out_unit,*) 'ERROR in ',name_sub
    write(out_unit,*) '  This subroutine needs psi with a grid representation'
    write(out_unit,*) '  but,  GridRep',psi%GridRep
    write(out_unit,*) '  and  BasisRep',psi%BasisRep
    write(out_unit,*) '  => use Channel_weight_SG4_basis'
    STOP 'ERROR in Channel_weight_SG4_grid: psi%GridRep=F'
  END IF

  nb_be = get_nb_be_FROM_psi(psi)
  nb_bi = get_nb_bi_FROM_psi(psi)

  d1=lbound(psi%BasisnD%para_SGType2%tab_nq_OF_SRep,dim=1)
  d2=ubound(psi%BasisnD%para_SGType2%tab_nq_OF_SRep,dim=1)

  IF(openmpi) THEN
    IF(MPI_scheme==1) THEN
      d1=iGs_MPI(1,MPI_id)
      d2=iGs_MPI(2,MPI_id)
    ELSEIF(MPI_scheme==3) THEN
      IF(MPI_nodes_p0) THEN
        d1=iGs_MPI(1,MPI_id)
        d2=iGs_MPI(2,MPI_sub_id(2))
      ELSE
        d1=0
        d2=-1
      ENDIF
    ENDIF
  ENDIF

  IF (psi%cplx) THEN
    ! psi almost in the right Smolyak representation.
    iq = 0
    DO iG=d1,d2
      nq_AT_iG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
      nb0      = psi%BasisnD%para_SGType2%nb0
!write(6,*) 'iG ',iG
!write(6,*) 'size Psi%CvecG(iq+1:iq+nb0*nq_AT_iG) ',size(Psi%CvecG(iq+1:iq+nb0*nq_AT_iG))
!write(6,*) 'nq_AT_iG,nb0 ',nq_AT_iG,psi%BasisnD%para_SGType2%nb0
!write(6,*) 'nq_AT_iG*nb0 ',nq_AT_iG*psi%BasisnD%para_SGType2%nb0
!flush(6)

      CVecG = reshape(Psi%CvecG(iq+1:iq+nb0*nq_AT_iG),shape=[nq_AT_iG,psi%BasisnD%para_SGType2%nb0])

      CALL Get_weight_FROM_OneDP(wrho,iG,                                   &
                    psi%BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,  &
                    psi%BasisnD%tab_basisPrimSG)

      !DO i=1,nq_AT_iG
      !  write(6,*) 'iq_a,wrho',iq+i,wrho(i),wrho(i)*Psi%BasisnD%WeightSG(iG)
      !END DO

      ib0 = 0
      DO i_be=1,nb_be
      DO i_bi=1,nb_bi
        ib0 = ib0 + 1
        tab_WeightChannels(i_bi,i_be) = tab_WeightChannels(i_bi,i_be) +        &
                          Psi%BasisnD%WeightSG(iG) * dot_product(CVecG(:,ib0), &
                                                            wrho * CVecG(:,ib0))
      END DO
      END DO
      iq = iq + nb0*nq_AT_iG
      IF (allocated(CVecG)) deallocate(CVecG)

    END DO

   ELSE
     ! psi almost in the right Smolyak representation.
     iq = 0
     DO iG=d1,d2
       nq_AT_iG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

       CALL Get_weight_FROM_OneDP(wrho,iG,                                   &
                    psi%BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,  &
                    psi%BasisnD%tab_basisPrimSG)

       RVecG = reshape(Psi%RvecG(iq+1:iq+nq_AT_iG),shape=[nq_AT_iG,psi%BasisnD%para_SGType2%nb0])

       ib0 = 0
       DO i_be=1,nb_be
       DO i_bi=1,nb_bi
         ib0 = ib0 + 1
         tab_WeightChannels(i_bi,i_be) = tab_WeightChannels(i_bi,i_be) +        &
                           Psi%BasisnD%WeightSG(iG) * dot_product(RVecG(:,ib0), &
                                                            wrho * RVecG(:,ib0))
       END DO
       END DO
       iq = iq + nq_AT_iG
       IF (allocated(RVecG)) deallocate(RVecG)

     END DO
  END IF

  IF(openmpi .AND. keep_MPI .AND. MPI_scheme/=2) THEN
    CALL MPI_Reduce_sum_matrix(tab_WeightChannels,1,nb_bi,1,nb_be,root_MPI,MS=MPI_scheme)
    CALL MPI_Bcast_matrix(tab_WeightChannels,1,nb_bi,1,nb_be,root_MPI,MS=MPI_scheme)
  ENDIF
  CALL dealloc_NParray(wrho,'wrho',name_sub)
!------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'tab_WeightChannels',tab_WeightChannels
    write(out_unit,*) 'END ',name_sub
  END IF
!------------------------------------------------------

END SUBROUTINE Channel_weight_SG4_grid
  SUBROUTINE Channel_weight_SG4_grid_old(tab_WeightChannels,psi)
  USE EVR_system_m
  USE mod_psi_set_alloc
  USE mod_basis_BtoG_GtoB_SGType4
  USE mod_MPI_aux
  IMPLICIT NONE

!- variables for the WP ----------------------------------------
  real (kind=Rkind), intent(inout), allocatable :: tab_WeightChannels(:,:)
  TYPE (param_psi),  intent(inout)              :: psi

!-- working variables ---------------------------------
  TYPE(Type_SmolyakRep)             :: WSRep ! smolyak rep for SparseGrid_type=4

  integer                           :: ib0,nb0,i_be,nb_be,i_bi,nb_bi
  integer                           :: iG,iq,nq_AT_iG
  complex (kind=Rkind), allocatable :: CVecG(:,:)
  real (kind=Rkind),    allocatable :: RVecG(:,:)
  Integer                           :: d1,d2

  !- for debuging --------------------------------------------------
  character(len=*), parameter :: name_sub='Channel_weight_SG4_grid_old'
  logical,parameter :: debug = .FALSE.
  !logical,parameter :: debug = .TRUE.
  !-------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'psi'
    CALL ecri_psi(psi=psi)
  END IF
!-------------------------------------------------------
  IF (.NOT. psi%GridRep) THEN
    write(out_unit,*) 'ERROR in ',name_sub
    write(out_unit,*) '  This subroutine needs psi with a grid representation'
    write(out_unit,*) '  but,  GridRep',psi%GridRep
    write(out_unit,*) '  and  BasisRep',psi%BasisRep
    write(out_unit,*) '  => use Channel_weight_SG4_basis'
    STOP 'ERROR in Channel_weight_SG4_grid: psi%GridRep=F'
  END IF

  CALL Set_weight_TO_SmolyakRep(WSRep,                                          &
                          psi%BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,  &
                          psi%BasisnD%tab_basisPrimSG)

   !write(out_unit,*) 'Weight'
   !CALL Write_SmolyakRep(WSRep)

  nb_be = get_nb_be_FROM_psi(psi)
  nb_bi = get_nb_bi_FROM_psi(psi)

  d1=lbound(WSRep%SmolyakRep,dim=1)
  d2=ubound(WSRep%SmolyakRep,dim=1)

  IF(openmpi) THEN
    IF(MPI_scheme==1) THEN
      d1=iGs_MPI(1,MPI_id)
      d2=iGs_MPI(2,MPI_id)
    ELSEIF(MPI_scheme==3) THEN
      IF(MPI_nodes_p0) THEN
        d1=iGs_MPI(1,MPI_id)
        d2=iGs_MPI(2,MPI_sub_id(2))
      ELSE
        d1=0
        d2=-1
      ENDIF
    ENDIF
  ENDIF

  IF (psi%cplx) THEN
    ! psi almost in the right Smolyak representation.
    iq = 0
    ! DO iG=lbound(WSRep%SmolyakRep,dim=1),ubound(WSRep%SmolyakRep,dim=1)
    DO iG=d1,d2
      nq_AT_iG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
      nb0      = psi%BasisnD%para_SGType2%nb0
!write(6,*) 'iG ',iG
!write(6,*) 'size Psi%CvecG(iq+1:iq+nb0*nq_AT_iG) ',size(Psi%CvecG(iq+1:iq+nb0*nq_AT_iG))
!write(6,*) 'nq_AT_iG,nb0 ',nq_AT_iG,psi%BasisnD%para_SGType2%nb0
!write(6,*) 'nq_AT_iG*nb0 ',nq_AT_iG*psi%BasisnD%para_SGType2%nb0
!flush(6)

      CVecG = reshape(Psi%CvecG(iq+1:iq+nb0*nq_AT_iG),shape=[nq_AT_iG,psi%BasisnD%para_SGType2%nb0])

      ib0 = 0
      DO i_be=1,nb_be
      DO i_bi=1,nb_bi
        ib0 = ib0 + 1
        tab_WeightChannels(i_bi,i_be) = tab_WeightChannels(i_bi,i_be) +        &
                          Psi%BasisnD%WeightSG(iG) * dot_product(CVecG(:,ib0), &
                                        WSRep%SmolyakRep(iG)%V * CVecG(:,ib0))
      END DO
      END DO
      iq = iq + nb0*nq_AT_iG
      IF (allocated(CVecG)) deallocate(CVecG)

    END DO

   ELSE
     ! psi almost in the right Smolyak representation.
     iq = 0
     DO iG=d1,d2
       nq_AT_iG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

       RVecG = reshape(Psi%RvecG(iq+1:iq+nq_AT_iG),shape=[nq_AT_iG,psi%BasisnD%para_SGType2%nb0])

       ib0 = 0
       DO i_be=1,nb_be
       DO i_bi=1,nb_bi
         ib0 = ib0 + 1
         tab_WeightChannels(i_bi,i_be) = tab_WeightChannels(i_bi,i_be) +        &
                           Psi%BasisnD%WeightSG(iG) * dot_product(RVecG(:,ib0), &
                                         WSRep%SmolyakRep(iG)%V * RVecG(:,ib0))
       END DO
       END DO
       iq = iq + nq_AT_iG
       IF (allocated(RVecG)) deallocate(RVecG)

     END DO
  END IF

  IF(openmpi .AND. keep_MPI .AND. MPI_scheme/=2) THEN
    CALL MPI_Reduce_sum_matrix(tab_WeightChannels,1,nb_bi,1,nb_be,root_MPI,MS=MPI_scheme)
    CALL MPI_Bcast_matrix(tab_WeightChannels,1,nb_bi,1,nb_be,root_MPI,MS=MPI_scheme)
  ENDIF

  CALL dealloc_SmolyakRep(WSRep)

!------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'tab_WeightChannels',tab_WeightChannels
    write(out_unit,*) 'END ',name_sub
  END IF
!------------------------------------------------------

END SUBROUTINE Channel_weight_SG4_grid_old
SUBROUTINE Channel_weight_SG4_basis(tab_WeightChannels,psi)
USE EVR_system_m
USE mod_basis
USE mod_psi_set_alloc
USE mod_basis_BtoG_GtoB_SGType4
USE mod_MPI_aux
IMPLICIT NONE

!- variables for the WP ----------------------------------------
real (kind=Rkind), intent(inout), allocatable :: tab_WeightChannels(:,:)
TYPE (param_psi),  intent(inout)              :: psi

!-- working variables ---------------------------------

TYPE (param_psi)          :: RCPsi(2)
TYPE(Type_SmolyakRep)     :: SRep ! smolyak rep for SparseGrid_type=4

integer           :: i_qa,i_qaie
integer           :: i_be,i_bi,i_ba,i_baie,ib0
integer           :: ii_baie,if_baie
real (kind=Rkind) :: WrhonD,temp,max_w
integer           :: nb_be,nb_bi
integer           :: iSG,iqSG,err_sub ! for SG4
TYPE (OldParam)   :: OldPara
real (kind=Rkind) :: WeightSG,Norm2



!- for debuging --------------------------------------------------
character(len=*), parameter :: name_sub='Channel_weight_SG4_basis'
logical,parameter :: debug = .FALSE.
!logical,parameter :: debug = .TRUE.
!-------------------------------------------------------
IF (debug) THEN
  write(out_unit,*) 'BEGINNING ',name_sub
  !write(out_unit,*) 'psi'
  !CALL ecri_psi(psi=psi)
END IF
!-------------------------------------------------------

IF (.NOT. psi%BasisRep) THEN
  write(out_unit,*) 'ERROR in ',name_sub
  write(out_unit,*) '  This subroutine needs psi with a basis representation'
  write(out_unit,*) '  but, BasisRep',psi%BasisRep
  write(out_unit,*) '  and   GridRep',psi%GridRep
  write(out_unit,*) '  => use Channel_weight_SG4_v0grid'
  STOP 'ERROR in Channel_weight_SG4_basis: psi%BasisRep=F'
END IF

nb_be = get_nb_be_FROM_psi(psi)
nb_bi = get_nb_bi_FROM_psi(psi)

IF (psi%cplx) THEN
  RCPsi = Psi

  !For the real part
  CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RCPsi(1)%RVecB,        &
                    psi%BasisnD%tab_basisPrimSG,psi%BasisnD%nDindB,  &
                    psi%BasisnD%para_SGType2)
  !write(out_unit,*) 'Real Part (Basis)'
  !CALL Write_SmolyakRep(SRep)

  ib0 = 0
  DO i_be=1,nb_be
  DO i_bi=1,nb_bi
    ib0 = ib0 + 1
    tab_WeightChannels(i_bi,i_be) =                                   &
      dot_product_SmolyakRep_Basis(SRep,SRep,psi%BasisnD%WeightSG,ib0)
  END DO
  END DO
  !write(out_unit,*) 'tab_WeightChannels',tab_WeightChannels

  !For the imaginary part
  CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RCPsi(2)%RVecB,        &
                    psi%BasisnD%tab_basisPrimSG,psi%BasisnD%nDindB,  &
                    psi%BasisnD%para_SGType2)
  !write(out_unit,*) 'Imaginary Part (Basis)'
  !CALL Write_SmolyakRep(SRep)

  ib0 = 0
  DO i_be=1,nb_be
  DO i_bi=1,nb_bi
    ib0 = ib0 + 1
    tab_WeightChannels(i_bi,i_be) = tab_WeightChannels(i_bi,i_be) +   &
      dot_product_SmolyakRep_Basis(SRep,SRep,psi%BasisnD%WeightSG,ib0)
  END DO
  END DO
  !CALL Write_SmolyakRep(SRep)

 ELSE
  CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,Psi%RVecB,              &
                    psi%BasisnD%tab_basisPrimSG,psi%BasisnD%nDindB,   &
                    psi%BasisnD%para_SGType2)
  ib0 = 0
  DO i_be=1,nb_be
  DO i_bi=1,nb_bi
    ib0 = ib0 + 1
    tab_WeightChannels(i_bi,i_be) =                                   &
      dot_product_SmolyakRep_Basis(SRep,SRep,psi%BasisnD%WeightSG,ib0)
  END DO
  END DO
  !CALL Write_SmolyakRep(SRep)
END IF

!write(out_unit,*) 'Norm2_SG4',Norm2

CALL dealloc_psi(RCPsi(1))
CALL dealloc_psi(RCPsi(2))

CALL dealloc_SmolyakRep(SRep)
!------------------------------------------------------
IF (debug) THEN
  write(out_unit,*) 'tab_WeightChannels',tab_WeightChannels
  write(out_unit,*) 'END ',name_sub
END IF
!------------------------------------------------------

END SUBROUTINE Channel_weight_SG4_basis
      SUBROUTINE Channel_weight_contracHADA(w_harm,psi)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

       TYPE (param_psi)   :: psi
       real (kind=Rkind) :: w_harm(psi%nb_bi)

       integer       :: ih,ihk,k


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING Channel_weight_contracHADA'
         write(out_unit,*) 'nb_bi,nb_ai',psi%nb_bi
       END IF
!-----------------------------------------------------------


!        weight on the harmonic states
         DO ih=1,psi%nb_bi
           ihk = sum(psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:ih-1))
           w_harm(ih) = ZERO
           DO k=1,psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(ih)
             ihk = ihk + 1
             IF (psi%cplx) THEN
               w_harm(ih) = w_harm(ih) + abs(psi%CvecB(ihk))**2
             ELSE
               w_harm(ih) = w_harm(ih) + psi%RvecB(ihk)**2
             END IF
           END DO
         END DO

!----------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*) 'END Channel_weight_contracHADA'
        END IF
!----------------------------------------------------------


      end subroutine Channel_weight_contracHADA

!=======================================================================================
!> MPI version of Channel_weight_contracADA
!=======================================================================================
!SUBROUTINE Channel_weight_contracADA_MPI(w_harm,psi)
!  USE EVR_system_m
!  USE mod_psi_set_alloc
!  IMPLICIT NONE
!
!  TYPE(param_psi)                                       :: psi
!  Real(kind=Rkind)                                      :: w_harm(psi%nb_bi)
!
!  Integer                                               :: ih
!  Integer                                               :: ihk
!  Integer                                               :: k
!
!  ! weight on the harmonic states
!  DO ih=1,psi%nb_bi
!    ihk=sum(psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:ih-1))
!    w_harm(ih)=ZERO
!    DO k=1,psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(ih)
!      ihk=ihk+1
!      write(*,*) 'ih,ihk check',ih,ihk
!      IF(psi%cplx) THEN
!        w_harm(ih)=w_harm(ih)+abs(psi%CvecB(ihk))**2
!      ELSE
!        w_harm(ih)=w_harm(ih)+psi%RvecB(ihk)**2
!      ENDIF
!    ENDDO
!  ENDDO
!
!END SUBROUTINE Channel_weight_contracADA_MPI

!==============================================================
!
!      calc_DM
!
!==============================================================
      SUBROUTINE calc_DM(psi,max_1D,T,info,print_DM)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: psi


      complex (kind=Rkind),allocatable :: DM(:,:),DM2(:,:)
      complex (kind=Rkind),allocatable :: DM_v1(:,:)
      complex (kind=Rkind),allocatable :: DM2_v1(:,:)
      real (kind=Rkind),allocatable    :: weight(:)
      complex (kind=Rkind) :: tr,tr2

      real (kind=Rkind)    :: T ! time
      logical          :: print_DM
      character (len=*)    :: info

      integer          :: nb_ba
      integer          :: ib,ibie
      integer          :: max_1D,max_print

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_DM'
      logical, parameter :: debug =.FALSE.
!     logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_tot',psi%nb_tot
        write(out_unit,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
        write(out_unit,*) 'nb_bie',psi%nb_bi*psi%nb_be
      END IF
!-----------------------------------------------------------

      IF (print_DM) THEN
        max_print = min(50,psi%nb_tot)
        CALL alloc_NParray(weight,[psi%nb_tot],"weight",name_sub)
        weight(:) = real(conjg(psi%CvecB)*psi%CvecB,kind=Rkind)
!       write(out_unit,*) 'T weight_tot',T,weight(1:max_print)
        write(out_unit,'(a,f12.1,1x,a,50(1x,f8.6))')                           &
                 'T weight_tot ',T,trim(info),weight(1:max_print)
        CALL dealloc_NParray(weight,"weight",name_sub)
      END IF
      IF (psi%nb_baie /= psi%nb_tot) THEN
        CALL alloc_NParray(DM ,[psi%nb_tot,psi%nb_tot],"DM" ,name_sub)
        CALL alloc_NParray(DM2,[psi%nb_tot,psi%nb_tot],"DM2",name_sub)

        DM(:,:) = CZERO
        DM(:,:) = matmul(                                             &
                conjg(reshape(psi%CvecB,[psi%nb_tot,1]) ),          &
                      reshape(psi%CvecB,[1,psi%nb_tot])  )

        DM2 = matmul(DM,DM)

        tr  = CZERO
        tr2 = CZERO
        DO  ibie=1,psi%nb_tot
          tr = tr + DM(ibie,ibie)
          tr2 = tr2 + DM2(ibie,ibie)
        END DO
        write(out_unit,*) 'T tr(DM),tr(DM2)',T,abs(tr),abs(tr2)
        CALL dealloc_NParray(DM ,"DM" ,name_sub)
        CALL dealloc_NParray(DM2,"DM2",name_sub)
      END IF

      IF (psi%nb_baie == psi%nb_tot .AND. psi%nb_bi*psi%nb_be > 1) THEN
        nb_ba = psi%nb_ba
        CALL alloc_NParray(DM_v1 ,[nb_ba,nb_ba],"DM_v1",name_sub)
        CALL alloc_NParray(DM2_v1,[nb_ba,nb_ba],"DM2_v1",name_sub)

        DM_v1(:,:) = matmul(                                            &
                conjg(reshape(psi%CvecB,[nb_ba,1]) ),                 &
                      reshape(psi%CvecB,[1,nb_ba])  )
!       DM_v1(:,:) = DM(1:psi%nb_ba,1:psi%nb_ba)
        DM2_v1 = matmul(DM_v1,DM_v1)
        tr  = CZERO
        tr2 = CZERO
        DO  ib=1,psi%nb_ba
          tr = tr + DM_v1(ib,ib)
          tr2 = tr2 + DM2_v1(ib,ib)
        END DO
        write(out_unit,*) 'T v1 tr(DM),tr(DM2)',T,abs(tr),abs(tr2)

        CALL dealloc_NParray(DM_v1 ,"DM_v1" ,name_sub)
        CALL dealloc_NParray(DM2_v1,"DM2_v1",name_sub)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------



      END SUBROUTINE calc_DM

!==============================================================
!
!      calc_nDTk(psi)
!
!==============================================================
      SUBROUTINE calc_nDTk(psi,T)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: psi


      real (kind=Rkind)    :: T ! time

      integer          :: k
      integer          :: iei,ib,ibie,nb_bie

      real (kind=Rkind),allocatable :: w_ei(:)

      complex (kind=Rkind),allocatable :: Ek_ei(:)
      complex (kind=Rkind),allocatable :: Ck_ei(:),Sk_ei(:)
      complex (kind=Rkind)    :: a

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_nDTk'
      integer :: err_mem,memory
      logical, parameter :: debug =.FALSE.
!     logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
      END IF
!-----------------------------------------------------------

      IF (psi%nb_baie /= psi%nb_tot) RETURN


      nb_bie = psi%nb_bi*psi%nb_be
      CALL alloc_NParray(w_ei,[nb_bie],"w_ei",name_sub)
      CALL alloc_NParray(Ek_ei,[nb_bie],"Ek_ei",name_sub)
      CALL alloc_NParray(Ck_ei,[nb_bie],"Ck_ei",name_sub)
      CALL alloc_NParray(Sk_ei,[nb_bie],"Sk_ei",name_sub)

      DO ib=1,psi%nb_ba
        k = int(ib/2)
        DO iei=1,nb_bie
          ibie = ib + (iei-1) * psi%nb_ba

          IF (psi%cplx) a = psi%CvecB(ibie)
          IF (.NOT. psi%cplx) a=cmplx(psi%RvecB(ibie),ZERO,kind=Rkind)
          IF (mod(ib-1,2) == 0) THEN
            Ck_ei(iei) = a
          ELSE
            Sk_ei(iei) = a
          END IF
        END DO
        IF (mod(ib-1,2) == 0) THEN

          Ek_ei = CHALF*(Ck_ei-EYE*Sk_ei)
          DO iei=1,nb_bie
            w_ei(iei) = abs(Ek_ei(iei))**2
          END DO
          write(998,*) 'T,+k',T,k,w_ei

          Ek_ei = CHALF*(Ck_ei+EYE*Sk_ei)
          DO iei=1,nb_bie
            w_ei(iei) = abs(Ek_ei(iei))**2
          END DO
          write(999,*) 'T,-k',T,k,w_ei

        END IF
      END DO
      write(998,*)
      write(999,*)
      flush(998)
      flush(999)

      CALL dealloc_NParray(w_ei,"w_ei",name_sub)
      CALL dealloc_NParray(Ek_ei,"Ek_ei",name_sub)
      CALL dealloc_NParray(Ck_ei,"Ck_ei",name_sub)
      CALL dealloc_NParray(Sk_ei,"Sk_ei",name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------



      END SUBROUTINE calc_nDTk

END MODULE mod_ana_psi
