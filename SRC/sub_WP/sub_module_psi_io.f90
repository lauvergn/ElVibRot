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
      MODULE mod_psi_io
      USE EVR_system_m
      USE mod_nDindex
      USE mod_basis
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: Read_psi_nDBasis,Write_Psi_nDBasis,sub_save_psi
      PUBLIC :: Write_header_saveFile_psi,Read_header_saveFile_psi
      PUBLIC :: lect_psiBasisRep,lect_psiBasisRepnotall,lect_psiBasisRepnotall_nD
      PUBLIC :: ecri_psiBasisRepnotall_nD

      CONTAINS

!=======================================================================================
!     Save vectors
!=======================================================================================
      SUBROUTINE sub_save_psi(psi,nb_save,file_WP)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (File_t) :: file_WP

      integer            :: nb_save
      TYPE (param_psi)   :: psi(nb_save)

!------ working parameters --------------------------------
      integer       :: j,nioWP,ioerr,Version_File
      !integer  :: Version_File,nb_psi,nb_tot

!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='sub_save_psi'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) ' nb_save',nb_save
        write(out_unit,*) ' nioWP',file_WP%unit
        write(out_unit,*) ' name',file_WP%name
        write(out_unit,*)
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      Version_File  = FilePsiVersion

      CALL file_open(file_WP,iunit=nioWP,lformatted=file_WP%formatted)

      CALL Write_header_saveFile_psi(psi,nb_save,file_WP)

      DO j=1,nb_save
          IF (psi(j)%BasisRep) THEN
            ! output
            CALL Write_Psi_nDBasis(psi(j),nioWP,j,ZERO,file_WP%formatted,Version_File)
          ELSE IF (psi(j)%GridRep) THEN
            IF (psi(j)%cplx) THEN
              write(nioWP,*) psi(j)%CvecG
            ELSE
              write(nioWP,*) psi(j)%RvecG
            END IF
          END IF
          flush(nioWP)
      END DO

      CALL file_close(file_WP)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'END ',name_sub
        flush(out_unit)
       END IF
!----------------------------------------------------------

      END SUBROUTINE sub_save_psi
!=======================================================================================

!=======================================================================================
!     reading BasisRep Wave packet
!---------------------------------------------------------------------------------------
      SUBROUTINE lect_psiBasisRep(psiBasisRep,WP0cplx,                  &
                                  nb_a,n_h,nb_elec,WP0n_h,WP0nb_elec)
      USE EVR_system_m
      IMPLICIT NONE

      integer          :: nb_a,n_h,nb_elec
      integer          :: WP0n_h,WP0nb_elec
      complex (kind=Rkind) :: psiBasisRep(nb_a,n_h,nb_elec)
      real (kind=Rkind)    :: a,b
      logical          :: WP0cplx

      integer          :: i_b,i_e,i_h,ini_h,ini_e,fin_h,fin_e
!----- for debuging --------------------------------------------------
!     logical, parameter :: debug =.FALSE.
      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING lect_psiBasisRep'
        write(out_unit,*) 'nb_a,n_h,nb_elec',nb_a,n_h,nb_elec
        write(out_unit,*) 'WP0n_h,WP0nb_elec',WP0n_h,WP0nb_elec
      END IF
!-----------------------------------------------------------

      ini_h = WP0n_h
      ini_e = WP0nb_elec
      fin_h = WP0n_h
      fin_e = WP0nb_elec
      IF (WP0n_h == 0) THEN
        ini_h = 1
        fin_h = n_h
      END IF
      IF (WP0nb_elec == 0) THEN
        ini_e = 1
        fin_e = nb_elec
      END IF

      IF ( WP0cplx ) THEN
         DO i_e=ini_e,fin_e
         DO i_h=ini_h,fin_h
           DO i_b=1,nb_a
             read(in_unit,*) a,b
             psiBasisRep(i_b,i_h,i_e) = cmplx(a,b,kind=Rkind)
           END DO
           read(in_unit,*)
         END DO
         END DO
      ELSE
         DO i_e=ini_e,fin_e
         DO i_h=ini_h,fin_h
           DO i_b=1,nb_a
             read(in_unit,*) a
             psiBasisRep(i_b,i_h,i_e) = cmplx(a,ZERO,kind=Rkind)
           END DO
           read(in_unit,*)
         END DO
         END DO
       END IF

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'psiBasisRep'
         DO i_e=ini_e,fin_e
         DO i_h=ini_h,fin_h
           DO i_b=1,nb_a
             write(out_unit,*) i_e,i_h,i_b,psiBasisRep(i_b,i_h,i_e)
           END DO
           write(out_unit,*)
         END DO
         END DO
         write(out_unit,*) 'END lect_psiBasisRep'
       END IF
!-----------------------------------------------------------

      END SUBROUTINE lect_psiBasisRep
!=======================================================================================

!==============================================================
!     reading BasisRep Wave packet
!==============================================================
      SUBROUTINE lect_psiBasisRepnotall(WP0,WP0cplx)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: WP0

      real (kind=Rkind)  :: a,b
      logical            :: WP0cplx

      integer            :: i_bhe,i_b,i_e,i_h
      integer            :: err
!----- for debuging --------------------------------------------------
      !logical, parameter :: debug =.FALSE.
      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING lect_psiBasisRepnotall'
      END IF
!-----------------------------------------------------------

      IF (WP0%cplx) THEN
        WP0%CvecB(:) = CZERO
      ELSE
        WP0%RvecB(:) = ZERO
      END IF

      DO
        a = ZERO
        b = ZERO
        IF (WP0cplx) read(in_unit,*,iostat=err) i_b,i_h,i_e,a,b
        IF (.NOT. WP0cplx) read(in_unit,*,iostat=err) i_b,i_h,i_e,a
!       write(out_unit,*) i_b,i_h,i_e,a,b,'err',err
        IF (err /= 0) EXIT
        i_bhe = i_b + ( (i_h-1)+ (i_e-1)*WP0%nb_bi ) * WP0%nb_ba
        IF (WP0%cplx) THEN
          WP0%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
        ELSE
          WP0%RvecB(i_bhe) = a
        END IF
      END DO

!-----------------------------------------------------------
       IF (debug) THEN

         write(out_unit,*) 'WP0BasisRep'
           CALL ecri_psi(T=ZERO,psi=WP0)
         write(out_unit,*) 'END lect_psiBasisRepnotall'
       END IF
!-----------------------------------------------------------

      END SUBROUTINE lect_psiBasisRepnotall
      SUBROUTINE lect_psiBasisRepnotall_nD(WP0,nioWP,WP0cplx,lformated)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: WP0

      real (kind=Rkind)  :: a,b
      logical            :: WP0cplx
      integer            :: nioWP
      logical            :: lformated
      logical            :: lformated_loc
      character (len=10) :: name_end

      integer            :: i_bhe,i_b,i_e,i_h,i_bguess
      integer            :: Rerr
      integer            :: i,ndim
      integer            :: ind_contractcHAC(3)
      integer            :: ind_ndim(WP0%BasisnD%ndim+2)
      logical, save      :: done_basis_is_smaller=.FALSE.
      logical            :: basis_is_OK
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='lect_psiBasisRepnotall_nD'
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'cplx',WP0cplx
        write(out_unit,*) 'WP0%cplx',WP0%cplx
        write(out_unit,*) 'nb_basis_act1',WP0%nb_basis_act1
        write(out_unit,*) 'nioWP',nioWP
        write(out_unit,*) 'lformated',lformated
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      WP0 = ZERO

      ndim = WP0%BasisnD%ndim

      lformated_loc = lformated
      IF (nioWP == 5) lformated_loc = .TRUE.

!-----------------------------------------------------------
!-----------------------------------------------------------
      i_bguess = 1
      IF (lformated_loc) THEN
        IF (WP0%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
          DO
            read(nioWP,'(a)',iostat=Rerr) name_end
            IF (index(name_end,"end") > 0) EXIT
            backspace(nioWP)


            a = ZERO
            b = ZERO
            IF (WP0cplx) THEN
             read(nioWP,*,iostat=Rerr) ind_contractcHAC(:),a,b
            ELSE
             read(nioWP,*,iostat=Rerr) ind_contractcHAC(:),a
            END IF
            IF (debug) write(out_unit,*) ind_contractcHAC(:),a,b,'Rerr',Rerr
            IF (Rerr /= 0) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) ' problem while reading the WP'
              write(out_unit,*) ind_contractcHAC(:),a,b,'Rerr',Rerr
              STOP
            END IF

            ! test the electronic and adiabatic channels
            i_h = ind_contractcHAC(2)
            i_e = ind_contractcHAC(3)
            basis_is_OK = i_h <= WP0%nb_bi .AND. i_e <= WP0%nb_be

            IF (.NOT. basis_is_OK) THEN
              IF (.NOT. done_basis_is_smaller)                          &
                     write(out_unit,*) ' WARNNING the basis is smaller'
              done_basis_is_smaller = .TRUE.
            END IF

            IF (.NOT. basis_is_OK) CYCLE

            i_b = ind_contractcHAC(1)


            basis_is_OK = i_b <= WP0%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
            IF (debug) write(out_unit,*) 'i_b,i_h,i_e',i_b,i_h,i_e

            IF (.NOT. basis_is_OK) THEN
              IF (.NOT. done_basis_is_smaller)                          &
                     write(out_unit,*) ' WARNNING the basis is smaller'
              done_basis_is_smaller = .TRUE.
            END IF
            IF (.NOT. basis_is_OK) CYCLE

            i_bhe = sum(WP0%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1)) + i_b
            IF (i_bhe > WP0%nb_tot) STOP 'ERROR i_bhe>nb_tot'

            IF (debug) write(out_unit,*) 'i_bhe,i_b,i_h,i_e,a,b',i_bhe,ind_contractcHAC(:),a,b
            IF (WP0%cplx) THEN
              WP0%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
            ELSE
              WP0%RvecB(i_bhe) = a
            END IF
         END DO

        ELSE IF (WP0%nb_tot == WP0%nb_baie) THEN
          DO
            a = ZERO
            b = ZERO
            IF (WP0cplx) THEN
             read(nioWP,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a,b
            ELSE
             read(nioWP,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a
            END IF
            IF (debug) write(out_unit,*) ind_ndim(:),a,b,'Rerr',Rerr
            IF (Rerr /= 0) EXIT

            ! test the electronic and adiabatic channels
            i_h = ind_ndim(ndim+1)
            i_e = ind_ndim(ndim+2)
            basis_is_OK = i_h <= WP0%nb_bi .AND. i_e <= WP0%nb_be

            IF (.NOT. basis_is_OK) THEN
              IF (.NOT. done_basis_is_smaller)                          &
                     write(out_unit,*) ' WARNNING the basis is smaller'
              done_basis_is_smaller = .TRUE.
            END IF

            IF (.NOT. basis_is_OK) CYCLE

            CALL calc_InD_FROM_ndim_index(WP0%BasisnD,ind_ndim(1:ndim),i_b,i_bguess)

            basis_is_OK = i_b <= WP0%nb_ba

            IF (.NOT. basis_is_OK) THEN
              IF (.NOT. done_basis_is_smaller)                          &
                         write(out_unit,*) ' WARNNING the basis is smaller'
              done_basis_is_smaller = .TRUE.
            END IF
            IF (.NOT. basis_is_OK) CYCLE
            i_bguess = i_b + 1


            i_bhe = i_b + ( (i_h-1)+ (i_e-1)*WP0%nb_bi ) * WP0%nb_ba
            IF (debug) write(out_unit,*) 'i_bhe,indnD,a,b',i_bhe,ind_ndim(:),a,b
            IF(keep_MPI) THEN
              IF (WP0%cplx) THEN
                WP0%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
              ELSE
                WP0%RvecB(i_bhe) = a
              END IF
            ENDIF
          END DO
        ELSE
          RETURN
        END IF
!-----------------------------------------------------------
!-----------------------------------------------------------
      ELSE
      DO
        a = ZERO
        b = ZERO
        IF (WP0cplx) THEN
          read(nioWP,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a,b
        ELSE
          read(nioWP,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a
        END IF
        IF (debug) write(out_unit,*) ind_ndim(:),a,b,'Rerr',Rerr
        IF (Rerr /= 0) EXIT

        ! test the electronic and adiabatic channels
        i_h = ind_ndim(ndim+1)
        i_e = ind_ndim(ndim+2)
        basis_is_OK = i_h <= WP0%nb_bi .AND. i_e <= WP0%nb_be

        IF (.NOT. basis_is_OK) THEN
          IF (.NOT. done_basis_is_smaller)                              &
                         write(out_unit,*) ' WARNNING the basis is smaller'
          done_basis_is_smaller = .TRUE.
        END IF

        IF (.NOT. basis_is_OK) CYCLE

        CALL calc_InD_FROM_ndim_index(WP0%BasisnD,ind_ndim(1:ndim),i_b,i_bguess)
        basis_is_OK = i_b <= WP0%nb_ba

        IF (.NOT. basis_is_OK) THEN
          IF (.NOT. done_basis_is_smaller)                              &
                         write(out_unit,*) ' WARNNING the basis is smaller'
          done_basis_is_smaller = .TRUE.
        END IF
        IF (.NOT. basis_is_OK) CYCLE
        i_bguess = i_b + 1

        i_bhe = i_b + ( (i_h-1)+ (i_e-1)*WP0%nb_bi ) * WP0%nb_ba
        IF (debug) write(out_unit,*) 'i_bhe,indnD,a,b',i_bhe,ind_ndim(:),a,b
        IF (WP0%cplx) THEN
          WP0%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
        ELSE
          WP0%RvecB(i_bhe) = a
        END IF
      END DO
      END IF
!-----------------------------------------------------------
!-----------------------------------------------------------

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'WP0BasisRep'
         CALL ecri_psi(T=ZERO,psi=WP0)
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
      END SUBROUTINE lect_psiBasisRepnotall_nD

!=======================================================================================
      SUBROUTINE ecri_psiBasisRepnotall_nD(WP0,nio,epsi,lformated,iwp)
      USE EVR_system_m
      USE mod_psi_set_alloc
      !USE mod_psi_Op
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: WP0
      real (kind=Rkind)  :: epsi
      integer            :: nio
      logical            :: lformated
      integer            :: iwp
      logical            :: lformated_loc

      real (kind=Rkind)  :: a,b

      integer            :: i_bhe,i_b,i_e,i_h
      integer            :: ind_ndim(WP0%BasisnD%ndim)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='ecri_psiBasisRepnotall_nD'
      logical, parameter :: debug =.FALSE.
!     logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'cplx',WP0%cplx
        write(out_unit,*) 'nb_basis_act1',WP0%nb_basis_act1
        write(out_unit,*) 'nb_ba',WP0%nb_ba
        write(out_unit,*) 'epsi',epsi
      END IF
!-----------------------------------------------------------

      lformated_loc = lformated
      IF (nio == 6) lformated_loc = .TRUE.

      IF (WP0%nb_tot == WP0%nb_baie) THEN
!-----------------------------------------------------------
        IF (lformated_loc) THEN
          i_bhe = 0
          DO i_e=1,WP0%nb_be
          DO i_h=1,WP0%nb_bi
          DO i_b=1,WP0%nb_ba
          i_bhe = i_bhe + 1

          CALL Rec_ndim_index(WP0%BasisnD,ind_ndim(:),i_b)

          IF (WP0%cplx) THEN
            a = real(WP0%CvecB(i_bhe),kind=Rkind)
            b = aimag(WP0%CvecB(i_bhe))
            IF (sqrt(a*a+b*b) >= epsi) write(nio,*) ind_ndim(:),i_h,i_e,a,b
          ELSE
            a = WP0%RvecB(i_bhe)
!           write(out_unit,*) indnD(:),a
            IF (abs(a) >= epsi) write(nio,*) ind_ndim(:),i_h,i_e,a
          END IF
          END DO
          END DO
          END DO
          write(nio,*) 'end wp ',iwp
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        ELSE
          i_bhe = 0
          DO i_e=1,WP0%nb_be
          DO i_h=1,WP0%nb_bi
          DO i_b=1,WP0%nb_ba
          i_bhe = i_bhe + 1

          CALL Rec_ndim_index(WP0%BasisnD,ind_ndim(:),i_b)

          IF (WP0%cplx) THEN
            a = real(WP0%CvecB(i_b),kind=Rkind)
            b = aimag(WP0%CvecB(i_b))
            IF (sqrt(a*a+b*b) >= epsi) write(nio) ind_ndim(:),i_h,i_e,a,b
          ELSE
            a = WP0%RvecB(i_b)
!           write(out_unit,*) indnD(:),a
            IF (abs(a) >= epsi) write(nio) ind_ndim(:),i_h,i_e,a
          END IF
          END DO
          END DO
          END DO
          write(nio) 'end wp ',iwp
        END IF
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 ELSE IF (WP0%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
   IF (lformated_loc) THEN
        DO i_e=1,WP0%nb_be
        DO i_h=1,WP0%nb_bi
           i_bhe = sum(WP0%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,WP0%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             IF (WP0%cplx) THEN
               a = real(WP0%CvecB(i_bhe),kind=Rkind)
               b = aimag(WP0%CvecB(i_bhe))
               IF (sqrt(a*a+b*b) >= epsi) write(nio,*) i_b,i_h,i_e,a,b
             ELSE
               a = WP0%RvecB(i_bhe)
!              write(out_unit,*) indnD(:),a
               IF (abs(a) >= epsi) write(nio,*) i_b,i_h,i_e,a
             END IF
           END DO
        END DO
        END DO
        write(nio,*) 'end wp ',iwp
    ELSE
        DO i_e=1,WP0%nb_be
        DO i_h=1,WP0%nb_bi
           i_bhe = sum(WP0%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,WP0%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             IF (WP0%cplx) THEN
               a = real(WP0%CvecB(i_bhe),kind=Rkind)
               b = aimag(WP0%CvecB(i_bhe))
               IF (sqrt(a*a+b*b) >= epsi) write(nio,*) i_b,i_h,i_e,a,b
             ELSE
               a = WP0%RvecB(i_bhe)
!              write(out_unit,*) indnD(:),a
               IF (abs(a) >= epsi) write(nio) i_b,i_h,i_e,a
             END IF
           END DO
        END DO
        END DO
        write(nio) 'end wp ',iwp
    END IF
 END IF
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      END SUBROUTINE ecri_psiBasisRepnotall_nD
!=======================================================================================

SUBROUTINE Read_psi_nDBasis(Psi,nioPsi,lformated,version,  &
                            list_nDindBasis1_TO_nDindBasis2,nb_tot)
USE EVR_system_m
USE mod_psi_set_alloc
IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
TYPE (param_psi), intent(inout)   :: Psi
integer,          intent(in)      :: nioPsi,nb_tot,version
logical,          intent(in)      :: lformated
integer,          intent(in)      :: list_nDindBasis1_TO_nDindBasis2(nb_tot)

real (kind=Rkind)  :: a,b,T
logical            :: lformated_loc
character (len=10) :: name_end

real (kind=Rkind),    allocatable  :: Rvec(:)
complex (kind=Rkind), allocatable  :: Cvec(:)


integer            :: i_bhe,i_bhe_read,i_b,i_e,i_h
integer            :: Rerr
integer            :: i,ndim
integer            :: ind_contractcHAC(3)
integer            :: ind_ndim(Psi%BasisnD%ndim+2)
logical, save      :: done_basis_is_smaller=.FALSE.
logical            :: basis_is_OK
!----- for debuging --------------------------------------------------
character (len=*), parameter :: name_sub='Read_psi_nDBasis'
logical, parameter :: debug =.FALSE.
!logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
IF (debug) THEN
  write(out_unit,*) 'BEGINNING ',name_sub
  write(out_unit,*) 'Psi%cplx',Psi%cplx
  write(out_unit,*) 'nb_basis_act1',Psi%nb_basis_act1
  write(out_unit,*) 'nioPsi',nioPsi
  write(out_unit,*) 'lformated',lformated
  write(out_unit,*) 'version: ',version
  flush(out_unit)
END IF
!-----------------------------------------------------------

Psi = ZERO

ndim = Psi%BasisnD%ndim


lformated_loc = lformated
IF (nioPsi == 5 .OR. nioPsi == in_unit) lformated_loc = .TRUE.
!-----------------------------------------------------------

SELECT CASE (version)
CASE(0)
  IF (lformated_loc) THEN
    IF (Psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
      DO
        read(nioPsi,'(a)',iostat=Rerr) name_end
        IF (index(name_end,"end") > 0) EXIT
        backspace(nioPsi)


        a = ZERO
        b = ZERO
        IF (Psi%cplx) THEN
         read(nioPsi,*,iostat=Rerr) ind_contractcHAC(:),a,b
        ELSE
         read(nioPsi,*,iostat=Rerr) ind_contractcHAC(:),a
        END IF
        IF (debug) write(out_unit,*) ind_contractcHAC(:),a,b,'Rerr',Rerr
        IF (Rerr /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' problem while reading the WP'
          write(out_unit,*) ind_contractcHAC(:),a,b,'Rerr',Rerr
          STOP
        END IF

        ! test the electronic and adiabatic channels
        i_h = ind_contractcHAC(2)
        i_e = ind_contractcHAC(3)
        basis_is_OK = i_h <= Psi%nb_bi .AND. i_e <= Psi%nb_be

        IF (.NOT. basis_is_OK) THEN
          IF (.NOT. done_basis_is_smaller)                          &
                 write(out_unit,*) ' WARNNING the basis is smaller'
          done_basis_is_smaller = .TRUE.
        END IF

        IF (.NOT. basis_is_OK) CYCLE

        i_b = ind_contractcHAC(1)


        basis_is_OK = i_b <= Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
        IF (debug) write(out_unit,*) 'i_b,i_h,i_e',i_b,i_h,i_e

        IF (.NOT. basis_is_OK) THEN
          IF (.NOT. done_basis_is_smaller)                          &
                 write(out_unit,*) ' WARNNING the basis is smaller'
          done_basis_is_smaller = .TRUE.
        END IF
        IF (.NOT. basis_is_OK) CYCLE

        i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1)) + i_b
        IF (i_bhe > Psi%nb_tot) STOP 'ERROR i_bhe>nb_tot'

        IF (debug) write(out_unit,*) 'i_bhe,i_b,i_h,i_e,a,b',i_bhe,ind_contractcHAC(:),a,b
        IF (Psi%cplx) THEN
          Psi%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
        ELSE
          Psi%RvecB(i_bhe) = a
        END IF
     END DO

    ELSE IF (Psi%nb_tot == Psi%nb_baie) THEN
      DO i_bhe_read=1,nb_tot
        a = ZERO
        b = ZERO
        IF (Psi%cplx) THEN
         read(nioPsi,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a,b
        ELSE
         read(nioPsi,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a
        END IF
        IF (debug) write(out_unit,*) ind_ndim(:),a,b,'Rerr',Rerr
        IF (Rerr /= 0) EXIT

        i_bhe = list_nDindBasis1_TO_nDindBasis2(i_bhe_read)

        IF (debug) write(out_unit,*) 'i_bhe,indnD,a,b',i_bhe,ind_ndim(:),a,b

        IF (i_bhe > 0 .AND. i_bhe <= Psi%nb_baie) THEN
          IF (Psi%cplx) THEN
            Psi%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
          ELSE
            Psi%RvecB(i_bhe) = a
          END IF
        END IF

      END DO
      read(nioPsi,*,iostat=Rerr)   ! for the end wp

    ELSE
      RETURN
    END IF
  ELSE
    DO i_bhe_read=1,nb_tot
      a = ZERO
      b = ZERO
      IF (Psi%cplx) THEN
       read(nioPsi,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a,b
      ELSE
       read(nioPsi,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a
      END IF
      IF (debug) write(out_unit,*) ind_ndim(:),a,b,'Rerr',Rerr
      IF (Rerr /= 0) EXIT

      i_bhe = list_nDindBasis1_TO_nDindBasis2(i_bhe_read)

      IF (debug) write(out_unit,*) 'i_bhe,indnD,a,b',i_bhe,ind_ndim(:),a,b

      IF (i_bhe > 0) THEN
        IF (Psi%cplx) THEN
          Psi%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
        ELSE
          Psi%RvecB(i_bhe) = a
        END IF
      END IF

    END DO
    !read(nioPsi,iostat=Rerr)   ! for the end wp

  END IF
CASE(1)
  IF (Psi%cplx) THEN
    CALL alloc_NParray(Cvec,shape(list_nDindBasis1_TO_nDindBasis2),'Cvec',name_sub)
    IF (lformated_loc) THEN
      read(nioPsi,*,iostat=Rerr) Cvec(:)
    ELSE
      read(nioPsi,iostat=Rerr) Cvec(:)
    END IF
    DO i_bhe_read=1,nb_tot
      Psi%CvecB((list_nDindBasis1_TO_nDindBasis2(i_bhe_read))) = Cvec(i_bhe_read)
    END DO
    CALL dealloc_NParray(Cvec,'Cvec',name_sub)

  ELSE

    CALL alloc_NParray(Rvec,shape(list_nDindBasis1_TO_nDindBasis2),'Rvec',name_sub)
    IF (lformated_loc) THEN
      read(nioPsi,*,iostat=Rerr) Rvec(:)
    ELSE
      read(nioPsi,iostat=Rerr) Rvec(:)
    END IF
    DO i_bhe_read=1,nb_tot
      Psi%RvecB((list_nDindBasis1_TO_nDindBasis2(i_bhe_read))) = Rvec(i_bhe_read)
    END DO

    CALL dealloc_NParray(Rvec,'Rvec',name_sub)
  END IF

  IF (Rerr /= 0) THEN
    write(out_unit,*) 'ERROR in ',name_sub
    write(out_unit,*) 'Problem while reading the WP with version=',version
    STOP
  END IF
CASE(2)
  IF (lformated_loc) THEN
    IF (Psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
      DO
        read(nioPsi,'(a)',iostat=Rerr) name_end
        IF (index(name_end,"end") > 0) EXIT
        backspace(nioPsi)


        a = ZERO
        b = ZERO
        IF (Psi%cplx) THEN
         read(nioPsi,*,iostat=Rerr) T,ind_contractcHAC(:),a,b
        ELSE
         read(nioPsi,*,iostat=Rerr) T,ind_contractcHAC(:),a
        END IF
        IF (debug) write(out_unit,*) T,ind_contractcHAC(:),a,b,'Rerr',Rerr
        IF (Rerr /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' problem while reading the WP'
          write(out_unit,*) ind_contractcHAC(:),a,b,'Rerr',Rerr
          STOP
        END IF

        ! test the electronic and adiabatic channels
        i_h = ind_contractcHAC(2)
        i_e = ind_contractcHAC(3)
        basis_is_OK = i_h <= Psi%nb_bi .AND. i_e <= Psi%nb_be

        IF (.NOT. basis_is_OK) THEN
          IF (.NOT. done_basis_is_smaller)                          &
                 write(out_unit,*) ' WARNNING the basis is smaller'
          done_basis_is_smaller = .TRUE.
        END IF

        IF (.NOT. basis_is_OK) CYCLE

        i_b = ind_contractcHAC(1)


        basis_is_OK = i_b <= Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
        IF (debug) write(out_unit,*) 'i_b,i_h,i_e',i_b,i_h,i_e

        IF (.NOT. basis_is_OK) THEN
          IF (.NOT. done_basis_is_smaller)                          &
                 write(out_unit,*) ' WARNNING the basis is smaller'
          done_basis_is_smaller = .TRUE.
        END IF
        IF (.NOT. basis_is_OK) CYCLE

        i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1)) + i_b
        IF (i_bhe > Psi%nb_tot) STOP 'ERROR i_bhe>nb_tot'

        IF (debug) write(out_unit,*) 'i_bhe,i_b,i_h,i_e,a,b',i_bhe,ind_contractcHAC(:),a,b
        IF (Psi%cplx) THEN
          Psi%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
        ELSE
          Psi%RvecB(i_bhe) = a
        END IF
     END DO

    ELSE IF (Psi%nb_tot == Psi%nb_baie) THEN
      DO i_bhe_read=1,nb_tot

        a = ZERO
        b = ZERO
        IF (Psi%cplx) THEN
         read(nioPsi,*,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a,b
        ELSE
         read(nioPsi,*,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a
        END IF
        IF (debug) write(out_unit,*) T,ind_ndim(:),a,b,'Rerr',Rerr
        IF (Rerr /= 0) EXIT

        i_bhe = list_nDindBasis1_TO_nDindBasis2(i_bhe_read)

        IF (debug) write(out_unit,*) 'T,i_bhe,indnD,a,b',T,i_bhe,ind_ndim(:),a,b

        IF (i_bhe > 0 .AND. i_bhe <= Psi%nb_baie) THEN
          IF (Psi%cplx) THEN
            Psi%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
          ELSE
            Psi%RvecB(i_bhe) = a
          END IF
        END IF

      END DO
      read(nioPsi,*,iostat=Rerr)   ! for the end wp

    ELSE
      RETURN
    END IF
  ELSE
    DO i_bhe_read=1,nb_tot
      a = ZERO
      b = ZERO
      IF (Psi%cplx) THEN
       read(nioPsi,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a,b
      ELSE
       read(nioPsi,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a
      END IF
      IF (debug) write(out_unit,*) T,ind_ndim(:),a,b,'Rerr',Rerr
      IF (Rerr /= 0) EXIT

      i_bhe = list_nDindBasis1_TO_nDindBasis2(i_bhe_read)

      IF (debug) write(out_unit,*) 'T,i_bhe,indnD,a,b',T,i_bhe,ind_ndim(:),a,b

      IF (i_bhe > 0) THEN
        IF (Psi%cplx) THEN
          Psi%CvecB(i_bhe) = cmplx(a,b,kind=Rkind)
        ELSE
          Psi%RvecB(i_bhe) = a
        END IF
      END IF

    END DO
    !read(nioPsi,iostat=Rerr)   ! for the end wp

  END IF

CASE Default
  STOP 'no default in Read_psi_nDBasis'
END SELECT
!-----------------------------------------------------------
!-----------------------------------------------------------

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'PsiBasisRep'
         CALL ecri_psi(T=ZERO,psi=Psi)
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
      END SUBROUTINE Read_psi_nDBasis

      SUBROUTINE Read_list_nDindBasis1_TO_nDindBasis2(Psi,              &
               list_nDindBasis1_TO_nDindBasis2,nb_tot,nioPsi,lformated, &
               version)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_psi), intent(in)   :: Psi

      integer, intent(in)    :: nb_tot,nioPsi,version
      integer, intent(inout) :: list_nDindBasis1_TO_nDindBasis2(nb_tot)
      logical, intent(in)    :: lformated


      logical            :: lformated_loc
      integer            :: i_bhe,i_bhe_read,i_b,i_e,i_h,i_bguess
      integer            :: Rerr
      integer            :: i,ndim
      integer            :: ind_ndim(Psi%BasisnD%ndim+2)
      logical            :: done_basis_is_smaller
      logical            :: basis_is_OK
      real (kind=Rkind)  :: a,b,T
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_list_nDindBasis1_TO_nDindBasis2'
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_basis_act1',Psi%nb_basis_act1
        write(out_unit,*) 'nioPsi',nioPsi
        write(out_unit,*) 'lformated',lformated
        write(out_unit,*) 'version: ',version
        flush(out_unit)
      END IF
!-----------------------------------------------------------

done_basis_is_smaller = .FALSE.

lformated_loc = lformated
IF (nioPsi == 5 .OR. nioPsi == in_unit) lformated_loc = .TRUE.

SELECT CASE (version)
CASE(0)
  i_bguess   = 1
  ndim = Psi%BasisnD%ndim
  list_nDindBasis1_TO_nDindBasis2(:) = 0

  i_bhe_read = 0
  IF (Psi%nb_tot == Psi%nb_baie) THEN
    DO
      i_bhe_read = i_bhe_read + 1

      IF (lformated_loc) THEN  ! formated file
        IF (psi%cplx) THEN
           ! for ifort, it is important to read the coefficients (a,b), &
           !  otherwise if a and b are in a second line, it can be interpreted
           !  as an integer of the next line (ndim=3).
          read(nioPsi,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a,b
        ELSE
          read(nioPsi,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a
        END IF
      ELSE
        IF (psi%cplx) THEN
          read(nioPsi,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a,b
        ELSE
          read(nioPsi,iostat=Rerr) (ind_ndim(i),i=1,ndim+2),a
        END IF
      END IF

      IF (Rerr /= 0) EXIT

      ! test the electronic and adiabatic channels
      i_h = ind_ndim(ndim+1)
      i_e = ind_ndim(ndim+2)

      basis_is_OK = (i_h <= Psi%nb_bi .AND. i_e <= Psi%nb_be)
      IF (.NOT. basis_is_OK) done_basis_is_smaller = .TRUE.
      IF (.NOT. basis_is_OK) CYCLE

      CALL calc_InD_FROM_ndim_index(Psi%BasisnD,ind_ndim(1:ndim),i_b,i_bguess)

      basis_is_OK = (i_b <= Psi%nb_ba)
      IF (.NOT. basis_is_OK) done_basis_is_smaller = .TRUE.
      IF (.NOT. basis_is_OK) CYCLE

      i_bguess = i_b + 1


      i_bhe = i_b + ( (i_h-1)+ (i_e-1)*Psi%nb_bi ) * Psi%nb_ba
      IF (debug) write(out_unit,*) 'i_bhe,indnD',i_bhe,ind_ndim(:)
      list_nDindBasis1_TO_nDindBasis2(i_bhe_read) = i_bhe
    END DO
  ELSE
    RETURN
  END IF
CASE(1)
  i_bguess   = 1
  ndim = Psi%BasisnD%ndim
  list_nDindBasis1_TO_nDindBasis2(:) = 0

  IF (Psi%nb_tot == Psi%nb_baie) THEN
    DO i_bhe_read =1,nb_tot

      IF (lformated_loc) THEN  ! formated file
        read(nioPsi,*,iostat=Rerr) (ind_ndim(i),i=1,ndim+2)
      ELSE
        read(nioPsi,iostat=Rerr) (ind_ndim(i),i=1,ndim+2)
      END IF

      IF (Rerr /= 0) STOP 'ERROR in Read_list_nDindBasis1_TO_nDindBasis2 version=1'

      ! test the electronic and adiabatic channels
      i_h = ind_ndim(ndim+1)
      i_e = ind_ndim(ndim+2)

      basis_is_OK = (i_h <= Psi%nb_bi .AND. i_e <= Psi%nb_be)
      IF (.NOT. basis_is_OK) done_basis_is_smaller = .TRUE.
      IF (.NOT. basis_is_OK) CYCLE

      CALL calc_InD_FROM_ndim_index(Psi%BasisnD,ind_ndim(1:ndim),i_b,i_bguess)

      basis_is_OK = (i_b <= Psi%nb_ba)
      IF (.NOT. basis_is_OK) done_basis_is_smaller = .TRUE.
      IF (.NOT. basis_is_OK) CYCLE
      i_bguess = i_b + 1

      i_bhe = i_b + ( (i_h-1)+ (i_e-1)*Psi%nb_bi ) * Psi%nb_ba
      IF (debug) write(out_unit,*) 'i_bhe,indnD',i_bhe,ind_ndim(:)
      list_nDindBasis1_TO_nDindBasis2(i_bhe_read) = i_bhe
    END DO
  ELSE
    RETURN
  END IF
CASE(2)
  i_bguess   = 1
  ndim = Psi%BasisnD%ndim
  list_nDindBasis1_TO_nDindBasis2(:) = 0

  i_bhe_read = 0
  IF (Psi%nb_tot == Psi%nb_baie) THEN
    DO
      i_bhe_read = i_bhe_read + 1

      IF (lformated_loc) THEN  ! formated file
        IF (psi%cplx) THEN
           ! for ifort, it is important to read the coefficients (a,b), &
           !  otherwise if a and b are in a second line, it can be interpreted
           !  as an integer of the next line (ndim=3).
          read(nioPsi,*,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a,b
        ELSE
          read(nioPsi,*,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a
        END IF
      ELSE
        IF (psi%cplx) THEN
          read(nioPsi,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a,b
        ELSE
          read(nioPsi,iostat=Rerr) T,(ind_ndim(i),i=1,ndim+2),a
        END IF
      END IF

      IF (Rerr /= 0) EXIT

      ! test the electronic and adiabatic channels
      i_h = ind_ndim(ndim+1)
      i_e = ind_ndim(ndim+2)

      basis_is_OK = (i_h <= Psi%nb_bi .AND. i_e <= Psi%nb_be)
      IF (.NOT. basis_is_OK) done_basis_is_smaller = .TRUE.
      IF (.NOT. basis_is_OK) CYCLE

      CALL calc_InD_FROM_ndim_index(Psi%BasisnD,ind_ndim(1:ndim),i_b,i_bguess)

      basis_is_OK = (i_b <= Psi%nb_ba)
      IF (.NOT. basis_is_OK) done_basis_is_smaller = .TRUE.
      IF (.NOT. basis_is_OK) CYCLE

      i_bguess = i_b + 1


      i_bhe = i_b + ( (i_h-1)+ (i_e-1)*Psi%nb_bi ) * Psi%nb_ba
      IF (debug) write(out_unit,*) 'i_bhe,indnD',i_bhe,ind_ndim(:)
      list_nDindBasis1_TO_nDindBasis2(i_bhe_read) = i_bhe
    END DO
  ELSE
    RETURN
  END IF
CASE Default
  STOP 'no default in Read_list_nDindBasis1_TO_nDindBasis2'
END SELECT

IF (done_basis_is_smaller) write(out_unit,*) ' WARNNING the basis is smaller'

!-----------------------------------------------------------
IF (debug) THEN
  write(out_unit,*) ' list_nDindBasis1_TO_nDindBasis2'
  write(out_unit,'(10(I0,1X))') list_nDindBasis1_TO_nDindBasis2(:)
  write(out_unit,*) 'END ',name_sub
END IF

END SUBROUTINE Read_list_nDindBasis1_TO_nDindBasis2

SUBROUTINE Read_header_saveFile_psi(psi,nb_read,list_nDindBasis1_TO_nDindBasis2, &
                                    file_WP,Version_File)
USE EVR_system_m
USE mod_psi_set_alloc
IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
integer,              intent(inout) :: nb_read,Version_File
TYPE (param_psi),     intent(inout) :: psi(:)
TYPE (File_t),        intent(inout) :: file_WP
integer, allocatable, intent(inout) :: list_nDindBasis1_TO_nDindBasis2(:)



!------ working parameters --------------------------------


integer                  :: nioWP
character (len=Line_len) :: line

integer  :: n1,n2

integer   :: nb_psi,nb_tot
character :: ch
namelist / headerFile / Version_File,nb_psi,nb_tot

!----- for debuging --------------------------------------------------
integer :: err_mem,memory,ioerr
character (len=*), parameter ::name_sub='Read_header_saveFile_psi'
logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
IF (debug) THEN
  write(out_unit,*) 'BEGINNING ',name_sub
  write(out_unit,*)
END IF

! first we read the header of the file. There are 3 options:
! - one integer:  nb_readWP_file,                 hence Version_File=0, formatted only
! - two integers: nb_readWP_file and nb_tot_file, hence Version_File=0, formatted only
! - the &headerFile namelist,                     hence Version_File=1, formatted and unformatted (not a namelist)

nioWP = file_WP%unit

IF (file_WP%formatted) THEN
  Version_File = 0 ! default for the option 1 and 2
  !read(nioWP,*,iostat=ioerr) nb_readWP,nb_tot
  read(nioWP,'(a)') line
  read(line,*,iostat=ioerr) nb_read,nb_tot ! here we test the option 2
  IF (ioerr /= 0) THEN ! it means that nb_tot is not present or we have the namelist (option 1 or 3)
    read(line,*,iostat=ioerr) nb_read
    nb_tot = 0
  END IF

  IF (ioerr /= 0) THEN ! it means that no integer is present => It must be the namelist (option 3)
    ! the file has to be reopened
    CALL file_close(file_WP)
    CALL file_open(file_WP,nioWP,lformatted=file_WP%formatted,old=.TRUE.)

    Version_File = FilePsiVersion ! default from the module "EVR_system_m"
    nb_psi       = 0
    nb_tot       = 0
    read(nioWP,headerFile,iostat=ioerr)
    nb_read      = nb_psi

    ! True error while reading the file
    IF (ioerr /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Error while reading the formatted file header (', &
                          trim(file_WP%name),').'
      STOP
    END IF

!    CALL alloc_NParray(list_nDindBasis1_TO_nDindBasis2,[nb_tot],&
!                      "list_nDindBasis1_TO_nDindBasis2",name_sub)
!    CALL Read_list_nDindBasis1_TO_nDindBasis2(psi(1),             &
!                    list_nDindBasis1_TO_nDindBasis2,nb_tot,nioWP, &
!                    file_WP%formatted,Version_File)

  END IF



ELSE  ! not formated file
  ! here only the option 3, however not with a namelist
  read(nioWP,iostat=ioerr) ch
  IF (ioerr == 0 .OR. ch == 'X') THEN
    read(nioWP,iostat=ioerr) nb_read
    read(nioWP,iostat=ioerr) nb_tot
    read(nioWP,iostat=ioerr) Version_File
    read(nioWP,iostat=ioerr) ch
  END IF

  ! True error while reading the file
  IF (ioerr /= 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' Error while reading the unformatted file header (', &
                        trim(file_WP%name),').'
    STOP
  END IF

END IF

IF (debug) THEN
  write(out_unit,*) ' nb_readWP_file: ',nb_read
  write(out_unit,*) ' nb_tot_file   : ',nb_tot
  write(out_unit,*) ' Version_File  : ',Version_File
  flush(out_unit)
END IF

! When nb_tot > 0, the list_nDindBasis1_TO_nDindBasis2(:) has to be read.
IF (nb_tot > 0) THEN
  CALL alloc_NParray(list_nDindBasis1_TO_nDindBasis2,[nb_tot],&
                    "list_nDindBasis1_TO_nDindBasis2",name_sub)
  CALL Read_list_nDindBasis1_TO_nDindBasis2(psi(1),             &
                  list_nDindBasis1_TO_nDindBasis2,nb_tot,nioWP, &
                  file_WP%formatted,Version_File)

  IF (debug) THEN
    write(out_unit,*) 'list_nDindBasis1_TO_nDindBasis2'
    write(out_unit,'(10(I0,1X))') list_nDindBasis1_TO_nDindBasis2(:)
  END IF

  IF (Version_File == 0 .OR. Version_File == 2) THEN ! Version_File=0, option=2
    ! the file has to be reoponed (only for Version_File=0, option=2)
    ! because one psi has been read to get the indexes
    CALL file_close(file_WP)
    CALL file_open(file_WP,nioWP,lformatted=file_WP%formatted,old=.TRUE.)
    IF (file_WP%formatted) THEN
      read(nioWP,*,iostat=ioerr) n1,n2
    ELSE
      read(nioWP,iostat=ioerr) n1,n2
    END IF
  END IF

END IF

! here, the file must be open and ready to read the 1st psi.

IF (debug) THEN
  write(out_unit,*)
  write(out_unit,*) 'END ',name_sub
END IF

END SUBROUTINE Read_header_saveFile_psi

      SUBROUTINE Write_header_saveFile_psi(psi,nb_save,file_WP)
      USE EVR_system_m
      USE mod_psi_set_alloc
      !USE mod_psi_Op
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (File_t), intent(in)  :: file_WP
      integer          , intent(in)  :: nb_save
      TYPE (param_psi) , intent(in)  :: psi(:)


!------ working parameters --------------------------------
      integer  ::   ioerr

      integer  :: Version_File,nb_psi,nb_tot
      namelist / headerFile / Version_File,nb_psi,nb_tot

!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='Write_header_saveFile_psi'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) ' nb_save,ndim',nb_save,shape(psi)
        write(out_unit,*)
      END IF
!-----------------------------------------------------------

      Version_File  = FilePsiVersion
      nb_psi        = nb_save
      nb_tot        = psi(1)%nb_tot

      IF (file_WP%formatted) THEN
        IF (Version_File == 0) THEN
          write(file_WP%unit,*) nb_psi,nb_tot
        ELSE ! version 1
          write(file_WP%unit,nml=headerFile)
          CALL Write_list_nDindBasis(psi(1),file_WP%unit,file_WP%formatted,Version_File)
        END IF
      ELSE ! version=1
        write(file_WP%unit,iostat=ioerr) 'X'
        write(file_WP%unit,iostat=ioerr) nb_psi
        write(file_WP%unit,iostat=ioerr) nb_tot
        write(file_WP%unit,iostat=ioerr) 1
        write(file_WP%unit,iostat=ioerr) 'X'
        CALL Write_list_nDindBasis(psi(1),file_WP%unit,file_WP%formatted,Version_File)
      END IF

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'END ',name_sub
       END IF
!----------------------------------------------------------

      END SUBROUTINE Write_header_saveFile_psi
      SUBROUTINE Write_list_nDindBasis(Psi,nioPsi,lformated,version)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_psi), intent(in)   :: Psi
      integer,          intent(in)   :: nioPsi,version
      logical,          intent(in)   :: lformated


      logical            :: lformated_loc
      integer            :: i_bhe,i_b,i_e,i_h
      integer            :: ind_ndim(Psi%BasisnD%ndim)
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Write_list_nDindBasis'
      logical, parameter :: debug =.FALSE.
!      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_basis_act1',Psi%nb_basis_act1
        write(out_unit,*) 'nioPsi',nioPsi
        write(out_unit,*) 'lformated',lformated
        write(out_unit,*) 'version: ',version
        flush(out_unit)
      END IF
!-----------------------------------------------------------


lformated_loc = lformated
IF (nioPsi == 6 .OR. nioPsi == out_unit) lformated_loc = .TRUE.

SELECT CASE (version)
CASE(0)
  CONTINUE  ! nothing to be done
  ! because the indices are already with the vectors
CASE(1)
  IF (Psi%nb_tot == Psi%nb_baie) THEN

    IF (lformated_loc) THEN
      i_bhe = 0
      DO i_e=1,Psi%nb_be
      DO i_h=1,Psi%nb_bi
      DO i_b=1,Psi%nb_ba
        i_bhe = i_bhe + 1

        CALL Rec_ndim_index(Psi%BasisnD,ind_ndim(:),i_b)
        write(nioPsi,*) ind_ndim(:),i_h,i_e

      END DO
      END DO
      END DO

    ELSE ! not formatted
      i_bhe = 0
      DO i_e=1,Psi%nb_be
      DO i_h=1,Psi%nb_bi
      DO i_b=1,Psi%nb_ba
        i_bhe = i_bhe + 1

        CALL Rec_ndim_index(Psi%BasisnD,ind_ndim(:),i_b)
        write(nioPsi) ind_ndim(:),i_h,i_e

      END DO
      END DO
      END DO
   END IF

  ELSE IF (Psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
    IF (lformated_loc) THEN
        DO i_e=1,Psi%nb_be
        DO i_h=1,Psi%nb_bi
           i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             write(nioPsi,*) i_b,i_h,i_e
           END DO
        END DO
        END DO
    ELSE
        DO i_e=1,Psi%nb_be
        DO i_h=1,Psi%nb_bi
           i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             write(nioPsi) i_b,i_h,i_e
           END DO
        END DO
        END DO
    END IF
  END IF

CASE Default
  STOP 'no default in Write_list_nDindBasis'
END SELECT
!-----------------------------------------------------------
!-----------------------------------------------------------

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
      END SUBROUTINE Write_list_nDindBasis

!=======================================================================================
  SUBROUTINE Write_Psi_nDBasis(Psi,nioPsi,iPsi,epsi,lformated,version,T)
    USE EVR_system_m
    USE mod_psi_set_alloc
    IMPLICIT NONE


    !----- variables for the WP propagation ----------------------------
    TYPE (param_psi),  intent(in)           :: Psi
    integer,           intent(in)           :: nioPsi,version,iPsi
    logical,           intent(in)           :: lformated
    real (kind=Rkind), intent(in)           :: epsi
    real (kind=Rkind), intent(in), optional :: T

    real (kind=Rkind)  :: a,b
    logical            :: lformated_loc

    integer            :: i_bhe,i_bhe_read,i_b,i_e,i_h
    integer            :: Rerr
    integer            :: i,ndim
    integer            :: ind_ndim(Psi%BasisnD%ndim)
    real (kind=Rkind)  :: T_loc

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_psi_nDBasis'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'Psi%cplx',Psi%cplx
      write(out_unit,*) 'nioPsi',nioPsi
      write(out_unit,*) 'lformated',lformated
      write(out_unit,*) 'version: ',version
      IF (present(T))  write(out_unit,*) 'T:       ',T
      flush(out_unit)
    END IF
    !-----------------------------------------------------------

    ndim = Psi%BasisnD%ndim

    lformated_loc = lformated
    IF (nioPsi == 6 .OR. nioPsi == out_unit) lformated_loc = .TRUE.

!-----------------------------------------------------------
!-----------------------------------------------------------

SELECT CASE (version)
CASE(0)

  IF (Psi%nb_tot == Psi%nb_baie) THEN

    IF (lformated_loc) THEN
      i_bhe = 0
      DO i_e=1,Psi%nb_be
      DO i_h=1,Psi%nb_bi
      DO i_b=1,Psi%nb_ba
      i_bhe = i_bhe + 1

        CALL Rec_ndim_index(Psi%BasisnD,ind_ndim(:),i_b)

        IF (Psi%cplx) THEN
          a = real(Psi%CvecB(i_bhe),kind=Rkind)
          b = aimag(Psi%CvecB(i_bhe))
          IF (sqrt(a*a+b*b) >= epsi) write(nioPsi,*) ind_ndim(:),i_h,i_e,a,b
        ELSE
          a = Psi%RvecB(i_bhe)
          IF (abs(a) >= epsi) write(nioPsi,*) ind_ndim(:),i_h,i_e,a
        END IF
      END DO
      END DO
      END DO
      write(nioPsi,*) 'end wp ',ipsi

    ELSE ! not formatted
      i_bhe = 0
      DO i_e=1,Psi%nb_be
      DO i_h=1,Psi%nb_bi
      DO i_b=1,Psi%nb_ba
      i_bhe = i_bhe + 1

        CALL Rec_ndim_index(Psi%BasisnD,ind_ndim(:),i_b)


        IF (Psi%cplx) THEN
          a = real(Psi%CvecB(i_b),kind=Rkind)
          b = aimag(Psi%CvecB(i_b))
          IF (sqrt(a*a+b*b) >= epsi) write(nioPsi) ind_ndim(:),i_h,i_e,a,b
        ELSE
          a = Psi%RvecB(i_b)
!         write(out_unit,*) indnD(:),a
          IF (abs(a) >= epsi) write(nioPsi) ind_ndim(:),i_h,i_e,a
        END IF
      END DO
      END DO
      END DO
      !write(nioPsi) 'end wp ',ipsi
   END IF

  ELSE IF (Psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
    IF (lformated_loc) THEN
        DO i_e=1,Psi%nb_be
        DO i_h=1,Psi%nb_bi
           i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             IF (Psi%cplx) THEN
               a = real(Psi%CvecB(i_bhe),kind=Rkind)
               b = aimag(Psi%CvecB(i_bhe))
               IF (sqrt(a*a+b*b) >= epsi) write(nioPsi,*) i_b,i_h,i_e,a,b
             ELSE
               a = Psi%RvecB(i_bhe)
!              write(out_unit,*) indnD(:),a
               IF (abs(a) >= epsi) write(nioPsi,*) i_b,i_h,i_e,a
             END IF
           END DO
        END DO
        END DO
        write(nioPsi,*) 'end wp ',ipsi
    ELSE
        DO i_e=1,Psi%nb_be
        DO i_h=1,Psi%nb_bi
           i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             IF (Psi%cplx) THEN
               a = real(Psi%CvecB(i_bhe),kind=Rkind)
               b = aimag(Psi%CvecB(i_bhe))
               IF (sqrt(a*a+b*b) >= epsi) write(nioPsi,*) i_b,i_h,i_e,a,b
             ELSE
               a = Psi%RvecB(i_bhe)
!              write(out_unit,*) indnD(:),a
               IF (abs(a) >= epsi) write(nioPsi) i_b,i_h,i_e,a
             END IF
           END DO
        END DO
        END DO
        !write(nioPsi) 'end wp ',ipsi
    END IF
  END IF

CASE(1)
  IF (Psi%cplx) THEN
    IF (lformated_loc) THEN
      write(nioPsi,*,iostat=Rerr) T_loc,Psi%CvecB(:)
    ELSE
      write(nioPsi,iostat=Rerr) T_loc,Psi%CvecB(:)
    END IF
  ELSE
    IF (lformated_loc) THEN
      write(nioPsi,*,iostat=Rerr) T_loc,Psi%RvecB(:)
    ELSE
      write(nioPsi,iostat=Rerr) T_loc,Psi%RvecB(:)
    END IF
  END IF

  IF (Rerr /= 0) THEN
    write(out_unit,*) 'ERROR in ',name_sub
    write(out_unit,*) 'Problem while writing the WP with version=',version
    STOP
  END IF

CASE(2)
  T_loc = ZERO
  IF (present(T)) T_loc = T

  IF (Psi%nb_tot == Psi%nb_baie) THEN

    IF (lformated_loc) THEN
      i_bhe = 0
      DO i_e=1,Psi%nb_be
      DO i_h=1,Psi%nb_bi
      DO i_b=1,Psi%nb_ba
      i_bhe = i_bhe + 1

        CALL Rec_ndim_index(Psi%BasisnD,ind_ndim(:),i_b)

        IF (Psi%cplx) THEN
          a = real(Psi%CvecB(i_bhe),kind=Rkind)
          b = aimag(Psi%CvecB(i_bhe))
          IF (sqrt(a*a+b*b) >= epsi) write(nioPsi,*) T_loc,ind_ndim(:),i_h,i_e,a,b
        ELSE
          a = Psi%RvecB(i_bhe)
          IF (abs(a) >= epsi) write(nioPsi,*) T_loc,ind_ndim(:),i_h,i_e,a
        END IF
      END DO
      END DO
      END DO
      write(nioPsi,*) 'end wp ',ipsi

    ELSE ! not formatted
      i_bhe = 0
      DO i_e=1,Psi%nb_be
      DO i_h=1,Psi%nb_bi
      DO i_b=1,Psi%nb_ba
      i_bhe = i_bhe + 1

        CALL Rec_ndim_index(Psi%BasisnD,ind_ndim(:),i_b)


        IF (Psi%cplx) THEN
          a = real(Psi%CvecB(i_b),kind=Rkind)
          b = aimag(Psi%CvecB(i_b))
          IF (sqrt(a*a+b*b) >= epsi) write(nioPsi) ind_ndim(:),i_h,i_e,a,b
        ELSE
          a = Psi%RvecB(i_b)
!         write(out_unit,*) indnD(:),a
          IF (abs(a) >= epsi) write(nioPsi) ind_ndim(:),i_h,i_e,a
        END IF
      END DO
      END DO
      END DO
      !write(nioPsi) 'end wp ',ipsi
   END IF

  ELSE IF (Psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
    IF (lformated_loc) THEN
        DO i_e=1,Psi%nb_be
        DO i_h=1,Psi%nb_bi
           i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             IF (Psi%cplx) THEN
               a = real(Psi%CvecB(i_bhe),kind=Rkind)
               b = aimag(Psi%CvecB(i_bhe))
               IF (sqrt(a*a+b*b) >= epsi) write(nioPsi,*) T_loc,i_b,i_h,i_e,a,b
             ELSE
               a = Psi%RvecB(i_bhe)
!              write(out_unit,*) indnD(:),a
               IF (abs(a) >= epsi) write(nioPsi,*) T_loc,i_b,i_h,i_e,a
             END IF
           END DO
        END DO
        END DO
        write(nioPsi,*) 'end wp ',ipsi
    ELSE
        DO i_e=1,Psi%nb_be
        DO i_h=1,Psi%nb_bi
           i_bhe = sum(Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(1:i_h-1))
           DO i_b=1,Psi%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i_h)
             i_bhe = i_bhe + 1
             IF (Psi%cplx) THEN
               a = real(Psi%CvecB(i_bhe),kind=Rkind)
               b = aimag(Psi%CvecB(i_bhe))
               IF (sqrt(a*a+b*b) >= epsi) write(nioPsi,*) T_loc,i_b,i_h,i_e,a,b
             ELSE
               a = Psi%RvecB(i_bhe)
!              write(out_unit,*) indnD(:),a
               IF (abs(a) >= epsi) write(nioPsi) T_loc,i_b,i_h,i_e,a
             END IF
           END DO
        END DO
        END DO
        !write(nioPsi) 'end wp ',ipsi
    END IF
  END IF
CASE Default
  STOP 'no default in Write_psi_nDBasis'
END SELECT
!-----------------------------------------------------------
!-----------------------------------------------------------

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'PsiBasisRep'
         CALL ecri_psi(T=ZERO,psi=Psi)
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
      END SUBROUTINE Write_psi_nDBasis
!=======================================================================================

      END MODULE mod_psi_io
