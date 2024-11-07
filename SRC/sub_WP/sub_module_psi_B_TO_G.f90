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
MODULE mod_psi_B_TO_G
USE mod_basis
IMPLICIT NONE

PRIVATE
PUBLIC :: sub_PsiBasisRep_TO_GridRep,sub_PsiGridRep_TO_BasisRep,        &
          sub_d0d1d2PsiBasisRep_TO_GridRep

CONTAINS


!================================================================
!
!     transformation BasisRep to GridRep
!
!================================================================

      SUBROUTINE sub_PsiBasisRep_TO_GridRep(psi)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi), intent(inout)   :: psi

!------ working variables ---------------------------------
      integer       :: ibaie0,iqaie0,ibaie1,iqaie1,ibie


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_PsiBasisRep_TO_GridRep'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_ba,nb_qa_WithNoGrid',psi%nb_ba,psi%nb_qa_WithNoGrid
        write(out_unit,*) 'nb_act1',psi%nb_act1
        write(out_unit,*) 'asso BasisnD ',associated(psi%BasisnD)
        write(out_unit,*) 'nb_basis',psi%BasisnD%nb_basis
        write(out_unit,*)
        !write(out_unit,*) 'psi BasisRep'
        !CALL ecri_psi(ZERO,psi,ecri_BasisRep=.TRUE.,ecri_GridRep=.FALSE.)
      END IF
!-----------------------------------------------------------

      IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' HADA contraction: ',                       &
                         psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC
        write(out_unit,*) ' and ',name_sub,' is not possible'
        STOP
      END IF

!------ initisalisation ----------------------------------
      CALL alloc_psi(psi,GridRep=.TRUE.)
!------ end initisalisation -------------------------------

      IF (psi%cplx) THEN
        psi%CvecG(:) = CZERO
        IF ( .NOT. allocated(psi%CvecB) ) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' psi%CvecB MUST be allocated !!'
          STOP
        END IF
        IF (NewBasisEl) THEN
          CALL RecCvecB_TO_CvecG(psi%CvecB,psi%CvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE IF (psi%BasisnD%SparseGrid_type == 4) THEN ! special case with SG4
          CALL RecCvecB_TO_CvecG(psi%CvecB,psi%CvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE
          iqaie0 = 1
          ibaie0 = 1
          iqaie1 = psi%nb_qa_WithNoGrid
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecCvecB_TO_CvecG(psi%CvecB(ibaie0:ibaie1),            &
                                   psi%CvecG(iqaie0:iqaie1),            &
                                   psi%nb_ba,psi%nb_qa_WithNoGrid,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa_WithNoGrid
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa_WithNoGrid
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF
      ELSE
        psi%RvecG(:) = ZERO
        IF ( .NOT. allocated(psi%RvecB) ) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' psi%RvecB MUST be allocated !!'
          STOP
        END IF
        IF (NewBasisEl) THEN
          CALL RecRvecB_TO_RvecG(psi%RvecB,psi%RvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE IF (psi%BasisnD%SparseGrid_type == 4) THEN ! specail case with SG4
          CALL RecRvecB_TO_RvecG(psi%RvecB,psi%RvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE
          iqaie0 = 1
          ibaie0 = 1
          iqaie1 = psi%nb_qa_WithNoGrid
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecRvecB_TO_RvecG(psi%RvecB(ibaie0:ibaie1),            &
                                   psi%RvecG(iqaie0:iqaie1),            &
                                   psi%nb_ba,psi%nb_qa_WithNoGrid,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa_WithNoGrid
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa_WithNoGrid
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF
      END IF

!----------------------------------------------------------
      IF (debug) THEN
         write(out_unit,*)
         !write(out_unit,*) 'psiGridRep'
         !CALL ecri_psi(ZERO,psi,out_unit,ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
         write(out_unit,*)
         write(out_unit,*) ' END in ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE sub_PsiBasisRep_TO_GridRep


!================================================================
!
!     transformation GridRep to BasisRep
!
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE sub_PsiGridRep_TO_BasisRep(psi)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi), intent(inout)   :: psi

!------ working variables ---------------------------------
      integer       :: ibaie0,iqaie0,ibaie1,iqaie1,ibie

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_PsiGridRep_TO_BasisRep'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_ba,nb_qa_WithNoGrid',psi%nb_ba,psi%nb_qa_WithNoGrid
        write(out_unit,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
        write(out_unit,*) 'nb_act1',psi%nb_act1
        write(out_unit,*)
        write(out_unit,*) 'nb_basis',psi%BasisnD%nb_basis
        write(out_unit,*)

        write(out_unit,*) 'psi GridRep'
        CALL ecri_psi(ZERO,psi,out_unit,.TRUE.,.FALSE.)
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' HADA contraction:',psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC
        write(out_unit,*) ' and ',name_sub,' is not possible'
        STOP
      END IF


!---- initisalisation ----------------------------------
      CALL alloc_psi(psi,BasisRep=.TRUE.)
!---- end initisalisation -------------------------------
      IF (psi%cplx) THEN
        IF ( .NOT. allocated(psi%CvecG) .AND. keep_MPI) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' psi%CvecG MUST be allocated !!'
          STOP
        END IF
        IF(keep_MPI)psi%CvecB(:) = CZERO
        IF (NewBasisEl) THEN
          CALL RecCvecG_TO_CvecB(psi%CvecG,psi%CvecB,psi%nb_qaie,psi%nb_baie,psi%BasisnD)
        ELSE IF (psi%BasisnD%SparseGrid_type == 4) THEN ! special case with SG4
          CALL RecCvecG_TO_CvecB(psi%CvecG,psi%CvecB,psi%nb_qaie,psi%nb_baie,psi%BasisnD)
        ELSE
          iqaie0 = 1
          ibaie0 = 1
          iqaie1 = psi%nb_qa_WithNoGrid
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecCvecG_TO_CvecB(psi%CvecG(iqaie0:iqaie1),            &
                                       psi%CvecB(ibaie0:ibaie1),        &
                                       psi%nb_qa_WithNoGrid,psi%nb_ba,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa_WithNoGrid
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa_WithNoGrid
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF

      ELSE
        IF ( .NOT. allocated(psi%RvecG) .AND. keep_MPI) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' psi%RvecG MUST be allocated !!'
          STOP
        END IF
        psi%RvecB(:) = ZERO
        IF (NewBasisEl) THEN
          CALL RecRvecG_TO_RvecB(psi%RvecG,psi%RvecB,psi%nb_qaie,psi%nb_baie,psi%BasisnD)
        ELSE IF (psi%BasisnD%SparseGrid_type == 4) THEN ! specail case with SG4
          CALL RecRvecG_TO_RvecB(psi%RvecG,psi%RvecB,psi%nb_qaie,psi%nb_baie,psi%BasisnD)
        ELSE
          iqaie0 = 1
          ibaie0 = 1
          iqaie1 = psi%nb_qa_WithNoGrid
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecRvecG_TO_RvecB(psi%RvecG(iqaie0:iqaie1),            &
                                     psi%RvecB(ibaie0:ibaie1),          &
                                     psi%nb_qa_WithNoGrid,psi%nb_ba,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa_WithNoGrid
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa_WithNoGrid
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF
      END IF

!----------------------------------------------------------
      IF (debug) THEN
         write(out_unit,*)
         write(out_unit,*) 'psiBasisRep'
         CALL ecri_psi(ZERO,psi)
         write(out_unit,*)
         write(out_unit,*) 'sub_PsiGridRep_TO_BasisRep'
      END IF
!----------------------------------------------------------

      END SUBROUTINE sub_PsiGridRep_TO_BasisRep


!================================================================
!
!     transformation BasisRep to GridRep with derivative of psiBasisRep
!
!================================================================
      !!@description: ransformation BasisRep to GridRep with derivative of psiBasisRep
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE sub_d0d1d2PsiBasisRep_TO_GridRep(psi,tab_derQdyn)
      USE EVR_system_m
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi), intent(inout)   :: psi
      integer,          intent(in)      :: tab_derQdyn(2)

!------ working variables ---------------------------------
      integer       :: ibaie0,iqaie0,ibaie1,iqaie1,ibie


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_d0d1d2PsiBasisRep_TO_GridRep'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_ba,nb_qa_WithNoGrid',psi%nb_ba,psi%nb_qa_WithNoGrid
        write(out_unit,*) 'nb_act1',psi%nb_act1
        write(out_unit,*)
        write(out_unit,*) 'nb_basis',psi%BasisnD%nb_basis
        write(out_unit,*) 'tab_derQdyn',tab_derQdyn
        write(out_unit,*)
        write(out_unit,*) 'psi BasisRep'
        CALL ecri_psi(Psi=psi)
      END IF
!-----------------------------------------------------------

      IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' HADA contraction: ',                       &
                         psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC
        write(out_unit,*) '  and ',name_sub,' is not possible'
        STOP
      END IF


!------ initisalisation ----------------------------------
      psi%GridRep = .TRUE.
      CALL alloc_psi(psi)
!------ end initisalisation -------------------------------

      IF (psi%cplx) THEN
        IF ( .NOT. allocated(psi%CvecB) ) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) ' psi%CvecB MUST be allocated !!'
           STOP
        END IF

        psi%CvecG(:) = CZERO
        iqaie0 = 1
        ibaie0 = 1
        iqaie1 = psi%nb_qa_WithNoGrid
        ibaie1 = psi%nb_ba
        DO ibie=1,psi%nb_bi*psi%nb_be
          CALL RecCvecB_TO_CvecG(psi%CvecB(ibaie0:ibaie1),              &
                                 psi%CvecG(iqaie0:iqaie1),              &
                                 psi%nb_ba,psi%nb_qa_WithNoGrid,psi%BasisnD,       &
                                 tab_derQdyn)
          iqaie0 = iqaie0 + psi%nb_qa_WithNoGrid
          ibaie0 = ibaie0 + psi%nb_ba
          iqaie1 = iqaie1 + psi%nb_qa_WithNoGrid
          ibaie1 = ibaie1 + psi%nb_ba
        END DO
      ELSE
        IF ( .NOT. allocated(psi%RvecB) ) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' psi%RvecB MUST be allocated !!'
          STOP
        END IF

        psi%RvecG(:) = ZERO
        iqaie0 = 1
        ibaie0 = 1
        iqaie1 = psi%nb_qa_WithNoGrid
        ibaie1 = psi%nb_ba
        DO ibie=1,psi%nb_bi*psi%nb_be
          CALL RecRvecB_TO_RvecG(psi%RvecB(ibaie0:ibaie1),              &
                                 psi%RvecG(iqaie0:iqaie1),              &
                                 psi%nb_ba,psi%nb_qa_WithNoGrid,psi%BasisnD,       &
                                 tab_derQdyn)
          iqaie0 = iqaie0 + psi%nb_qa_WithNoGrid
          ibaie0 = ibaie0 + psi%nb_ba
          iqaie1 = iqaie1 + psi%nb_qa_WithNoGrid
          ibaie1 = ibaie1 + psi%nb_ba
        END DO
      END IF

!----------------------------------------------------------
      IF (debug) THEN
         write(out_unit,*)
         write(out_unit,*) 'psiGridRep'
         CALL ecri_psi(Psi=psi)
         write(out_unit,*)
         write(out_unit,*) 'END in ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE sub_d0d1d2PsiBasisRep_TO_GridRep

END MODULE mod_psi_B_TO_G
