!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
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
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
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
      USE mod_system
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
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_ba,nb_qa',psi%nb_ba,psi%nb_qa
        write(out_unitp,*) 'nb_act1',psi%nb_act1
        write(out_unitp,*) 'asso BasisnD ',associated(psi%BasisnD)
        write(out_unitp,*) 'nb_basis',psi%BasisnD%nb_basis
        write(out_unitp,*)
        !write(out_unitp,*) 'psi BasisRep'
        !CALL ecri_psi(ZERO,psi,ecri_BasisRep=.TRUE.,ecri_GridRep=.FALSE.)
      END IF
!-----------------------------------------------------------

      IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' HADA contraction: ',                       &
                         psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC
        write(out_unitp,*) ' and ',name_sub,' is not possible'
        STOP
      END IF

!------ initisalisation ----------------------------------
      CALL alloc_psi(psi,GridRep=.TRUE.)
!------ end initisalisation -------------------------------

      IF (psi%cplx) THEN
        psi%CvecG(:) = CZERO
        IF ( .NOT. allocated(psi%CvecB) ) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' psi%CvecB MUST be allocated !!'
          STOP
        END IF
        IF (NewBasisEl) THEN
          CALL RecCvecB_TO_CvecG(psi%CvecB,psi%CvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE IF (psi%BasisnD%SparseGrid_type == 4) THEN ! specail case with SG4
          CALL RecCvecB_TO_CvecG(psi%CvecB,psi%CvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE
          iqaie0 = 1
          ibaie0 = 1
          iqaie1 = psi%nb_qa
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecCvecB_TO_CvecG(psi%CvecB(ibaie0:ibaie1),            &
                                   psi%CvecG(iqaie0:iqaie1),            &
                                   psi%nb_ba,psi%nb_qa,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF
      ELSE
        psi%RvecG(:) = ZERO
        IF ( .NOT. allocated(psi%RvecB) ) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' psi%RvecB MUST be allocated !!'
          STOP
        END IF
        IF (NewBasisEl) THEN
          CALL RecRvecB_TO_RvecG(psi%RvecB,psi%RvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE IF (psi%BasisnD%SparseGrid_type == 4) THEN ! specail case with SG4
          CALL RecRvecB_TO_RvecG(psi%RvecB,psi%RvecG,psi%nb_ba,psi%nb_qa,psi%BasisnD)
        ELSE
          iqaie0 = 1
          ibaie0 = 1
          iqaie1 = psi%nb_qa
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecRvecB_TO_RvecG(psi%RvecB(ibaie0:ibaie1),            &
                                   psi%RvecG(iqaie0:iqaie1),            &
                                   psi%nb_ba,psi%nb_qa,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF
      END IF

!----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*)
         !write(out_unitp,*) 'psiGridRep'
         !CALL ecri_psi(ZERO,psi,out_unitp,ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
         write(out_unitp,*)
         write(out_unitp,*) ' END in ',name_sub
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
      USE mod_system
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
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_ba,nb_qa',psi%nb_ba,psi%nb_qa
        write(out_unitp,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
        write(out_unitp,*) 'nb_act1',psi%nb_act1
        write(out_unitp,*)
        write(out_unitp,*) 'nb_basis',psi%BasisnD%nb_basis
        write(out_unitp,*)

        write(out_unitp,*) 'psi GridRep'
        CALL ecri_psi(ZERO,psi,out_unitp,.TRUE.,.FALSE.)
        flush(out_unitp)
      END IF
!-----------------------------------------------------------

      IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' HADA contraction:',psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC
        write(out_unitp,*) ' and ',name_sub,' is not possible'
        STOP
      END IF


!---- initisalisation ----------------------------------
      CALL alloc_psi(psi,BasisRep=.TRUE.)
!---- end initisalisation -------------------------------
      IF (psi%cplx) THEN
        IF ( .NOT. allocated(psi%CvecG) .AND. keep_MPI) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' psi%CvecG MUST be allocated !!'
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
          iqaie1 = psi%nb_qa
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecCvecG_TO_CvecB(psi%CvecG(iqaie0:iqaie1),            &
                                       psi%CvecB(ibaie0:ibaie1),        &
                                       psi%nb_qa,psi%nb_ba,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF

      ELSE
        IF ( .NOT. allocated(psi%RvecG) .AND. keep_MPI) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' psi%RvecG MUST be allocated !!'
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
          iqaie1 = psi%nb_qa
          ibaie1 = psi%nb_ba
          DO ibie=1,psi%nb_bi*psi%nb_be
            CALL RecRvecG_TO_RvecB(psi%RvecG(iqaie0:iqaie1),            &
                                     psi%RvecB(ibaie0:ibaie1),          &
                                     psi%nb_qa,psi%nb_ba,psi%BasisnD)
            iqaie0 = iqaie0 + psi%nb_qa
            ibaie0 = ibaie0 + psi%nb_ba
            iqaie1 = iqaie1 + psi%nb_qa
            ibaie1 = ibaie1 + psi%nb_ba
          END DO
        END IF
      END IF

!----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'psiBasisRep'
         CALL ecri_psi(ZERO,psi)
         write(out_unitp,*)
         write(out_unitp,*) 'sub_PsiGridRep_TO_BasisRep'
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
      USE mod_system
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
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_ba,nb_qa',psi%nb_ba,psi%nb_qa
        write(out_unitp,*) 'nb_act1',psi%nb_act1
        write(out_unitp,*)
        write(out_unitp,*) 'nb_basis',psi%BasisnD%nb_basis
        write(out_unitp,*) 'tab_derQdyn',tab_derQdyn
        write(out_unitp,*)
        write(out_unitp,*) 'psi BasisRep'
        CALL ecri_psi(Psi=psi)
      END IF
!-----------------------------------------------------------

      IF (psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' HADA contraction: ',                       &
                         psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC
        write(out_unitp,*) '  and ',name_sub,' is not possible'
        STOP
      END IF


!------ initisalisation ----------------------------------
      psi%GridRep = .TRUE.
      CALL alloc_psi(psi)
!------ end initisalisation -------------------------------

      IF (psi%cplx) THEN
        IF ( .NOT. allocated(psi%CvecB) ) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' psi%CvecB MUST be allocated !!'
           STOP
        END IF

        psi%CvecG(:) = CZERO
        iqaie0 = 1
        ibaie0 = 1
        iqaie1 = psi%nb_qa
        ibaie1 = psi%nb_ba
        DO ibie=1,psi%nb_bi*psi%nb_be
          CALL RecCvecB_TO_CvecG(psi%CvecB(ibaie0:ibaie1),              &
                                 psi%CvecG(iqaie0:iqaie1),              &
                                 psi%nb_ba,psi%nb_qa,psi%BasisnD,       &
                                 tab_derQdyn)
          iqaie0 = iqaie0 + psi%nb_qa
          ibaie0 = ibaie0 + psi%nb_ba
          iqaie1 = iqaie1 + psi%nb_qa
          ibaie1 = ibaie1 + psi%nb_ba
        END DO
      ELSE
        IF ( .NOT. allocated(psi%RvecB) ) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' psi%RvecB MUST be allocated !!'
          STOP
        END IF

        psi%RvecG(:) = ZERO
        iqaie0 = 1
        ibaie0 = 1
        iqaie1 = psi%nb_qa
        ibaie1 = psi%nb_ba
        DO ibie=1,psi%nb_bi*psi%nb_be
          CALL RecRvecB_TO_RvecG(psi%RvecB(ibaie0:ibaie1),              &
                                 psi%RvecG(iqaie0:iqaie1),              &
                                 psi%nb_ba,psi%nb_qa,psi%BasisnD,       &
                                 tab_derQdyn)
          iqaie0 = iqaie0 + psi%nb_qa
          ibaie0 = ibaie0 + psi%nb_ba
          iqaie1 = iqaie1 + psi%nb_qa
          ibaie1 = ibaie1 + psi%nb_ba
        END DO
      END IF

!----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'psiGridRep'
         CALL ecri_psi(Psi=psi)
         write(out_unitp,*)
         write(out_unitp,*) 'END in ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE sub_d0d1d2PsiBasisRep_TO_GridRep

END MODULE mod_psi_B_TO_G
