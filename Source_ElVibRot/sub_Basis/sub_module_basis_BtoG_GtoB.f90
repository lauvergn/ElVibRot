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
      MODULE mod_basis_BtoG_GtoB
      USE mod_system
      USE mod_basis_set_alloc
      USE mod_basis_BtoG_GtoB_SGType2
      USE mod_param_SGType2
      USE mod_basis_BtoG_GtoB_MPI,ONLY:CVecB_TO_CVecG_R_MPI,CVecB_TO_CVecG_C_MPI
      IMPLICIT NONE

      CONTAINS
!      ==========================================================
!       VecGridRep <=> VecBasisRep
!      ==========================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
       RECURSIVE SUBROUTINE RecRVecG_TO_RvecB(RVecG,RvecB,nq,nb,basis_set)
       USE mod_basis_BtoG_GtoB_SGType4
       USE mod_basis_RCVec_SGType4
        IMPLICIT NONE


        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq,nb
        real (kind=Rkind), intent(inout) :: RVecG(:)
        real (kind=Rkind), intent(inout) :: RvecB(:)

        real (kind=Rkind), allocatable   :: RTempB(:,:)
        real (kind=Rkind), allocatable   :: RTempG(:,:,:)

        integer                          :: nbb,ibb1,ibb2

        real (kind=Rkind), allocatable   :: RB(:)
        real (kind=Rkind), allocatable   :: RG(:)

        integer                          :: nnb,nb2,ib,ib2,newnb2

        integer                          :: nnq,nq2,iq,iq2,nqi,nbi
        integer                          :: ibasis
        integer                          :: i_SG,iq0_SG,iq1_SG,nq_SG
        integer                          :: nb_thread
        integer                          :: der00(2) = [0,0]

        TYPE(Type_SmolyakRep)            :: SRep ! smolyak rep for SparseGrid_type=4
        integer                          :: itabR,iG,nR

!----- for debuging --------------------------------------------------
        integer :: err_mem,memory
        character (len=*), parameter :: name_sub='RecRVecG_TO_RvecB'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'RVecG(:)',RVecG(:)
        END IF

        IF (basis_set%packed_done) THEN
          CALL RG_TO_RB_basis(RVecG,RvecB,basis_set)
        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product

            nnq  = nq
            nnb  = 1
            nb2  = 1
            nbb  = 1
            nq2  = 1

            CALL alloc_NParray(RTempB,[nnq,nbb],"RTempB",name_sub)
            RTempB(:,:) = reshape(RVecG,shape=[nnq,nbb])


            DO ibasis=1,basis_set%nb_basis

              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnq = nnq / nq2

              CALL alloc_NParray(RTempG,[nnq,nq2,nnb],"RTempG",name_sub)
              RTempG(:,:,:) = reshape(RTempB,shape=[nnq,nq2,nnb])
              !write(out_unitp,*) 'ibasis shape G',ibasis,shape(RTempG)
              nbb = sum(basis_set%Tab_OF_Tabnb2(ibasis)%vec)

!              write(out_unitp,*) 'G=>B ibasis,nnb,Tab_OF_Tabnb2',ibasis,nnb,   &
!                                        basis_set%Tab_OF_Tabnb2(ibasis)%vec

              CALL dealloc_NParray(RTempB,"RTempB",name_sub)
              CALL alloc_NParray(RTempB,[nnq,nbb],"RTempB",name_sub)
              !write(out_unitp,*) 'ibasis shape B',ibasis,shape(RTempB)


              CALL alloc_NParray(RG,[nq2],"RG",name_sub)
              CALL alloc_NParray(RB,[nb2],"RB",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq
                    CALL RG_TO_RB_basis(RTempG(iq,:,ib),RTempB(iq,ibb1:ibb2),  &
                                        basis_set%tab_Pbasis(ibasis)%Pbasis)
                  END DO
                  ibb1 = ibb1 + newnb2
                END DO

              ELSE

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq

                  RG(:) = RTempG(iq,:,ib)
                  CALL RecRVecG_TO_RvecB(RG,RB,nq2,newnb2,basis_set%tab_Pbasis(ibasis)%Pbasis)

                  RTempB(iq,ibb1:ibb2) = RB(1:newnb2)

                  END DO

                  ibb1 = ibb1 + newnb2
                END DO
              END IF

              CALL dealloc_NParray(RB,"RB",name_sub)
              CALL dealloc_NParray(RG,"RG",name_sub)

              CALL dealloc_NParray(RTempG,"RTempG",name_sub)

              nnb = nbb

            END DO

            RvecB(:) = reshape(RTempB, shape=[nnq*nnb] )
            CALL dealloc_NParray(RTempB,"RTempB",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            RvecB(:) = ZERO
            CALL alloc_NParray(RB,[nb],"RB",name_sub)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
              nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq1_SG = iq1_SG + nq_SG
              CALL RecRVecG_TO_RvecB(RVecG(iq0_SG:iq1_SG),RB,nq_SG,nb,  &
                                     basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq0_SG = iq0_SG + nq_SG
              RvecB(:) = RvecB(:) + RB(:) * basis_set%WeightSG(i_SG)

            END DO
            CALL dealloc_NParray(RB,"RB",name_sub)

          CASE (2) ! Sparse basis (Smolyak 2d implementation)

            CALL sub_G_TO_B(RVecG,RVecB,basis_set%WeightSG,             &
                      basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)

            !write(out_unitp,*) ' ERROR in ',name_sub
            !STOP 'SparseGrid_type=2'

          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            !!  RVecG TO SRep
            CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                   nb0=basis_set%para_SGType2%nb0)
            CALL tabR2bis_TO_SmolyakRep1(SRep,RVecG) ! on the grid

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL GSmolyakRep_TO_BSmolyakRep(SRep,                     &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,  &
                    basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL GSmolyakRep_TO3_BSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF

            CALL SmolyakRepBasis_TO_tabPackedBasis(SRep,RVecB,          &
                             basis_set%tab_basisPrimSG,basis_set%nDindB,&
                             basis_set%para_SGType2,basis_set%WeightSG)

            CALL dealloc_SmolyakRep(SRep)

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'RvecB(:)',RvecB(:)
         write(out_unitp,*) 'END ',name_sub
       END IF
      END SUBROUTINE RecRVecG_TO_RvecB


      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecCVecG_TO_CVecB(CVecG,CVecB,nq,nb,basis_set)
       USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE

        TYPE (basis), intent(in)            :: basis_set
        integer, intent(in)                 :: nq,nb
        complex (kind=Rkind), intent(inout) :: CVecG(:)
        complex (kind=Rkind), intent(inout) :: CVecB(:)

        complex (kind=Rkind), allocatable   :: CTempG(:,:,:)
        complex (kind=Rkind), allocatable   :: CTempB(:,:)

        real (kind=Rkind), allocatable :: RVecG(:)
        real (kind=Rkind), allocatable :: RVecB(:)

        TYPE(Type_SmolyakRep)            :: SRep ! smolyak rep for SparseGrid_type=4
        integer                          :: itabR,iG,nR

        integer                             :: nbb,ibb1,ibb2

        complex (kind=Rkind), allocatable       :: CB(:)
        complex (kind=Rkind), allocatable       :: CG(:)
        complex (kind=Rkind), allocatable       :: d0cb(:,:)
        real    (kind=Rkind), allocatable       :: d0b(:,:)
        real    (kind=Rkind), allocatable       :: w(:)

        integer                             :: nnb,nb2,ib,ib2,newnb2
        integer                             :: nnq,nq2,iq,iq2

        integer                             :: ibasis
        integer                             :: i_SG,iq0_SG,iq1_SG,nq_SG
        integer                             :: nb_thread
        integer                             :: der00(2) = [0,0]

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecCVecG_TO_CVecB'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' BEGINNING in ',name_sub
        write(out_unitp,*) 'nq,nb',nq,nb
        write(out_unitp,*) 'shape CVecG',shape(CVecG)
        write(out_unitp,*) 'shape CVecB',shape(CVecB)
        flush(out_unitp)
      END IF

        IF (basis_set%packed_done) THEN
          CALL CG_TO_CB_basis(CVecG,CvecB,basis_set)
        ELSE
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product
            nnq  = nq
            nnb  = 1
            nb2  = 1
            nbb  = 1
            nq2  = 1

            CALL alloc_NParray(CTempB,[nnq,nbb],"CTempB",name_sub)
            CTempB(:,:) = reshape(CVecG,shape=[nnq,nbb])

            DO ibasis=1,basis_set%nb_basis
              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)

              nnq = nnq / nq2

              CALL alloc_NParray(CTempG,[nnq,nq2,nnb],"CTempG",name_sub)

              CTempG(:,:,:) = reshape(CTempB,shape=[nnq,nq2,nnb])

              nbb = sum(basis_set%Tab_OF_Tabnb2(ibasis)%vec)

              !write(out_unitp,*) 'G=>B ibasis,nnb,Tab_OF_Tabnb2',ibasis,nnb,      &
              !                          basis_set%Tab_OF_Tabnb2(ibasis)%vec

              CALL dealloc_NParray(CTempB,"CTempB",name_sub)
              CALL alloc_NParray(CTempB,[nnq,nbb],"CTempB",name_sub)


              CALL alloc_NParray(CG,[nq2],"CG",name_sub)
              CALL alloc_NParray(CB,[nb2],"CB",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN
                  ibb1 = 1
                  ibb2 = 0
                  DO ib=1,nnb
                    newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                    ibb2 = ibb2 + newnb2

                    DO iq=1,nnq
                      CALL CG_TO_CB_basis(CTempG(iq,:,ib),CTempB(iq,ibb1:ibb2),basis_set%tab_Pbasis(ibasis)%Pbasis)
                    END DO

                    ibb1 = ibb1 + newnb2
                  END DO
              ELSE

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq

                    CG(:) = CTempG(iq,:,ib)
                    CALL RecCVecG_TO_CVecB(CG,CB,nq2,newnb2,          &
                                          basis_set%tab_Pbasis(ibasis)%Pbasis)

                    CTempB(iq,ibb1:ibb2) = CB(1:newnb2)

                  END DO

                  ibb1 = ibb1 + newnb2
                END DO
              END IF

              CALL dealloc_NParray(CB,"CB",name_sub)
              CALL dealloc_NParray(CG,"CG",name_sub)
              CALL dealloc_NParray(CTempG,"CTempG",name_sub)

              nnb = nbb
            END DO

            CVecB(:) = reshape(CTempB, shape=[nnq*nnb] )
            CALL dealloc_NParray(CTempB,"CTempB",name_sub)
          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            CVecB(:) = CZERO
            CALL alloc_NParray(CB,[nb],"CB",name_sub)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
              nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq1_SG = iq1_SG + nq_SG
              CALL RecCVecG_TO_CVecB(CVecG(iq0_SG:iq1_SG),CB,           &
                                           nq_SG,nb,                    &
                                           basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq0_SG = iq0_SG + nq_SG
              CVecB(:) = CVecB(:) + CB(:) * cmplx(basis_set%WeightSG(i_SG),kind=Rkind)
            END DO
            CALL dealloc_NParray(CB,"CB",name_sub)

          CASE (2) ! Sparse basis (Smolyak 2d  implementation)

            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)
            CALL alloc_NParray(RVecG,shape(CVecG),'RVecG',name_sub)

            RVecG(:) = real(CVecG,kind=Rkind)

            CALL sub_G_TO_B(RVecG,RVecB,basis_set%WeightSG,             &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecB(:) = cmplx(RVecB,kind=Rkind)

            RVecG(:) = aimag(CVecG)
            CALL sub_G_TO_B(RVecG,RVecB,basis_set%WeightSG,             &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecB(:) = CVecB + EYE*cmplx(RVecB,kind=Rkind)

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)
            CALL dealloc_NParray(RVecG,'RVecG',name_sub)

          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)

            !!  RVecG TO SRep
            CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                   nb0=basis_set%para_SGType2%nb0)

            IF(keep_MPI) CALL tabR2bis_TO_SmolyakRep1(SRep,real(CVecG,kind=Rkind)) ! on the grid

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL GSmolyakRep_TO_BSmolyakRep(SRep,                        &
                        basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval, &
                        basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)

            ELSE
              CALL GSmolyakRep_TO3_BSmolyakRep(SRep,basis_set%para_SGType2,&
                                               basis_set%tab_basisPrimSG)

            END IF

            CALL SmolyakRepBasis_TO_tabPackedBasis(SRep,RVecB,          &
                             basis_set%tab_basisPrimSG,basis_set%nDindB,&
                             basis_set%para_SGType2,basis_set%WeightSG)

            CALL dealloc_SmolyakRep(SRep)

            IF(keep_MPI) CVecB(:) = cmplx(RVecB,kind=Rkind)


            !!  RVecG TO SRep
            CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                   nb0=basis_set%para_SGType2%nb0)
            IF(keep_MPI) CALL tabR2bis_TO_SmolyakRep1(SRep,aimag(CVecG)) ! on the grid
            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL GSmolyakRep_TO_BSmolyakRep(SRep,                     &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,  &
                    basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL GSmolyakRep_TO3_BSmolyakRep(SRep,basis_set%para_SGType2,&
                                               basis_set%tab_basisPrimSG)
            END IF

            CALL SmolyakRepBasis_TO_tabPackedBasis(SRep,RVecB,          &
                             basis_set%tab_basisPrimSG,basis_set%nDindB,&
                             basis_set%para_SGType2,basis_set%WeightSG)

            CALL dealloc_SmolyakRep(SRep)

            IF(keep_MPI) CVecB(:) = CVecB + EYE*cmplx(RVecB,kind=Rkind)

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT
        END IF

      IF (debug) write(out_unitp,*) ' END in ',name_sub


      END SUBROUTINE RecCVecG_TO_CVecB

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecRvecB_TO_RVecG(RvecB,RVecG,nb,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE
        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq,nb
        real (kind=Rkind), intent(inout) :: RVecG(:)
        real (kind=Rkind), intent(in)    :: RvecB(:)
        integer, optional                :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        real (kind=Rkind), allocatable       :: RTempG(:,:,:)
        real (kind=Rkind), allocatable       :: RTempB(:,:)
        TYPE(Type_SmolyakRep)                :: SRep ! smolyak rep for SparseGrid_type=4

        real (kind=Rkind), allocatable       :: RG(:)
        real (kind=Rkind), allocatable       :: RB(:)
        real (kind=Rkind), allocatable       :: dnb(:,:)

        integer                          :: nbb,ibb1,ibb2
        integer                          :: nnb,nb2,ib,ib2,newnb2
        integer                          :: nnq,nq2,iq,iq2
        integer                          :: ibasis
        integer                          :: i_SG,iq0_SG,iq1_SG,nq_SG,ii_SG,nq0_SG

        integer                          :: itabR,iG,nR


        integer                          :: nb_thread


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecRvecB_TO_RVecG'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!---------------------------------------------------------------------

        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'size RVecG',size(RVecG)
          write(out_unitp,*) 'size RVecB',size(RVecB)
          write(out_unitp,*) 'RvecB(:)',RvecB(:)
          !CALL RecWrite_basis(basis_set)
          flush(out_unitp)
        END IF

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        WHERE (tab_der_loc < 0) tab_der_loc = 0

        IF (BasisTOGrid_omp == 0) THEN
          nb_thread = 1
        ELSE
          nb_thread = BasisTOGrid_maxth
        END IF
        !write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        IF (basis_set%packed_done) THEN ! Packed basis
          CALL RB_TO_RG_basis(RvecB,RVecG,basis_set,tab_der_loc)
        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed_done !!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product
            nnb = nb
            nbb = nb
            nnq = 1
            nb2 = 1
            nq2 = 1

            CALL alloc_NParray(RTempG,[nnq,nq2,nnb],"RTempG",name_sub)

            RTempG(:,:,:) = reshape(RvecB,shape=[nnq,nq2,nnb])

            DO ibasis=basis_set%nb_basis,1,-1

              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnb = get_size(basis_set%Tab_OF_Tabnb2(ibasis))

              CALL alloc_NParray(RTempB, [nnq,nbb],"RTempB",name_sub)

              RTempB(:,:) = reshape(RTempG,shape=[nnq,nbb])

              CALL dealloc_NParray(RTempG,"RTempG",name_sub)
              CALL alloc_NParray(RTempG,[nnq,nq2,nnb],"RTempG",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2
                  DO iq=1,nnq
                    CALL RB_TO_RG_basis(RTempB(iq,ibb1:ibb2),RTempG(iq,:,ib),basis_set%tab_Pbasis(ibasis)%Pbasis,tab_der_loc)
                  END DO
                  ibb1 = ibb1 + newnb2
                END DO

              ELSE
                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq
                    CALL RecRvecB_TO_RVecG(RTempB(iq,ibb1:ibb2),RTempG(iq,:,ib),&
                                           newnb2,nq2,                          &
                                           basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                           tab_der_loc)
                  END DO
                  ibb1 = ibb1 + newnb2
                END DO

              END IF

              nbb = nnb
              nnq = nnq * nq2
              CALL dealloc_NParray(RTempB,"RTempB",name_sub)
            END DO

            RvecG(:) = reshape(RTempG, shape=[nnq*nnb] )
            CALL dealloc_NParray(RTempG,"RTempG",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               !write(out_unitp,*) 'i_SG,nq_SG',i_SG,nq_SG
!               iq1_SG = iq1_SG + nq_SG
               iq0_SG = 1
               DO ii_SG=1,i_SG-1
                 nq0_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(ii_SG)%Pbasis)
                 iq0_SG = iq0_SG + nq0_SG
               END DO
               iq1_SG = iq0_SG-1 + nq_SG
               !write(out_unitp,*) 'iq0_SG,iq1_SG',iq0_SG,iq1_SG
               CALL RecRvecB_TO_RVecG(RvecB,RVecG(iq0_SG:iq1_SG),       &
                                      nb,nq_SG,                         &
                                     basis_set%tab_PbasisSG(i_SG)%Pbasis,&
                                      tab_der_loc)
!               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d  implementation)

            CALL sub_B_TO_G(RVecB,RVecG,                                &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)

          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RVecB,          &
                             basis_set%tab_basisPrimSG,basis_set%nDindB,&
                             basis_set%para_SGType2)

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL BSmolyakRep_TO_GSmolyakRep(SRep,                         &
                 basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,         &
                 basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL BSmolyakRep_TO3_GSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF

            CALL SmolyakRep2_TO_tabR1bis(RVecG,SRep)
            CALL dealloc_SmolyakRep(SRep)


          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'RvecG(:)',RvecG(:)
         write(out_unitp,*) 'END ',name_sub
       END IF

      END SUBROUTINE RecRvecB_TO_RVecG

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecCVecB_TO_CVecG(CVecB,CVecG,nb,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        USE mod_MPI_aux
        IMPLICIT NONE
        TYPE (basis), intent(in)            :: basis_set
        integer, intent(in)                 :: nq,nb
        complex (kind=Rkind), intent(inout) :: CVecG(:)
        complex (kind=Rkind), intent(in)    :: CVecB(:)
        integer, optional                   :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        complex (kind=Rkind), allocatable   :: CTempG(:,:,:)
        complex (kind=Rkind), allocatable   :: CTempB(:,:)
        TYPE(Type_SmolyakRep)               :: SRep ! smolyak rep for SparseGrid_type=4

        complex (kind=Rkind), allocatable   :: CG(:)
        complex (kind=Rkind), allocatable   :: CB(:)
        real (kind=Rkind), allocatable      :: dnb(:,:)
        complex (kind=Rkind), allocatable   :: dncb(:,:)
        complex (kind=Rkind), allocatable   :: Cvec(:)
        complex (kind=Rkind), allocatable   :: Cvec_temp(:)

        real (kind=Rkind), allocatable      :: RVecG(:)
        real (kind=Rkind), allocatable      :: RVecB(:)

        integer                             :: nbb,ibb1,ibb2
        integer                             :: nnb,nb2,ib,ib2,newnb2
        integer                             :: nnq,nq2,iq,iq2
        integer                             :: ibasis
        integer                             :: i_SG,iq0_SG,iq1_SG,nq_SG

        integer                             :: itabR,iG,nR

        integer                             :: nb_thread
        integer                             :: d1,d2
        integer                             :: Cvec_length(0:MPI_np-1)

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecCVecB_TO_CVecG'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        WHERE (tab_der_loc < 0) tab_der_loc = 0
 
        IF (basis_set%packed_done) THEN
          CALL CB_TO_CG_basis(CVecB,CVecG,basis_set,tab_der_loc)
        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product

            nnb = nb
            nbb = nb
            nnq = 1
            nb2 = 1
            nq2 = 1

            CALL alloc_NParray(CTempG,[nnq,nq2,nnb],"CTempG",name_sub)
            CTempG(:,:,:) = reshape(CVecB,shape=[nnq,nq2,nnb])

            DO ibasis=basis_set%nb_basis,1,-1

              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnb = get_size(basis_set%Tab_OF_Tabnb2(ibasis))

!             write(out_unitp,*) 'B=>G ibasis,Tab_OF_Tabnb2',ibasis,              &
!                                        basis_set%Tab_OF_Tabnb2(ibasis)%vec

              CALL alloc_NParray(CTempB,[nnq,nbb],"CTempB",name_sub)

              CTempB(:,:) = reshape(CTempG,shape=[nnq,nbb])
              CALL dealloc_NParray(CTempG,"CTempG",name_sub)
              CALL alloc_NParray(CTempG,[nnq,nq2,nnb],"CTempG",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN
                  ibb1 = 1
                  ibb2 = 0
                  DO ib=1,nnb
                    newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                    ibb2 = ibb2 + newnb2
                    DO iq=1,nnq
                      CALL CB_TO_CG_basis(CTempB(iq,ibb1:ibb2),CTempG(iq,:,ib),basis_set%tab_Pbasis(ibasis)%Pbasis,tab_der_loc)
                    END DO
                    ibb1 = ibb1 + newnb2
                  END DO
              ELSE
                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2
                  DO iq=1,nnq
                    CALL RecCVecB_TO_CVecG(CTempB(iq,ibb1:ibb2),        &
                                           CTempG(iq,:,ib),newnb2,nq2,  &
                                   basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                                             tab_der_loc)
                  END DO
                  ibb1 = ibb1 + newnb2
                END DO
              END IF

              nbb = nnb
              nnq = nnq * nq2
              CALL dealloc_NParray(CTempB,'CTempB',name_sub)
            END DO

            CvecG(:) = reshape(CTempG, shape=[nnq*nnb] )
            CALL dealloc_NParray(CTempG,"CTempG",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG  = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               iq1_SG = iq1_SG + nq_SG
               CALL RecCVecB_TO_CVecG(CVecB,CVecG(iq0_SG:iq1_SG),       &
                                      nb,nq_SG,                         &
                                   basis_set%tab_PbasisSG(i_SG)%Pbasis, &
                                      tab_der_loc)
               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d  implementation)
            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)
            CALL alloc_NParray(RVecG,shape(CVecG),'RVecG',name_sub)

            RVecB(:) = real(CVecB,kind=Rkind)
            CALL sub_B_TO_G(RVecB,RVecG,                                &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecG(:) = cmplx(RVecG,kind=Rkind)

            RVecB(:) = aimag(CVecB)
            CALL sub_B_TO_G(RVecB,RVecG,                                &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecG(:) = CVecG + EYE*cmplx(RVecG,kind=Rkind)

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)
            CALL dealloc_NParray(RVecG,'RVecG',name_sub)
            !write(out_unitp,*) ' ERROR in ',name_sub
            !STOP 'SparseGrid_type=2'

          CASE (4) ! Sparse basis (Smolyak 4th implementation)
            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)

            RVecB(:) = real(CVecB,kind=Rkind)
            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RVecB,                 &
              basis_set%tab_basisPrimSG,basis_set%nDindB,basis_set%para_SGType2)
            !write(6,*) 'real(CVecB)' ; CALL Write_SmolyakRep(SRep)

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL BSmolyakRep_TO_GSmolyakRep(SRep,                         &
                 basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,         &
                 basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL BSmolyakRep_TO3_GSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF
            !write(6,*) 'real(CVecG)' ; CALL Write_SmolyakRep(SRep)

            IF(openmpi) THEN
              CALL CVecB_TO_CVecG_R_MPI(SRep,CVecG)
            ELSE
              itabR = 0
              DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
                nR = size(SRep%SmolyakRep(iG)%V)
                CVecG(itabR+1:itabR+nR) = cmplx(SRep%SmolyakRep(iG)%V,kind=Rkind)
                itabR = itabR + nR
              END DO
            ENDIF
            CALL dealloc_SmolyakRep(SRep)
            !write(6,*) 'real(CVecG)?',CVecG

            RVecB(:) = aimag(CVecB)
            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RVecB,              &
              basis_set%tab_basisPrimSG,basis_set%nDindB,basis_set%para_SGType2)
            !write(6,*) 'im(CVecB)' ; CALL Write_SmolyakRep(SRep)

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL BSmolyakRep_TO_GSmolyakRep(SRep,                         &
                 basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,         &
                 basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL BSmolyakRep_TO3_GSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF
            !write(6,*) 'im(CVecG)' ; CALL Write_SmolyakRep(SRep)

            IF(openmpi) THEN
              CALL CVecB_TO_CVecG_C_MPI(SRep,CVecG)
            ELSE
              itabR = 0
              DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
                nR = size(SRep%SmolyakRep(iG)%V)
                CVecG(itabR+1:itabR+nR) = CVecG(itabR+1:itabR+nR) +       &
                               EYE*cmplx(SRep%SmolyakRep(iG)%V,kind=Rkind)
                itabR = itabR + nR
              END DO
            ENDIF
            CALL dealloc_SmolyakRep(SRep)
            !write(6,*) 'CVecG?',CVecG

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

      END SUBROUTINE RecCVecB_TO_CVecG

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO

      RECURSIVE SUBROUTINE DerivOp_TO_RVecG(RVecG,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE
        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq
        real (kind=Rkind), intent(inout) :: RVecG(:)
        integer, optional                :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        real (kind=Rkind), allocatable       :: RG1(:,:,:)
        real (kind=Rkind), allocatable       :: RG2(:,:,:)
        real (kind=Rkind), pointer           :: B3GG(:,:)

        integer                          :: nnq1,nnq3,nq2,iq1,iq2,iq3
        integer                          :: iqi,iqe
        logical                          :: skip
        integer                          :: ibasis
        integer                          :: i_SG,iq0_SG,iq1_SG,nq_SG,ii_SG,nq0_SG
        integer                          :: nb_thread
        TYPE(Type_SmolyakRep)            :: SRep ! smolyak rep for SparseGrid_type=4
        integer                          :: itabR,iG,nR

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='DerivOp_TO_RVecG'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!---------------------------------------------------------------------

        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'RvecG(:)',RvecG(:)
          write(out_unitp,*) 'packed_done',basis_set%packed_done
          write(out_unitp,*) 'nb_basis',basis_set%nb_basis
          write(out_unitp,*) 'SparseGrid_type',basis_set%SparseGrid_type

          !CALL RecWrite_basis(basis_set)
          flush(out_unitp)
        END IF

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        IF (debug) write(out_unitp,*) 'tab_der_loc',tab_der_loc
        WHERE (tab_der_loc < 0) tab_der_loc = 0

        IF (BasisTOGrid_omp == 0) THEN
          nb_thread = 1
        ELSE
          nb_thread = BasisTOGrid_maxth
        END IF
        !write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        IF (basis_set%packed_done) THEN ! Packed basis

          dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

          flush(out_unitp)
          IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN ! dnba_ind(:)=0 => no derivative
            CONTINUE ! notthing to do
          ELSE
            nq2 = get_nq_FROM_basis(basis_set)
            CALL Get3_MatdnRGG(basis_set,B3GG,dnba_ind)
            RVecG(:) = matmul(B3GG,RVecG(:))
          END IF

        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed_done !!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product
            nnq3 = nq
            nnq1 = 1
            nq2  = 1

            CALL alloc_NParray(RG1,[nnq1,nq2,nnq3],"RG1",name_sub)
            RG1(:,:,:) = reshape(RvecG,shape=[nnq1,nq2,nnq3])

            DO ibasis=basis_set%nb_basis,1,-1

              nq2  = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnq3 = nnq3 / nq2



              dnba_ind(:) = basis_set%tab_Pbasis(ibasis)%Pbasis%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

              IF (dnba_ind(1) /= 0 .OR. dnba_ind(2) /= 0) THEN
                CALL alloc_NParray(RG2,[nnq1,nq2,nnq3],"RG2",name_sub)
                RG2(:,:,:) = reshape(RG1,shape=[nnq1,nq2,nnq3])

                IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN

                  CALL Get3_MatdnRGG(basis_set%tab_Pbasis(ibasis)%Pbasis,B3GG,dnba_ind)

                 !$OMP parallel do default(none)                        &
                 !$OMP shared(B3GG,RG2,nnq3,nnq1)                        &
                 !$OMP private(iq1,iq3)                                 &
                 !$OMP num_threads(nb_thread)
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                     RG2(iq1,:,iq3) = matmul(B3GG,RG2(iq1,:,iq3))
                  END DO
                  END DO
                 !$OMP end parallel do

                ELSE
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                    CALL DerivOp_TO_RVecG(RG2(iq1,:,iq3),nq2,           &
                                   basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                                            tab_der_loc)
                  END DO
                  END DO
                END IF

                CALL dealloc_NParray(RG1,"RG1",name_sub)
                CALL alloc_NParray(RG1,[nnq1,nq2,nnq3],"RG1",name_sub)
                RG1 = RG2
                CALL dealloc_NParray(RG2,"RG2",name_sub)
              END IF

              nnq1 = nnq1 * nq2

            END DO

            RvecG(:) = reshape(RG1, shape=[nq] )
            CALL dealloc_NParray(RG1,"RG1",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' it does not work with SparseGrid_type=1'
            STOP 'SparseGrid_type=1'

            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               !write(out_unitp,*) 'i_SG,nq_SG',i_SG,nq_SG
!               iq1_SG = iq1_SG + nq_SG
               iq0_SG = 1
               DO ii_SG=1,i_SG-1
                 nq0_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(ii_SG)%Pbasis)
                 iq0_SG = iq0_SG + nq0_SG
               END DO
               iq1_SG = iq0_SG-1 + nq_SG
               !write(out_unitp,*) 'iq0_SG,iq1_SG',iq0_SG,iq1_SG
               CALL DerivOp_TO_RVecG(RVecG(iq0_SG:iq1_SG),nq_SG,        &
                                    basis_set%tab_PbasisSG(i_SG)%Pbasis,&
                                                            tab_der_loc)
!               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d implementation)
            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
              CONTINUE ! nothing to do RvecG is unchanged
            ELSE
              CALL DerivOp_TO_RVecG_SGType2(RVecG,nq,                   &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis,        &
                                     tab_der=tab_der_loc)
            END IF

          CASE (4) ! Sparse basis (Smolyak 4th implementation)
            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
              CONTINUE ! nothing to do RvecG is unchanged
            ELSE
              !!  RVecG TO SRep
              CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,         &
                                   nb0=basis_set%para_SGType2%nb0)
              CALL tabR2bis_TO_SmolyakRep1(SRep,RVecG)

              IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
                CALL DerivOp_TO3_GSmolyakRep(SRep,basis_set%para_SGType2,&
                           basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              ELSE
                CALL DerivOp_TO3_GSmolyakRep(SRep,basis_set%para_SGType2,&
                      basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              END IF

              CALL SmolyakRep2_TO_tabR1bis(RVecG,SRep)

              CALL dealloc_SmolyakRep(SRep)

            END IF

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'RvecG(:)',RvecG(:)
         write(out_unitp,*) 'END ',name_sub
       END IF

      END SUBROUTINE DerivOp_TO_RVecG


      RECURSIVE SUBROUTINE DerivOp_TO_CVecG(CVecG,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE
        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq
        complex (kind=Rkind), intent(inout) :: CVecG(:)
        integer, optional                :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        complex (kind=Rkind), allocatable       :: CG1(:,:,:)
        complex (kind=Rkind), allocatable       :: CG2(:,:,:)
        real (kind=Rkind), pointer              :: B3GG(:,:)

        real (kind=Rkind), allocatable :: RVecG(:)
        !TYPE(Type_SmolyakRep)               :: SRep ! smolyak rep for SparseGrid_type=4
        TYPE(Type_SmolyakRepC)               :: SRep ! smolyak rep for SparseGrid_type=4


        integer                          :: nnq1,nnq3,nq2,iq1,iq2,iq3
        logical                          :: skip

        integer                          :: ibasis
        integer                          :: i,i_SG,iq0_SG,iq1_SG,nq_SG,ii_SG,nq0_SG
        integer                          :: nb_thread
        real (kind=Rkind) :: a

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='DerivOp_TO_CVecG'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!---------------------------------------------------------------------

        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'CvecG(:)',CvecG(:)
          !CALL RecWrite_basis(basis_set)
          flush(out_unitp)
        END IF

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        WHERE (tab_der_loc < 0) tab_der_loc = 0

        IF (BasisTOGrid_omp == 0) THEN
          nb_thread = 1
        ELSE
          nb_thread = BasisTOGrid_maxth
        END IF
        !write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        IF (basis_set%packed_done) THEN ! Packed basis

          dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

          IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN ! dnba_ind(:)=0 => no derivative
            CONTINUE ! notthing to do
          ELSE
            CALL Get3_MatdnRGG(basis_set,B3GG,dnba_ind)
            CVecG(:) = matmul(B3GG,CVecG)
          END IF

        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed_done !!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product

            nnq3 = nq
            nnq1 = 1
            nq2  = 1

            CALL alloc_NParray(CG1,[nnq1,nq2,nnq3],"CG1",name_sub)
            CG1(:,:,:) = reshape(CvecG,shape=[nnq1,nq2,nnq3])

            DO ibasis=basis_set%nb_basis,1,-1

              nq2  = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnq3 = nnq3 / nq2


              dnba_ind(:) = basis_set%tab_Pbasis(ibasis)%Pbasis%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

              IF (dnba_ind(1) /= 0 .OR. dnba_ind(2) /= 0) THEN
                CALL alloc_NParray(CG2,[nnq1,nq2,nnq3],"CG2",name_sub)
                CG2(:,:,:) = reshape(CG1,shape=[nnq1,nq2,nnq3])

                IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN



                  CALL Get3_MatdnRGG(basis_set%tab_Pbasis(ibasis)%Pbasis,B3GG,dnba_ind)

                 !$OMP parallel do default(none)                        &
                 !$OMP shared(B3GG,CG2,nnq3,nnq1)                        &
                 !$OMP private(iq1,iq3)                                 &
                 !$OMP num_threads(nb_thread)
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                     CG2(iq1,:,iq3) = matmul(B3GG,CG2(iq1,:,iq3))
                  END DO
                  END DO
                 !$OMP end parallel do

                ELSE
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                    CALL DerivOp_TO_CVecG(CG2(iq1,:,iq3),nq2,           &
                                   basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                                            tab_der_loc)
                  END DO
                  END DO
                END IF

                CALL dealloc_NParray(CG1,"CG1",name_sub)
                CALL alloc_NParray(CG1,[nnq1,nq2,nnq3],"CG1",name_sub)
                CG1 = CG2
                CALL dealloc_NParray(CG2,"CG2",name_sub)
              END IF

              nnq1 = nnq1 * nq2

            END DO

            CvecG(:) = reshape(CG1, shape=[nq] )
            CALL dealloc_NParray(CG1,"CG1",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' it does not work with SparseGrid_type=1'
            STOP 'SparseGrid_type=1'

            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               !write(out_unitp,*) 'i_SG,nq_SG',i_SG,nq_SG
!               iq1_SG = iq1_SG + nq_SG
               iq0_SG = 1
               DO ii_SG=1,i_SG-1
                 nq0_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(ii_SG)%Pbasis)
                 iq0_SG = iq0_SG + nq0_SG
               END DO
               iq1_SG = iq0_SG-1 + nq_SG
               !write(out_unitp,*) 'iq0_SG,iq1_SG',iq0_SG,iq1_SG
               CALL DerivOp_TO_CVecG(CVecG(iq0_SG:iq1_SG),nq_SG,        &
                                    basis_set%tab_PbasisSG(i_SG)%Pbasis,&
                                                            tab_der_loc)
!               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d and 4th implementation)

            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
              CONTINUE ! nothing to do CvecG is unchanged
            ELSE

              CALL alloc_NParray(RVecG,shape(CVecG),'RVecG',name_sub)

              RVecG(:) = real(CVecG,kind=Rkind)
              CALL DerivOp_TO_RVecG_SGType2(RVecG,nq,                   &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis,        &
                                     tab_der=tab_der_loc)

              DO i=1,size(CVecG)
                a = RVecG(i)
                RVecG(i) = aimag(CVecG(i))
                CVecG(i) = cmplx(a,kind=Rkind)
              END DO
              CALL DerivOp_TO_RVecG_SGType2(RVecG,nq,                   &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis,        &
                                     tab_der=tab_der_loc)
              CVecG(:) = CVecG(:) + EYE*cmplx(RVecG(:),kind=Rkind)

              CALL dealloc_NParray(RVecG,'RVecG',name_sub)

            END IF

          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (all(dnba_ind == 0)) THEN
              CONTINUE ! nothing to do RvecG is unchanged
            ELSE
              !!  RVecG TO SRep
              CALL alloc2_SmolyakRepC(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                     basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                     nb0=basis_set%para_SGType2%nb0)
              CALL tabC2bis_TO_SmolyakRepC1(SRep,CVecG)

              IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
                CALL DerivOp_TO3_GSmolyakRepC(SRep,basis_set%para_SGType2,&
                      basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              ELSE
                CALL DerivOp_TO3_GSmolyakRepC(SRep,basis_set%para_SGType2,&
                      basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              END IF

              CALL SmolyakRepC2_TO_tabC1bis(CVecG,SRep)

              CALL dealloc_SmolyakRepC(SRep)

            END IF

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'CvecG(:)',CvecG(:)
         write(out_unitp,*) 'END ',name_sub
       END IF

      END SUBROUTINE DerivOp_TO_CVecG

SUBROUTINE RVecBC_TO_RvecB(RVecBC,RvecB,nbc,nb,basis_set)
  IMPLICIT NONE


  TYPE (basis),       intent(in)    :: basis_set
  integer,            intent(in)    :: nbc,nb
  real (kind=Rkind),  intent(in)    :: RVecBC(:)
  real (kind=Rkind),  intent(inout) :: RvecB(:)

  real (kind=Rkind), allocatable   :: RB(:,:,:)
  real (kind=Rkind), allocatable   :: RBC(:,:,:)
  integer :: i,i1,i2,i3,N1,NC3,n2,nc2
  integer, allocatable :: tab_nb(:),tab_nbc(:)

!----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='RVecBC_TO_RvecB'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'size RVecBC(:)',size(RVecBC)
  END IF

  IF (basis_set%packed_done) THEN

    RvecB(:) = matmul(basis_set%Rvec(:,1:nbc),RVecBC)
    !RvecB(:) = ZERO
    !DO ibc=1,nbc
    !  RvecB(:) = RvecB(:) + basis_set%Rvec(:,ibc)*RVecBC(ibc)
    !END DO

  ELSE ! basis_set%nb_basis MUST BE > 0
    IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'
    IF (basis_set%SparseGrid_type /= 0) STOP 'ERROR SparseGrid_type must be 0'

    ! tab_nb and tab_nbc calculations
    CALL alloc_NParray(tab_nb, [basis_set%nb_basis], 'tab_nb',  name_sub)
    CALL alloc_NParray(tab_nbc,[basis_set%nb_basis], 'tab_nbc', name_sub)
    DO i=1,basis_set%nb_basis
      tab_nb(i)  = get_nb_FROM_basis(basis_set%tab_Pbasis(i)%Pbasis)
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        tab_nbc(i) = basis_set%tab_Pbasis(i)%Pbasis%nbc
      ELSE
        tab_nbc(i) = tab_nb(i)
      END IF
    END DO

    NC3 = product(tab_nbc)
    N1  = 1
    n2  = 1
    nc2 = 1

    CALL alloc_NParray(RB,[N1,n2,NC3],'RC',name_sub)
    RB(:,:,:) = reshape(RVecBC,shape=[N1,n2,NC3])

    DO i=basis_set%nb_basis,1,-1
      n2  = tab_nb(i)
      nc2 = tab_nbc(i)
      NC3 = NC3 / nc2
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        CALL alloc_NParray(RBC,[N1,nc2,NC3],'RBC',name_sub)
        RBC(:,:,:) = reshape(RB,shape=[N1,nc2,NC3])

        CALL dealloc_NParray(RB,'RB',name_sub)
        CALL alloc_NParray(RB,[N1,n2,NC3],'RC',name_sub)

        DO i3=1,NC3
        DO i1=1,N1
          RB(i1,:,i3) = matmul(basis_set%tab_Pbasis(i)%Pbasis%Rvec(:,1:nc2),RBC(i1,:,i3))
        END DO
        END DO

        CALL dealloc_NParray(RBC,'RBC',name_sub)
      END IF

      N1 = N1 * n2

    END DO

    RvecB(:) = reshape(RB, shape=[N1])
    CALL dealloc_NParray(RB,"RB",name_sub)

    CALL dealloc_NParray(tab_nb, 'tab_nb',  name_sub)
    CALL dealloc_NParray(tab_nbc,'tab_nbc', name_sub)

  END IF

  IF (debug) THEN
    write(out_unitp,*) 'RvecB(:)',RvecB(:)
    write(out_unitp,*) 'END ',name_sub
  END IF
END SUBROUTINE RVecBC_TO_RvecB
SUBROUTINE RVecB_TO_RvecBC(RvecB,RVecBC,nb,nbc,basis_set)
  IMPLICIT NONE


  TYPE (basis),       intent(in)    :: basis_set
  integer,            intent(in)    :: nbc,nb
  real (kind=Rkind),  intent(inout) :: RVecBC(:)
  real (kind=Rkind),  intent(in)    :: RvecB(:)

  real (kind=Rkind), allocatable   :: RB(:,:,:)
  real (kind=Rkind), allocatable   :: RBC(:,:,:)
  integer :: i,i1,i2,i3,NC1,N3,n2,nc2
  integer, allocatable :: tab_nb(:),tab_nbc(:)

!----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='RVecB_TO_RvecBC'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'size RVecB(:)',size(RVecB)
  END IF

  IF (basis_set%packed_done) THEN

    RvecBC(:) = matmul(transpose(basis_set%Rvec(:,1:nbc)),RVecB)
    !RvecB(:) = ZERO
    !DO ibc=1,nbc
    !  RvecB(:) = RvecB(:) + basis_set%Rvec(:,ibc)*RVecBC(ibc)
    !END DO

  ELSE ! basis_set%nb_basis MUST BE > 0
    IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'
    IF (basis_set%SparseGrid_type /= 0) STOP 'ERROR SparseGrid_type must be 0'

    ! tab_nb and tab_nbc calculations
    CALL alloc_NParray(tab_nb, [basis_set%nb_basis], 'tab_nb',  name_sub)
    CALL alloc_NParray(tab_nbc,[basis_set%nb_basis], 'tab_nbc', name_sub)
    DO i=1,basis_set%nb_basis
      tab_nb(i)  = get_nb_FROM_basis(basis_set%tab_Pbasis(i)%Pbasis)
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        tab_nbc(i) = basis_set%tab_Pbasis(i)%Pbasis%nbc
      ELSE
        tab_nbc(i) = tab_nb(i)
      END IF
    END DO

    N3  = product(tab_nb)
    NC1 = 1
    n2  = 1
    nc2 = 1

    CALL alloc_NParray(RBC,[NC1,n2,N3],'RC',name_sub)
    RBC(:,:,:) = reshape(RVecB,shape=[NC1,n2,N3])

    DO i=basis_set%nb_basis,1,-1

      n2  = tab_nb(i)
      nc2 = tab_nbc(i)
      N3  = N3 / n2
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        CALL alloc_NParray(RB,[NC1,n2,N3],'RB',name_sub)
        RB(:,:,:) = reshape(RBC,shape=[NC1,n2,N3])

        CALL dealloc_NParray(RBC,'RBC',name_sub)
        CALL alloc_NParray(RBC,[NC1,nc2,N3],'RBC',name_sub)

        DO i3=1,N3
        DO i1=1,NC1
          RBC(i1,:,i3) = matmul(transpose(basis_set%tab_Pbasis(i)%Pbasis%Rvec(:,1:nc2)),RB(i1,:,i3))
        END DO
        END DO

        CALL dealloc_NParray(RB,'RB',name_sub)
      END IF

      NC1 = NC1 * nc2

    END DO

    RvecBC(:) = reshape(RBC, shape=[NC1])
    CALL dealloc_NParray(RBC,"RBC",name_sub)

    CALL dealloc_NParray(tab_nb, 'tab_nb',  name_sub)
    CALL dealloc_NParray(tab_nbc,'tab_nbc', name_sub)

  END IF

  IF (debug) THEN
    write(out_unitp,*) 'RvecBC(:)',RvecBC(:)
    write(out_unitp,*) 'END ',name_sub
  END IF
END SUBROUTINE RVecB_TO_RvecBC


SUBROUTINE CVecBC_TO_CvecB(CVecBC,CvecB,nbc,nb,basis_set)
  IMPLICIT NONE


  TYPE (basis),       intent(in)        :: basis_set
  integer,            intent(in)        :: nbc,nb
  complex (kind=Rkind),  intent(in)     :: CvecBC(:)
  complex (kind=Rkind),  intent(inout)  :: CvecB(:)

  complex (kind=Rkind), allocatable   :: CB(:,:,:)
  complex (kind=Rkind), allocatable   :: CBC(:,:,:)
  integer :: i,i1,i2,i3,N1,NC3,n2,nc2
  integer, allocatable :: tab_nb(:),tab_nbc(:)

!----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='CVecBC_TO_CvecB'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'size CvecBC(:)',size(CvecBC)
  END IF

  IF (basis_set%packed_done) THEN

    CvecB(:) = matmul(basis_set%Rvec(:,1:nbc),CvecBC)
    !CvecB(:) = ZERO
    !DO ibc=1,nbc
    !  CvecB(:) = CvecB(:) + basis_set%Rvec(:,ibc)*CvecBC(ibc)
    !END DO

  ELSE ! basis_set%nb_basis MUST BE > 0
    IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'
    IF (basis_set%SparseGrid_type /= 0) STOP 'ERROR SparseGrid_type must be 0'

    ! tab_nb and tab_nbc calculations
    CALL alloc_NParray(tab_nb, [basis_set%nb_basis], 'tab_nb',  name_sub)
    CALL alloc_NParray(tab_nbc,[basis_set%nb_basis], 'tab_nbc', name_sub)
    DO i=1,basis_set%nb_basis
      tab_nb(i)  = get_nb_FROM_basis(basis_set%tab_Pbasis(i)%Pbasis)
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        tab_nbc(i) = basis_set%tab_Pbasis(i)%Pbasis%nbc
      ELSE
        tab_nbc(i) = tab_nb(i)
      END IF
    END DO

    NC3 = product(tab_nbc)
    N1  = 1
    n2  = 1
    nc2 = 1

    CALL alloc_NParray(CB,[N1,n2,NC3],'RC',name_sub)
    CB(:,:,:) = reshape(CvecBC,shape=[N1,n2,NC3])

    DO i=basis_set%nb_basis,1,-1
      n2  = tab_nb(i)
      nc2 = tab_nbc(i)
      NC3 = NC3 / nc2
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        CALL alloc_NParray(CBC,[N1,nc2,NC3],'CBC',name_sub)
        CBC(:,:,:) = reshape(CB,shape=[N1,nc2,NC3])

        CALL dealloc_NParray(CB,'CB',name_sub)
        CALL alloc_NParray(CB,[N1,n2,NC3],'RC',name_sub)

        DO i3=1,NC3
        DO i1=1,N1
          CB(i1,:,i3) = matmul(basis_set%tab_Pbasis(i)%Pbasis%Rvec(:,1:nc2),CBC(i1,:,i3))
        END DO
        END DO

        CALL dealloc_NParray(CBC,'CBC',name_sub)
      END IF

      N1 = N1 * n2

    END DO

    CvecB(:) = reshape(CB, shape=[N1])
    CALL dealloc_NParray(CB,"CB",name_sub)

    CALL dealloc_NParray(tab_nb, 'tab_nb',  name_sub)
    CALL dealloc_NParray(tab_nbc,'tab_nbc', name_sub)

  END IF

  IF (debug) THEN
    write(out_unitp,*) 'CvecB(:)',CvecB(:)
    write(out_unitp,*) 'END ',name_sub
  END IF
END SUBROUTINE CVecBC_TO_CvecB
SUBROUTINE CVecB_TO_CvecBC(CvecB,CvecBC,nb,nbc,basis_set)
  IMPLICIT NONE


  TYPE (basis),       intent(in)    :: basis_set
  integer,            intent(in)    :: nbc,nb
  complex (kind=Rkind),  intent(inout) :: CvecBC(:)
  complex (kind=Rkind),  intent(in)    :: CvecB(:)

  complex (kind=Rkind), allocatable   :: CB(:,:,:)
  complex (kind=Rkind), allocatable   :: CBC(:,:,:)
  integer :: i,i1,i2,i3,NC1,N3,n2,nc2
  integer, allocatable :: tab_nb(:),tab_nbc(:)

!----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='CVecB_TO_CvecBC'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'size CvecB(:)',size(CvecB)
  END IF

  IF (basis_set%packed_done) THEN

    CvecBC(:) = matmul(transpose(basis_set%Rvec(:,1:nbc)),CvecB)
    !CvecB(:) = ZERO
    !DO ibc=1,nbc
    !  CvecB(:) = CvecB(:) + basis_set%Rvec(:,ibc)*CvecBC(ibc)
    !END DO

  ELSE ! basis_set%nb_basis MUST BE > 0
    IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'
    IF (basis_set%SparseGrid_type /= 0) STOP 'ERROR SparseGrid_type must be 0'

    ! tab_nb and tab_nbc calculations
    CALL alloc_NParray(tab_nb, [basis_set%nb_basis], 'tab_nb',  name_sub)
    CALL alloc_NParray(tab_nbc,[basis_set%nb_basis], 'tab_nbc', name_sub)
    DO i=1,basis_set%nb_basis
      tab_nb(i)  = get_nb_FROM_basis(basis_set%tab_Pbasis(i)%Pbasis)
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        tab_nbc(i) = basis_set%tab_Pbasis(i)%Pbasis%nbc
      ELSE
        tab_nbc(i) = tab_nb(i)
      END IF
    END DO

    N3  = product(tab_nb)
    NC1 = 1
    n2  = 1
    nc2 = 1

    CALL alloc_NParray(CBC,[NC1,n2,N3],'RC',name_sub)
    CBC(:,:,:) = reshape(CvecB,shape=[NC1,n2,N3])

    DO i=basis_set%nb_basis,1,-1

      n2  = tab_nb(i)
      nc2 = tab_nbc(i)
      N3  = N3 / n2
      IF (basis_set%tab_Pbasis(i)%Pbasis%contrac_RVecOnly) THEN
        CALL alloc_NParray(CB,[NC1,n2,N3],'CB',name_sub)
        CB(:,:,:) = reshape(CBC,shape=[NC1,n2,N3])

        CALL dealloc_NParray(CBC,'CBC',name_sub)
        CALL alloc_NParray(CBC,[NC1,nc2,N3],'CBC',name_sub)

        DO i3=1,N3
        DO i1=1,NC1
          CBC(i1,:,i3) = matmul(transpose(basis_set%tab_Pbasis(i)%Pbasis%Rvec(:,1:nc2)),CB(i1,:,i3))
        END DO
        END DO

        CALL dealloc_NParray(CB,'CB',name_sub)
      END IF

      NC1 = NC1 * nc2

    END DO

    CvecBC(:) = reshape(CBC, shape=[NC1])
    CALL dealloc_NParray(CBC,"CBC",name_sub)

    CALL dealloc_NParray(tab_nb, 'tab_nb',  name_sub)
    CALL dealloc_NParray(tab_nbc,'tab_nbc', name_sub)

  END IF

  IF (debug) THEN
    write(out_unitp,*) 'CvecBC(:)',CvecBC(:)
    write(out_unitp,*) 'END ',name_sub
  END IF
END SUBROUTINE CVecB_TO_CvecBC

END MODULE mod_basis_BtoG_GtoB
