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
     MODULE mod_Auto_Basis
     IMPLICIT NONE

     PRIVATE
     PUBLIC :: Auto_basis, RecAuto_basis, RecSet_EneH0, AutoParam_basis
     PUBLIC :: basis_TO_AllBasis, All_param_TO_para_H, sub_MatOp_HADA

     CONTAINS

!=======================================================================================
!
!        Automatic calculation of a contracted basis set (in nD)
!
!=======================================================================================
      SUBROUTINE Auto_basis(para_Tnum,mole,para_AllBasis,para_ReadOp)
      use mod_Coord_KEO
      use mod_PrimOp
      USE mod_nDindex
      use mod_basis, only: param_allbasis, sgtype, get_nq_from_basis,  &
                           get_nqa_from_basis,                         &
                           get_nb_bi_from_allbasis, recwrite_basis,    &
                           basis, sub_dngb_to_dnbb, clean_basis,       &
                           sub_dngb_to_dngg, construct_primitive_basis,&
                           sub_contraction_basis, check_ortho_basis,   &
                           dealloc_allbasis, alloc_allbasis,           &
                           basis2tobasis1,dealloc_basis_ext2n
      USE mod_Op
      USE mod_MPI
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), target :: mole,mole_loc
      TYPE (Tnum),      target :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis), target :: para_AllBasis

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in) :: para_ReadOp


!----- local variables
      integer :: ib
      integer, allocatable :: tab_nbc(:)
      integer :: nqa,nb_bi
      logical :: With_BGG,contrac_RVecOnly
      real(kind=Rkind), allocatable :: TDParam(:)
      integer :: nb_TDParam

!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Auto_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------
      With_BGG = .TRUE. ! default
      IF (SGtype == 2) With_BGG = .TRUE.
      IF (SGtype == 1) With_BGG = .FALSE.
      IF (debug) write(out_unit,*) 'SGtype,With_BGG',SGtype,With_BGG
      CALL Set_dnGGRep(para_AllBasis%BasisnD,With_BGG)

      mole_loc = mole

      !> primary memory here
      CALL RecAuto_basis(para_Tnum,mole_loc,para_AllBasis%BasisnD,para_ReadOp)

      contrac_RVecOnly = .FALSE.
      IF (para_AllBasis%BasisnD%SparseGrid_type == 0) THEN ! not sparse basis
        DO ib=1,para_AllBasis%BasisnD%nb_basis
          contrac_RVecOnly = contrac_RVecOnly .OR. para_AllBasis%BasisnD%tab_Pbasis(ib)%Pbasis%contrac_RVecOnly
        END DO
      END IF
      IF (contrac_RVecOnly) THEN
        write(out_unit,*) 'Special traitment when contrac_RVecOnly=t'
        IF (allocated(para_AllBasis%BasisnD%nDindB_contracted)) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'nDindB_contracted is allocated'
          STOP 'ERROR in Auto_basis: nDindB_contracted is allocated'
        END IF
        CALL alloc_NParray(para_AllBasis%BasisnD%nDindB_contracted,'nDindB_contracted',name_sub)

        CALL alloc_NParray(tab_nbc,[para_AllBasis%BasisnD%nb_basis],'tab_nbc',name_sub)
        DO ib=1,para_AllBasis%BasisnD%nb_basis
          IF (para_AllBasis%BasisnD%tab_Pbasis(ib)%Pbasis%contrac_RVecOnly) THEN
            tab_nbc(ib) = para_AllBasis%BasisnD%tab_Pbasis(ib)%Pbasis%nbc
          ELSE
            tab_nbc(ib) = get_nb_FROM_basis(para_AllBasis%BasisnD%tab_Pbasis(ib)%Pbasis)
          END IF
        END DO
        write(out_unit,*) 'tab_nbc(:)',tab_nbc
        IF (product(tab_nbc) <= 0) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'Wrong values in tab_nbc(:)',tab_nbc
          STOP 'ERROR in Auto_basis: Wrong values in tab_nbc(:)'
        END IF
        CALL init_nDindexPrim(para_AllBasis%BasisnD%nDindB_contracted,          &
                              para_AllBasis%BasisnD%nb_basis,tab_nbc,           &
                              type_OF_nDindex=1) ! standard direct-product
        CALL dealloc_NParray(tab_nbc,'tab_nbc',name_sub)
        !CALL Write_nDindex(para_AllBasis%BasisnD%nDindB_contracted)
      END IF

      CALL dealloc_CoordType(mole_loc)
      !CALL Write_SymAbelian(para_AllBasis%BasisnD%P_SymAbelian)

      nqa = get_nqa_FROM_basis(para_AllBasis%BasisnD)

      IF (print_level > -1 .AND. MPI_id==0) THEN
        write(out_unit,*) '==================================================='
        write(out_unit,*) '=== NEW BASIS ====================================='
        write(out_unit,*) 'packed',para_AllBasis%BasisnD%packed
        write(out_unit,*) 'nb_ba,nb_qa',para_AllBasis%BasisnD%nb,nqa
        write(out_unit,*) 'nb_bi',get_nb_bi_FROM_AllBasis(para_AllBasis)
        write(out_unit,*) 'nb_elec',para_ReadOp%nb_elec
        write(out_unit,*) 'nb_inact2n',mole%nb_inact2n
        flush(out_unit)
      END IF
      IF (nqa < 1) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'The number of grid points is < 1'
        write(out_unit,*) ' nb_qa',nqa
        STOP 'ERROR in Auto_basis: The number of grid points is < 1'
      END IF
      IF(keep_MPI) THEN
        IF (allocated(para_AllBasis%BasisnD%nDindB%nDsize))             &
          write(out_unit,*) 'nDindB%nDsize',para_AllBasis%BasisnD%nDindB%nDsize
        !CALL Write_nDindex(para_AllBasis%BasisnD%nDindB,"BasisnD%nDinB ")

         IF (allocated(para_AllBasis%BasisnD%nDindG%nDsize))             &
           write(out_unit,*) 'nDindG%nDsize',para_AllBasis%BasisnD%nDindG%nDsize
         !CALL Write_nDindex(para_AllBasis%BasisnD%nDindG,"BasisnD%nDinG ")
      END IF

      nb_TDParam = get_nb_TDParam_FROM_basis(para_AllBasis%BasisnD)
      IF (nb_TDParam > 0) then
        CALL alloc_NParray(TDParam,[nb_TDParam],'nb_TDParam',name_sub)
        CALL Get_TDParam_FROM_basis(para_AllBasis%BasisnD,TDParam)
        write(out_unit,*) 'TDParam ',TDParam(:)
        CALL dealloc_NParray(TDParam,'nb_TDParam',name_sub)
      ELSE
        write(out_unit,*) 'no TDParam'
      END IF

      IF (debug) THEN
        write(out_unit,*) '==== NEW BASIS ======================================='
        CALL RecWrite_basis(para_AllBasis%BasisnD,write_all=.TRUE.)
        write(out_unit,*) '==== END NEW BASIS ==================================='
        write(out_unit,*) 'END ',name_sub
      END IF
      IF (print_level > 1 ) write(out_unit,*) 'nrho in ',name_sub,para_AllBasis%BasisnD%nrho(:)
      IF (print_level > -1) write(out_unit,*) '==================================================='
      CALL RecWriteMiniMini_basis(para_AllBasis%BasisnD)

      END SUBROUTINE Auto_basis
!=======================================================================================


!=======================================================================================
! Auto contraction
!=======================================================================================
      RECURSIVE SUBROUTINE RecAuto_basis(para_Tnum,mole,BasisnD,para_ReadOp)
      USE EVR_system_m
      USE mod_PrimOp
      USE mod_basis
      USE BasisMakeGrid
      USE mod_Op
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), target :: mole
      TYPE (Tnum),      target :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis) :: basisnD

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in) :: para_ReadOp


!----- local variables
      TYPE (basis) :: basisnD_contrac
      logical :: Rec_call
      integer :: i,j,iact,isym,ibasis,jbasis,nb_ba,nb0,nbc0
      integer :: nb_elec_save,n_h_save
      real(kind=Rkind) :: ene0,Effi

      integer :: ib,ibi,val,ii,nq,i_SG,nb_SG
      logical :: Print_basis

       integer, save :: rec = 0
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'RecAuto_basis'
!---------------------------------------------------------------------
      rec = rec + 1
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub,rec
      END IF
!---------------------------------------------------------------------
       Print_basis = BasisnD%print_info_OF_basisDP .AND. print_level > -1 .OR. debug


      IF (Print_basis .AND. MPI_id==0) THEN
        write(out_unit,*) '==============================================='
        write(out_unit,*) '= Rec Auto Basis, layer: ',rec
        write(out_unit,*) '==============================================='

        write(out_unit,*) 'packed_done:                   ',BasisnD%packed_done
        write(out_unit,*) 'SparseGrid_type:               ',BasisnD%SparseGrid_type
        write(out_unit,*) 'nb_basis:                      ',BasisnD%nb_basis
        flush(out_unit)
      END IF

      IF (BasisnD%nb_basis > 0 .AND. .NOT. BasisnD%packed_done .AND. .NOT. BasisnD%BuildBasis_done) THEN

        IF (BasisnD%SparseGrid_type == 1) THEN
          ! For Sparse Grid
          CALL RecSparseGrid_ForDP_type1(BasisnD,para_Tnum,mole,para_ReadOp)

          !- d1b => d1BasisRep and  d2b => d2BasisRep ------------
          CALL sub_dnGB_TO_dnBB(BasisnD)

          CALL clean_basis(BasisnD)


          IF (Print_basis) write(out_unit,*) 'Sparse Grid type1 done. Layer: ',rec

        ELSE IF (BasisnD%SparseGrid_type == 2) THEN
          CALL RecSparseGrid_ForDP_type2(BasisnD,para_Tnum,mole,para_ReadOp)

          !- d1b => dnBGG%d1 and  d2b => dnBGG%d2 ------------
          CALL sub_dnGB_TO_dnGG(BasisnD)

          !CALL clean_basis(BasisnD)

          IF (BasisnD%print_info_OF_basisDP .AND. print_level > -1)     &
                write(out_unit,*) 'Sparse Grid type2 done. Layer: ',rec

        ELSE IF (BasisnD%SparseGrid_type == 4) THEN

          CALL RecSparseGrid_ForDP_type4(BasisnD,para_Tnum,mole,para_ReadOp)

          !- d1b => dnBGG%d1 and  d2b => dnBGG%d2 ------------
          CALL sub_dnGB_TO_dnGG(BasisnD)

          !CALL clean_basis(BasisnD)

          IF (Print_basis .AND. MPI_id==0) write(out_unit,*) 'Sparse Grid type4 done. Layer: ',rec

        ELSE
          ! For direct product basis (recursive)

          IF (Print_basis) write(out_unit,*) 'direct_product%...%nb:         ',  &
                  (BasisnD%tab_Pbasis(i)%Pbasis%nb,i=1,BasisnD%nb_basis)


          IF (.NOT. BasisnD%tab_basis_done .OR. .NOT. BasisnD%BuildBasis_done) THEN
            DO ibasis=1,BasisnD%nb_basis


              IF (Print_basis) write(out_unit,*)                       &
                                'direct_product%tab_Pbasis%(i): ',ibasis

              IF (BasisnD%Type_OF_nDindB == 0)                          &
                       BasisnD%tab_Pbasis(ibasis)%Pbasis%packed = .TRUE.

              IF (BasisnD%tab_Pbasis(ibasis)%Pbasis%L_SparseGrid == -1) &
                      BasisnD%tab_Pbasis(ibasis)%Pbasis%L_SparseGrid =  &
                                                   BasisnD%L_SparseGrid

              IF (BasisnD%tab_Pbasis(ibasis)%Pbasis%L_SparseBasis == -1) &
                      BasisnD%tab_Pbasis(ibasis)%Pbasis%L_SparseBasis = &
                                                   BasisnD%L_SparseBasis

              IF (BasisnD%tab_Pbasis(ibasis)%Pbasis%Norm_OF_nDindB == huge(ONE)) &
                      BasisnD%tab_Pbasis(ibasis)%Pbasis%Norm_OF_nDindB = &
                                                   BasisnD%Norm_OF_nDindB

              BasisnD%tab_Pbasis(ibasis)%Pbasis%print_info_OF_basisDP = &
                                           BasisnD%print_info_OF_basisDP

              IF (.NOT. BasisnD%check_nq_OF_basis)                      &
              BasisnD%tab_Pbasis(ibasis)%Pbasis%check_nq_OF_basis = .FALSE.
              IF (.NOT. BasisnD%check_basis)                            &
                 BasisnD%tab_Pbasis(ibasis)%Pbasis%check_basis = .FALSE.

              CALL RecAuto_basis(para_Tnum,mole,                        &
                                 BasisnD%tab_Pbasis(ibasis)%Pbasis,     &
                                 para_ReadOp)

              IF (Print_basis) write(out_unit,*)                       &
                    'direct_product%tab_Pbasis(i): ',ibasis,' done. Layer: ',rec
              flush(out_unit)
            END DO
          ELSE
            IF (Print_basis) write(out_unit,*)                         &
                'direct_product%tab_basis_done is already done. Layer: ',rec

          END IF
          BasisnD%tab_basis_done = .TRUE.

          IF (Print_basis) write(out_unit,*) 'direct_product%...%nb',  &
                 (BasisnD%tab_Pbasis(i)%Pbasis%nb,i=1,BasisnD%nb_basis)
          IF (Print_basis) write(out_unit,*) 'direct_product%...%nq',  &
           (get_nq_FROM_basis(BasisnD%tab_Pbasis(i)%Pbasis),i=1,BasisnD%nb_basis)
          flush(out_unit)

          ! direct product construction
          CALL sub_DirProd_basis(BasisnD)
          IF (Print_basis) write(out_unit,*) 'direct_product done. Layer:',rec
        END IF
      ELSE ! primitive basis
        IF (Print_basis) THEN
          write(out_unit,*) 'Active basis:                  ',BasisnD%active
          write(out_unit,*) 'check_nq_OF_basis:             ',BasisnD%check_nq_OF_basis
          flush(out_unit)
        END IF

        CALL construct_primitive_basis(BasisnD)

        IF (BasisnD%auto_basis) THEN
          CALL AutoParam_basis(BasisnD,para_Tnum,mole,para_ReadOp)
          CALL construct_primitive_basis(BasisnD)
        END IF

        IF (Print_basis) write(out_unit,*) 'Primitive basis done. Layer:      ',rec
        flush(out_unit)

      END IF
      flush(out_unit)

      ! For the recursivity ....

      ! To optimize cubature grid (experimental) .....
      IF (BasisnD%make_cubature .OR. BasisnD%Restart_make_cubature) THEN
        IF (Print_basis) write(out_unit,*) 'Cubature basis. Layer:         ',rec
        CALL Make_grid_basis(BasisnD)
        IF (Print_basis) write(out_unit,*) 'Cubature basis done. Layer:    ',rec
        STOP
      END IF

      ! Now the contraction .....
      IF (BasisnD%contrac .AND. .NOT. BasisnD%BuildBasis_done) THEN
        IF (BasisnD%auto_contrac) THEN

          !PODVR
          IF (BasisnD%ndim == 1 .AND. BasisnD%POGridRep) THEN
            nb0  = BasisnD%nb  ! save the value before the contraction
            nbc0 = BasisnD%nbc ! save the value before the contraction

            CALL Autocontract_basis(BasisnD,para_Tnum,mole,para_ReadOp)

            IF (Print_basis) write(out_unit,*) 'Autocontract_basis (POGridRep_poly) done. Layer: ',rec
            flush(out_unit)

            BasisnD%nqc =  BasisnD%nbc
            !CALL POGridRep2_basis(BasisnD,nb0)
            CALL POGridRep_basis(BasisnD,nb0)

            IF (Print_basis) write(out_unit,*) 'POGridRep_basis done. Layer:   ',rec
            flush(out_unit)

            !- d1b => d1BasisRep and  d2b => d2BasisRep ------------
            CALL sub_dnGB_TO_dnBB(BasisnD)

            !- d1b => d1GRep and  d2b => d2GRep ------------
            CALL sub_dnGB_TO_dnGG(BasisnD)

            BasisnD%nbc = nbc0
          END IF

          BasisnD%POGridRep_polyortho = .FALSE.

          IF (BasisnD%contrac_RVecOnly .OR. BasisnD%contrac_analysis) THEN
            CALL basis2TObasis1(BasisnD_contrac,BasisnD)
            CALL Autocontract_basis(BasisnD_contrac,para_Tnum,mole,para_ReadOp)
            !CALL alloc_NParray(BasisnD%Rvec,shape(BasisnD_contrac%Rvec),&
            !                  'BasisnD%Rvec',name_sub)
            !BasisnD%Rvec(:,:) = BasisnD_contrac%Rvec
            CALL move_alloc(FROM=BasisnD_contrac%Rvec,TO=BasisnD%Rvec)
            BasisnD%nbc       = BasisnD_contrac%nbc
            CALL move_alloc(FROM=BasisnD_contrac%EneH0, TO=BasisnD%EneH0)
            CALL dealloc_basis(BasisnD_contrac)
          ELSE
            CALL Autocontract_basis(BasisnD,para_Tnum,mole,para_ReadOp)
          END IF

          IF (Print_basis) write(out_unit,*) 'Autocontract_basis done. Layer:',rec
          flush(out_unit)

        ELSE ! just contraction (not the automatic procedure)

          IF (BasisnD%contrac_RVecOnly .OR. BasisnD%contrac_analysis) THEN
            CALL basis2TObasis1(BasisnD_contrac,BasisnD)
            CALL sub_contraction_basis(BasisnD_contrac,.FALSE.)
            !CALL alloc_NParray(BasisnD%Rvec,shape(BasisnD_contrac%Rvec),        &
            !                  'BasisnD%Rvec',name_sub)
            !BasisnD%Rvec(:,:) = BasisnD_contrac%Rvec
            CALL move_alloc(FROM=BasisnD_contrac%Rvec,TO=BasisnD%Rvec)
            BasisnD%nbc       = BasisnD_contrac%nbc
            CALL dealloc_basis(BasisnD_contrac)
          ELSE
            CALL sub_contraction_basis(BasisnD,.FALSE.)
          END IF

          IF (Print_basis) write(out_unit,*) 'Contract_basis done. Layer:    ',rec

        END IF
        !- d1b => d1BasisRep and  d2b => d2BasisRep ------------
        CALL sub_dnGB_TO_dnBB(BasisnD)

        !---------------------------------------------------------------
        !- check the overlap matrix -----------------------------
        CALL check_ortho_basis(BasisnD)
      END IF
      BasisnD%BuildBasis_done = .TRUE.

      IF (Print_basis .AND. MPI_id==0) write(out_unit,*) 'Auto_basis done. Layer:        ',rec

      IF (debug) THEN
        CALL RecWrite_basis(BasisnD)
      END IF
      IF (Print_basis .AND. MPI_id==0) THEN
        write(out_unit,*) '==============================================='
        write(out_unit,*) '= END Rec Auto Basis for Layer:',rec
        write(out_unit,*) '==================================================='
        flush(out_unit)
      END IF
      rec = rec - 1

      END SUBROUTINE RecAuto_basis
!=======================================================================================

!=======================================================================================
      RECURSIVE SUBROUTINE RecSet_EneH0(para_Tnum,mole,BasisnD,para_ReadOp)

      USE EVR_system_m
      USE mod_nDindex
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),    intent(in)    :: mole
      TYPE (Tnum),         intent(in)    :: para_Tnum
!----- for the basis set ----------------------------------------------
      TYPE (basis),        intent(inout) :: basisnD
!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)    :: para_ReadOp

!----- local variables
      integer             :: i,ib,L,nb
      integer             :: nDval(basisnD%nb_basis)
      TYPE(REAL_WU)       :: RWU_E

      integer, save :: rec = 0
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'RecSet_EneH0'
!---------------------------------------------------------------------
      rec = rec + 1
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub,rec
        IF (allocated(BasisnD%EneH0)) write(out_unit,*) 'BasisnD%EneH0',BasisnD%EneH0
        DO ib=1,BasisnD%nb_basis
          IF (allocated(BasisnD%tab_Pbasis(ib)%Pbasis%EneH0)) THEN
            write(out_unit,*) ib,'tab_Pbasis(i)....%EneH0',BasisnD%tab_Pbasis(ib)%Pbasis%EneH0
          ELSE
            write(out_unit,*) ib,'tab_Pbasis(i)....%EneH0: not allocated'
          END IF
        END DO

      END IF
!---------------------------------------------------------------------

      IF (BasisnD%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unit,*) '==============================================='
        write(out_unit,*) '= RecSet_EneH0, layer: ',rec
        write(out_unit,*) '==============================================='

        write(out_unit,*) 'packed_done:                   ',BasisnD%packed_done
        write(out_unit,*) 'SparseGrid_type:               ',BasisnD%SparseGrid_type
        write(out_unit,*) 'nb_basis:                      ',BasisnD%nb_basis
        flush(out_unit)
      END IF

      IF (BasisnD%nb_basis > 0 .AND. .NOT. BasisnD%packed_done) THEN

        IF (BasisnD%contrac_RVecOnly) THEN
          IF (.NOT. BasisnD%auto_contrac) THEN
            CALL Set_EneH0_OF_ContracBasis(BasisnD,para_Tnum,mole,para_ReadOp)
          END IF
        ELSE

          IF (allocated(BasisnD%EneH0))    THEN
            CALL dealloc_NParray(BasisnD%EneH0,"BasisnD%EneH0",name_sub)
          END IF
          CALL alloc_NParray(BasisnD%EneH0,[BasisnD%nb],"BasisnD%EneH0",name_sub)

          SELECT CASE (BasisnD%SparseGrid_type)
          CASE (0) ! Direct product

            nb = 1
            DO i=1,BasisnD%nb_basis
              CALL RecSet_EneH0(para_Tnum,mole,                           &
                                BasisnD%tab_Pbasis(i)%Pbasis,para_ReadOp)
              nb = nb * size(BasisnD%tab_Pbasis(i)%Pbasis%EneH0)
            END DO

            IF (allocated(BasisnD%nDindB_contracted)) THEN

              IF (nb /= BasisnD%nDindB_contracted%Max_nDI) STOP 'ERROR wrong EneH0 size'
              CALL dealloc_NParray(BasisnD%EneH0,"BasisnD%EneH0",name_sub)
              CALL alloc_NParray(BasisnD%EneH0,[BasisnD%nDindB_contracted%Max_nDI],"BasisnD%EneH0",name_sub)

              DO ib=1,BasisnD%nDindB_contracted%Max_nDI
                CALL calc_nDindex(BasisnD%nDindB_contracted,ib,nDval)
                BasisnD%EneH0(ib) = ZERO
                DO i=1,BasisnD%nb_basis
                  BasisnD%EneH0(ib) = BasisnD%EneH0(ib) +                       &
                                BasisnD%tab_Pbasis(i)%Pbasis%EneH0(nDval(i))
                END DO
              END DO
            ELSE
              DO ib=1,BasisnD%nb
                CALL calc_nDindex(BasisnD%nDindB,ib,nDval)
                BasisnD%EneH0(ib) = ZERO
                DO i=1,BasisnD%nb_basis
                  BasisnD%EneH0(ib) = BasisnD%EneH0(ib) +                       &
                                BasisnD%tab_Pbasis(i)%Pbasis%EneH0(nDval(i))
                END DO
              END DO
            END IF

          CASE (1) ! Sparse basis (Smolyak 1st implementation)

            L = BasisnD%L_SparseBasis

            DO i=1,BasisnD%nb_basis
              CALL RecSet_EneH0(para_Tnum,mole,                           &
                                BasisnD%tab_basisPrimSG(i,L),para_ReadOp)
            END DO

            DO ib=1,BasisnD%nb
              CALL calc_nDindex(BasisnD%nDindB,ib,nDval)

              BasisnD%EneH0(ib) = ZERO
              DO i=1,BasisnD%nb_basis
                BasisnD%EneH0(ib) = BasisnD%EneH0(ib) +                 &
                              BasisnD%tab_basisPrimSG(i,L)%EneH0(nDval(i))
              END DO
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d implementation)

            L = BasisnD%L_SparseBasis
            DO i=1,BasisnD%nb_basis
              CALL RecSet_EneH0(para_Tnum,mole,                           &
                                BasisnD%tab_basisPrimSG(L,i),para_ReadOp)
            END DO

            DO ib=1,BasisnD%nb
              CALL calc_nDindex(BasisnD%nDindB,ib,nDval)

              BasisnD%EneH0(ib) = ZERO
              DO i=1,BasisnD%nb_basis
                BasisnD%EneH0(ib) = BasisnD%EneH0(ib) +                 &
                              BasisnD%tab_basisPrimSG(L,i)%EneH0(nDval(i))
              END DO
            END DO


          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            L = BasisnD%L_SparseBasis
            DO i=1,BasisnD%nb_basis
              CALL RecSet_EneH0(para_Tnum,mole,                           &
                                BasisnD%tab_basisPrimSG(L,i),para_ReadOp)
            END DO

            CALL init_nDval_OF_nDindex(BasisnD%nDindB,nDval)
            DO ib=1,BasisnD%nb
              CALL ADD_ONE_TO_nDindex(BasisnD%nDindB,nDval)
              !CALL calc_nDindex(BasisnD%nDindB,ib,nDval)

              BasisnD%EneH0(ib) = ZERO
              DO i=1,BasisnD%nb_basis
                BasisnD%EneH0(ib) = BasisnD%EneH0(ib) +                 &
                              BasisnD%tab_basisPrimSG(L,i)%EneH0(nDval(i))
              END DO
            END DO

          CASE DEFAULT
            write(out_unit,*) ' ERROR in',name_sub
            write(out_unit,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
            write(out_unit,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT
        END IF

        IF (debug) THEN
          write(out_unit,*) '<d0b(:,ib) I H0 I d0b(:,ib)>'
          DO i=1,size(BasisnD%EneH0)
            RWU_E  = REAL_WU(BasisnD%EneH0(i),'au','E')
            write(out_unit,*) i,RWU_Write(RWU_E,WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
          END DO
          flush(out_unit)
        END IF

      ELSE ! packed basis
        IF (debug) write(out_unit,*) 'RecSet_EneH0: packed basis'

        IF (.NOT. BasisnD%auto_contrac) THEN ! Because, it is done with an auto_contracted basis
          CALL Set_EneH0_OF_PackedBasis(BasisnD,para_Tnum,mole,para_ReadOp)
        ELSE
          IF (debug) THEN
            write(out_unit,*) 'auto_contrac basis <d0b(:,ib) I H0 I d0b(:,ib)>'
            DO i=1,size(BasisnD%EneH0)
              RWU_E  = REAL_WU(BasisnD%EneH0(i),'au','E')
              write(out_unit,*) i,RWU_Write(RWU_E,WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
            END DO
            flush(out_unit)
          END IF
        END IF
      END IF ! for BasisnD%nb_basis > 0 .AND. .NOT. BasisnD%packed_done

      IF (BasisnD%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unit,*) '==============================================='
        write(out_unit,*) '= END RecSet_EneH0 for Layer:',rec
        write(out_unit,*) '==================================================='
        flush(out_unit)
      END IF
      rec = rec - 1

      END SUBROUTINE RecSet_EneH0
!=======================================================================================

!=======================================================================================
      SUBROUTINE Set_EneH0_OF_PackedBasis(basis_Set,para_Tnum,mole,para_ReadOp)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

      !----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),    intent(in)    :: mole
      TYPE (Tnum),         intent(in)    :: para_Tnum

      !----- for the basis set ----------------------------------------------
      TYPE (basis),        intent(inout) :: basis_Set

      !----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)    :: para_ReadOp

!----- local variables -----------------------------------------------
      real(kind=Rkind), allocatable :: Hmat(:,:)
      real (kind=Rkind)             :: pot_Qref
      integer                       :: i
      TYPE(REAL_WU)                 :: RWU_E
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Set_EneH0_OF_PackedBasis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.TRUE.)
      END IF
!---------------------------------------------------------------------
      IF (allocated(basis_Set%EneH0))    THEN
        CALL dealloc_NParray(basis_Set%EneH0,"basis_Set%EneH0",name_sub)
      END IF
      CALL alloc_NParray(basis_Set%EneH0,[basis_Set%nb],"basis_Set%EneH0",name_sub)



      CALL Get_Hmat_FROM_basis(Hmat,basis_Set,para_Tnum,mole,para_ReadOp,pot_Qref)

      !----- diagonal elements of Rmat ---------------------------------
      ! Warning (DML 12/01/2021), pot0 remove from EneH0, otherwise ...
      ! ... %pot_Qref is present several times in the mutidimentional basis (DP, Smolyak)
      DO i=1,size(basis_set%EneH0)
        basis_set%EneH0(i) = Hmat(i,i) - pot_Qref
      END DO

      IF (print_level > -1 .OR. debug) THEN
        IF(MPI_id==0) write(out_unit,*) '<d0b(:,ib) I H0 I d0b(:,ib)>'
        DO i=1,size(basis_set%EneH0)
          RWU_E  = REAL_WU(basis_set%EneH0(i),'au','E')
          IF(MPI_id==0) write(out_unit,*) i,RWU_Write(RWU_E,WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
        END DO
        flush(out_unit)
      END IF

      !-----------------------------------------------------------------
      ! deallocation ....
      IF (allocated(Hmat))  CALL dealloc_NParray(Hmat, 'Hmat', name_sub)
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.TRUE.)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)

  END SUBROUTINE Set_EneH0_OF_PackedBasis
  SUBROUTINE Set_EneH0_OF_ContracBasis(basis_Set,para_Tnum,mole,para_ReadOp)

      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

      !----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),    intent(in)    :: mole
      TYPE (Tnum),         intent(in)    :: para_Tnum
      !----- for the basis set ----------------------------------------------
      TYPE (basis),        intent(inout) :: basis_Set
      !----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)    :: para_ReadOp

!----- local variables -----------------------------------------------
      real(kind=Rkind), allocatable :: Hmat(:,:)
      real (kind=Rkind)             :: pot_Qref
      integer                       :: i
      TYPE(REAL_WU)                 :: RWU_E
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Set_EneH0_OF_ContracBasis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.TRUE.)
      END IF
!---------------------------------------------------------------------


      CALL Get_Hmat_FROM_basis(Hmat,basis_Set,para_Tnum,mole,para_ReadOp,pot_Qref)

      IF (.NOT. allocated(basis_set%Rvec)) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' basis_set%Rvec is not allocated!'
        write(out_unit,*) ' CHECK the source!!'
        STOP
      END IF

      IF (allocated(basis_Set%EneH0))    THEN
        CALL dealloc_NParray(basis_Set%EneH0,"basis_Set%EneH0",name_sub)
      END IF
      CALL alloc_NParray(basis_Set%EneH0,[basis_Set%nbc],"basis_Set%EneH0",name_sub)

      !----- diagonal elements of Rmat ---------------------------------
      ! Warning (DML 12/01/2021), pot0 remove from EneH0, otherwise ...
      ! ... %pot_Qref is present several times in the mutidimentional basis (DP, Smolyak)
      DO i=1,size(basis_set%EneH0)
        basis_set%EneH0(i) = dot_product(basis_set%Rvec(:,i),                   &
                                    matmul(Hmat,basis_set%Rvec(:,i))) - pot_Qref
      END DO

      IF (print_level > -1 .OR. debug) THEN
        IF(MPI_id==0) write(out_unit,*) '<d0b(:,ib) I H0 I d0b(:,ib)>'
        DO i=1,size(basis_set%EneH0)
          RWU_E  = REAL_WU(basis_set%EneH0(i),'au','E')
          IF(MPI_id==0) write(out_unit,*) i,RWU_Write(RWU_E,WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
        END DO
        flush(out_unit)
      END IF

      !-----------------------------------------------------------------
      ! deallocation ....
      IF (allocated(Hmat))  CALL dealloc_NParray(Hmat, 'Hmat', name_sub)
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.TRUE.)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)

  END SUBROUTINE Set_EneH0_OF_ContracBasis

  ! basis parameters: for HO basis set (scaleQ)
  SUBROUTINE AutoParam_basis(basis_Set,para_Tnum,mole,para_ReadOp)

      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(in) :: mole
      TYPE (Tnum),      intent(in) :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: basis_Set

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)         :: para_ReadOp



      logical       :: HObasis


      TYPE (basis)      :: basis_temp
      integer           :: print_level_loc

!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'AutoParam_basis'
!---------------------------------------------------------------------

      ! test for HO basis
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        CALL RecWrite_basis(basis_Set,write_all=.FALSE.)
      END IF
!---------------------------------------------------------------------
      print_level_loc = print_level

      IF (debug) THEN
        CALL set_print_level(2,force=.TRUE.)
      ELSE
        CALL set_print_level(-1,force=.TRUE.)
      END IF

      CALL basis2TObasis1(basis_temp,basis_Set)

      basis_temp%nb            = 2
      CALL Set_nq_OF_basis(basis_temp,nq=10)
      basis_temp%With_L        = .FALSE.
      basis_temp%L_SparseGrid  = -1
      basis_temp%L_SparseBasis = -1
      basis_temp%Type_OF_nDindB = 1

      CALL construct_primitive_basis(basis_temp)


      !write(out_unit,*) 'In ',name_sub,' basis_temp%type',basis_temp%type,HObasis
      HObasis = (basis_temp%type == 20 .OR. basis_temp%type == 21 .OR.     &
                 basis_temp%type == 200 .OR. basis_temp%type == 201)

      IF (HObasis .AND. basis_temp%auto_basis) THEN
        write(out_unit,*) 'Initial Q0 and scaleQ',basis_Set%Q0(1),basis_Set%scaleQ(1)
        CALL AutoParam_basis_scaleQ(basis_temp,para_Tnum,mole,para_ReadOp)
        CALL AutoParam_basis_Q0(basis_temp,para_Tnum,mole,para_ReadOp)
        CALL AutoParam_basis_scaleQ(basis_temp,para_Tnum,mole,para_ReadOp)

        basis_Set%Q0(1)      = basis_temp%Q0(1)
        basis_Set%scaleQ(1)  = basis_temp%scaleQ(1)
      END IF

      basis_Set%auto_basis = .FALSE.
      CALL set_print_level(print_level_loc,force=.TRUE.)

      CALL dealloc_basis(basis_temp)
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.TRUE.)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)

      END SUBROUTINE AutoParam_basis
      ! basis parameters: for HO basis set (scaleQ)
      SUBROUTINE AutoParam_basis_scaleQ(basis_temp,para_Tnum,mole,para_ReadOp)

      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(in) :: mole
      TYPE (Tnum),      intent(in) :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: basis_temp

!----- variables pour la namelist minimum ----------------------------
      integer :: nb_scalar_Op,type_HamilOp
      logical :: calc_scalar_Op,direct_KEO

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)         :: para_ReadOp


       integer       :: i

      real (kind=Rkind) :: Ene1_10SQrange(29)      !  range [1.0, 1.5 ....14.5 15.0]
      real (kind=Rkind) :: Enep1_p9SQrange(9)      !  range [0.1 0.2            0.9] around the optimal value of the previous range
      real (kind=Rkind) :: Enep01_p09SQrange(9)    !  range [0.01 0.02         0.09]
      real (kind=Rkind) :: Enep001_p009SQrange(19) !  range [0.001 0.002      0.009] around the optimal value of the previous range

      real (kind=Rkind) :: E,Emin,SQmin
      real (kind=Rkind) :: auTOcm_inv

!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'AutoParam_basis_scaleQ'
!---------------------------------------------------------------------

      ! test for HO basis
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        !CALL RecWrite_basis(basis_temp,write_all=.FALSE.)
      END IF
!---------------------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      basis_temp%scaleQ(1) = ONE


      Emin = HUGE(ONE)
      SQMin = basis_temp%scaleQ(1)

      DO i=1,size(Ene1_10SQrange)

        CALL dealloc_dnb_OF_basis(basis_temp)
        CALL dealloc_xw_OF_basis(basis_temp)
        CALL dealloc_nDindex(basis_temp%nDindG)

        CALL construct_primitive_basis(basis_temp)

        E = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)
        IF (debug) write(out_unit,*) 'SQ,E',basis_temp%scaleQ(1),E*auTOcm_inv
        flush(out_unit)

        IF (E < Emin) THEN
          Emin  = E
          SQMin = basis_temp%scaleQ(1)
        END IF
        Ene1_10SQrange(i) = E

        basis_temp%scaleQ(1) = basis_temp%scaleQ(1) + HALF

      END DO

      !short range (evry 0.1)
      IF (SQMin == ONE) THEN
        basis_temp%scaleQ(1) = ONETENTH
      ELSE IF (SQMin == 15._Rkind) THEN
        basis_temp%scaleQ(1) = 15._Rkind
      ELSE
        basis_temp%scaleQ(1) = SQMin-0.4_Rkind
      END IF
      DO i=1,size(Enep1_p9SQrange)

        CALL dealloc_dnb_OF_basis(basis_temp)
        CALL dealloc_xw_OF_basis(basis_temp)
        CALL construct_primitive_basis(basis_temp)

        E = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)
        IF (debug) write(out_unit,*) 'SQ,E',basis_temp%scaleQ(1),E*auTOcm_inv

        IF (E < Emin) THEN
          Emin  = E
          SQMin = basis_temp%scaleQ(1)
        END IF
        Enep1_p9SQrange(i) = E

        basis_temp%scaleQ(1) = basis_temp%scaleQ(1) + ONETENTH

      END DO


      IF (SQMin == ONETENTH) THEN
        basis_temp%scaleQ(1) = ONETENTH**2
        DO i=1,size(Enep01_p09SQrange)

        CALL dealloc_dnb_OF_basis(basis_temp)
        CALL dealloc_xw_OF_basis(basis_temp)
        CALL construct_primitive_basis(basis_temp)

        E = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)
        IF (debug) write(out_unit,*) 'SQ,E',basis_temp%scaleQ(1),E*auTOcm_inv

        IF (E < Emin) THEN
          Emin  = E
          SQMin = basis_temp%scaleQ(1)
        END IF
        Enep01_p09SQrange(i) = E

        basis_temp%scaleQ(1) = basis_temp%scaleQ(1) + ONETENTH**2

        END DO

        !short range (evry 0.001)
        IF (SQMin == ONETENTH**2) THEN
          basis_temp%scaleQ(1) = ONETENTH**3
        ELSE
          basis_temp%scaleQ(1) = SQMin-0.009_Rkind
        END IF
        DO i=1,size(Enep001_p009SQrange)

          CALL dealloc_dnb_OF_basis(basis_temp)
          CALL dealloc_xw_OF_basis(basis_temp)
          CALL construct_primitive_basis(basis_temp)

          E = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)
          IF (debug) write(out_unit,*) 'SQ,E',basis_temp%scaleQ(1),E*auTOcm_inv

          IF (E < Emin) THEN
            Emin  = E
            SQMin = basis_temp%scaleQ(1)
          END IF
          Enep001_p009SQrange(i) = E

          basis_temp%scaleQ(1) = basis_temp%scaleQ(1) + ONETENTH**3

        END DO


      END IF
      write(out_unit,*) 'Optimal scaleQ, E(cm-1)',SQMin,Emin*auTOcm_inv


      basis_temp%scaleQ(1)  = SQMin

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        !CALL RecWrite_basis(basis_temp,write_all=.TRUE.)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)

      END SUBROUTINE AutoParam_basis_scaleQ
      SUBROUTINE AutoParam_basis_Q0(basis_temp,para_Tnum,mole,para_ReadOp)

      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(in) :: mole
      TYPE (Tnum),      intent(in) :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: basis_temp

!----- variables pour la namelist minimum ----------------------------
      integer :: nb_scalar_Op,type_HamilOp
      logical :: calc_scalar_Op,direct_KEO

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)         :: para_ReadOp


      integer           :: i,type_Q
      real (kind=Rkind) :: Ene1_Q0range(21)      !  range Q0 +/- scaleQ/10*i
      real (kind=Rkind) :: Ene1_Q0range2(21)      !  range Q0 +/- scaleQ/100*i
      real (kind=Rkind) :: E,Emin,Q0min,Q0ini
      real (kind=Rkind) :: auTOcm_inv

!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'AutoParam_basis_Q0'
!---------------------------------------------------------------------

      ! test for HO basis
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        !CALL RecWrite_basis(basis_temp,write_all=.FALSE.)
      END IF
!---------------------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')


      !write(out_unit,*) 'type_Qin',mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin(:)
      !write(out_unit,*) 'iQdyn',basis_temp%iQdyn(1)
      !write(out_unit,*) 'type_Qin',mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin(basis_temp%iQdyn(1))
      type_Q = mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin(basis_temp%iQdyn(1))


      Q0ini = basis_temp%Q0(1)
      CALL dealloc_dnb_OF_basis(basis_temp)
      CALL dealloc_xw_OF_basis(basis_temp)
      CALL dealloc_nDindex(basis_temp%nDindG)

      CALL construct_primitive_basis(basis_temp)

      Q0Min = Q0ini
      Emin  = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)

      !write(out_unit,*) 'Initial Q0,E',Q0ini,Emin*auTOcm_inv

      basis_temp%Q0(1) = Q0ini - basis_temp%scaleQ(1)
      IF (check_OutOfRange(basis_temp%Q0(1),type_Q) .OR. (type_Q==2 .AND. basis_temp%Q0(1) < Q0ini/TWO)) THEN
        IF (type_Q==2) basis_temp%Q0(1) = Q0ini/TWO
        CALL Set_InitRange(basis_temp%Q0(1),type_Q)
      END IF

      DO i=1,size(Ene1_Q0range)

        CALL dealloc_dnb_OF_basis(basis_temp)
        CALL dealloc_xw_OF_basis(basis_temp)
        CALL dealloc_nDindex(basis_temp%nDindG)

        CALL construct_primitive_basis(basis_temp)

        E = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)
        IF (debug) write(out_unit,*) 'Q0,E',basis_temp%Q0(1),E*auTOcm_inv
        flush(out_unit)

        IF (E < Emin) THEN
          Emin  = E
          Q0Min = basis_temp%Q0(1)
        END IF
        Ene1_Q0range(i) = E

        basis_temp%Q0(1) = basis_temp%Q0(1) + basis_temp%scaleQ(1)/TEN

        IF (check_OutOfRange(basis_temp%Q0(1),type_Q)) EXIT

      END DO


      basis_temp%Q0(1) = Q0Min - basis_temp%scaleQ(1)/TEN
      IF (check_OutOfRange(basis_temp%Q0(1),type_Q) .OR. (type_Q==2 .AND. basis_temp%Q0(1) < Q0ini/TWO)) THEN
        IF (type_Q==2) basis_temp%Q0(1) = Q0ini/TWO
        CALL Set_InitRange(basis_temp%Q0(1),type_Q)
      END IF

      DO i=1,size(Ene1_Q0range2)

        CALL dealloc_dnb_OF_basis(basis_temp)
        CALL dealloc_xw_OF_basis(basis_temp)
        CALL dealloc_nDindex(basis_temp%nDindG)

        CALL construct_primitive_basis(basis_temp)

        E = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)
        IF (debug) write(out_unit,*) 'Q0,E',basis_temp%Q0(1),E*auTOcm_inv
        flush(out_unit)

        IF (E < Emin) THEN
          Emin  = E
          Q0Min = basis_temp%Q0(1)
        END IF
        Ene1_Q0range2(i) = E

        basis_temp%Q0(1) = basis_temp%Q0(1) + basis_temp%scaleQ(1)/TEN**2

        IF (check_OutOfRange(basis_temp%Q0(1),type_Q)) EXIT


      END DO

      basis_temp%Q0(1) = Q0Min - basis_temp%scaleQ(1)/TEN**2
      IF (check_OutOfRange(basis_temp%Q0(1),type_Q) .OR. (type_Q==2 .AND. basis_temp%Q0(1) < Q0ini/TWO)) THEN
        IF (type_Q==2) basis_temp%Q0(1) = Q0ini/TWO
        CALL Set_InitRange(basis_temp%Q0(1),type_Q)
      END IF

      DO i=1,size(Ene1_Q0range2)

        CALL dealloc_dnb_OF_basis(basis_temp)
        CALL dealloc_xw_OF_basis(basis_temp)
        CALL dealloc_nDindex(basis_temp%nDindG)

        CALL construct_primitive_basis(basis_temp)

        E = Ene_FROM_basis(basis_temp,para_Tnum,mole,para_ReadOp)
        IF (debug) write(out_unit,*) 'Q0,E',basis_temp%Q0(1),E*auTOcm_inv
        flush(out_unit)

        IF (E < Emin) THEN
          Emin  = E
          Q0Min = basis_temp%Q0(1)
        END IF
        Ene1_Q0range2(i) = E

        basis_temp%Q0(1) = basis_temp%Q0(1) + basis_temp%scaleQ(1)/TEN**3

        IF (check_OutOfRange(basis_temp%Q0(1),type_Q)) EXIT


      END DO

      write(out_unit,*) 'Optimal Q0, E (cm-1)',Q0Min,Emin*auTOcm_inv


      basis_temp%Q0(1)  = Q0Min

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        !CALL RecWrite_basis(basis_temp,write_all=.TRUE.)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)

      END SUBROUTINE AutoParam_basis_Q0
      FUNCTION check_OutOfRange(Q,type_Q)
      USE EVR_system_m
      IMPLICIT NONE

      logical                          :: check_OutOfRange
      integer,           intent(in)    :: type_Q
      real (kind=Rkind), intent(inout) :: Q


      check_OutOfRange = (type_Q ==  3 .AND. Q <= ZERO)  .OR.   &
                         (type_Q ==  3 .AND. Q >= Pi)    .OR.   &
                         (type_Q == -3 .AND. Q <= -ONE)  .OR.   &
                         (type_Q == -3 .AND. Q >=  ONE)  .OR.   &
                         (type_Q ==  2 .AND. Q <= ZERO)

      END FUNCTION check_OutOfRange
      SUBROUTINE Set_InitRange(Q,type_Q)
      USE EVR_system_m
      IMPLICIT NONE

      integer,           intent(in)    :: type_Q
      real (kind=Rkind), intent(inout) :: Q


      IF ( type_Q == 3 .AND. Q <= ZERO) Q = ONETENTH
      IF ( type_Q == 3 .AND. Q >= Pi)   Q = PI/TWO

      IF ( type_Q == -3 .AND. Q < -ONE) Q = -ONE + ONETENTH
      IF ( type_Q == -3 .AND. Q >  ONE) Q = ZERO

      IF ( type_Q == 2 .AND. Q <= ZERO) Q = HALF

      END SUBROUTINE Set_InitRange
  FUNCTION Ene_FROM_basis(basis_Set,para_Tnum,mole,para_ReadOp)

      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

      real(kind=Rkind) :: Ene_FROM_basis

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),    intent(in)  :: mole
      TYPE (Tnum),         intent(in)  :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis),        intent(in)  :: basis_Set

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)  :: para_ReadOp

!----- local variables -----------------------------------------------
      real(kind=Rkind), allocatable :: Hmat(:,:),Rdiag(:),Rvp(:,:)
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Ene_FROM_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.FALSE.)
      END IF
!---------------------------------------------------------------------

      CALL Get_Hmat_FROM_basis(Hmat,basis_Set,para_Tnum,mole,para_ReadOp)

      !---------------------------------------------------------------
      ! digonalization of Hmat
      CALL alloc_NParray(Rdiag,[basis_Set%nb],'Rdiag',name_sub)
      CALL alloc_NParray(Rvp,shape(Hmat),'Rvp',name_sub)

      CALL sub_diago_H(Hmat,Rdiag,Rvp,basis_Set%nb,.TRUE.)

      Ene_FROM_basis = Rdiag(1)

      !-----------------------------------------------------------------
      ! deallocation ....
      IF (allocated(Rvp))   CALL dealloc_NParray(Rvp,  'Rvp',  name_sub)
      IF (allocated(Rdiag)) CALL dealloc_NParray(Rdiag,'Rdiag',name_sub)
      IF (allocated(Hmat))  CALL dealloc_NParray(Hmat, 'Hmat', name_sub)
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.TRUE.)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)

  END FUNCTION Ene_FROM_basis

  SUBROUTINE Get_Hmat_FROM_basis(Hmat,basis_Set,para_Tnum,mole,para_ReadOp,pot_Qref)

      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- Hamiltonian matrix -----------------------------------------------
      real(kind=Rkind), allocatable,  intent(inout)           :: Hmat(:,:)
      real(kind=Rkind),               intent(inout), optional :: pot_Qref

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),               intent(in)  :: mole
      TYPE (Tnum),                    intent(in)  :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis),                   intent(in)  :: basis_Set

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t),            intent(in)  :: para_ReadOp


!----- local variables -----------------------------------------------
!----- Operators and para_ReadOp ----------------
      TYPE (ReadOp_t)         :: ReadOp_AutoBasis
      TYPE (param_AllOp)          :: para_AllOp_loc
!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis)       :: para_AllBasis_loc
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)            :: mole_loc
      TYPE (Tnum)                 :: para_Tnum_loc

      integer                     :: i,nq
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Get_Hmat_FROM_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.FALSE.)
      END IF
!---------------------------------------------------------------------

      ! modification of mole => mole_loc and Tnum => Tnum_loc
      para_Tnum_loc                   = para_Tnum
      para_Tnum_loc%JJ                = 0
      para_Tnum_loc%NonGcteRange(:)   = 0

      mole_loc                        = mole
      mole_loc%Cart_transfo           = .FALSE.
      ! If needed, change RPH transfo in flexible transfo
      CALL CoordTypeRPH_TO_CoordTypeFlex(mole_loc)

      CALL basis_TO_Allbasis(basis_Set,para_AllBasis_loc,mole_loc)
      !CALL RecWrite_basis(para_AllBasis_loc%BasisnD,write_all=.FALSE.)

      nq = get_nq_FROM_basis(para_AllBasis_loc%BasisnD)
      CALL ReadOp2_TO_ReadOp1_FOR_AutoBasis(ReadOp_AutoBasis,para_ReadOp,nq)

      !---------------------------------------------------------------
      ! make Operators: H and S
      ! allocation of tab_Op
      para_AllOp_loc%nb_Op = 2 ! just H and S
      CALL alloc_array(para_AllOp_loc%tab_Op,[para_AllOp_loc%nb_Op],            &
                      'para_AllOp_loc%tab_Op',name_sub)
      !i=1 => for H
      CALL All_param_TO_para_H(para_Tnum_loc,mole_loc,para_AllBasis_loc,        &
                               ReadOp_AutoBasis,para_AllOp_loc%tab_Op(1))


      ! old direct=2 with a matrix
      para_AllOp_loc%tab_Op(1)%make_Mat                                = .TRUE.
      para_AllOp_loc%tab_Op(1)%para_ReadOp%para_FileGrid%Save_MemGrid  = .TRUE.
      para_AllOp_loc%tab_Op(1)%para_ReadOp%para_FileGrid%Save_FileGrid = .FALSE.

      !i=2 => for S
      i=2
      CALL param_Op1TOparam_Op2(para_AllOp_loc%tab_Op(1),                       &
                                para_AllOp_loc%tab_Op(i))
      para_AllOp_loc%tab_Op(i)%name_Op = 'S'
      para_AllOp_loc%tab_Op(i)%n_Op    = -1

      CALL Init_TypeOp(para_AllOp_loc%tab_Op(i)%param_TypeOp,                   &
                       type_Op=0,nb_Qact=mole_loc%nb_act1,cplx=.FALSE.,         &
                       JRot=para_Tnum_loc%JJ,direct_KEO=.FALSE.,direct_ScalOp=.FALSE.)
      CALL derive_termQact_TO_derive_termQdyn(                                  &
                            para_AllOp_loc%tab_Op(i)%derive_termQdyn,           &
                            para_AllOp_loc%tab_Op(i)%derive_termQact,           &
                              mole_loc%ActiveTransfo%list_QactTOQdyn)

      !---------------------------------------------------------------
      ! make the Grid
      CALL sub_qa_bhe(para_AllOp_loc)

      IF (present(pot_Qref)) pot_Qref = para_AllOp_loc%tab_Op(1)%para_ReadOp%pot_Qref

      !---------------------------------------------------------------
      ! make the matrix of H
      CALL sub_MatOp(para_AllOp_loc%tab_Op(1),debug)

      CALL move_alloc(FROM=para_AllOp_loc%tab_Op(1)%Rmat,TO=Hmat)

      !-----------------------------------------------------------------
      ! deallocation ....
      CALL dealloc_AllBasis(para_AllBasis_loc)
      CALL dealloc_para_AllOp(para_AllOp_loc)
      CALL dealloc_CoordType(mole_loc)
      CALL dealloc_Tnum(para_Tnum_loc)
      CALL dealloc_ReadOp(ReadOp_AutoBasis)


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        !CALL RecWrite_basis(basis_Set,write_all=.TRUE.)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)

  END SUBROUTINE Get_Hmat_FROM_basis

  SUBROUTINE Autocontract_basis(basis_AutoContract,para_Tnum,mole,para_ReadOp)

      USE EVR_system_m
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      USE BasisMakeGrid
      USE mod_Op
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (basis),        intent(inout) :: basis_AutoContract

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),    intent(in)    :: mole
      TYPE (Tnum),         intent(in)    :: para_Tnum

!----- variables for the construction of H ---------------------------
      TYPE (ReadOp_t), intent(in)    :: para_ReadOp


!----- local variables -----------------------------------------------
      real(kind=Rkind), allocatable :: Hmat(:,:),Rdiag(:),Rvp(:,:)
      integer                       :: i,nbc
      real(kind=Rkind)              :: pot_Qref,ene0
      TYPE(REAL_WU)                 :: RWU_E,RWU_DE
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Autocontract_basis'
!---------------------------------------------------------------------
      IF (print_level > -1) write(out_unit,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unit,*)
        !write(out_unit,*) 'BEGINNING ',name_sub
        !CALL Write_CoordType(mole_loc)
        write(out_unit,*) 'basis before contraction'
        write(out_unit,*)
        CALL RecWrite_basis(basis_AutoContract)
      END IF
!---------------------------------------------------------------------
      IF (basis_AutoContract%max_nbc < 1) basis_AutoContract%max_nbc =  &
                                                 basis_AutoContract%nb
      IF (basis_AutoContract%min_nbc < 1) basis_AutoContract%min_nbc = 2
      IF (basis_AutoContract%min_nbc > basis_AutoContract%max_nbc)      &
                basis_AutoContract%min_nbc = basis_AutoContract%max_nbc

      IF (print_level > -1) THEN
        write(out_unit,*) 'min_nbc,max_nbc,max_ene_contrac (ua)',      &
                 basis_AutoContract%min_nbc,basis_AutoContract%max_nbc, &
                                    basis_AutoContract%max_ene_contrac
      END IF

      CALL Get_Hmat_FROM_basis(Hmat,basis_AutoContract,para_Tnum,mole,para_ReadOp,pot_Qref)

      !---------------------------------------------------------------
      ! digonalization of Hmat
      CALL alloc_NParray(Rdiag,[basis_AutoContract%nb],'Rdiag',name_sub)
      CALL alloc_NParray(Rvp,shape(Hmat),'Rvp',name_sub)

      CALL sub_diago_H(Hmat,Rdiag,Rvp, basis_AutoContract%nb,.TRUE.)

      IF (debug) THEN
        write(out_unit,*) 'Eigenvalues for the contraction'
        CALL Write_Vec_MPI(Rdiag,out_unit,5)
        write(out_unit,*) 'Eigenvectors for the contraction'
        CALL Write_Mat_MPI(Rvp,out_unit,5)
      END IF
      !---------------------------------------------------------------
      ! Energy levels + the new nbc
      ene0 = Rdiag(1)

      IF (basis_AutoContract%nbc > 0) THEN
        nbc = basis_AutoContract%nbc
      ELSE
        nbc = count(Rdiag(:) <= basis_AutoContract%max_ene_contrac+ene0)
        IF (nbc < basis_AutoContract%min_nbc)                           &
                                        nbc = basis_AutoContract%min_nbc
        IF (nbc > basis_AutoContract%max_nbc)                           &
                                        nbc = basis_AutoContract%max_nbc
      END IF
      IF (nbc > basis_AutoContract%nb) nbc = basis_AutoContract%nb

      basis_AutoContract%nqc = nbc + basis_AutoContract%nqPLUSnbc_TO_nqc
      basis_AutoContract%nqPLUSnbc_TO_nqc = 0

      IF (print_level > -1) THEN
        write(out_unit,*) 'nb levels:',nbc

        write(out_unit,*) 'levels: '
        DO i=1,nbc
          RWU_E  = REAL_WU(Rdiag(i),     'au','E')
          RWU_DE = REAL_WU(Rdiag(i)-ene0,'au','E')
          write(out_unit,*) i,RWU_Write(RWU_E ,WithUnit=.FALSE.,WorkingUnit=.FALSE.),&
                           " ",RWU_Write(RWU_DE,WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
        END DO
        flush(out_unit)
      END IF
      basis_AutoContract%nbc = nbc

      ! Energy levels
      IF (allocated(basis_AutoContract%EneH0))    THEN
        CALL dealloc_NParray(basis_AutoContract%EneH0,"basis_set%EneH0",name_sub)
      END IF
      CALL alloc_NParray(basis_AutoContract%EneH0,[nbc],                    &
                        "basis_AutoContract%EneH0",name_sub)
      basis_AutoContract%EneH0(:) = Rdiag(1:nbc) - pot_Qref

      !---------------------------------------------------------------
      ! contraction => nb=nbc
      IF (print_level > -1) write(out_unit,*) 'alloc Rvec',allocated(basis_AutoContract%Rvec)
      IF (allocated(basis_AutoContract%Rvec))  THEN
        CALL dealloc_NParray(basis_AutoContract%Rvec,                     &
                            "basis_AutoContract%Rvec",name_sub)
      END IF
      CALL alloc_NParray(basis_AutoContract%Rvec,                         &
                     [basis_AutoContract%nb,basis_AutoContract%nb], &
                                     "basis_AutoContract%Rvec",name_sub)

      IF (basis_AutoContract%POGridRep_polyortho .AND. basis_AutoContract%ndim == 1) THEN
        CALL sub_make_polyorthobasis(basis_AutoContract,Rvp)
      ELSE
        basis_AutoContract%Rvec = Rvp
        CALL sub_contraction_basis(basis_AutoContract,.TRUE.)
      END IF
      !basis_AutoContract%POGridRep_polyortho = .FALSE.
      IF (debug) THEN
        write(out_unit,*) 'Eigenvectors on the grid'
        DO i=1,get_nq_FROM_basis(basis_AutoContract)
          write(out_unit,*) i,basis_AutoContract%x(:,i),basis_AutoContract%dnRGB%d0(i,:)
        END DO
        flush(out_unit)
      END IF

      !-----------------------------------------------------------------
      ! deallocation ....
      IF (allocated(Rvp))   CALL dealloc_NParray(Rvp,  'Rvp',  name_sub)
      IF (allocated(Rdiag)) CALL dealloc_NParray(Rdiag,'Rdiag',name_sub)
      IF (allocated(Hmat))  CALL dealloc_NParray(Hmat, 'Hmat', name_sub)
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'basis after contraction'
        write(out_unit,*)
        CALL RecWrite_basis(basis_AutoContract)
        write(out_unit,*)
        !write(out_unit,*) 'END ',name_sub
      END IF
      IF (print_level > -1) write(out_unit,*) 'END ',name_sub
      flush(out_unit)

  END SUBROUTINE Autocontract_basis
  SUBROUTINE basis_TO_AllBasis(basis_temp,Allbasis,mole)
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: basis_temp
      TYPE (param_AllBasis) :: AllBasis


       integer :: i,iact,isym,i_Q

!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'basis_TO_Allbasis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub

        write(out_unit,*) 'basis'
        write(out_unit,*)
        CALL RecWrite_basis(basis_temp)

      END IF
!---------------------------------------------------------------------
      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 1)                &
            mole%ActiveTransfo%list_act_OF_Qdyn(i) =                    &
                                        basis_temp%auto_contrac_type1_TO
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 21)               &
            mole%ActiveTransfo%list_act_OF_Qdyn(i) =                    &
                                       basis_temp%auto_contrac_type21_TO
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 22)               &
            mole%ActiveTransfo%list_act_OF_Qdyn(i) =                    &
                                       basis_temp%auto_contrac_type21_TO
      END DO

        DO i=1,basis_temp%ndim
          isym = basis_temp%iQdyn(i)
          mole%ActiveTransfo%list_act_OF_Qdyn(isym) = 1
        END DO

        CALL type_var_analysis_OF_CoordType(mole)

        IF (debug) THEN
          write(out_unit,*) 'mole:'
          CALL Write_CoordType(mole)
        END IF
        !---------------------------------------------------------------
        ! make AllBasis with ONE basis set

        CALL alloc_AllBasis(AllBasis)

        AllBasis%nb_be = 1

        CALL basis2TObasis1(AllBasis%BasisnD,basis_temp)

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'param of Allbasis '
        write(out_unit,*)
        write(out_unit,*) 'nb_ba,nb_qa',AllBasis%BasisnD%nb,           &
                                     get_nq_FROM_basis(AllBasis%BasisnD)
        write(out_unit,*) 'nDindB%nDsize',AllBasis%BasisnD%nDindB%nDsize
        write(out_unit,*) 'nDindG%nDsize',AllBasis%BasisnD%nDindG%nDsize
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

  END SUBROUTINE basis_TO_AllBasis

  SUBROUTINE All_param_TO_para_H(para_Tnum,mole,para_AllBasis,para_ReadOp,para_H)

      USE EVR_system_m
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_analysis
      USE mod_propa
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),target :: mole
      TYPE (Tnum),     target :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis),target :: para_AllBasis

!----- variables for the construction of H ---------------------------
      TYPE (param_Op)     :: para_H
      TYPE (ReadOp_t), intent(in) :: para_ReadOp


!----- working variables ---------------------------------------------
      integer :: i,j,i_term,rk,rl
      integer :: err_mem,memory
!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'All_param_TO_para_H'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
      END IF


!---------------------------------------------------------------------
!     - copy parameters in para_H --------
      para_H%alloc         = .FALSE.
      para_H%init_var      = .TRUE.

!     ------------------------------
      para_H%mole      => mole
      para_H%para_Tnum => para_Tnum
!     ------------------------------

!     ------------------------------
      para_H%para_AllBasis => para_AllBasis
      para_H%BasisnD       => para_H%para_AllBasis%BasisnD
      para_H%Basis2n       => para_H%para_AllBasis%Basis2n
!     ------------------------------

      para_H%para_ReadOp      = para_ReadOp

      para_H%n_Op             = 0
      para_H%name_Op          = 'H'


      para_H%nb_OpPsi         = 0

      para_H%nb_ba            = get_nb_FROM_basis(para_AllBasis%BasisnD)
      IF (para_ReadOp%Op_WithContracRVec) THEN
        para_H%nbc_ba         = para_AllBasis%BasisnD%nDindB_contracted%Max_nDI
      ELSE
        para_H%nbc_ba         = para_H%nb_ba
      END IF
      para_H%nb_qa            = get_nqa_FROM_basis(para_AllBasis%BasisnD)

      para_H%nb_qa_WithNoGrid = get_nqa_FROM_basis_withNoGrid(para_AllBasis%BasisnD)

      IF (mole%nb_inact2n == 0) THEN
        para_H%nb_bi          = 1
      ELSE
        para_H%nb_bi          = get_nb_bi_FROM_AllBasis(para_AllBasis)
      END IF
      para_H%nb_be            = para_ReadOp%nb_elec


      para_H%nb_bie           = para_H%nb_bi * para_H%nb_be

      para_H%nb_bai           = para_H%nb_ba * para_H%nb_bi
      para_H%nb_qai           = para_H%nb_qa * para_H%nb_bi

      para_H%nb_baie          = para_H%nb_bai * para_H%nb_be
      para_H%nb_qaie          = para_H%nb_qai * para_H%nb_be
      para_H%nq_tot           = para_H%nb_qa_WithNoGrid * para_H%nb_bie


      para_H%nb_bRot          = para_ReadOp%nb_bRot

      para_H%nb_tot           = para_H%nbc_ba * para_H%nb_bie * para_H%nb_bRot
      para_H%nb_tot_ini       = para_H%nb_baie * para_H%nb_bRot

      para_H%Make_Mat         = para_ReadOp%Make_Mat
      para_H%pack_Op          = para_ReadOp%pack_Op
      para_H%read_Op          = para_ReadOp%read_Op
      para_H%tol_pack         = para_ReadOp%tol_pack
      para_H%tol_nopack       = para_ReadOp%tol_nopack
      para_H%spectral         = para_ReadOp%spectral
      para_H%spectral_Op      = para_ReadOp%spectral_Op


      para_H%nb_act1          = mole%nb_act1

      CALL Init_TypeOp(para_H%param_TypeOp,                             &
                       type_Op=para_ReadOp%type_HamilOp,                &
                       nb_Qact=mole%nb_act1,cplx=para_ReadOp%pot_cplx,  &
                       JRot=para_H%Para_Tnum%JJ,                        &
                       direct_KEO=para_ReadOp%direct_KEO,               &
                       direct_ScalOp=para_ReadOp%direct_ScalOp)

      CALL derive_termQact_TO_derive_termQdyn(para_H%derive_termQdyn,   &
              para_H%derive_termQact,mole%ActiveTransfo%list_QactTOQdyn)

      para_H%sym_Hamil = .NOT. (Para_Tnum%nrho == 0 .OR. para_H%type_Op == 10)

      para_H%scaled          = .FALSE.
      para_H%E0              = ZERO
      para_H%Esc             = ONE
      para_H%Hmin            = huge(ONE)
      para_H%Hmax            = -huge(ONE)
      para_H%pot0            = para_ReadOp%pot0
      para_H%pot_only        = para_ReadOp%pot_only
      para_H%T_only          = para_ReadOp%T_only

      IF (debug) THEN
        CALL write_param_Op(para_H)
        write(out_unit,*) 'END ',name_sub
      END IF

  END SUBROUTINE All_param_TO_para_H

!=======================================================================================
!     adiabatic contraction
!=======================================================================================
  SUBROUTINE sub_MatOp_HADA(para_H,para_ana,para_intensity,para_AllOp,const_phys)
      USE EVR_system_m
      USE mod_Constant
      USE mod_nDindex
      USE mod_Op
      USE mod_psi,          ONLY : param_psi,alloc_psi,alloc_NParray,dealloc_NParray
      USE mod_propa
      USE mod_analysis
      USE mod_fullanalysis
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_AllOp)          :: para_AllOp
      TYPE (param_Op)             :: para_H


!----- physical and mathematical constants ---------------------------
      TYPE (constant)            :: const_phys

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)       :: para_ana
      TYPE (param_intensity) :: para_intensity


      ! local variables
      TYPE (param_AllBasis), target        :: para_AllBasis_loc
      TYPE (param_psi),      allocatable   :: Tab_Psi(:)
      TYPE (param_Op)                      :: para_H_HADA
      integer :: max_ana_save

!------ quadrature points and weight -----------------------------

      real (kind=Rkind) :: WnD


!----- for the gestion of the memory ---------------------------------
      integer   :: nreste,nplus


      real (kind=Rkind) :: Qdyn(para_H%mole%nb_var)
      real (kind=Rkind) :: Qact(para_H%mole%nb_act1)
      real (kind=Rkind), allocatable :: H_HADA(:,:,:)

!------ for td0b ...         -------------------------------------
      integer                           :: kmem
      real (kind=Rkind), allocatable    :: td0b(:,:)
      TYPE (param_d0MatOp), allocatable :: d0MatOpd0bWrho(:,:)
      TYPE (param_d0MatOp) :: d0MatOp
      integer                           :: type_Op

!----- divers ----------------------------------------------------
      integer  :: nb_ba,nb_bie,print_psi_save
      real (kind=Rkind) :: ZPE

      integer  :: i,k,nb_blocks

      integer  :: iterm_Op,i1_h,ib1,ib2

      integer :: nDGridI,nb_bie_read

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='sub_MatOp_HADA'
!-----------------------------------------------------------

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub,para_H%n_Op
        write(out_unit,*) 'nb_act1,nb_var',para_H%mole%nb_act1,               &
                                   para_H%mole%nb_var
        write(out_unit,*) 'nb_ba,nb_bie',para_H%nb_ba,para_H%nb_bie
        write(out_unit,*) 'max_nb_ba_ON_HAC',para_H%para_AllBasis%basis_ext2n%max_nb_ba_ON_HAC
        write(out_unit,*) 'max_ene_ON_HAC  ',para_H%para_AllBasis%basis_ext2n%max_ene_ON_HAC
        CALL write_param_Op(para_H)
      END IF
!-----------------------------------------------------------

!     ----------------------------------------------------------------
!     - test on nb_bi : nb_bi>0
      IF (para_H%nb_bi <= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_bi MUST be >0'
        write(out_unit,*) ' nb_bi',para_H%nb_bi
        write(out_unit,*)
        CALL write_param_Op(para_H)
        STOP
      END IF
!     ----------------------------------------------------------------

!-------------------------------------------------------------
!-     memories allocation: td0b, Opd0bWrho mat1, mat2, mat3
!-------------------------------------------------------------
      nb_ba   = para_H%nb_ba
      nb_bie  = para_H%nb_bie
      type_Op = para_H%para_ReadOp%Type_HamilOp ! H
      IF (type_Op /= 1) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '    Type_HamilOp',type_Op
        write(out_unit,*) '    Type_HamilOp MUST be equal to 1 for HADA or cHAC'
        write(out_unit,*) '    CHECK your data!!'
        STOP
      END IF

      memory = nb_ba**2 * nb_bie
      CALL alloc_NParray(H_HADA, [nb_ba,nb_ba,nb_bie],'H_HADA',name_sub)
      H_HADA(:,:,:) = ZERO

!-------------------------------------------------------------
!-     memories allocation: td0b, Opd0bWrho matRV, mat2, mat3
!-------------------------------------------------------------
      ! selected the optimal value of kmem, as function of max_mem and mem_tot
      kmem = min(10,para_H%nb_qa)
      !kmem = para_H%nb_qa
      allocate(d0MatOpd0bWrho(kmem,nb_ba),stat=err_mem)
      memory = kmem*nb_ba
      CALL error_memo_allo(err_mem,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')
      DO i=1,nb_ba
      DO k=1,kmem
        CALL Init_d0MatOp(d0MatOpd0bWrho(k,i),type_Op,0,nb_bie,iQact=0,JRot=0,cplx=para_H%cplx)
      END DO
      END DO

      nb_blocks = para_H%nb_qa/kmem
      IF (mod(para_H%nb_qa,kmem) /= 0) nb_blocks = nb_blocks + 1
      IF (print_level>0) write(out_unit,*) 'number of blocks',nb_blocks

      CALL alloc_NParray(td0b,[nb_ba,kmem],'td0b',name_sub)

      CALL Init_d0MatOp(d0MatOp,para_H%param_TypeOp,nb_bie)
      !--------------------------------------------------------


!-------------------------------------------------------------
!     - Real part of the Operator ---------------------------
!-------------------------------------------------------------

!     --------------------------------------------------------
!     - loop of kmem block -----------------------------------
      nreste = para_H%nb_qa
      nDGridI = 0

 98   CONTINUE

      IF (print_level>-1) write(out_unit,'(a)',ADVANCE='no') 'Block of ReMatOp(:,:) (%): '
      flush(out_unit)

      nplus = min(kmem,nreste)
      DO k=1,nplus

        nDGridI = nDGridI + 1
        IF (mod(nDGridI,max(1,int(para_H%nb_qa/10))) == 0 .AND.         &
            print_level>-1) THEN
           write(out_unit,'(a,i3)',ADVANCE='no') ' -',                 &
                      int(real(nDGridI,kind=Rkind)*HUNDRED/para_H%nb_qa)
           flush(out_unit)
        END IF

        CALL sub_reading_Op(nDGridI,para_H%nb_qa,d0MatOp,para_H%n_Op,   &
                           Qdyn,para_H%mole%nb_var,para_H%mole%nb_act1, &
                                Qact,WnD,para_H%file_grid)

        iterm_Op = d0MatOp%derive_term_TO_iterm(0,0)
        DO i1_h=1,para_H%nb_bie
          para_H%Hmin = min(para_H%Hmin,d0MatOp%ReVal(i1_h,i1_h,iterm_Op))
          para_H%Hmax = max(para_H%Hmax,d0MatOp%ReVal(i1_h,i1_h,iterm_Op))
        END DO

        CALL calc_td0b_OpRVd0bW(nDGridI,k,td0b,d0MatOpd0bWrho,          &
                                WnD,kmem,d0MatOp,para_H,                &
                                para_H%BasisnD)

      END DO

      DO iterm_Op=1,d0MatOpd0bWrho(1,1)%nb_term

        !================================================================
        ! loop on i1_h
        !================================================================
        DO i1_h=1,nb_bie
          DO ib2=1,para_H%nb_ba
          DO ib1=1,para_H%nb_ba
            DO k=1,nplus
             H_HADA(ib2,ib1,i1_h) = H_HADA(ib2,ib1,i1_h) +              &
                td0b(ib2,k) *  d0MatOpd0bWrho(k,ib1)%ReVal(i1_h,i1_h,iterm_Op)
            END DO
          END DO
          END DO
        END DO
      END DO

      nreste = nreste - nplus
      IF (nreste > 0) GOTO 98
!     - END of loop of kmem block ----------------------------
!     --------------------------------------------------------
      IF (print_level>-1) write(out_unit,'(a)',ADVANCE='yes') ' - End'
      flush(out_unit)

      CALL dealloc_NParray(td0b,     'td0b',     name_sub)

      DO i=1,nb_ba
      DO k=1,kmem
        CALL dealloc_d0MatOp(d0MatOpd0bWrho(k,i))
      END DO
      END DO
      deallocate(d0MatOpd0bWrho,stat=err_mem)
      memory = kmem*nb_ba
      CALL error_memo_allo(err_mem,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')

      CALL dealloc_d0MatOp(d0MatOp)


!     --------------------------------------------------------
!     Diagonalization of HADA channels
      DO i1_h=1,nb_bie
        CALL sub_diago_H(H_HADA(:,:,i1_h),                              &
                para_H%para_AllBasis%basis_ext2n%Eneba_ON_HAC(:,i1_h),  &
                para_H%para_AllBasis%basis_ext2n%d0Cba_ON_HAC(:,:,i1_h),&
                nb_ba,para_H%sym_Hamil)
      END DO


      max_ana_save       = para_ana%max_ana
      print_psi_save     = para_ana%print_psi
      para_ana%print_psi = 0



      ! now, we are working on part of para_H : para_H_HADA. We need to set a local para_AllBasis
      para_AllBasis_loc%alloc   = .TRUE.
      para_AllBasis_loc%nb_be   = 1
      para_AllBasis_loc%nb_bi   = 1
      para_AllBasis_loc%BasisnD => para_H%para_AllBasis%BasisnD ! we don't need to copy
      para_AllBasis_loc%Basis2n => para_H%para_AllBasis%Basis2n ! we don't need to copy

      nb_ba  = get_nb_FROM_Basis(para_H%para_AllBasis%BasisnD)
      para_AllBasis_loc%basis_ext2n%contrac_ba_ON_HAC = .TRUE.
      para_AllBasis_loc%basis_ext2n%max_nb_ba_ON_HAC  =                 &
                        para_H%para_AllBasis%basis_ext2n%max_nb_ba_ON_HAC
      CALL alloc_basis_ext2n(para_AllBasis_loc%basis_ext2n,nb_ba,nb_bie=1)

      para_H_HADA%para_AllBasis => para_AllBasis_loc

      CALL param_HTOparam_H_HADA(1,para_H,para_H_HADA)


!     --------------------------------------------------------
!     Analysis of HADA channels
      ! but first the ZPE from all the channels
      ZPE = minval(para_H%para_AllBasis%basis_ext2n%Eneba_ON_HAC)
      CALL Set_ZPE_OF_Op(para_H_HADA,ZPE=ZPE,forced=.TRUE.)

      ! Analysis of the HADA channels
      CALL alloc_NParray(Tab_Psi,[nb_ba],'Tab_Psi',name_sub)
      DO i=1,nb_ba
        CALL init_psi(Tab_psi(i),para_H_HADA,para_H_HADA%cplx)
        CALL alloc_psi(Tab_Psi(i))
      END DO


      DO i1_h=1,nb_bie

        DO i=1,nb_ba
            Tab_Psi(i)%RvecB(:)   = para_H%para_AllBasis%basis_ext2n%d0Cba_ON_HAC(:,i,i1_h)
            Tab_psi(i)%CAvOp      = para_H%para_AllBasis%basis_ext2n%Eneba_ON_HAC(i,i1_h)
            Tab_psi(i)%IndAvOp    = para_H%n_Op  ! it should be 0
        END DO


        para_ana%max_ana = count(                                       &
           (para_H%para_AllBasis%basis_ext2n%Eneba_ON_HAC(:,i1_h) -     &
            para_H_HADA%ZPE) < para_H%para_AllBasis%basis_ext2n%max_ene_ON_HAC)
        para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i1_h) =           &
            min(para_ana%max_ana,para_H%para_AllBasis%basis_ext2n%max_nb_ba_ON_HAC)

        write(out_unit,*) 'Number of level(s) on the HAC channel: ',   &
                i1_h,para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(i1_h)

        CALL sub_analyse(Tab_Psi,nb_ba,para_H_HADA,                     &
                           para_ana,para_intensity,para_AllOp,const_phys)

      END DO
      CALL dealloc_NParray(Tab_Psi,'Tab_Psi',name_sub)


!---------------------------------------------------------------
      IF (sum(para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC) == 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Sum of nb_ba_ON_HAC(:) == 0 !!'
        STOP
      END IF
      write(out_unit,*) ' Sum of levels on HAC channels: ',            &
                     sum(para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC)

      para_ana%print_psi = print_psi_save
      para_ana%max_ana   = max_ana_save

      CALL dealloc_para_Op(para_H_HADA)
      CALL dealloc_NParray(H_HADA,'H_HADA',name_sub)

      ! para_AllBasis_loc cannot be deallocated by dealloc_AllBasis ....
      ! otherwise, BasisnD qnd Basis2n will be deallocated.
      nullify(para_AllBasis_loc%BasisnD)
      nullify(para_AllBasis_loc%Basis2n)
      CALL dealloc_basis_ext2n(para_AllBasis_loc%basis_ext2n)


!---------------------------------------------------------------

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*)
         write(out_unit,*)
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

  END SUBROUTINE sub_MatOp_HADA
!=======================================================================================

END MODULE mod_Auto_Basis
