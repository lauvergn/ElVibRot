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
MODULE mod_Set_paraRPH
IMPLICIT NONE

PRIVATE
PUBLIC Set_paraPRH
CONTAINS

      SUBROUTINE Set_paraPRH(mole,para_Tnum,BasisnD)
      USE EVR_system_m
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (basis),     intent(in)    :: BasisnD
      TYPE (CoordType), intent(inout) :: mole
      TYPE (Tnum),      intent(in)    :: para_Tnum


!------ working variables ---------------------------------
      integer :: ib,nq_part,iq,iq_list,nb_act1_RPH,nb_inact21_RPH,it,comp,iq_list_small


      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var)

      logical :: Find_in_List,tab_skip_transfo(mole%nb_Qtransfo),RPHCoord_IN_OneBasis



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_paraPRH'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (mole%tab_Qtransfo(mole%itRPH)%skip_transfo) RETURN


      nb_act1_RPH    = mole%RPHTransfo%nb_act1
      nb_inact21_RPH = mole%RPHTransfo%nb_inact21

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub

        write(out_unit,*) 'nb_qa',get_nq_FROM_basis(BasisnD)

        !CALL RecWrite_basis(BasisnD)

        write(out_unit,*) 'nb_act1_RPH',nb_act1_RPH
        write(out_unit,*) 'nb_inact21_RPH',nb_inact21_RPH

        CALL Write_RPHTransfo(mole%RPHTransfo)

        flush(out_unit)
      END IF

      DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
        tab_skip_transfo(it) = mole%tab_Qtransfo(it)%skip_transfo
        mole%tab_Qtransfo(it)%skip_transfo = .TRUE.
      END DO

      ! for tab_RPHpara_AT_Qact1
      IF (.NOT. associated(mole%RPHTransfo%tab_RPHpara_AT_Qact1)) THEN
        IF (debug) write(out_unit,*) ' tab_RPHpara_AT_Qact1'
        CALL alloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,[0],      &
                        'mole%RPHTransfo%tab_RPHpara_AT_Qact1',name_sub,[0])

        CALL get_Qact0(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
        CALL Set_RPHpara_AT_Qact1(mole%RPHTransfo%tab_RPHpara_AT_Qact1(0), &
                                  Qact,para_Tnum,mole)
        mole%RPHTransfo%init_Qref = .TRUE.

        write(out_unit,*) ' Frequencies, normal modes at the reference geometry'

        write(out_unit,11) Qact(1:nb_act1_RPH), &
               mole%RPHTransfo%tab_RPHpara_AT_Qact1(0)%dnEHess%d0(:)*get_Conv_au_TO_unit('E','cm-1')
 11     format(' frequencies : ',30f10.4)
        write(out_unit,*) 'dnQopt'
        CALL Write_dnVec(mole%RPHTransfo%tab_RPHpara_AT_Qact1(0)%dnQopt)
        write(out_unit,*) 'dnC_inv'
        CALL Write_dnMat(mole%RPHTransfo%tab_RPHpara_AT_Qact1(0)%dnC_inv)
        flush(out_unit)

      END IF

      ! Check if the nb_act1_RPH coordinates belong to one basis set (primitive ?)
      ! 1) RPHTransfo MUST be the 2d transformation after the active one.
      !write(out_unit,*) 'asso RPH, itRPH,nb_Qtransfo',associated(mole%RPHTransfo),mole%itRPH,mole%nb_Qtransfo
      RPHCoord_IN_OneBasis = associated(mole%RPHTransfo) .AND.          &
                                      (mole%itRPH == mole%nb_Qtransfo-1)

      RPHCoord_IN_OneBasis = RPHCoord_IN_OneBasis .AND.                 &
        (count(mole%RPHTransfo%list_act_OF_Qdyn(1:nb_act1_RPH) == 1) == nb_act1_RPH)
      !write(out_unit,*) 'list_act_OF_Qdyn',mole%RPHTransfo%list_act_OF_Qdyn


      ! 2) basis functions of BasisnD are defined as a product (BasisnD%nb_basis > 0)
      !  => if true, RPHCoord_IN_OneBasis CAN be true
      RPHCoord_IN_OneBasis = RPHCoord_IN_OneBasis .AND. (BasisnD%nb_basis > 0)

      ! 3) Check nb_act1_RPH coordinates belong to one primitive basis set
      IF (RPHCoord_IN_OneBasis) THEN
        DO ib=1,BasisnD%nb_basis
          !write(out_unit,*) 'ib,iQdyn',ib,':',BasisnD%tab_Pbasis(ib)%Pbasis%iQdyn(:)
          IF (BasisnD%tab_Pbasis(ib)%Pbasis%ndim == nb_act1_RPH) THEN
            IF (all(BasisnD%tab_Pbasis(ib)%Pbasis%iQdyn ==              &
                  mole%RPHTransfo%list_QactTOQdyn(1:nb_act1_RPH)) ) EXIT
          END IF
        END DO
        RPHCoord_IN_OneBasis = RPHCoord_IN_OneBasis .AND. (ib <= BasisnD%nb_basis)
      END IF

      write(out_unit,*) 'RPHCoord_IN_OneBasis',RPHCoord_IN_OneBasis
      flush(out_unit)

      IF (RPHCoord_IN_OneBasis) THEN
        CALL Set_paraPRH_ONEBasis(mole,para_Tnum,BasisnD,ib)
      ELSE
        CALL Set_paraPRH_gene(mole,para_Tnum,BasisnD)
      END IF

      flush(out_unit)
      mole%RPHTransfo%init = .TRUE.

      DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
        mole%tab_Qtransfo(it)%skip_transfo = tab_skip_transfo(it)
      END DO
!     -------------------------------------------------------
      IF (debug) THEN
        CALL Write_RPHTransfo(mole%RPHTransfo)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------
      END SUBROUTINE Set_paraPRH
      SUBROUTINE Set_paraPRH_OneBasis(mole,para_Tnum,BasisnD,ib)
      USE EVR_system_m
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (basis)     :: BasisnD
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum
      integer, intent(in) :: ib ! index of the basis in tab_Pbasis(:) or tab_basisPrimSG(:,:)


!------ working variables ---------------------------------
      integer :: L,iq,iq_list,nb_act1_RPH,nb_inact21_RPH,it,comp,iq_list_small,nq

      real (kind=Rkind), allocatable :: Qact1_fromBasisnD(:)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var)
      real (kind=Rkind), allocatable :: List_Qact1(:,:),List_tmp_Qact1(:,:)

      logical :: Find_in_List,iqLoop_end
      TYPE (Type_RPHpara_AT_Qact1) :: Type_RPHpara_AT_Qref


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_paraPRH_OneBasis'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

      nb_act1_RPH    = mole%RPHTransfo%nb_act1
      nb_inact21_RPH = mole%RPHTransfo%nb_inact21

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'ib',ib

        !CALL RecWrite_basis(BasisnD)

        write(out_unit,*) 'nb_act1_RPH',nb_act1_RPH
        write(out_unit,*) 'nb_inact21_RPH',nb_inact21_RPH

        CALL Write_RPHTransfo(mole%RPHTransfo)

      END IF

      ! for tab_RPHpara_AT_Qact1
      !----------------------------------------------------------------
      !--- First the number of grid points ----------------------------
      CALL time_perso('Grid RPH')
      write(out_unit,*) 'Grid RPH'

      CALL alloc_NParray(Qact1_fromBasisnD,[nb_act1_RPH],'Qact1_fromBasisnD',name_sub)

      iqLoop_end = .FALSE.
      iq = 1
      L  = 0
      SELECT CASE (BasisnD%SparseGrid_type)
      CASE (0) ! normal direct product basis/grid
        nq = get_nq_FROM_basis(BasisnD%tab_Pbasis(ib)%Pbasis) ! it's used one for SparseGrid_type=0
      CASE (1) ! First Smolyak
           STOP ' SG1 not yet'
      CASE (2,4) ! First Smolyak
        nq = get_nq_FROM_basis(BasisnD%tab_basisPrimSG(L,ib))
      CASE Default
         STOP ' no default'
      END SELECT
      write(out_unit,*) 'L,ib,nq',L,ib,nq ; flush(out_unit)

      DO

        SELECT CASE (BasisnD%SparseGrid_type)
        CASE (0) ! normal direct product basis/grid
          Qact1_fromBasisnD = BasisnD%tab_Pbasis(ib)%Pbasis%x(:,iq)
          iq = iq + 1
          iqLoop_end = (iq > nq)
        CASE (1) ! First Smolyak
           STOP ' SG1 not yet'
        CASE (2,4) ! Second and fourth Smolyak

          Qact1_fromBasisnD = BasisnD%tab_basisPrimSG(L,ib)%x(:,iq)
          iq = iq + 1

          IF (iq > nq) THEN
            iq = 1
            L  = L + 1
            IF (L <= BasisnD%L_SparseGrid)                              &
                   nq = get_nq_FROM_basis(BasisnD%tab_basisPrimSG(L,ib))
          END IF

          iqLoop_end = (L > BasisnD%L_SparseGrid)

        CASE Default
           STOP ' no default'
        END SELECT

        !write(out_unit,*) 'L,iq,Qact1_fromBasisnD',L,iq-1,':',Qact1_fromBasisnD


        IF (.NOT. allocated(List_Qact1)) THEN ! first point
          CALL alloc_NParray(List_Qact1,[nb_act1_RPH, 1],'List_tmp_Qact1',name_sub)
          List_Qact1(:,1) = Qact1_fromBasisnD(:)
        ELSE

          Find_in_List = .FALSE.
          iq_list_small = 0
          DO iq_list=1,size(List_Qact1,dim=2)
            comp = compar_GridPoint(List_Qact1(:,iq_list),Qact1_fromBasisnD,nb_act1_RPH)
            IF (comp == -1) iq_list_small = iq_list
            Find_in_List = (comp == 0)

            IF (Find_in_List) EXIT
          END DO

          IF (.NOT. Find_in_List) THEN ! add the new point in the list

            CALL alloc_NParray(List_tmp_Qact1,[nb_act1_RPH, size(List_Qact1,dim=2)+1],'List_tmp_Qact1',name_sub)

            ! find the position to add the point
            IF (iq_list_small == 0) THEN ! add in the first point
              List_tmp_Qact1(:,1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,2:size(List_Qact1,dim=2)+1) = List_Qact1(:,:)
            ELSE IF (iq_list_small == size(List_Qact1,dim=2)) THEN ! add in the last point
              List_tmp_Qact1(:,1:iq_list_small)   = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
            ELSE
              List_tmp_Qact1(:,1:iq_list_small)   = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,  iq_list_small+2:size(List_Qact1,dim=2)+1) = List_Qact1(:,iq_list_small+1:size(List_Qact1,dim=2))
            END IF
            CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)

            CALL alloc_NParray(List_Qact1,shape(List_tmp_Qact1),'List_Qact1',name_sub)
            List_Qact1(:,:) = List_tmp_Qact1
            CALL dealloc_NParray(List_tmp_Qact1,'List_tmp_Qact1',name_sub)

          END IF


        END IF

        IF (iqLoop_end) EXIT


      END DO
      CALL time_perso('Grid RPH')
      write(out_unit,*) 'nb of Qact1 grid points',size(List_Qact1,dim=2)
      flush(out_unit)

      !----------------------------------------------------------------
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      !---- allocation of tab_RPHpara_AT_Qact1 ------------------------
      IF (associated(mole%RPHTransfo%tab_RPHpara_AT_Qact1)) THEN
        CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(                    &
                               mole%RPHTransfo%tab_RPHpara_AT_Qact1(0), &
                               Type_RPHpara_AT_Qref)
        CALL dealloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,        &
                          'mole%RPHTransfo%tab_RPHpara_AT_Qact1',name_sub)
      END IF
      mole%RPHTransfo%nb_Qa = size(List_Qact1,dim=2)
      CALL alloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,[mole%RPHTransfo%nb_Qa],&
                      'mole%RPHTransfo%tab_RPHpara_AT_Qact1',name_sub,[0])
      IF (Type_RPHpara_AT_Qref%init_done > 0) THEN
        CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(Type_RPHpara_AT_Qref,&
                               mole%RPHTransfo%tab_RPHpara_AT_Qact1(0))
        !mole%RPHTransfo%tab_RPHpara_AT_Qact1(0) = Type_RPHpara_AT_Qref
        CALL dealloc_RPHpara_AT_Qact1(Type_RPHpara_AT_Qref)
      END IF
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !---- Multidimensional loop -------------------------------------
        DO iq_list=1,size(List_Qact1,dim=2)

          Qact(:) = ZERO ! initialization to zero, otherwise some values are not initialized (later they will)
          Qact(1:nb_act1_RPH) = List_Qact1(:,iq_list)

          CALL Adding_InactiveCoord_TO_Qact(Qact,mole%ActiveTransfo) ! add rigid, flexible coordinates
          write(out_unit,*) 'new RPH point',iq_list
          !write(out_unit,*) 'new RPH point',Qact(:)

          flush(out_unit)

          CALL Set_RPHpara_AT_Qact1(                                    &
                          mole%RPHTransfo%tab_RPHpara_AT_Qact1(iq_list),&
                                                    Qact,para_Tnum,mole)

        END DO

      write(out_unit,*) 'END Grid RPH'
      flush(out_unit)

      CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)
      CALL dealloc_NParray(Qact1_fromBasisnD,'Qact1_fromBasisnD',name_sub)

!     -------------------------------------------------------
      IF (debug) THEN
        CALL Write_RPHTransfo(mole%RPHTransfo)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------
      END SUBROUTINE Set_paraPRH_OneBasis
      SUBROUTINE Set_paraPRH_gene(mole,para_Tnum,BasisnD)
      USE EVR_system_m
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      IMPLICIT NONE

!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (basis)     :: BasisnD
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum


!------ working variables ---------------------------------
      integer :: nq_part,iq,iq_list,nb_act1_RPH,nb_inact21_RPH,it,comp,iq_list_small

      real (kind=Rkind), allocatable :: Qact1_fromBasisnD(:)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var),auTOcm_inv
      real (kind=Rkind), allocatable :: List_Qact1(:,:),List_tmp_Qact1(:,:)

      logical :: Find_in_List

      TYPE (OldParam) :: OldPara ! to iq, ...



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_paraPRH_gene'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

      nb_act1_RPH    = mole%RPHTransfo%nb_act1
      nb_inact21_RPH = mole%RPHTransfo%nb_inact21

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub

        write(out_unit,*) 'nb_qa',get_nq_FROM_basis(BasisnD)

        !CALL RecWrite_basis(BasisnD)

        write(out_unit,*) 'nb_act1_RPH',nb_act1_RPH
        write(out_unit,*) 'nb_inact21_RPH',nb_inact21_RPH

        CALL Write_RPHTransfo(mole%RPHTransfo)

      END IF

      ! for tab_RPHpara_AT_Qact1
      !----------------------------------------------------------------
      !--- First the number of grid points ----------------------------
      CALL time_perso('Grid RPH')
      write(out_unit,*) 'Grid RPH'

      CALL alloc_NParray(Qact1_fromBasisnD,[nb_act1_RPH],'Qact1_fromBasisnD',name_sub)

      nq_part = get_nq_FROM_basis(BasisnD)/100
      DO iq=1,get_nq_FROM_basis(BasisnD)

        IF (debug .AND. (mod(iq,nq_part)==0)) THEN
          write(out_unit,*) 'iq,nq',iq,get_nq_FROM_basis(BasisnD)
          flush(out_unit)
        END IF

        Qact(:) = ZERO
        CALL Rec_Qact(Qact,BasisnD,iq,mole,OldPara)
        !write(out_unit,*) 'iq,size(List_Qact1,dim=2),Qact',iq,size(List_Qact1,dim=2),Qact
        !flush(out_unit)
        CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
        !write(out_unit,*) 'Qdyn',Qdyn
        !flush(out_unit)

        Qact1_fromBasisnD(:) = Qdyn(mole%RPHTransfo%list_QactTOQdyn(1:nb_act1_RPH))
        !write(out_unit,*) 'Qact1_fromBasisnD',Qact1_fromBasisnD


        IF (.NOT. allocated(List_Qact1)) THEN ! first point
          CALL alloc_NParray(List_Qact1,[nb_act1_RPH, 1],'List_tmp_Qact1',name_sub)
          List_Qact1(:,1) = Qact1_fromBasisnD(:)
        ELSE

          Find_in_List = .FALSE.
          iq_list_small = 0
          DO iq_list=1,size(List_Qact1,dim=2)
            comp = compar_GridPoint(List_Qact1(:,iq_list),Qact1_fromBasisnD,nb_act1_RPH)
            IF (comp == -1) iq_list_small = iq_list
            Find_in_List = (comp == 0)
            IF (Find_in_List) EXIT
          END DO

          IF (.NOT. Find_in_List) THEN ! add the new point in the list

            CALL alloc_NParray(List_tmp_Qact1,[nb_act1_RPH, size(List_Qact1,dim=2)+1],'List_tmp_Qact1',name_sub)

            ! find the position to add the point
            IF (iq_list_small == 0) THEN ! add in the first point
              List_tmp_Qact1(:,1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,2:size(List_Qact1,dim=2)+1) = List_Qact1(:,:)
            ELSE IF (iq_list_small == size(List_Qact1,dim=2)) THEN ! add in the last point
              List_tmp_Qact1(:,1:iq_list_small) = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
            ELSE
              List_tmp_Qact1(:,1:iq_list_small) = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,  iq_list_small+2:size(List_Qact1,dim=2)+1) = List_Qact1(:,iq_list_small+1:size(List_Qact1,dim=2))
            END IF
            CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)

            CALL alloc_NParray(List_Qact1,shape(List_tmp_Qact1),'List_Qact1',name_sub)
            List_Qact1 = List_tmp_Qact1
            CALL dealloc_NParray(List_tmp_Qact1,'List_tmp_Qact1',name_sub)

          END IF


        END IF

      END DO
      CALL time_perso('Grid RPH')
      write(out_unit,*) 'nb of Qact1 grid points',size(List_Qact1,dim=2)
      flush(out_unit)
      !----------------------------------------------------------------
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      !---- allocation of tab_RPHpara_AT_Qact1 ------------------------
      mole%RPHTransfo%nb_Qa = size(List_Qact1,dim=2)
      CALL alloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,[mole%RPHTransfo%nb_Qa],&
                      'mole%RPHTransfo%tab_RPHpara_AT_Qact1',name_sub)
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !---- Multidimensional loop -------------------------------------
        DO iq_list=1,size(List_Qact1,dim=2)

          Qact(:)             = ZERO
          Qact(1:nb_act1_RPH) = List_Qact1(:,iq_list)
          CALL Adding_InactiveCoord_TO_Qact(Qact,mole%ActiveTransfo) ! add rigid, flexible coordinates

          write(out_unit,*) 'new RPH point',iq_list
          flush(out_unit)

          CALL Set_RPHpara_AT_Qact1(                                    &
                          mole%RPHTransfo%tab_RPHpara_AT_Qact1(iq_list),&
                                                    Qact,para_Tnum,mole)

        END DO

      write(out_unit,*) 'END Grid RPH'
      flush(out_unit)

      CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)
      CALL dealloc_NParray(Qact1_fromBasisnD,'Qact1_fromBasisnD',name_sub)

!     -------------------------------------------------------
      IF (debug) THEN
        CALL Write_RPHTransfo(mole%RPHTransfo)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------
      END SUBROUTINE Set_paraPRH_gene

      FUNCTION compar_GridPoint(Q1,Q2,n)
      USE EVR_system_m
      IMPLICIT NONE
        integer :: compar_GridPoint
        integer, intent(in) :: n
        real(kind=Rkind), intent(in) :: Q1(n),Q2(n)
        integer :: i

        DO i=1,n
          IF (abs(Q1(i)-Q2(i)) < ONETENTH**10) THEN
            compar_GridPoint = 0
          ELSE
            IF (Q1(i) < Q2(i)) THEN
              compar_GridPoint = -1
            ELSE
              compar_GridPoint =  1
            END IF
            EXIT
          END IF
        END DO
      END FUNCTION compar_GridPoint

END MODULE mod_Set_paraRPH
