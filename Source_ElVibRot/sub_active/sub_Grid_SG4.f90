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
      SUBROUTINE Make_Grid_SG4(para_AllOp)
      USE mod_system
!$    USE omp_lib, only : OMP_GET_THREAD_NUM
      USE mod_nDindex
      USE mod_Op
      IMPLICIT NONE

!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
  TYPE (param_AllOp), intent(inout) :: para_AllOp


!------ working variables ---------------------------------
  TYPE (TypeRVec)    :: TabAllOp_Rvec(para_AllOp%nb_Op)
  logical, allocatable           :: Grid_cte(:)

  integer,            allocatable       :: tab_l(:)
  integer                               :: iG,ith,nb_thread,iOp,iterm,i_e,id1,id2


  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  !----- for debuging --------------------------------------------------
  !integer :: err_mem,memory
  character (len=*), parameter :: name_sub='Make_Grid_SG4'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------
  mole    => para_AllOp%tab_Op(1)%mole
  BasisnD => para_AllOp%tab_Op(1)%BasisnD

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
  END IF
  !-----------------------------------------------------------

  !----- built tables --------------------------------------------------
  !Define the volume element (nrho of Basis => nrho of Tnum
  CALL nrho_Basis_TO_nhro_Tnum(para_AllOp%tab_Op(1)%para_AllBasis,mole)

  DO iOp=1,para_AllOp%nb_Op
    IF (para_AllOp%tab_Op(iOp)%n_op == -1) CYCLE ! for S

    IF (.NOT. para_AllOp%tab_Op(iOp)%alloc_Grid) THEN

      CALL alloc_NParray(Grid_cte,[para_AllOp%tab_Op(iOp)%nb_term],&
                        "Grid_cte",name_sub)
      Grid_cte(:) = .FALSE.

      IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%n_op == 0) THEN
        Grid_cte(:) = .TRUE.  ! KEO
        Grid_cte(1) = .FALSE. ! potential
      END IF
      IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%name_Op(1:3) == "Mu_") THEN
        Grid_cte(:) = .TRUE.  ! KEO, rotation, Mu_xx,Mu_xy...
      END IF
      IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%name_Op(1:5) == "Corr_") THEN
        Grid_cte(:) = .TRUE.  ! KEO, Coriolis
      END IF

      CALL dealloc_NParray(Grid_cte,"Grid_cte",name_sub)
    END IF
  END DO

  !----- Transfert the constant KEO to Mate_cte -----------------
  IF (para_AllOp%tab_Op(1)%para_Tnum%Gcte) THEN
    DO iOp=1,para_AllOp%nb_Op

      IF (para_AllOp%tab_Op(iOp)%n_op == 0) THEN

        DO iterm=1,para_AllOp%tab_Op(iOp)%nb_term

          IF (.NOT. para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Grid_cte) CYCLE

          id1 = para_AllOp%tab_Op(iOp)%derive_termQact(1,iterm)
          id2 = para_AllOp%tab_Op(iOp)%derive_termQact(2,iterm)
          IF (id1 /= 0 .AND. id2 /= 0) THEN ! f2 for G
            IF (id1 == id2) THEN
              DO i_e=1,para_AllOp%tab_Op(iOp)%para_ReadOp%nb_elec
                para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Mat_cte(i_e,i_e) = &
                  -HALF* para_AllOp%tab_Op(iOp)%para_Tnum%Gref(id1,id2)
              END DO
            ELSE
             DO i_e=1,para_AllOp%tab_Op(iOp)%para_ReadOp%nb_elec
                para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Mat_cte(i_e,i_e) = &
                  - para_AllOp%tab_Op(iOp)%para_Tnum%Gref(id1,id2)
              END DO
            END IF
          END IF

        END DO

      END IF
      IF (para_AllOp%tab_Op(iOp)%name_Op(1:3) == "Mu_") THEN
        STOP 'Rot cte not yet!'
      END IF
      IF (para_AllOp%tab_Op(iOp)%name_Op(1:5) == "Corr_") THEN
        STOP 'Corriolis cte not yet!'
      END IF

    END DO
  END IF
  !-----------------------------------------------------
  !-----------------------------------------------------
  !-- Multidimensional loop ----------------------------
  IF (Grid_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = BasisnD%para_SGType2%nb_threads
  END IF
  IF (print_level > 1) write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

  IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
      mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
    write(out_unitp,'(a)',ADVANCE='no') '---'
    flush(out_unitp)
  END IF

  !$OMP parallel                                                &
  !$OMP default(none)                                           &
  !$OMP shared(para_AllOp,BasisnD,print_level,out_unitp,MPI_id) &
  !$OMP private(iG,tab_l,ith)                                   &
  !$OMP num_threads(nb_thread)

  !--------------------------------------------------------------
  !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
  ith = 0
  !$ ith = OMP_GET_THREAD_NUM()
  CALL alloc_NParray(tab_l,[BasisnD%nb_basis],'tab_l',name_sub)
  tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
  !--------------------------------------------------------------

  ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
  DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)

    CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

    !CALL Make_Grid_OF_ONEGDP_SG4(TabAllOp_Rvec,tab_l,iG,para_AllOp)

    IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
        mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
      write(out_unitp,'(a)',ADVANCE='no') '---'
      flush(out_unitp)
    END IF

    !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
  END DO
  CALL dealloc_NParray(tab_l,'tab_l',name_sub)


  !$OMP   END PARALLEL

  IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
    write(out_unitp,'(a)',ADVANCE='yes') '----]'
  END IF
  flush(out_unitp)

  !-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF
  !-------------------------------------------------------

END SUBROUTINE Make_Grid_SG4

