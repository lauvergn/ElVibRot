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
MODULE mod_Davidson
USE mod_Constant
USE mod_MPI
IMPLICIT NONE

PRIVATE
PUBLIC :: sub_propagation_Davidson

CONTAINS

!===============================================================================
      SUBROUTINE sub_propagation_Davidson(psi,Ene,nb_diago,max_diago,   &
                                          para_H,para_Davidson,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,dealloc_psi,norm2_psi, &
                             alloc_NParray,dealloc_NParray,ecri_psi,    &
                             sub_LCpsi_TO_psi,sub_save_psi
      USE mod_Op
      USE mod_propa,       ONLY:param_propa,param_Davidson
      USE mod_propa_MPI,   ONLY:MPI_Bcast_param_Davidson
      USE mod_Davidson_MPI,ONLY:MakeResidual_Davidson_MPI4,             &
                                MakeResidual_Davidson_j_MPI3
      USE mod_ana_psi_MPI, ONLY:share_psi_nodes_MPI
      USE mod_psi_Op_MPI,  ONLY:sub_LCpsi_TO_psi_MPI
      USE mod_MPI_aux
      IMPLICIT NONE

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op)   :: para_H

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa)        :: para_propa
      TYPE (param_Davidson)     :: para_Davidson

      integer                   :: nb_diago,nb_diagoR,max_diago
      TYPE (param_psi)          :: psi(max_diago)



      TYPE (param_psi), allocatable  :: psi0(:)
      TYPE (param_psi), allocatable  :: Hpsi(:)
      TYPE (param_psi)          :: g
      real (kind=Rkind)         :: Ene(max_diago),Ene0(max_diago),Ene1(max_diago)
      real (kind=Rkind)         :: DEne(max_diago),EneRef(max_diago)
      real (kind=Rkind)         :: conv_Ene,rms_Ene,Di,ZPE,DE
      logical                   :: save_WP
      real (kind=Rkind),allocatable :: H(:,:)
      real (kind=Rkind),allocatable :: H_overlap(:,:) !< <psi|Hpsi> overlap
      real (kind=Rkind),allocatable :: S_overlap(:,:) !< <psi|psi> overlap

      real (kind=Rkind),allocatable :: Vec0(:,:)
      real (kind=Rkind),allocatable :: Vec(:,:)

      logical                   :: VecToBeIncluded(max_diago)
      logical                   :: converge(max_diago)
      logical                   :: convergeEne(max_diago),convergeResi(max_diago)
      real (kind=Rkind)         :: non_hermitic
      integer                   :: ierr,ndim,ndim0,iresidu,fresidu,ndim_Vec0,ndim_ini0
      complex (kind=Rkind)      :: CS
      real (kind=Rkind)         :: RS,a,max_Sii,max_Sij
      real (kind=Rkind)         :: tab_norm2g(max_diago)
      character (len=Name_len)  :: info
      real (kind=Rkind)         :: tab_a(max_diago)

      !------ working parameters --------------------------------
      TYPE(File_t)  :: Log_file
      integer           :: iunit

      complex (kind=Rkind) :: Overlap

      integer           :: it,it_all,no,nb_min_Ene
      integer           :: i,k,ii,j,j_ini,j_end,iqa
      integer           :: i1,i2,i3,ib
      integer           :: err
      integer           :: max_ecri
      real (kind=Rkind) :: DeltaT,T      ! time
      real (kind=Rkind) :: DeltaE,Deltapsi,epsi,norm2g,th,Hcv
      real (kind=Rkind) :: E1,RE0,auTOcm_inv,auTOene
      character (len=Name_len) :: WriteUnit

      real (kind=Rkind) :: RealTime,min_Ene
      TYPE (Time_t) :: DavidsonTime

      logical           :: conv
      integer           :: nb_conv_states,nb_unconv_states,nb_added_states

      integer           :: nb_thread
      integer, pointer  :: list_readWP(:)

      logical           :: With_Grid,With_Basis,Hmin_OR_Hmax
      logical           :: exit_Davidson
      logical           :: Forced_sym_Hamil = .FALSE.
      logical           :: sym_Hamil

      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_propagation_Davidson'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      Hmin_OR_Hmax = para_Davidson%Hmax_propa .OR.           &
                     para_Davidson%Hmin_propa
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' propagation Davidson'
        write(out_unitp,*) ' nb_diago',nb_diago
        write(out_unitp,*) ' max_diago',max_diago
        write(out_unitp,*) ' Hmin_OR_Hmax',Hmin_OR_Hmax

        write(out_unitp,*) ' para_Davidson',para_Davidson
        write(out_unitp,*)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
      auTOene    = get_Conv_au_TO_WriteUnit('E',WriteUnit)
      sym_Hamil  = (Forced_sym_Hamil .OR. para_H%sym_Hamil)

      IF(MPI_id==0) THEN
        write(out_unitp,*) 'all_lower_states',para_Davidson%all_lower_states
        write(out_unitp,*) 'lower_states    ',para_Davidson%lower_states
        write(out_unitp,*) 'project_WP0     ',para_Davidson%project_WP0
        write(out_unitp,*) 'NewVec_type     ',para_Davidson%NewVec_type
        write(out_unitp,*) 'sym_Hamil       ',sym_Hamil
      ENDIF


      IF (para_Davidson%Op_Transfo .AND. para_H%para_ReadOp%Op_Transfo) THEN

        DE = para_Davidson%max_ene
        IF(MPI_id==0) THEN
          write(out_unitp,*) 'DE:             ',DE
          write(out_unitp,*) 'degree_Transfo: ',para_H%para_ReadOp%degree_Transfo
          write(out_unitp,*) 'alloc Poly_Transfo:   ',allocated(para_H%para_ReadOp%Poly_Transfo)
          write(out_unitp,*) 'Poly_Transfo:   ',para_H%para_ReadOp%Poly_Transfo
        ENDIF

        para_Davidson%max_ene = para_H%para_ReadOp%Poly_Transfo(0)
        DO i=1,para_H%para_ReadOp%degree_Transfo
          para_Davidson%max_ene = para_Davidson%max_ene + &
            para_H%para_ReadOp%Poly_Transfo(i) * DE**i
        END DO
      END IF
      para_Davidson%save_max_ene = para_Davidson%max_ene * &
                                              para_Davidson%scaled_max_ene

      !------ initialization -------------------------------------
      Log_file%name='Davidson.log'
      IF(MPI_id==0) CALL file_open(Log_file,iunit)

      With_Grid  = para_Davidson%With_Grid
      With_Basis = .NOT. para_Davidson%With_Grid

      !para_propa%file_WP%formatted = .TRUE.
      para_propa%file_WP%name         = para_Davidson%name_file_saveWP
      para_Davidson%formatted_file_WP = para_propa%file_WP%formatted

      CALL alloc_NParray(psi0,shape(psi),'psi0',name_sub)
      CALL alloc_NParray(Hpsi,shape(psi),'Hpsi',name_sub)

      DO i=1,max_diago
        CALL init_psi( psi(i),para_H,para_H%cplx)
        CALL init_psi(Hpsi(i),para_H,para_H%cplx)
      END DO

      IF (With_Grid) THEN
        DO i=1,max_diago
          Hpsi(i)%BasisRep = .FALSE.
          Hpsi(i)%GridRep  = .TRUE.
        END DO
      END IF

      CALL init_psi(g,para_H,para_H%cplx)
      CALL alloc_psi(g,      BasisRep=With_Basis,GridRep=With_Grid)

      IF(keep_MPI) THEN
        ! read guess on master
        CALL ReadWP0_Davidson(psi,psi0,Vec0,nb_diago,max_diago,   &
                              para_Davidson,para_H%cplx)
      ENDIF
      ! save the nb_diago wp
      IF(MPI_id==0) CALL sub_save_psi(psi,nb_diago,para_propa%file_WP)
      write(out_unitp,*) '  sub_save_psi: psi done'
      RealTime = Delta_RealTime(DavidsonTime)

      flush(out_unitp)
      IF(openmpi .AND. MPI_scheme/=1) THEN
        CALL MPI_Bcast_(nb_diago,size1_MPI,root_MPI)
        CALL MPI_Bcast_param_Davidson(para_Davidson)
      ENDIF
      !CALL time_perso('Davidson psi0')
!para_mem%mem_debug = .TRUE.
      !===================================================================
      !===================================================================
      ! LOOP
      !===================================================================
      !===================================================================
      nb_diago        = max(1,nb_diago) ! number of Eigenvalues
      epsi            = para_Davidson%conv_resi
      norm2g          = HUNDRED * epsi
      conv_Ene        = HUNDRED * epsi
      tab_norm2g(:)   = norm2g
      convergeEne(:)  = .FALSE.
      convergeResi(:) = .FALSE.
      converge(:)     = .FALSE.

      it      = 0
      ndim    = nb_diago
      ndim0   = 0
      iresidu = 0
      conv    = .FALSE.
      exit_Davidson=.FALSE.

      CALL alloc_NParray(vec,[ndim,ndim],"vec",name_sub)
      IF (MatOp_omp /= 2) THEN
        nb_thread = 1
      ELSE
        nb_thread = MatOp_maxth
      END IF

      write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread
      write(out_unitp,*) 'Beginning Davidson iteration'
      flush(out_unitp)

      IF(MPI_id==0) THEN
        write(iunit,*) 'Beginning Davidson iteration' ; flush(iunit)
      ENDIF

      !--------------------------------------------------------------------------------
      ! loop for davidson with maximum iter number para_Davidson%max_it
      !--------------------------------------------------------------------------------
      DO it=0,para_Davidson%max_it

        IF(MPI_id==0) THEN
          write(out_unitp,*) '--------------------------------------------------'
          write(iunit,*) 'Davidson iteration',it ; flush(iunit)
        ENDIF

        !CALL time_perso('Beginining it')

        !---------------------------------------------------------------
        !- Hpsi(:) -----------------------------------------------------
        IF (debug) write(out_unitp,*) 'Hpsi(:)',it,ndim,ndim0
        IF (debug) flush(out_unitp)
        CALL sub_MakeHPsi_Davidson(it,psi(1:ndim),Hpsi(1:ndim),Ene,ndim0,para_H,       &
                                   para_Davidson,iunit,H_overlap,S_overlap,save_WP)
        save_WP = .FALSE.
        !CALL time_perso('MakeHPsi done') 

        IF (debug) write(out_unitp,*) 'Hpsi(:) done',it,ndim,ndim0
        IF (debug) flush(out_unitp)
        !- Hpsi(:) -----------------------------------------------------
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        !- build H
        IF (debug) write(out_unitp,*) 'H mat',it,ndim,ndim0
        IF (debug) flush(out_unitp)

        ! H built from psi and Hpsi
        ! MPI: H is now assigned from H_overlap directly on all the threads
        CALL sub_MakeH_Davidson(it,psi(1:ndim),Hpsi(1:ndim),H,H_overlap,S_overlap,     &
                                para_Davidson)
        ndim0 = ndim
        !CALL time_perso('MakeH done')

        ! if symmetric
        CALL sub_hermitic_H(H,ndim,non_hermitic,sym_Hamil)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)
        IF (debug) write(out_unitp,*) 'shape H',shape(H)

        IF (non_hermitic > FOUR*ONETENTH**4) THEN
          write(out_unitp,*) 'WARNING: non_hermitic is BIG'
          write(out_unitp,31) non_hermitic
31          format(' Hamiltonien: ',f16.12,' au')
        ELSE
          write(out_unitp,51) non_hermitic*auTOcm_inv
51          format(' Hamiltonien: ',f16.12,' cm-1')
        END IF

        IF (sym_Hamil) THEN
          epsi = max(para_Davidson%conv_resi,                &
                   TEN**para_Davidson%conv_hermitian *       &
                                               non_hermitic)
        ELSE
          epsi = para_Davidson%conv_resi
        END IF
        !CALL time_perso('MakeH done')

        IF (debug) write(out_unitp,*) 'H mat',it,ndim,ndim0
        flush(out_unitp)
        !- build H
        !----------------------------------------------------------

        !----------------------------------------------------------
        !- diagonalization
        IF (debug) write(out_unitp,*) 'diago',it,ndim,ndim0
        IF (debug) flush(out_unitp)

        CALL dealloc_NParray(vec,"vec",name_sub)
        CALL alloc_NParray(vec,[ndim,ndim],"vec",name_sub)
        Ene(:) = ZERO

        ! write(out_unitp,*) 'ndim',ndim
        ! write(out_unitp,*) 'shape ..',shape(H),shape(Vec),shape(Ene),shape(trav)
        IF (sym_Hamil) THEN
          ! CALL diagonalization(H,Ene(1:ndim),Vec,ndim,3,1,.FALSE.)
          ! consider the MPI of diagonalization
          CALL diagonalization(H,Ene(1:ndim),Vec,ndim,3,1,.True.)
          ! CALL diagonalization(H,Ene(1:ndim),Vec,ndim,2,1,.FALSE.)
        ELSE
          ! CALL diagonalization(H,Ene(1:ndim),Vec,ndim,4,1,.FALSE.)
          CALL diagonalization(H,Ene(1:ndim),Vec,ndim,4,1,.True.)
        END IF
        IF (debug) write(out_unitp,*) 'Ene',Ene(1:ndim)

        IF (it == 0 .OR. (it > 1 .AND.                                    &
              mod(it-1,para_Davidson%num_resetH) == 0) ) THEN
          EneRef(:) = Ene(:)
        END IF

        IF (it == 0) Ene0(:) = Ene(:) + conv_Ene
        !CALL time_perso('diago done')
        IF (debug) write(out_unitp,*) 'diago',it,ndim,ndim0
        IF (debug) flush(out_unitp)
        !- diagonalization
        !----------------------------------------------------------

        !----------------------------------------------------------
        ! Save vec(:) on vec0(:)
        IF (debug) write(out_unitp,*) 'selec',it,ndim,ndim0
        IF (debug) flush(out_unitp)
        ! find the physical min_ene (to avoid holes)
        IF (para_H%para_ReadOp%min_pot == huge(ONE)) THEN
          min_Ene = -huge(ONE)
        ELSE
          min_Ene = para_H%para_ReadOp%min_pot
        END IF
        ZPE = huge(ONE)
        IF (para_H%para_ReadOp%Op_Transfo) THEN
          CALL Set_ZPE_OF_Op(para_H,Ene(1:ndim),forced=.TRUE.)
        ELSE
          CALL Set_ZPE_OF_Op(para_H,Ene(1:ndim),Ene_min=min_Ene,forced=.TRUE.)
        END IF
        ZPE     = para_H%ZPE
        IF (debug) write(out_unitp,*) 'ZPE',it,ZPE
        !write(out_unitp,*) 'min_Ene,min_pot,ZPE',it,min_Ene,para_H%para_ReadOp%min_pot,ZPE
        CALL sub_projec_Davidson(Ene,VecToBeIncluded,nb_diago,min_Ene,  &
                                 para_H%para_ReadOp%min_pot,ZPE,        &
                                 psi,psi0,Vec,Vec0,para_Davidson,it,.TRUE.)
        !CALL time_perso('projec done')
          IF (debug) write(out_unitp,*) 'selec',it,ndim,ndim0
          IF (debug) flush(out_unitp)
          ! CALL Write_Mat(Vec,out_unitp,5)
          ! Save vec(:) on vec0(:)
          !----------------------------------------------------------

          IF (para_Davidson%all_lower_states) THEN
            nb_diago = count((Ene(1:count(VecToBeIncluded))-ZPE) <= para_Davidson%max_Ene)
            VecToBeIncluded = .FALSE.
            VecToBeIncluded(1:nb_diago) = .TRUE.
          END IF

          IF (debug) write(out_unitp,*) 'nb_diago,ZPE,min_Ene',nb_diago,ZPE,min_Ene
          IF (debug) write(out_unitp,*) 'it,Ene(:)',it,Ene(1:ndim)*auTOene

          DEne(1:nb_diago) = Ene(1:nb_diago)-Ene0(1:nb_diago)
          conv_Ene = maxval(abs(DEne(1:nb_diago)))
          write(out_unitp,41) 'convergence (it, norm2g/epsi, conv_Ene): ',&
                                                   it,norm2g/epsi,conv_Ene
          IF (para_Davidson%Hmax_propa) THEN
            write(out_unitp,21) it,ndim,norm2g/epsi,iresidu,              &
                             -Ene(1:nb_diago)*auTOene
          ELSE
            write(out_unitp,21) it,ndim,norm2g/epsi,iresidu,              &
                              Ene(1:nb_diago)*auTOene
          END IF

          DO j=1,nb_diago
            convergeEne(j) = abs(DEne(j)) < para_Davidson%conv_Ene
          END DO
          write(out_unitp,41) 'it Diff Ene (' // trim(WriteUnit) // '):    ',it, &
                           DEne(1:nb_diago) * auTOene
          write(out_unitp,42) 'it convergenceEne(:):  ',it,               &
                                                    convergeEne(1:nb_diago)
          flush(out_unitp)
          IF(MPI_id==0) THEN
            write(iunit,21) it,ndim,norm2g/epsi,iresidu,Ene(1:ndim)*auTOene
            flush(iunit)
          ENDIF
          !----------------------------------------------------------
          !- residual vector ---------------------------
          !-  and convergence --------------------------
          IF (debug) write(out_unitp,*) 'residual',it,ndim,ndim0
          IF (debug) flush(out_unitp)
          norm2g    = -ONE
          fresidu   = 0

        IF(openmpi) THEN
          CALL MakeResidual_Davidson_MPI4(ndim,g,psi,Hpsi,Ene,Vec,conv,converge,       &
                                       VecToBeIncluded,tab_norm2g,norm2g,convergeResi, &
                                       convergeEne,fresidu,iresidu,nb_diago,epsi)
        ELSE
          DO j=1,ndim
            IF (.NOT. converge(j) .AND. VecToBeIncluded(j)) THEN
              CALL MakeResidual_Davidson(j,g,psi,Hpsi,Ene,Vec)

              CALL norm2_psi(g)
              tab_norm2g(j) = sqrt(g%norm2)
              IF (fresidu == 0) fresidu = j
              IF (tab_norm2g(j) > norm2g) THEN
                iresidu = j
                norm2g = tab_norm2g(iresidu)
              END IF

              convergeResi(j) = tab_norm2g(j) < epsi
            END IF
            converge(j) = (convergeEne(j) .AND. convergeResi(j))

          END DO
          conv = all(converge(1:nb_diago))
        ENDIF

!#if(run_MPI)
!        CALL MPI_BCAST(conv,size1_MPI,MPI_LOGICAL,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif

          Ene0(1:nb_diago) = Ene(1:nb_diago)
          write(out_unitp,41) 'it tab_norm2g          ',it,tab_norm2g(1:nb_diago)
          write(out_unitp,42) 'it convergenceResi(:): ',it,convergeResi(1:nb_diago)
41        format(a,i3,100(1x,e9.2))
42        format(a,i3,100(1x,l9))

          !CALL time_perso('residual done')

          IF (debug) write(out_unitp,*) 'residual',it,ndim,ndim0
          IF (debug) flush(out_unitp)
          !- residual vector and convergence ------------------------
          !----------------------------------------------------------

          !----------------------------------------------------------
          !- convergence ? ------------------------------------------
          nb_conv_states   = count( converge(1:nb_diago) )
          nb_unconv_states = nb_diago-nb_conv_states
          write(out_unitp,*) 'it, conv',it,conv
          write(out_unitp,*) 'it, nb   converged state(s)',it,nb_conv_states
          write(out_unitp,*) 'it, nb unconverged state(s)',it,nb_unconv_states
          flush(out_unitp)
          !- convergence ? ------------------------------------------
          !----------------------------------------------------------

          CALL sub_NewVec_Davidson(it,psi,Hpsi,Ene,Ene0,EneRef,Vec,                    &
                                   converge,VecToBeIncluded,nb_diago,max_diago,        &
                                   para_Davidson,fresidu,ndim,                         &
                                   para_H%para_ReadOp%Op_Transfo,                      &
                                   para_H%para_ReadOp%E0_Transfo,S_overlap)

          nb_added_states = ndim-ndim0
          IF (nb_added_states == 0) THEN
            conv = .TRUE.
            write(out_unitp,*) ' WARNING: The conv is forced to true, because nb_added_states=0.'
          END IF
          save_WP = (ndim == max_diago) .OR. conv .OR.                  &
                     it == para_Davidson%max_it .OR.                    &
                     (it > 0 .AND. mod(it,para_Davidson%num_resetH) == 0)
          Save_WP = Save_WP .AND. .NOT. Hmin_OR_Hmax
          !- new vectors --------------------------------------------
          !----------------------------------------------------------

          IF(MPI_id==0) CALL exit_Davidson_external(exit_Davidson,save_WP,it)
          IF(openmpi) THEN
            CALL MPI_BCAST_(exit_Davidson,size1_MPI,root_MPI)
            CALL MPI_BCAST_(save_WP,size1_MPI,root_MPI)
          ENDIF 

          write(out_unitp,*) 'ndim,max_diago',ndim,max_diago
          write(out_unitp,*) 'save_WP,conv',save_WP,conv
          IF (debug) THEN
            write(out_unitp,*) 'it,max_it',it,para_Davidson%max_it, &
                                     (it == para_Davidson%max_it)
            write(out_unitp,*) 'mod(it,para_Davidson%num_resetH)',  &
                                mod(it,para_Davidson%num_resetH)
          END IF
          flush(out_unitp)

          !----------------------------------------------------------
          !- save psi(:) on file
          IF (save_WP) THEN
            IF (para_H%para_ReadOp%Op_Transfo) THEN
              CALL Set_ZPE_OF_Op(para_H,Ene(1:count(VecToBeIncluded)),forced=.TRUE.)
            ELSE
              CALL Set_ZPE_OF_Op(para_H,Ene(1:count(VecToBeIncluded)),Ene_min=min_Ene,forced=.TRUE.)
            END IF
            ZPE = para_H%ZPE
            CALL sub_projec_Davidson(Ene,VecToBeIncluded,nb_diago,           &
                                     min_Ene,para_H%para_ReadOp%min_pot,ZPE, &
                                     psi,psi0,Vec,Vec0,para_Davidson,it,.TRUE.)

            !CALL time_perso('projec done')

            write(out_unitp,*) 'save psi(:)',it,ndim,ndim0
            flush(out_unitp)

            !- check the orthogonality ------------------------
            CALL sub_MakeS_Davidson(it,psi(1:ndim),With_Grid,debug,S_overlap)
            !CALL time_perso('MakeS done')

            IF(openmpi) THEN
              CALL sub_LCpsi_TO_psi_MPI(psi,Vec,ndim0,nb_diago)
            ELSE
              CALL sub_LCpsi_TO_psi(psi,Vec,ndim0,nb_diago)
            ENDIF
            write(out_unitp,*) '  sub_LCpsi_TO_psi: psi done',ndim0,nb_diago
            flush(out_unitp)

            IF(keep_MPI) THEN
            ! move the new vectors (nb_added_states), after the nb_diago ones
            DO i=1,nb_added_states
              psi(nb_diago+i) = psi(ndim0+i)
            END DO
            write(out_unitp,*) '  move psi done'
            flush(out_unitp)
            ENDIF ! for keep_MPI

            ! deallocation
            DO i=nb_diago+nb_added_states+1,ndim
              !write(out_unitp,*) 'dealloc psi(i), Hpsi(i)',i
              CALL dealloc_psi(psi(i))
            END DO
            write(out_unitp,*) '  deallocation psi done'
            flush(out_unitp)

            ! save the nb_diago wp
            IF(MPI_id==0) CALL sub_save_psi(psi,nb_diago,para_propa%file_WP)
            write(out_unitp,*) '  sub_save_psi: psi done'
            flush(out_unitp)

            IF(openmpi) THEN
              CALL sub_LCpsi_TO_psi_MPI(Hpsi,Vec,ndim0,nb_diago)
            ELSE
              CALL sub_LCpsi_TO_psi(Hpsi,Vec,ndim0,nb_diago)
            ENDIF
            write(out_unitp,*) '  sub_LCpsi_TO_psi: Hpsi done',ndim0,nb_diago
            flush(out_unitp)

            DO i=nb_diago+1,ndim
              !write(out_unitp,*) 'dealloc psi(i), Hpsi(i)',i
              CALL dealloc_psi(Hpsi(i))
            END DO

            write(out_unitp,*) '  deallocation Hpsi done'
            flush(out_unitp)

            Vec = Identity_Mat(ndim0)

            ndim0 = nb_diago
            ndim  = nb_diago+nb_added_states

            IF (allocated(Vec0))  THEN
              CALL dealloc_NParray(Vec0,"Vec0",name_sub)
            END IF
            CALL alloc_NParray(Vec0,[ndim0,ndim0],"Vec0",name_sub)
            Vec0 = Identity_Mat(ndim0)

            IF (allocated(H))  THEN
              CALL dealloc_NParray(H,"H",name_sub)
            END IF
          END IF ! for save_WP
          !- save psi(:) on file
          !----------------------------------------------------------

          RealTime = Delta_RealTime(DavidsonTime)
          IF (RealTime < TEN) THEN
            write(out_unitp,'(a,i0,a,i0)') 'At Davidson iteration: ',it,', Delta Real time (ms): ',int(10**3*RealTime)
          ELSE
            write(out_unitp,'(a,i0,a,i0)') 'At Davidson iteration: ',it,', Delta Real time (s): ',int(RealTime)
          END IF
          flush(out_unitp)

        IF (conv) EXIT

        !IF(openmpi .AND. exit_Davidson) EXIT
        IF(exit_Davidson) EXIT

!#if(run_MPI)
!        ! boardcast new ndim value
!        CALL MPI_BCAST(ndim,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!        CALL MPI_BCAST(ndim0,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
        ! boardcast new ndim value
        IF(openmpi .AND. MPI_scheme/=1) THEN
          CALL MPI_BCAST_(ndim ,size1_MPI,root_MPI)
          CALL MPI_BCAST_(ndim0,size1_MPI,root_MPI)
        ENDIF

      END DO ! for it=0,para_Davidson%max_it
      write(out_unitp,*) '--------------------------------------------------'

      !===================================================================
      !===================================================================
      ! END LOOP
      !===================================================================
      !===================================================================
      IF(MPI_id==0) THEN
        write(iunit,*) 'End Davidson ' ; flush(iunit)
      ENDIF

 21   format(' Davidson: ',2(i5,1x),e10.3,i5,1x,50(1x,f18.4))

      write(out_unitp,*)
      write(out_unitp,*) '==========================================='
      write(out_unitp,*) '==========================================='
      IF (conv) THEN
        IF(MPI_id==0) write(out_unitp,*) ' Davidson has converged after ',it,' iterations'
      ELSE
        write(out_unitp,*) ' WARNNING: Davidson has NOT converged after',it,' iterations'
      END IF
      flush(out_unitp)

      IF (para_H%para_ReadOp%Op_Transfo) THEN
        ! The energies have to be recalculated without T(Op)
        para_H%para_ReadOp%Op_Transfo = .FALSE.
        CALL sub_MakeHPsi_Davidson(it,psi(1:nb_diago),Hpsi(1:nb_diago),Ene,0,para_H,   &
                                   para_Davidson,iunit,H_overlap,S_overlap,.FALSE.)
        DO j=1,nb_diago
          g = Hpsi(j) - psi(j) * Ene(j)
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          write(out_unitp,*) 'lev:',j,Ene(j)*auTOene,tab_norm2g(j)
        END DO
      END IF
      IF(MPI_id==0) CALL file_close(Log_file)

      IF (.NOT. Hmin_OR_Hmax) THEN
        CALL Set_ZPE_OF_Op(para_H,Ene(1:nb_diago),                      &
                                          Ene_min=min_Ene,forced=.TRUE.)
        ZPE = para_H%ZPE
      END IF

      DO j=1,nb_diago
        psi(j)%CAvOp    = Ene(j)
        psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
        psi(j)%convAvOp = convergeEne(j) .AND. convergeResi(j)
      END DO

      IF(MPI_id==0) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'Number of Hamiltonian operations (H I psi >)',para_H%nb_OpPsi
        write(out_unitp,*)
      ENDIF

      IF (para_Davidson%Hmax_propa) THEN
        para_propa%Hmax = -Ene(1)
        para_H%Hmax     = -Ene(1)
        write(out_unitp,*) 'Hmax (ua)  : ',para_propa%Hmax
        write(out_unitp,*) 'Hmax (cm-1): ',para_propa%Hmax*auTOene
      END IF
      IF (para_Davidson%Hmin_propa) THEN
        para_propa%Hmin = Ene(1)
        para_H%Hmin     = Ene(1)
        write(out_unitp,*) 'Hmin (ua)  : ',para_propa%Hmin
        write(out_unitp,*) 'Hmin (cm-1): ',para_propa%Hmin*auTOene
      END IF

      write(out_unitp,*)
      write(out_unitp,*) '==========================================='
      write(out_unitp,*) '==========================================='
      flush(out_unitp)

      !----------------------------------------------------------
      IF (allocated(Vec))  THEN
        CALL dealloc_NParray(Vec,"Vec",name_sub)
      END IF
      IF (allocated(H))  THEN
        CALL dealloc_NParray(H,"H",name_sub)
      END IF
      IF (allocated(Vec0))  THEN
        CALL dealloc_NParray(Vec0,"Vec0",name_sub)
      END IF

      IF(keep_MPI) THEN
        DO i=1,size(Hpsi,dim=1)
          CALL dealloc_psi(Hpsi(i),delete_all=.TRUE.)
        END DO
        CALL dealloc_psi(g,delete_all=.TRUE.)

        DO i=nb_diago+1,max_diago
          CALL dealloc_psi(psi(i),delete_all=.TRUE.)
        END DO
        DO i=1,max_diago
          CALL dealloc_psi(psi0(i),delete_all=.TRUE.)
        END DO

        CALL dealloc_NParray(psi0,'psi0',name_sub)
        CALL dealloc_NParray(Hpsi,'Hpsi',name_sub)

        write(out_unitp,*) 'total memory',para_mem%mem_tot
      ENDIF

      !----------------------------------------------------------
      IF (debug) THEN
        DO i=1,nb_diago
          write(out_unitp,*) 'psi(i), i:',i
          CALL ecri_psi(psi=psi(i))
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
      !----------------------------------------------------------

 END SUBROUTINE sub_propagation_Davidson
!=======================================================================================

!=======================================================================================
 SUBROUTINE ReadWP0_Davidson(psi,psi0,Vec0,nb_diago,max_diago,          &
                             para_Davidson,cplx)
 USE mod_system
 USE mod_psi,     ONLY : param_psi,alloc_psi,norm2_psi,renorm_psi,      &
                   Set_symab_OF_psiBasisRep,sub_PsiBasisRep_TO_GridRep, &
                   param_WP0,Set_Random_psi,ecri_init_psi

 USE mod_WP0
 USE mod_propa,   ONLY : param_Davidson
 IMPLICIT NONE

 TYPE (param_Davidson) :: para_Davidson
 logical               :: cplx

 !----- WP, energy ... -----------------------------------
 TYPE (param_WP0)          :: para_WP0

 integer                   :: nb_diago,max_diago
 TYPE (param_psi)          :: psi(max_diago)
 TYPE (param_psi)          :: psi0(max_diago)

 real (kind=Rkind),allocatable :: Vec0(:,:)
 integer :: i
 real (kind=Rkind) :: a

 !----- for debuging --------------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='ReadWP0_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------

 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) ' nb_diago',nb_diago
   write(out_unitp,*) ' max_diago',max_diago

   write(out_unitp,*) ' para_Davidson',para_Davidson
   write(out_unitp,*)
   flush(out_unitp)
 END IF
 !-----------------------------------------------------------

 !------ read guess vectors ---------------------------------
 IF (nb_diago > 0 .OR. para_Davidson%read_WP) THEN
   para_WP0%nb_WP0              = nb_diago
   para_WP0%read_file           = para_Davidson%read_WP
   IF (para_Davidson%nb_readWP_OF_List > 0) THEN
     para_WP0%read_listWP0        = .FALSE.
   ELSE
     para_WP0%read_listWP0        = para_Davidson%read_listWP
   END IF
   para_WP0%WP0cplx             = cplx
   para_WP0%file_WP0%name       = para_Davidson%name_file_readWP
   para_WP0%file_WP0%formatted  = para_Davidson%formatted_file_readWP

   CALL sub_read_psi0(psi,para_WP0,max_diago,                      &
                      symab=para_Davidson%symab,ortho=.TRUE.)

   nb_diago = para_WP0%nb_WP0
   para_Davidson%nb_WP0 = para_WP0%nb_WP0

   IF (nb_diago < 1) THEN
     write(out_unitp,*) ' ERROR while reading the vector(s)'
     write(out_unitp,*) ' ... the number is 0'
     write(out_unitp,*) ' Probably, in the namelist (davidson)...'
     write(out_unitp,*) '   ... you select a wrong symmetry',para_Davidson%symab
     STOP
   END IF

   DO i=1,nb_diago
     CALL norm2_psi(psi(i))
     IF (debug) write(out_unitp,*) '   norm^2 of psi(i)',i,psi(i)%norm2
     IF ( abs(psi(i)%norm2-ONE) > ONETENTH**8) THEN
       write(out_unitp,*) ' ERROR while reading the vector(s)'
       write(out_unitp,*) ' ... the norm^2 of psi(i) is /= 1',i,psi(i)%norm2
       STOP
     END IF
   END DO

 ELSE
   para_Davidson%nb_WP0 = 0
   nb_diago = 1
   CALL alloc_psi(psi(1))
   CALL Set_Random_psi(psi(1))
   CALL Set_symab_OF_psiBasisRep(psi(1),para_Davidson%symab)
   CALL renorm_psi(psi(1),BasisRep=.TRUE.)
 END IF

 IF (para_Davidson%With_Grid) THEN
   DO i=1,max(1,nb_diago)
     CALL sub_PsiBasisRep_TO_GridRep(psi(i))
     CALL alloc_psi(psi(i),BasisRep=.FALSE.,GridRep=.TRUE.)
     !CALL Overlap_psi1_psi2(Overlap,psi(i),psi(i),With_Grid=.TRUE.)
     !write(out_unitp,*) 'Norm DVR',Overlap
   END DO
   DO i=max(1,nb_diago),max_diago
     psi(i)%BasisRep  = .FALSE.
     psi(i)%GridRep   = .TRUE.
   END DO
 END IF

 IF (para_Davidson%project_WP0) THEN
   CALL alloc_NParray(vec0,[nb_diago,nb_diago],"vec0",name_sub)
   vec0 = Identity_Mat(nb_diago)

   write(out_unitp,*) ' copy psi(:) to psi0(:)',para_Davidson%nb_WP0
   DO i=1,nb_diago
     psi0(i) = psi(i)
     !write(out_unitp,*) ' psi0(i)',i
     !CALL ecri_init_psi(psi0(i))
   END DO
 END IF

 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'psi(:), nb_diago',nb_diago
   DO i=1,nb_diago
     write(out_unitp,*) ' psi(i)',i
     CALL ecri_init_psi(psi(i))
   END DO
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE ReadWP0_Davidson
!=======================================================================================

!=======================================================================================
 SUBROUTINE sub_MakeHPsi_Davidson(it,psi,Hpsi,Ene,ndim0,                &
                                  para_H,para_Davidson,iunit,           &
                                  H_overlap,S_overlap,save_WP)
 USE mod_system
 USE mod_psi,         ONLY : param_psi,Overlap_psi1_psi2
 USE mod_psi_Op_MPI,  ONLY : Overlap_HS_matrix_MPI3,Overlap_H_matrix_MPI4

 USE mod_Op,          ONLY : param_Op, sub_TabOpPsi,sub_scaledOpPsi
 USE mod_propa,       ONLY : param_Davidson
 USE mod_MPI_aux
 IMPLICIT NONE


 integer           :: it,ndim0,iunit

 !----- WP, energy ... -----------------------------------
 TYPE (param_Davidson)     :: para_Davidson
 TYPE (param_psi)          :: psi(:)
 TYPE (param_psi)          :: Hpsi(:)

 real (kind=Rkind)         :: Ene(:)

 !----- Operator: Hamiltonian ----------------------------
 TYPE (param_Op)   :: para_H

 !------ working parameters --------------------------------
 integer              :: i,nb_thread,ndim
 logical              :: With_Grid,Op_Transfo
 complex (kind=Rkind) :: Overlap
 real (kind=Rkind)    :: auTOene
 
 Real (kind=Rkind),allocatable,intent(inout) :: H_overlap(:,:)
 Real (kind=Rkind),allocatable,intent(inout) :: S_overlap(:,:)
 Logical,                      intent(in)    :: save_WP

!----- for debuging --------------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='sub_MakeHPsi_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   !write(out_unitp,*) ' para_Davidson',para_Davidson
   write(out_unitp,*) 'para_Davidson%With_Grid',para_Davidson%With_Grid
   write(out_unitp,*) 'para_Davidson%Op_Transfo',para_Davidson%Op_Transfo
   write(out_unitp,*) 'MatOp_omp',MatOp_omp
   flush(out_unitp)
 END IF
 !-----------------------------------------------------------

 auTOene = get_Conv_au_TO_WriteUnit('E')
 ndim = size(psi)

!#if(run_MPI) 
!  CALL MPI_BCAST(ndim,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
 IF(openmpi .AND. MPI_scheme/=1) CALL MPI_BCAST_(ndim,size1_MPI,root_MPI)

 With_Grid  = para_Davidson%With_Grid
 Op_Transfo = para_Davidson%Op_Transfo

 IF (MatOp_omp /= 2) THEN
   nb_thread = 1
 ELSE
   nb_thread = MatOp_maxth
 END IF


 !----------------------------------------------------------
 !- Hpsi(:) ------------------------------------------
 IF (debug) write(out_unitp,*) 'Hpsi(:)',it,ndim,ndim0
 IF (debug) flush(out_unitp)

 CALL sub_TabOpPsi(psi(ndim0+1:ndim),Hpsi(ndim0+1:ndim),para_H,With_Grid=With_Grid,TransfoOp=Op_Transfo)

  ! DO scaling first
  IF (para_Davidson%Hmax_propa) THEN
    DO i=ndim0+1,ndim
      CALL sub_scaledOpPsi(psi(i),Hpsi(i),ZERO,-ONE)
    ENDDO
  ENDIF

  IF(openmpi) THEN
    !> calculate matrix H_overlap(i,j) for <psi(i)|Hpsi(j)>  (ndim0<=i,j<=ndim)
    !> shared by all threads
    IF((.NOT. allocated(H_overlap)) .OR. ndim0<ndim .OR. (save_WP)) THEN
      ! first iteration or restart iteration
      IF(it==0 .OR. (save_WP)) THEN
        CALL Overlap_HS_matrix_MPI3(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
      ELSE
        ! S_overlap is calculated in sub_NewVec_Davidson
        CALL Overlap_H_matrix_MPI4(H_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
      ENDIF
    ENDIF
  ENDIF

 DO i=ndim0+1,ndim

   !CALL sub_OpPsi(psi(i),Hpsi(i),para_H,With_Grid=With_Grid,TransfoOp=Op_Transfo)

   IF(openmpi) THEN
     Ene(i)=H_overlap(i,i)
   ELSE
     CALL Overlap_psi1_psi2(Overlap,psi(i),Hpsi(i),With_Grid=With_Grid)
     Ene(i) = real(Overlap,kind=Rkind)
   ENDIF

   IF (debug) write(out_unitp,*) 'Davidson Hpsi done',i,                &
                    Ene(i)*auTOene,(Ene(i)-Ene(1))*auTOene
   IF(MPI_id==0) THEN
     write(iunit,*) 'Davidson Hpsi done',i,                               &
                      Ene(i)*auTOene,(Ene(i)-Ene(1))*auTOene
     flush(iunit)
   ENDIF
 END DO

 ! combined with previous lines
! IF (para_Davidson%Hmax_propa) THEN
!
!   DO i=ndim0+1,ndim
!
!     CALL sub_scaledOpPsi(psi(i),Hpsi(i),ZERO,-ONE)    ! scaling
!
!     CALL Overlap_psi1_psi2(Overlap,psi(i),Hpsi(i),With_Grid=With_Grid)
!     Ene(i) = real(Overlap,kind=Rkind)
!
!     write(iunit,*) 'Davidson Hpsi done',i,                      &
!                    Ene(i)*auTOene,(Ene(i)-Ene(1))*auTOene
!     flush(iunit)
!   END DO
! END IF

 IF (debug) write(out_unitp,*) 'Hpsi(:) done ',it,ndim,ndim0
 IF (debug) flush(out_unitp)
 !- Hpsi(:) ------------------------------------------
 !----------------------------------------------------------

 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE sub_MakeHPsi_Davidson
!=======================================================================================

!=======================================================================================
 SUBROUTINE sub_MakeH_Davidson(it,psi,Hpsi,H,H_overlap,S_overlap,para_Davidson)
 USE mod_system
 USE mod_psi,    ONLY : param_psi,Overlap_psi1_psi2
 USE mod_propa,  ONLY : param_Davidson
 USE mod_MPI_aux
 IMPLICIT NONE

 integer,                  intent(in)     :: it

 !----- WP ... -----------------------------------
 TYPE (param_Davidson),    intent(in)     :: para_Davidson
 TYPE (param_psi),         intent(in)     :: psi(:)
 TYPE (param_psi),         intent(in)     :: Hpsi(:)

 !----- Operator: Hamiltonian ----------------------------
 real (kind=Rkind), allocatable, intent(inout) :: H(:,:)
 Real (kind=Rkind),              intent(in)    :: H_overlap(:,:)
 Real (kind=Rkind),              intent(in)    :: S_overlap(:,:)

 !------ working parameters --------------------------------
 integer              :: i,j,ndim,ndim0
 real (kind=Rkind), allocatable :: H0(:,:)
 complex (kind=Rkind) :: Overlap

!----- for debuging --------------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='sub_MakeH_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   CALL sub_MakeS_Davidson(it,psi,para_Davidson%With_Grid,.FALSE.,S_overlap)
   flush(out_unitp)
 END IF
 !-----------------------------------------------------------

 IF ( size(psi) /= size(Hpsi) .AND. MPI_id==0) THEN
   write(out_unitp,*) 'ERROR in ',name_sub
   write(out_unitp,*) 'size(psi) /= size(Hpsi)',size(psi),size(Hpsi)
   STOP
 END IF
 ndim = size(psi)

 IF (allocated(H)) THEN
   ndim0 = size(H,dim=1)
   !deallocate(H) ; ndim0 = 0
 ELSE
   ndim0 = 0
 END IF
 
 IF(openmpi) THEN
!   CALL MPI_BCAST(ndim, size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!   CALL MPI_BCAST(ndim0,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
   IF(MPI_scheme/=1) THEN
     CALL MPI_BCAST_(ndim, size1_MPI,root_MPI)
     CALL MPI_BCAST_(ndim0,size1_MPI,root_MPI)
   ENDIF

   CALL allocate_array(H,1,ndim,1,ndim)
!   IF(allocated(H)) CALL dealloc_NParray(H,"H",name_sub)
!   CALL alloc_NParray(H,[ndim,ndim],"H",name_sub)
   H(1:ndim,1:ndim)=H_overlap(1:ndim,1:ndim)
 ELSE

 !----------------------------------------------------------
 !- First save H in H0
 IF (debug) write(out_unitp,*) 'H mat',it,ndim,ndim0
 IF (debug) write(out_unitp,*) 'shape H mat',it,shape(H)

 IF (debug) flush(out_unitp)

 !block 1,1: ndim0*ndim0
 IF (ndim0 > 0) THEN
   CALL alloc_NParray(H0,[ndim0,ndim0],"H0",name_sub)
   H0(:,:) = H(:,:)


   CALL dealloc_NParray(H,"H",name_sub)
   CALL alloc_NParray(H,[ndim,ndim],"H",name_sub)
   H(:,:) = ZERO

   H(1:ndim0,1:ndim0) = H0(:,:)

   CALL dealloc_NParray(H0,"H0",name_sub)
 ELSE
   CALL alloc_NParray(H,[ndim,ndim],"H",name_sub)
   H(:,:) = ZERO
 END IF

 !block: 2,1 (ndim0-ndim)*ndim0
 DO i=1,ndim0
 DO j=ndim0+1,ndim
   CALL Overlap_psi1_psi2(Overlap,psi(j),Hpsi(i),With_Grid=para_Davidson%With_Grid)
   H(j,i) = real(Overlap,kind=Rkind)
 END DO
 END DO

 !blocks: 1,2 ndim0*(ndim0-ndim) + 2,2: (ndim0-ndim)*(ndim0-ndim)
 DO i=ndim0+1,ndim
 DO j=1,ndim
   CALL Overlap_psi1_psi2(Overlap,psi(j),Hpsi(i),With_Grid=para_Davidson%With_Grid)
   !write(out_unitp,*) 'H,i,j',i,j,Overlap
   H(j,i) = real(Overlap,kind=Rkind)
 END DO
 END DO

 ENDIF ! for openmpi

 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE sub_MakeH_Davidson
!=======================================================================================

 SUBROUTINE sub_MakeS_Davidson(it,psi,With_Grid,Print_Mat,S_overlap)

 USE mod_system
 USE mod_psi,         ONLY : param_psi,Overlap_psi1_psi2
 USE mod_propa,       ONLY : param_Davidson
 IMPLICIT NONE


 integer,          intent(in)         :: it
 TYPE (param_psi), intent(in)         :: psi(:)
 logical,          intent(in)         :: Print_Mat,With_Grid
 Real(kind=Rkind), intent(in)         :: S_overlap(:,:)

 !------ working parameters ------------------------------
 !----- Operator: Hamiltonian ----------------------------
 real (kind=Rkind), allocatable :: S(:,:)
 real (kind=Rkind) :: max_Sii,max_Sij

 integer              :: i,j,ndim
 complex (kind=Rkind) :: Overlap

!----- for debuging --------------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='sub_MakeS_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*)
   flush(out_unitp)
 END IF
 !-----------------------------------------------------------
ndim = size(psi)

 !- check the orthogonality ------------------------

IF(openmpi) THEN
  CALL sub_ana_S(S_overlap(1:ndim,1:ndim),ndim,max_Sii,max_Sij,.TRUE.)
  IF (Print_Mat ) CALL Write_Mat(S_overlap(1:ndim,1:ndim),out_unitp,5)
  flush(out_unitp)
ELSE
 CALL alloc_NParray(S,[ndim,ndim],"S",name_sub)

 DO j=1,ndim
 DO i=1,ndim
   CALL Overlap_psi1_psi2(Overlap,psi(i),psi(j),With_Grid)
   S(i,j) = real(Overlap,kind=Rkind)
 END DO
 END DO

 CALL sub_ana_S(S,ndim,max_Sii,max_Sij,.TRUE.)
 IF (Print_Mat) CALL Write_Mat(S,out_unitp,5)
 flush(out_unitp)
 CALL dealloc_NParray(S,"S",name_sub)
ENDIF ! openmpi

 !- check the orthogonality ------------------------

 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
   flush(out_unitp)
 END IF
 !----------------------------------------------------------

 END SUBROUTINE sub_MakeS_Davidson

 SUBROUTINE MakeResidual_Davidson(j,g,psi,Hpsi,Ene,Vec)
 USE mod_system
 USE mod_psi,    ONLY : param_psi,Set_symab_OF_psiBasisRep,             &
                        Overlap_psi1_psi2,dealloc_psi
 USE mod_propa,  ONLY : param_Davidson
 IMPLICIT NONE


 !----- WP, energy ... -----------------------------------
 integer, intent(in)   :: j
 TYPE (param_psi), intent(inout)       :: g
 TYPE (param_psi), intent(in)          :: psi(:),Hpsi(:)
 real (kind=Rkind), intent(in)         :: Ene(:)
 real (kind=Rkind), intent(in)         :: Vec(:,:)


 !------ working parameters --------------------------------
 integer              :: i,ndim,isym

 !----- for debuging --------------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='MakeResidual_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) '  ndim',size(Vec,dim=1)
   flush(out_unitp)
 END IF
 !-----------------------------------------------------------

 ndim = size(Vec,dim=1)
 isym = maxloc(abs(Vec(:,j)),dim=1) ! to find the rigth symmetry

 IF (allocated(g%RvecB)) THEN
   g%RvecB = ZERO
   DO i=1,ndim
     g%RvecB = g%RvecB + Hpsi(i)%RvecB*Vec(i,j) - psi(i)%RvecB * (Ene(j) * Vec(i,j))
   END DO
 ELSE IF (allocated(g%CvecB)) THEN
   g%CvecB = ZERO
   DO i=1,ndim
     g%CvecB = g%CvecB + Hpsi(i)%CvecB*Vec(i,j) - psi(i)%CvecB * (Ene(j) * Vec(i,j))
   END DO
 ELSE IF (allocated(g%RvecG)) THEN
   g%RvecG = ZERO
   DO i=1,ndim
     g%RvecG = g%RvecG + Hpsi(i)%RvecG*Vec(i,j) - psi(i)%RvecG * (Ene(j) * Vec(i,j))
   END DO
 ELSE IF (allocated(g%CvecG)) THEN
   g%CvecG = ZERO
   DO i=1,ndim
     g%CvecG = g%CvecG + Hpsi(i)%CvecG*Vec(i,j) - psi(i)%CvecG * (Ene(j) * Vec(i,j))
   END DO
 END IF

 CALL Set_symab_OF_psiBasisRep(g,symab=psi(isym)%symab)


! g = ZERO
! DO i=1,ndim
!   psiTemp = Hpsi(i) - psi(i) * Ene(j)
!   g = g + psiTemp*Vec(i,j)
! END DO
! CALL dealloc_psi(psiTemp,delete_all=.TRUE.)


 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

END SUBROUTINE MakeResidual_Davidson

 SUBROUTINE sub_NewVec_Davidson(it,psi,Hpsi,Ene,Ene0,EneRef,Vec,   &
                                converge,VecToBeIncluded,          &
                                nb_diago,max_diago,                &
                                para_Davidson,fresidu,ndim,        &
                                Op_Transfo,E0_Transfo,S_overlap)
 USE mod_system
 USE mod_basis
 USE mod_ana_psi_MPI, ONLY:norm_psi_MPI
 USE mod_Davidson_MPI,ONLY:Schmidt_process_MPI,MakeResidual_Davidson_j_MPI3
 USE mod_psi_Op_MPI,  ONLY:Set_symab_OF_psiBasisRep_MPI,calculate_overlap_MPI,         &
                          calculate_overlap1D_MPI
 USE mod_psi,         ONLY:param_psi,dealloc_psi,norm2_psi,renorm_psi,   &
                          Set_symab_OF_psiBasisRep,Overlap_psi1_psi2
 USE mod_propa,       ONLY:param_Davidson
 USE mod_MPI_aux
 IMPLICIT NONE

 integer           :: it,fresidu,ndim,isym
 integer           :: nb_diago,max_diago
 logical           :: Op_Transfo
 real (kind=Rkind) :: E0_Transfo



 !----- WP, energy ... -----------------------------------
 TYPE (param_Davidson)     :: para_Davidson
 TYPE (param_psi)          :: psi(max_diago)
 TYPE (param_psi)          :: Hpsi(max_diago)

 real (kind=Rkind)         :: Ene(max_diago),Ene0(max_diago)
 real (kind=Rkind)         :: EneRef(max_diago)
 real (kind=Rkind)         :: Vec(:,:)

 logical                   :: converge(max_diago)
 logical                   :: VecToBeIncluded(max_diago)

  Real(kind=Rkind),intent(inout),allocatable  :: S_overlap(:,:)
  Real(kind=Rkind)                            :: S_overlap_add(2*ndim,2*ndim)
  Real(kind=Rkind)                            :: Norm_Const
  Logical                                     :: sym_change

 !------ working parameters --------------------------------
 integer              :: i,j,j_ini,j_end,nb_added_states,ib,ndim0,iresidual
 complex (kind=Rkind) :: Overlap
 real (kind=Rkind)    :: a,Di,RS,Ene0_At_ib
 TYPE (param_psi)     :: psiTemp

 !----- for debuging --------------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='sub_NewVec_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) '  nb_diago   ',nb_diago
   write(out_unitp,*) '  max_diago  ',max_diago
   write(out_unitp,*) '  NewVec_type',para_Davidson%NewVec_type
   write(out_unitp,*) '  VecToBeIncluded',VecToBeIncluded(1:ndim)
   flush(out_unitp)
 END IF
 !-----------------------------------------------------------

 ndim0 = size(Vec(:,1))

 !write(out_unitp,*) ' symab: psi(:)',psi(1:ndim0)%symab


 !----------------------------------------------------------
 !- new vectors ------------------------------
 IF (debug) write(out_unitp,*) '  new psi',it,ndim,ndim0
 IF (debug) flush(out_unitp)
 IF (para_Davidson%residual_max_nb == 1) THEN
   j_ini = fresidu
   j_end = fresidu
 ELSE
   j_ini = 1
   j_end = ndim
 END IF
 IF (debug) write(out_unitp,*) '  j_ini,j_end',it,j_ini,j_end
 iresidual = 0
 DO j=j_ini,j_end
   IF (ndim == max_diago) EXIT
   IF (.NOT. converge(j) .AND. VecToBeIncluded(j)) THEN

     iresidual = iresidual + 1
     IF (iresidual > para_Davidson%residual_max_nb) EXIT

     isym = maxloc(abs(Vec(:,j)),dim=1) ! to find the rigth symmetry
     IF(keep_MPI) psi(ndim+1) = psi(isym) ! to be allocated correctly
     IF(keep_MPI) psi(ndim+1) = ZERO

     SELECT CASE (para_Davidson%NewVec_type)
     CASE (1) ! just the residual: equivament to Lanczos

       IF(openmpi) THEN
         CALL MakeResidual_Davidson_j_MPI3(j,psi(ndim+1),psi,Hpsi,Ene,Vec)
       ELSE
         CALL MakeResidual_Davidson(j,psi(ndim+1),psi,Hpsi,Ene,Vec)
       ENDIF

     !CASE 2 default
     CASE (2) ! Davidson
       !write(out_unitp,*) ' symab: psi(isym)',isym,psi(isym)%symab
       !write(out_unitp,*) ' vec(:,j)',j,Vec(:,j)

       DO i=1,ndim0
         Di = EneRef(i)
         a = Vec(i,j)
         IF (abs(Di - Ene(j)) > para_Davidson%conv_resi) THEN
           a = a / (Di - Ene(j))
         ELSE
           a = a / (Di - 0.999_Rkind*Ene0(j))
         END IF
         !write(out_unitp,*) 'i,j',i,j,a
         IF (abs(a) > ONETENTH**10) THEN
           !write(out_unitp,*) ' symab: psi(i)',i,psi(i)%symab

           IF(keep_MPI) psiTemp = Hpsi(i) - psi(i) * Ene(j)
           !write(out_unitp,*) ' symab: psiTemp',psi(i)%symab

           IF(keep_MPI) psi(ndim+1) = psi(ndim+1) + psiTemp * a
         END IF
       END DO

     CASE (3) ! Davidson

       DO i=1,ndim0
         Di = EneRef(i)
         a = Vec(i,j)
         IF (abs(Di - Ene(j)) > para_Davidson%conv_resi) THEN
           a = a / (Di - Ene(j))
         ELSE
           a = ZERO
         END IF
         !write(out_unitp,*) 'i,j',i,j,a
         IF (abs(a) > ONETENTH**10) THEN
          IF(keep_MPI) psiTemp = Hpsi(i) - psi(i) * Ene(j)
          IF(keep_MPI) psi(ndim+1) = psi(ndim+1) + psiTemp * a
         END IF
       END DO

     CASE (4) ! Davidson+precondioner
       ! first the residual
       IF(openmpi) THEN
         ! note, part of psi(ndim+1)%RvecB are kept on each threads
         CALL MakeResidual_Davidson_j_MPI3(j,psi(ndim+1),psi,Hpsi,Ene,Vec)
       ELSE
         CALL MakeResidual_Davidson(j,psi(ndim+1),psi,Hpsi,Ene,Vec)
       ENDIF

       !@chen psiTemp%nb_tot always 0? 
       ! then the scaling with respect to 1/(H0-Ene(j))
       DO ib=1,psiTemp%nb_tot
         Ene0_At_ib = get_Ene0_AT_ib_FROM_Basis(ib,psi(ndim+1)%BasisnD, &
                                                   psi(ndim+1)%Basis2n)
         IF (Op_transfo) THEN
           Di = Ene(j) - (Ene0_At_ib-E0_Transfo)**2
         ELSE
           Di = Ene(j) -  Ene0_At_ib
         END IF
         IF (abs(Di) > para_Davidson%conv_resi) THEN
           a = ONE / Di
         ELSE
           !a = ZERO
           a = ONE / (Di +ONETENTH**3)
         END IF
         IF(keep_MPI) psi(ndim+1)%RvecB(ib) = psi(ndim+1)%RvecB(ib) * a
       END DO

     CASE DEFAULT
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' This NewVec_type is unknown: ',para_Davidson%NewVec_type
       STOP 'Unkonwn NewVec_type!'
     END SELECT
     !write(out_unitp,*) ' symab: psi(isym), new vec',psi(isym)%symab,psi(ndim+1)%symab

     ! use Set_symab_OF_psiBasisRep_MPI later
     IF(keep_MPI) CALL Set_symab_OF_psiBasisRep(psi(ndim+1),psi(isym)%symab)
     !write(out_unitp,*) ' symab: psi(isym), new vec set sym',psi(isym)%symab,psi(ndim+1)%symab

     !- new vectors -------------------------------
     !- Schmidt ortho ------------------------------------
     !write(out_unitp,*) 'Schmidt ortho',it
     IF(openmpi) THEN
       CALL Schmidt_process_MPI(S_overlap_add(:,ndim+1-ndim0),psi,ndim,isym,           &
                                With_Grid=para_Davidson%With_Grid) 
       ! check and normalization
       CALL norm_psi_MPI(psi(ndim+1),1) ! get normalization constant
       IF(psi(ndim+1)%norm2<ONETENTH**10) CYCLE ! error case
       Norm_Const=psi(ndim+1)%norm2
       
       CALL Set_symab_OF_psiBasisRep_MPI(psi(ndim+1),psi(isym)%symab,sym_change)

       CALL norm_psi_MPI(psi(ndim+1),2) ! normalization

       ! synchronism S_overlap_add
       IF(sym_change) THEN ! Rvec changed in Set_symab_OF_psiBasisRep_MPI
         CALL calculate_overlap1D_MPI(psi,ndim+1,With_Grid=para_Davidson%With_Grid,    &
                                      S_overlap1D=S_overlap_add(:,ndim+1-ndim0))
       ELSE
         S_overlap_add(:,ndim+1-ndim0)=S_overlap_add(:,ndim+1-ndim0)/Norm_Const
       ENDIF
     ELSE

       IF (para_Davidson%With_Grid) THEN
         CALL renorm_psi(psi(ndim+1))
         DO i=1,ndim
           CALL Overlap_psi1_psi2(Overlap,psi(ndim+1),psi(i),      &
                                  With_Grid=para_Davidson%With_Grid)
           RS = real(Overlap,kind=Rkind)
           IF (RS == ZERO) CYCLE
           psi(ndim+1)%RvecG = (psi(ndim+1)%RvecG - psi(i)%RvecG * RS)/sqrt(ONE-RS**2)
         END DO
         CALL renorm_psi(psi(ndim+1))
         DO i=1,ndim
           CALL Overlap_psi1_psi2(Overlap,psi(ndim+1),psi(i),      &
                                  With_Grid=para_Davidson%With_Grid)
           RS = real(Overlap,kind=Rkind)
           IF (RS == ZERO) CYCLE
           psi(ndim+1)%RvecG = (psi(ndim+1)%RvecG - psi(i)%RvecG * RS)/sqrt(ONE-RS**2)
         END DO

       ELSE

         RS = dot_product(psi(ndim+1)%RvecB,psi(ndim+1)%RvecB)
         psi(ndim+1)%RvecB = psi(ndim+1)%RvecB / sqrt(RS)       
         DO i=1,ndim
           RS = dot_product(psi(ndim+1)%RvecB,psi(i)%RvecB)
           IF (RS == ZERO) CYCLE
           psi(ndim+1)%RvecB = psi(ndim+1)%RvecB - psi(i)%RvecB * RS
         END DO

         RS = dot_product(psi(ndim+1)%RvecB,psi(ndim+1)%RvecB)
         psi(ndim+1)%RvecB = psi(ndim+1)%RvecB / sqrt(RS)
         DO i=1,ndim
           RS = dot_product(psi(ndim+1)%RvecB,psi(i)%RvecB)
           IF (RS == ZERO) CYCLE
           psi(ndim+1)%RvecB = psi(ndim+1)%RvecB - psi(i)%RvecB * RS
         END DO

       END IF

       CALL norm2_psi(psi(ndim+1))
       IF (psi(ndim+1)%norm2 < ONETENTH**10) CYCLE ! otherwise dependent vector
       
       !write(out_unitp,*) ' symab: psi(isym), new vec ortho',psi(isym)%symab,psi(ndim+1)%symab

       IF (debug) THEN
          write(out_unitp,*) '  symab new vectors',psi(ndim+1)%symab
          CALL norm2_psi(psi(ndim+1))
          write(out_unitp,*) '  norm^2 of new vectors',psi(ndim+1)%norm2
       END IF

       CALL Set_symab_OF_psiBasisRep(psi(ndim+1),psi(isym)%symab)
       !write(out_unitp,*) ' symab: psi(isym), new vec set sym ',psi(isym)%symab,psi(ndim+1)%symab

       CALL renorm_psi(psi(ndim+1))
       !write(out_unitp,*) ' symab: psi(isym), new renormalized vector ',psi(isym)%symab,psi(ndim+1)%symab

       !- Schmidt ortho ------------------------------------
       !write(out_unitp,*) 'n+1, vec',ndim+1,psi(ndim+1)%RvecB
       !write(out_unitp,*) ' new vec symab, bits(symab)',WriteTOstring_symab(psi(ndim+1)%symab)
     ENDIF
     ndim = ndim + 1

   END IF

   IF (debug) write(out_unitp,*) '  new psi',it,ndim,ndim0
   IF (debug) flush(out_unitp)
 END DO

 nb_added_states = ndim-ndim0
 write(out_unitp,*) 'it, nb         new state(s)',it,nb_added_states
 !- new vectors -------------------------------
 !----------------------------------------------------------

  IF(openmpi) THEN
    ! get new S_overlap
    IF(ndim>ndim0) CALL increase_martix(S_overlap,name_sub,ndim0,ndim)
    DO i=ndim0+1,ndim
      DO j=1,i
        S_overlap(i,j)=S_overlap_add(j,i-ndim0)
        S_overlap(j,i)=S_overlap(i,j)
      ENDDO
    ENDDO
  ENDIF

 CALL dealloc_psi(psiTemp,delete_all=.TRUE.)

 !----------------------------------------------------------
 IF (debug) THEN
   CALL sub_MakeS_Davidson(it,psi(1:ndim),With_Grid=para_Davidson%With_Grid,           &
                           Print_Mat=.TRUE.,S_overlap=S_overlap)
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

END SUBROUTINE sub_NewVec_Davidson

 SUBROUTINE Sort_VecToBeIncluded_Davidson(Ene,Vec,VecToBeIncluded)
 USE mod_system
 IMPLICIT NONE

 real (kind=Rkind), intent(inout) :: Vec(:,:)
 real (kind=Rkind), intent(inout) :: Ene(:)
 logical, intent(inout)           :: VecToBeIncluded(:)

 !------ working parameters --------------------------------
 integer              :: i,i0
 integer              :: ndim,nb_included
 real (kind=Rkind)    :: E1
 real (kind=Rkind)    :: trav(size(Vec(:,1)))

 !----- for debuging --------------------------------------------------
 character (len=*), parameter ::name_sub='Sort_VecToBeIncluded_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 ndim  = size(Vec(:,1))


 ! sort the vectors to be included
 i0 = 0
 DO i=1,ndim
   IF (VecToBeIncluded(i)) THEN
     i0 = i0 + 1
     IF (i0 /= i) THEN ! permutation
       trav(:)   = vec(:,i)
       vec(:,i)  = vec(:,i0)
       vec(:,i0) = trav(:)

       E1        = Ene(i)
       Ene(i)    = Ene(i0)
       Ene(i0)   = E1
     END IF
   END IF

 END DO

 nb_included = count(VecToBeIncluded)
 VecToBeIncluded(:) = .FALSE.
 VecToBeIncluded(1:nb_included) = .TRUE.


 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE Sort_VecToBeIncluded_Davidson

!=======================================================================================
!     Sort Davidson
!=======================================================================================
 SUBROUTINE sub_projec_Davidson(Ene,VecToBeIncluded,nb_diago,min_Ene,min_pot,ZPE,  &
                                psi,psi0,Vec,Vec0,para_Davidson,it,print_project)
 USE mod_system
 USE mod_psi,     ONLY : param_psi,Overlap_psi1_psi2
 USE mod_propa,   ONLY : param_Davidson
 USE mod_MPI_aux
 IMPLICIT NONE

 TYPE (param_Davidson) :: para_Davidson
 logical               :: print_project

 logical           :: VecToBeIncluded(:)
 real (kind=Rkind) :: min_Ene,min_pot,ZPE,Ene(:)
 integer           :: it,nb_diago,kmin
 TYPE (param_psi),intent(in)  :: psi(:),psi0(:) ! size max_WP


 integer           :: k,i,ndim,ndim0
 real (kind=Rkind), allocatable :: Vec0(:,:),Vec(:,:)
 real (kind=Rkind) :: S0itmax,S0it
 integer :: klowestWP
 logical :: test_vec0

 complex (kind=Rkind) :: Overlap
 real (kind=Rkind)    :: Spsi_psi0(size(Vec(:,1)))

 !----- for debuging --------------------------------------------------
 character (len=*), parameter ::name_sub='sub_projec_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) ' lower_states',para_Davidson%lower_states
   write(out_unitp,*) ' ndim ',size(Vec,dim=1)
   IF (allocated(Vec0)) write(out_unitp,*) ' ndim0',size(Vec0,dim=1)
   write(out_unitp,*) ' nb_diago',nb_diago
   IF (allocated(Vec)) THEN
     write(out_unitp,*) 'Vec',shape(Vec)
     !CALL Write_Mat(Vec,out_unitp,5)
   END IF
   IF (allocated(Vec0)) THEN
     write(out_unitp,*) 'Vec0',shape(Vec0)
     !CALL Write_Mat(Vec0,out_unitp,5)
   END IF
   write(out_unitp,*) 'Ene',Ene(1:ndim)
   write(out_unitp,*) 'min_Ene',min_Ene
   write(out_unitp,*) 'min_pot',min_pot
   write(out_unitp,*) 'ZPE',ZPE
   write(out_unitp,*) 'Max_ene',para_Davidson%Max_ene
 END IF

 ndim  = size(Vec,dim=1)
 IF (debug .OR. print_project) write(out_unitp,*) 'Vector projections'

 ! find the physical min_ene (to avoid holes)
 IF (min_pot == huge(ONE)) THEN
   min_Ene = -huge(ONE)
 ELSE
   min_Ene = min_pot
 END IF
 klowestWP = 0
 IF (ndim > 0 .AND. para_Davidson%num_LowestWP > 0 .AND. para_Davidson%num_LowestWP <= ndim) THEN

   test_vec0 = allocated(vec0)
   IF (test_vec0) test_vec0 = (para_Davidson%num_LowestWP <= size(Vec0,dim=1))

   IF ( test_vec0 ) THEN

     IF(keep_MPI) THEN
     ndim0 = size(Vec0,dim=1)
     klowestWP = 0
     S0itmax   = ZERO
     DO k=1,ndim
       S0it = abs(dot_product(Vec0(:,para_Davidson%num_LowestWP),Vec(1:ndim0,k)))
       IF (S0it > S0itmax) THEN
         klowestWP = k
         S0itmax   = S0it
         min_Ene   = ene(k)
       END IF
     END DO
     ENDIF ! keep_MPI
!#if(run_MPI) 
!     CALL MPI_Bcast(min_Ene,size1_MPI,Real_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
     IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(min_Ene,size1_MPI,root_MPI)

   ELSE IF (para_Davidson%project_WP0 .AND. para_Davidson%num_LowestWP <= para_Davidson%nb_WP0) THEN

     IF(keep_MPI) THEN
     DO k=1,ndim
       CALL Overlap_psi1_psi2(Overlap,psi(k),psi0(para_Davidson%num_LowestWP))
       Spsi_psi0(k) = real(Overlap,kind=Rkind)
     END DO

     klowestWP = 0
     S0itmax   = ZERO
     DO k=1,ndim
       S0it = abs(dot_product(Spsi_psi0,Vec(:,k)))
       IF (S0it > S0itmax) THEN
         klowestWP = k
         S0itmax   = S0it
         min_Ene   = ene(k)
       END IF
     END DO
     ENDIF ! keep_MPI
!#if(run_MPI) 
!     CALL MPI_Bcast(min_Ene,size1_MPI,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
     IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(min_Ene,size1_MPI,root_MPI)

   END IF
 END IF

 IF (para_Davidson%lower_states) THEN
   IF (para_Davidson%Hmax_propa) THEN
     VecToBeIncluded(:) = .FALSE.
     VecToBeIncluded(1) = .TRUE.
   ELSE
     VecToBeIncluded(:) = .FALSE.
     IF (para_Davidson%all_lower_states) THEN
      DO i=1,ndim
        VecToBeIncluded(i) =(Ene(i) >= min_Ene .AND. Ene(i)-ZPE <= para_Davidson%Max_ene)
      END DO
     ELSE
      DO i=1,ndim
        VecToBeIncluded(i) =(Ene(i) >= min_Ene .AND. count(VecToBeIncluded) < nb_diago )
      END DO
      IF (count(VecToBeIncluded) /= nb_diago) THEN
        write(out_unitp,*) ' ERROR in ',name_sub,' from ', MPI_id
        write(out_unitp,*) 'nb_diago ',nb_diago
        write(out_unitp,*) 'VecToBeIncluded ',VecToBeIncluded
        write(out_unitp,*) 'Ene(:) ',Ene
        write(out_unitp,*) 'min_Ene ',min_Ene
 
        write(out_unitp,*) ' number of VecToBeIncluded /= nb_diago'
        STOP
      END IF
     END IF
   END IF

   CALL Sort_VecToBeIncluded_Davidson(Ene,Vec,VecToBeIncluded)

 ELSE
   IF (para_Davidson%project_WP0) THEN
     ndim0 = para_Davidson%nb_WP0
     CALL sub_projec3_Davidson(Vec,psi(1:ndim),psi0(1:ndim0),           &
                               VecToBeIncluded,para_Davidson%thresh_project,   &
                               Ene,min_Ene,print_project)

     nb_diago = count(VecToBeIncluded)

   ELSE

     ! need futher improvement for MPI here
     CALL sub_projec2_Davidson(Vec,Vec0,VecToBeIncluded,para_Davidson%thresh_project,Ene,min_Ene,print_project)
     nb_diago = count(VecToBeIncluded)

     IF(allocated(Vec0)) CALL dealloc_NParray(Vec0,"Vec0",name_sub)
     CALL alloc_NParray(Vec0,[ndim,ndim],"Vec0",name_sub)
     Vec0(:,:) = Vec(:,:)
   END IF
 END IF

 IF(MPI_id==0) Then
  IF (debug .OR. print_project) THEN 
     write(out_unitp,*) 'number of VecToBeIncluded',it,count(VecToBeIncluded)
     write(out_unitp,*) 'VecToBeIncluded',it,VecToBeIncluded(1:ndim)
     write(out_unitp,*) 'End Vector projections'
  END IF
 ENDIF

 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE sub_projec_Davidson
!=======================================================================================

 SUBROUTINE sub_projec2_Davidson(Vec,Vec0,VecToBeIncluded,thresh,Ene,min_Ene,print_project)
 USE mod_system
 USE mod_MPI_aux
 IMPLICIT NONE

 real (kind=Rkind) :: Vec0(:,:),Vec(:,:)
 real (kind=Rkind) :: Ene(:),min_Ene
 logical           :: VecToBeIncluded(:)
 real (kind=Rkind) :: thresh
 logical           :: print_project

 !------ working parameters --------------------------------
 real (kind=Rkind)    :: p,max_p,auTOene
 integer              :: i,j,i0,i1
 integer              :: ndim,ndim0,nb_included
 real (kind=Rkind)    :: E1
 real (kind=Rkind)    :: trav(size(Vec(:,1)))

 !----- for debuging --------------------------------------------------
 character (len=*), parameter ::name_sub='sub_projec2_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 ndim  = size(Vec(:,1))
 ndim0 = size(Vec0(:,1))
 auTOene = get_Conv_au_TO_WriteUnit('E')

 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) ' ndim,ndim0',ndim,ndim0
   write(out_unitp,*) ' shape Vec,Vec0',shape(Vec),shape(Vec0)

   write(out_unitp,21) Ene(1:ndim)*auTOene
 21     format(' Ene: ',50(1x,f18.4))

   !write(out_unitp,*)
   !write(out_unitp,*) 'Vec',ndim
   !CALL Write_Mat(Vec,out_unitp,5)
   !write(out_unitp,*) 'Vec0',ndim0
   !CALL Write_Mat(Vec0,out_unitp,5)
 END IF


 IF (debug .OR. print_project) write(out_unitp,*) 'Projection <psi|psi0>'
 VecToBeIncluded(:) = .FALSE.
 
 IF(keep_MPI) THEN
 DO j=1,ndim0
   max_p = ZERO
   DO i=1,ndim
      p = abs(dot_product(Vec(1:ndim0,i),Vec0(:,j)))
      max_p = max(p,max_p)
   END DO
   IF (debug .OR. print_project) write(out_unitp,*) ' New j of Vec0',j,'max_p',max_p
   DO i=1,ndim
      !write(out_unitp,*) i,Vec(1:nb_Vec0,i)
      p = abs(dot_product(Vec(1:ndim0,i),Vec0(:,j)))
      IF (p >= thresh*max_p) THEN
        VecToBeIncluded(i) = .TRUE.
        IF (debug .OR. print_project) write(out_unitp,*) 'i of Vec, p,Ene',i,p,Ene(i)*auTOene
      ELSE
        IF (debug) write(out_unitp,*) 'i of Vec, p,Ene',i,p,Ene(i)*auTOene
      END IF
   END DO
 END DO

 !Set to F, all vectors such Ene(i) < min_Ene
 DO i=1,ndim
   IF (VecToBeIncluded(i)) VecToBeIncluded(i) = (Ene(i) >= min_Ene)
 END DO
 ENDIF ! keep_MPI

!#if(run_MPI)
! CALL MPI_Bcast(VecToBeIncluded,ndim,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err) 
!#endif
 IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(VecToBeIncluded,ndim,root_MPI)

 CALL Sort_VecToBeIncluded_Davidson(Ene,Vec,VecToBeIncluded)

 IF (debug .OR. print_project) write(out_unitp,*) 'End projection <psi|psi0>'


 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE sub_projec2_Davidson

 SUBROUTINE sub_projec3_Davidson(Vec,psi,psi0,VecToBeIncluded,thresh,   &
                                 Ene,min_Ene,print_project)
 USE mod_system
 USE mod_psi,     ONLY : param_psi,Overlap_psi1_psi2
 USE mod_MPI_aux
 IMPLICIT NONE

 TYPE (param_psi)  :: psi(:),psi0(:)
 real (kind=Rkind) :: Vec(:,:)


 real (kind=Rkind) :: Ene(:),min_Ene
 logical           :: VecToBeIncluded(:)
 real (kind=Rkind) :: thresh
 logical           :: print_project

 !------ working parameters --------------------------------
 real (kind=Rkind)    :: Ene_save,max_p,auTOene
 integer              :: i,j,k,i0
 integer              :: ndim,ndim0,nb_included
 complex (kind=Rkind) :: Overlap
 real (kind=Rkind)    :: Spsi_psi0(size(Vec(:,1)),size(psi0))

 real (kind=Rkind)    :: Spsi_psi0_j(size(Vec(:,1)))
  real (kind=Rkind)   :: trav(size(Vec(:,1)))


 !----- for debuging --------------------------------------------------
 character (len=*), parameter ::name_sub='sub_projec3_Davidson'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------
 ndim  = size(Vec(:,1))
 ndim0 = size(psi0)
 IF (ndim < 1 .OR. ndim0 < 1) RETURN

 auTOene = get_Conv_au_TO_WriteUnit('E')

 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) ' ndim,ndim0',ndim,ndim0
   write(out_unitp,*) ' shape(vec)',shape(vec)

   write(out_unitp,21) Ene(1:ndim)*auTOene
 21     format(' Ene: ',50(1x,f18.4))

!   write(out_unitp,*) '================================='
!   write(out_unitp,*) '================================='
!   write(out_unitp,*) '================================='
!   write(out_unitp,*) ' psi0',ndim0
!   DO i=1,ndim0
!     CALL ecri_psi(ZERO,psi0(i))
!     write(out_unitp,*) '================================='
!   END DO
!   write(out_unitp,*) '================================='
!   write(out_unitp,*) '================================='
!   write(out_unitp,*) '================================='
!   write(out_unitp,*) ' psi',ndim
!   DO i=1,ndim
!     CALL ecri_psi(ZERO,psi(i))
!     write(out_unitp,*) '================================='
!   END DO
!   write(out_unitp,*) '================================='
!   write(out_unitp,*) '================================='
 END IF


 IF (debug .OR. print_project) write(out_unitp,*) 'Projection <psi|psi0>'
 VecToBeIncluded(:) = .FALSE.

 IF(keep_MPI) THEN
 DO j=1,ndim0
 DO k=1,ndim
   CALL Overlap_psi1_psi2(Overlap,psi(k),psi0(j))
   Spsi_psi0(k,j) = real(Overlap,kind=Rkind)
 END DO
 END DO

 DO j=1,ndim0
   DO i=1,ndim
     Spsi_psi0_j(i) = abs(dot_product(Vec(:,i),Spsi_psi0(:,j)))
   END DO
   max_p = maxval(Spsi_psi0_j)
   IF (debug .OR. print_project) write(out_unitp,*) ' New j of psi0',j,'max_p',max_p
   DO i=1,ndim
      !write(out_unitp,*) i,Vec(1:nb_Vec0,i)
      IF (Spsi_psi0_j(i) >= thresh*max_p) THEN
        VecToBeIncluded(i) = .TRUE.
        IF (debug .OR. print_project) write(out_unitp,*) 'i of psi, p,Ene',i,Spsi_psi0_j(i),Ene(i)*auTOene
      ELSE
        IF (debug) write(out_unitp,*) 'i of psi, p,Ene',i,Spsi_psi0_j(i),Ene(i)*auTOene
      END IF
   END DO
 END DO

 !Set to F, all vectors such Ene(i) < min_Ene
 DO i=1,ndim
   IF (VecToBeIncluded(i)) VecToBeIncluded(i) = (Ene(i) >= min_Ene)
 END DO
 ENDIF ! keep_MPI

!#if(run_MPI)
! CALL MPI_Bcast(VecToBeIncluded,ndim,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err) 
!#endif
 IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(VecToBeIncluded,ndim,root_MPI)

 CALL Sort_VecToBeIncluded_Davidson(Ene,Vec,VecToBeIncluded)


 IF (debug .OR. print_project) write(out_unitp,*) 'End projection <psi|psi0>'


 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE sub_projec3_Davidson


!=======================================================================================
!@brief set save_WP=.true. with external control
!
! create a file with name "Davidson_exit", the program wil exit
! and change the file name to "done_Davidson_exit".
!=======================================================================================
 SUBROUTINE exit_Davidson_external(exit_Davidson,save_WP,it)
  IMPLICIT NONE

  Logical,                       intent(inout) :: exit_Davidson
  Logical,                       intent(inout) :: save_WP
  Integer,                       intent(in)    :: it
  Logical                                      :: exist
  Integer                                      :: stat,iunit

  IF(it>2) THEN
    INQUIRE(FILE='Davidson_exit',EXIST=exist)
    IF(exist) THEN
      save_WP=.TRUE.
      exit_Davidson=.TRUE.
      open(NEWUNIT=iunit, FILE='Davidson_exit')
      close(iunit, status='delete') ! CALL file_close cannot be used
      open(NEWUNIT=iunit, FILE='done_Davidson_exit') 
      close(iunit) ! CALL file_close cannot be used
    ENDIF
  ENDIF
END SUBROUTINE exit_Davidson_external
END MODULE mod_Davidson
