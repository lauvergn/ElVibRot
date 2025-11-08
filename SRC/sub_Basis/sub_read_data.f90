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
! ++    basis set reading (1D+nD)
!
!       type_basis_set,n_contrac,A,B ...
!================================================================
      SUBROUTINE read_basis5(BasisnD,mole)
      USE EVR_system_m
      USE mod_basis
      use mod_Coord_KEO, only: CoordType
      IMPLICIT NONE

      !----- for the active basis set ------------------------------------
      TYPE (basis)               :: BasisnD

      !-- for Tnum and mole ----------------------------------------------
      TYPE (CoordType), intent(in) :: mole


      integer       :: i,i0,i1,nb_act1_test,i_Qdyn,i_Qbasis
      integer       :: liste_var_basis(mole%nb_var)
      logical       :: check_ba
      TYPE (basis)  :: BasisnD_loc

      character (len=Name_len)   :: name

!---------------------------------------------------------------------
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'read_basis5'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------
      CALL alloc_array(BasisnD_loc%tab_Pbasis,[mole%nb_act1+1],     &
                      'BasisnD_loc%tab_Pbasis',name_sub)
      DO i=1,size(BasisnD_loc%tab_Pbasis)
        CALL alloc_array(BasisnD_loc%tab_Pbasis(i)%Pbasis,              &
                        'BasisnD_loc%tab_Pbasis(i)%Pbasis',name_sub)
      END DO

      IF(MPI_id==0) THEN
        write(out_unit,*) '---------------------------------------'
        write(out_unit,*) '----------- BasisnD -------------------'
        write(out_unit,*) '---------------------------------------'
        write(out_unit,*)
      ENDIF

      ! parameter for BasisnD
      BasisnD_loc%type           = 1
      BasisnD_loc%name           = 'direct_prod'
      BasisnD_loc%contrac        = .FALSE.

!---------------------------------------------------------------
!    - read the parameters of all basis set --------------------
      !read(in_unit,*) name
      name='basisnd'
      CALL string_uppercase_TO_lowercase(name)
      !write(out_unit,*) 'name ',name

      IF (name .EQ. "nd" .OR. name .EQ. "basisnd") THEN
!        - read the parameters of all nD-basis set ---------------
         nb_act1_test = 0
         i = 1
         BasisnD_loc%nb_basis = 0
         BasisnD_loc%opt_param = 0
         DO WHILE (nb_act1_test < mole%nb_act1)

           CALL RecRead5_Basis(BasisnD_loc%tab_Pbasis(i)%Pbasis,mole)

           IF (BasisnD_loc%tab_Pbasis(i)%Pbasis%active) THEN
             BasisnD_loc%opt_param = BasisnD_loc%opt_param +            &
                              BasisnD_loc%tab_Pbasis(i)%Pbasis%opt_param

             nb_act1_test = nb_act1_test + BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
             i = i + 1
           END IF
           !write(out_unit,*) ' nb_act1_test: ',nb_act1_test,mole%nb_act1
         END DO
         BasisnD_loc%nb_basis = i-1

         IF (nb_act1_test /= mole%nb_act1) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) 'nb_act1_test /= nb_act1 !!',nb_act1_test,mole%nb_act1
           STOP
         END IF
      ELSE
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '   The old way to define the basis-set is not possible'
        STOP
      END IF

!     - check if the active coordinates are associated with a basis
      liste_var_basis(:) = 0
      DO i=1,BasisnD_loc%nb_basis
        DO i_Qbasis=1,BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
          i_Qdyn = BasisnD_loc%tab_Pbasis(i)%Pbasis%iQdyn(i_Qbasis)
          IF (liste_var_basis(i_Qdyn) == 1) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' The coordinate (iQdyn',i_Qdyn,                 &
                        ') is associated with TWO basis sets !!!'
            write(out_unit,*) ' CHECK your DATA'
            STOP
          END IF
          liste_var_basis(i_Qdyn) = 1
        END DO
      END DO

      check_ba = .TRUE.
      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 1 .AND.           &
                                           liste_var_basis(i) /= 1) THEN
          check_ba = .FALSE.
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' NO basis set for the ACTIVE coordinate',i
          write(out_unit,*) '   list_act_OF_Qdyn',mole%ActiveTransfo%list_act_OF_Qdyn(:)
          write(out_unit,*) '   liste_var_basis',liste_var_basis
          write(out_unit,*) ' CHECK your DATA'
        END IF
      END DO
      IF (.NOT. check_ba) STOP

      ! set ndim, iQdyn, nrho, ....
      BasisnD_loc%ndim = 0
      DO i=1,BasisnD_loc%nb_basis
        BasisnD_loc%ndim = BasisnD_loc%ndim +                           &
                            BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
      END DO

      CALL alloc_init_basis(BasisnD_loc)

      ! transfert of iQdyn, nrho --------------------------------
      i0 = 0
      DO i=1,BasisnD_loc%nb_basis
        IF (BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim == 0) CYCLE
        i1 = i0 + BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
        BasisnD_loc%iQdyn(i0+1:i1) = BasisnD_loc%tab_Pbasis(i)%Pbasis%iQdyn(:)
        BasisnD_loc%nrho(i0+1:i1)  = BasisnD_loc%tab_Pbasis(i)%Pbasis%nrho(:)
        i0 = i1
      END DO

      ! Set Tabder_Qdyn_TO_Qbasis(:)
      CALL alloc_NParray(BasisnD_loc%Tabder_Qdyn_TO_Qbasis,               &
                                                   [mole%nb_var],   &
                      'BasisnD_loc%Tabder_Qdyn_TO_Qbasis',name_sub,[0])
      BasisnD_loc%Tabder_Qdyn_TO_Qbasis(:) = 0
      DO i=1,BasisnD_loc%ndim
        BasisnD_loc%Tabder_Qdyn_TO_Qbasis(BasisnD_loc%iQdyn(i)) = i
      END DO

      BasisnD_loc%active            = .TRUE.
      IF(MPI_id==0) THEN
        write(out_unit,*)
        write(out_unit,*)
        write(out_unit,*) 'Number of active basis sets:',BasisnD_loc%nb_basis
      ENDIF

      IF (BasisnD_loc%nb_basis == 1) THEN
        IF(MPI_id==0) write(out_unit,*) 'WARNING: ONE layer of basis has been removed!!'
        CALL basis2TObasis1(BasisnD,BasisnD_loc%tab_Pbasis(1)%Pbasis)
      ELSE
        CALL basis2TObasis1(BasisnD,BasisnD_loc)
      END IF

      IF (BasisnD%MaxCoupling_OF_nDindB < 1) BasisnD%MaxCoupling_OF_nDindB = BasisnD%nb_basis

      IF(MPI_id==0) THEN
        write(out_unit,*) 'Tabder_Qdyn_TO_Qbasis',BasisnD%Tabder_Qdyn_TO_Qbasis(1:mole%nb_var)
        write(out_unit,*)
        write(out_unit,*) 'BasisnD%opt_param',BasisnD%opt_param
        write(out_unit,*)
      ENDIF

      CALL dealloc_basis(BasisnD_loc)
      !CALL RecWriteMini_basis(BasisnD)

      IF(MPI_id==0) THEN
        write(out_unit,*) '---------------------------------------'
        write(out_unit,*) '---------------------------------------'
      ENDIF
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BasisnD'
        CALL RecWrite_basis(BasisnD,.TRUE.)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      end subroutine read_basis5

      RECURSIVE SUBROUTINE RecRead5_Basis(basis_temp,mole)

      USE EVR_system_m
      USE mod_basis
      USE mod_Coord_KEO, only: CoordType
      IMPLICIT NONE

!----- for the active basis set ---------------------------------------
      TYPE (basis)  :: basis_temp

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole


      integer       :: i,i0,i1

!---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='RecRead5_Basis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

!     - read the namelist -------------------------------
      CALL read5_basis_nD(basis_temp,mole)

      IF (basis_temp%nb_basis > 0) THEN ! Direct_product, SparseBasis
        basis_temp%active = .TRUE.
        CALL alloc_array(basis_temp%tab_Pbasis,[basis_temp%nb_basis],&
                        'basis_temp%tab_Pbasis',name_sub)
        basis_temp%opt_param = 0
        DO i=1,basis_temp%nb_basis
          CALL alloc_array(basis_temp%tab_Pbasis(i)%Pbasis,             &
                          'basis_temp%tab_Pbasis(i)%Pbasis',name_sub)
          IF (basis_temp%packed) basis_temp%tab_Pbasis(i)%Pbasis%packed = .TRUE.
          IF (basis_temp%With_L) basis_temp%tab_Pbasis(i)%Pbasis%With_L = .TRUE.
          CALL RecRead5_Basis(basis_temp%tab_Pbasis(i)%Pbasis,mole)

          basis_temp%opt_param = basis_temp%opt_param +                 &
                              basis_temp%tab_Pbasis(i)%Pbasis%opt_param

          basis_temp%active = basis_temp%active .AND. basis_temp%tab_Pbasis(i)%Pbasis%active
          IF (.NOT. basis_temp%active) EXIT
        END DO

        IF (basis_temp%active) THEN
          basis_temp%ndim = 0
          DO i=1,basis_temp%nb_basis
            basis_temp%ndim = basis_temp%ndim +                         &
                            basis_temp%tab_Pbasis(i)%Pbasis%ndim
          END DO


          CALL alloc_init_basis(basis_temp)

          ! transfert of iQdyn, nrho, --------------------------------
          i0 = 0
          DO i=1,basis_temp%nb_basis
            IF (basis_temp%tab_Pbasis(i)%Pbasis%ndim == 0) CYCLE
            i1 = i0 + basis_temp%tab_Pbasis(i)%Pbasis%ndim
            basis_temp%iQdyn(i0+1:i1) = basis_temp%tab_Pbasis(i)%Pbasis%iQdyn(:)
            basis_temp%nrho(i0+1:i1)  = basis_temp%tab_Pbasis(i)%Pbasis%nrho(:)
            i0 = i1
          END DO
        END IF

      END IF

      IF (basis_temp%active) THEN

        CALL alloc_NParray(basis_temp%Tabder_Qdyn_TO_Qbasis,[mole%nb_var], &
                        'basis_temp%Tabder_Qdyn_TO_Qbasis',name_sub,[0])
        basis_temp%Tabder_Qdyn_TO_Qbasis(:) = 0
        ! now basis_temp%iQdyn(:) and Tabder_Qdyn_TO_Qbasis(:) are set-up
        DO i=1,basis_temp%ndim
          basis_temp%Tabder_Qdyn_TO_Qbasis(basis_temp%iQdyn(i)) = i
        END DO
        IF (print_level > 1) write(out_unit,*) 'Tabder_Qdyn_TO_Qbasis',&
                        basis_temp%Tabder_Qdyn_TO_Qbasis(1:mole%nb_var)

      ELSE
        CALL dealloc_basis(basis_temp)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE RecRead5_Basis
!================================================================
! ++    read basis_nD
!
!================================================================
      SUBROUTINE read5_basis_nD(basis_temp,mole)
      USE EVR_system_m
      use mod_Constant
      use mod_Coord_KEO, only: CoordType
      USE mod_basis
      IMPLICIT NONE

      TYPE (basis),     intent(inout)  :: basis_temp
      TYPE (CoordType), intent(inout)  :: mole

!----- for the active basis set ---------------------------------------
      integer              :: nb,nq,nq_extra,nbc,nqc
      integer              :: Nested,nq_max_Nested
      integer              :: SparseGrid_type
      logical              :: SparseGrid,With_L
      logical              :: SparseGrid_With_Cuba    ! When 2 or more are true, the program choses the optimal one
      logical              :: SparseGrid_With_Smolyak ! When only one is true, the program tries to use only one
      logical              :: SparseGrid_With_DP      ! Remark: when only SparseGrid_With_Cuba=T, and the grid does not exit the program stops
      integer              :: L_TO_n_type,max_nb,max_nq
      integer              :: L_SparseGrid,L_TO_nq_A,L_TO_nq_B,L_TO_nq_C,Lexpo_TO_nq
      integer              :: Li_SparseGrid(mole%nb_act1),Li_SparseBasis(mole%nb_act1)
      integer              :: L1_SparseGrid,L2_SparseGrid,Num_OF_Lmax
      integer              :: Type_OF_nDindB,MaxCoupling_OF_nDindB,nDinit_OF_nDindB
      real (kind=Rkind)    :: Norm_OF_nDindB,weight_OF_nDindB
      integer              :: L_SparseBasis,L_TO_nb_A,L_TO_nb_B,Lexpo_TO_nb
      integer              :: L1_SparseBasis,L2_SparseBasis
      integer, allocatable :: Tab_L_TO_n(:)
      logical              :: read_L_TO_n
      logical              :: dnBBRep
      integer              :: ndim,nb_basis

      logical              :: packed

      logical              :: contrac,contrac_analysis,contrac_RVecOnly
      logical              :: read_contrac_file
      logical              :: contrac_WITH_nDindB
      logical              :: auto_basis,auto_contrac,POGridRep,POGridRep_polyortho
      TYPE (REAL_WU)       :: max_ene_contrac
      integer              :: max_nbc,min_nbc,nqPLUSnbc_TO_nqc
      integer              :: auto_contrac_type1_TO,auto_contrac_type21_TO
      character (len=Line_len) :: name_contrac_file

      logical              :: restart_make_cubature,make_cubature,read_para_cubature
      logical              :: cplx,BasisEl

      integer              :: iQact(mole%nb_act1)
      integer              :: iQdyn(mole%nb_act1)
      integer              :: symab(3),index_symab(3)

      character (len=Name_len)     :: name
      character (len=Name_longlen) :: dummy_name

      real (kind=Rkind) :: cte(20,mole%nb_act1)
      real (kind=Rkind) :: A(mole%nb_act1),B(mole%nb_act1)
      real (kind=Rkind) :: Q0(mole%nb_act1),scaleQ(mole%nb_act1)
      real (kind=Rkind) :: k_HO(mole%nb_act1),m_HO(mole%nb_act1),G_HO(mole%nb_act1)
      integer           :: opt_A(mole%nb_act1),opt_B(mole%nb_act1)
      integer           :: opt_Q0(mole%nb_act1),opt_scaleQ(mole%nb_act1)
      logical           :: TD_Q0(mole%nb_act1),TD_scaleQ(mole%nb_act1)

      integer       :: i,nbLi

      NAMELIST /basis_nD/iQact,iQdyn,name,                                      &
                         nb,nq,nq_extra,                                        &
                         nbc,nqc,contrac,contrac_analysis,contrac_RVecOnly,     &
                         cte,cplx,                                              &
                         auto_basis,A,B,opt_A,opt_B,                            &
                         Q0,scaleQ,opt_Q0,opt_scaleQ,k_HO,m_HO,G_HO,            &
                         TD_Q0,TD_scaleQ,                                       &
                         symab,index_symab,                                     &
                         L_TO_n_type,                                           &
                         L_SparseGrid,L_TO_nq_A,L_TO_nq_B,L_TO_nq_C,            &
                         L1_SparseGrid,L2_SparseGrid,Num_OF_Lmax,               &
                         Lexpo_TO_nq,Lexpo_TO_nb,max_nb,max_nq,                 &
                         L_SparseBasis,L_TO_nb_A,L_TO_nb_B,read_L_TO_n,         &
                         L1_SparseBasis,L2_SparseBasis,                         &
                         Li_SparseGrid,Li_SparseBasis,                          &
                         SparseGrid,SparseGrid_type,With_L,                     &
                         SparseGrid_With_Cuba,SparseGrid_With_Smolyak,          &
                         SparseGrid_With_DP,                                    &
                         Nested,nq_max_Nested,                                  &
                         Type_OF_nDindB,Norm_OF_nDindB,weight_OF_nDindB,        &
                         MaxCoupling_OF_nDindB,nDinit_OF_nDindB,contrac_WITH_nDindB,   &
                         packed,dnBBRep,                                        &
                         name_contrac_file,auto_contrac,                        &
                         make_cubature,restart_make_cubature,read_para_cubature,&
                         POGridRep_polyortho,POGridRep,nb_basis,                &
                         max_ene_contrac,max_nbc,min_nbc,nqPLUSnbc_TO_nqc,      &
                         auto_contrac_type1_TO,auto_contrac_type21_TO

!----- for debuging --------------------------------------------------
      integer :: err_io
      character (len=*), parameter :: name_sub='read5_basis_nD'
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!------- read the namelist --------------------------------------
      name                     = "0"
      packed                   = .FALSE.
      dnBBRep                  = .FALSE.

      contrac                  = .FALSE.
      contrac_analysis         = .FALSE.
      contrac_RVecOnly         = .FALSE.

      make_cubature            = .FALSE.
      restart_make_cubature    = .FALSE.
      read_para_cubature       = .FALSE.
      auto_contrac             = .FALSE.
      POGridRep                = .FALSE.
      POGridRep_polyortho      = .FALSE.
      max_nbc                  = 0
      min_nbc                  = 0
      auto_contrac_type1_TO    = 100
      auto_contrac_type21_TO   = 200
      name_contrac_file        = " "
      max_ene_contrac          = REAL_WU(TEN**4,'cm-1','E') ! 10000 cm-1
      cplx                     = .FALSE.

      cte(:,:)                 = ZERO
      A(:)                     = ZERO
      B(:)                     = ZERO
      Q0(:)                    = ZERO
      scaleQ(:)                = ZERO
      k_HO(:)                  = ZERO
      m_HO(:)                  = ZERO
      G_HO(:)                  = ZERO

      opt_A(:)                 = 0
      opt_B(:)                 = 0
      opt_Q0(:)                = 0
      opt_scaleQ(:)            = 0
      auto_basis               = .FALSE.

      TD_Q0(:)                 = .FALSE.
      TD_scaleQ(:)             = .FALSE.

      iQact(:)                 = 0
      iQdyn(:)                 = 0
      nb                       = 0
      nq                       = 0
      nq_extra                 = 0
      Nested                   = 0
      nq_max_Nested            = -1
      nbc                      = -1
      nqc                      = -1
      symab(:)                 = -1
      index_symab(:)           = -1
      nqPLUSnbc_TO_nqc         = 2
      nb_basis                 = 0

      L_TO_n_type              = 0   ! for both nq and nb
      max_nq                   = huge(1)
      L_TO_nq_A                = 1   ! nq(L) = L_TO_nq_A + L_TO_nq_B * L**expo
      L_TO_nq_B                = 1   ! nq(L) = L_TO_nq_A + L_TO_nq_B * L**expo
      Lexpo_TO_nq              = 1
      L_TO_nq_C                = 0
      L_SparseGrid             = -1
      L1_SparseGrid            = huge(1)
      L2_SparseGrid            = huge(1)
      Li_SparseGrid(:)         = huge(1)
      Li_SparseBasis(:)        = huge(1)
      Num_OF_Lmax              = 0

      max_nb                   = huge(1)
      L_TO_nb_A                = -1   ! nb(L) = L_TO_nb_A + L_TO_nb_B * L**expo
      L_TO_nb_B                = -1   ! nb(L) = L_TO_nb_A + L_TO_nb_B * L**expo
      Lexpo_TO_nb              = -1
      read_L_TO_n              = .FALSE.
      L_SparseBasis            = -1
      L1_SparseBasis           = huge(1)
      L2_SparseBasis           = huge(1)

      SparseGrid               = .FALSE.
      SparseGrid_type          = -1
      With_L                   = basis_temp%With_L  ! because, it can be set up in RecRead5_Basis

      SparseGrid_With_Cuba     = .TRUE.
      SparseGrid_With_Smolyak  = .TRUE.
      SparseGrid_With_DP       = .TRUE.
      Type_OF_nDindB           = 1
      Norm_OF_nDindB           = huge(ONE)
      weight_OF_nDindB         = ONE
      nDinit_OF_nDindB         = 1
      MaxCoupling_OF_nDindB    = -1
      contrac_WITH_nDindB      = .FALSE.

      read(in_unit,basis_nD,IOSTAT=err_io)
      IF (err_io < 0) THEN
        write(out_unit,basis_nD)
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the namelist "basis_nD"'
        write(out_unit,*) ' end of file or end of record'
        write(out_unit,*) ' Probably, you forget a basis set ...'
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF
      IF (err_io > 0) THEN
        write(out_unit,basis_nD)
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the namelist "basis_nD"'
        write(out_unit,*) ' Probably, some arguments of namelist are wrong.'
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF

      IF (print_level > 1 .OR. debug) THEN
         IF(MPI_id==0) write(out_unit,basis_nD)
      ELSE
         IF(MPI_id==0) write(out_unit,*) 'Basis name : ',name
      END IF

      IF (count(cte /= ZERO) > 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Do not use cte(:,:) to define the scaling or the range parameters'
        write(out_unit,*) ' Instead, use A and B or Q0 and scaleQ'
        write(out_unit,*) ' CHECK your data'
        STOP
      END IF
      packed = packed .OR. basis_temp%packed .OR. (nb_basis < 1)
      packed = packed .OR. ((contrac .OR. auto_contrac) .AND. .NOT. (contrac_RVecOnly .OR. contrac_analysis))

      ! Here only iQact(:) will be set-up, although iQdyn are read from the data
      ndim = count(iQact(:) > 0)
      IF (ndim == 0) THEN
        ndim = count(iQdyn(:) > 0)
        DO i=1,ndim
          iQact(i) = mole%ActiveTransfo%list_QdynTOQact(iQdyn(i))
        END DO
      END IF
      ! Here only iQact(:) is set-up
      basis_temp%active = .TRUE.
      DO i=1,ndim
        basis_temp%active = basis_temp%active .AND. (iQact(i) <= mole%nb_act1)
      END DO

      IF (.NOT. basis_temp%active) RETURN

      basis_temp%name                   = name

      CALL string_uppercase_TO_lowercase(basis_temp%name)

      IF (trim(adjustl(basis_temp%name)) == 'el') THEN
        BasisEl = .TRUE.
      ELSE
        BasisEl = .FALSE.
      END IF


      IF (ndim <= 0 .AND. nb_basis == 0 .AND. .NOT. BasisEl) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' The primitive basis has no coordinates !'
        write(out_unit,*) ' Specified the active coordinates with iQact(:) or iQdyn(:)'
        write(out_unit,*) ' or "nb_basis" for direct-product basis'
        write(out_unit,*) ' STOP in ',name_sub
        write(out_unit,basis_nD)
        STOP
      END IF

      IF (auto_contrac_type1_TO /= 0 .AND. auto_contrac_type1_TO /= 100) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' The possible values for "auto_contrac_type1_TO" are: '
         write(out_unit,*) ' "0" or "100"'
         write(out_unit,*) ' Your value is:',auto_contrac_type1_TO
         STOP
      END IF
      IF (auto_contrac_type21_TO /= 0 .AND. auto_contrac_type21_TO /= 100 .AND. &
          auto_contrac_type21_TO /= 20 .AND. auto_contrac_type21_TO /= 200) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' the possible values for "auto_contrac_type21_TO" are:'
         write(out_unit,*) ' "0", "100", "20" or "200"'
         write(out_unit,*) ' Your value is:',auto_contrac_type21_TO
         STOP
      END IF

      read_contrac_file = (name_contrac_file /= " ")

      IF (POGridRep_polyortho) POGridRep = .TRUE.
      IF (.NOT. POGridRep) nqPLUSnbc_TO_nqc = 0
      IF (nqPLUSnbc_TO_nqc < 0) nqPLUSnbc_TO_nqc = 0

      basis_temp%ndim                   = ndim
      basis_temp%nb                     = nb
      basis_temp%nbc                    = nbc
      basis_temp%nqc                    = nqc
      basis_temp%nqPLUSnbc_TO_nqc       = nqPLUSnbc_TO_nqc
      basis_temp%nq_max_Nested          = nq_max_Nested
      basis_temp%Nested                 = Nested
      CALL Set_nq_OF_basis(basis_temp,nq)

      IF (L_SparseBasis /= -1 .AND. Norm_OF_nDindB == huge(ONE)) THEN
        Type_OF_nDindB = 0
        Norm_OF_nDindB = real(L_SparseBasis,kind=Rkind)
      END IF
      IF (Norm_OF_nDindB /= huge(ONE) .AND. L_SparseBasis == -1) THEN
        L_SparseBasis = int(Norm_OF_nDindB)
      END IF

      IF (SparseGrid .AND. SparseGrid_type == -1) THEN
        SparseGrid_type = 1 !with SG (old way)
      ELSE IF (SparseGrid .AND. SparseGrid_type == 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  SparseGrid=T (with SG) and SparseGrid_type=0 (without SG)'
        write(out_unit,*) '  Do not use the "SparseGrid" variable, ...'
        write(out_unit,*) '  ... use only the "SparseGrid_type"'
        write(out_unit,*) '  Check your data!!'
        STOP
      ELSE IF (.NOT. SparseGrid .AND. SparseGrid_type == -1) THEN
        SparseGrid_type = 0 ! without SG
      END IF
      SparseGrid = (SparseGrid_type > 0)

      basis_temp%SparseGrid_type             = SparseGrid_type
      basis_temp%L_SparseGrid                = L_SparseGrid
      basis_temp%L_SparseBasis               = L_SparseBasis

      IF (L1_SparseGrid  == huge(1) .AND. L1_SparseBasis < huge(1)) L1_SparseGrid  = L1_SparseBasis
      IF (L1_SparseBasis == huge(1) .AND. L1_SparseGrid  < huge(1)) L1_SparseBasis = L1_SparseGrid
      IF (L1_SparseBasis > L1_SparseGrid) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  L1_SparseBasis > L1_SparseGrid : it is impossible'
          write(out_unit,*) '  L1_SparseBasis,L1_SparseGrid',L1_SparseBasis,L1_SparseGrid
          write(out_unit,*) '  Check your data'
          STOP
      END IF
      basis_temp%para_SGType2%L1_SparseGrid  = L1_SparseGrid
      basis_temp%para_SGType2%L1_SparseBasis = L1_SparseBasis

      IF (L2_SparseGrid  == huge(1) .AND. L2_SparseBasis < huge(1)) L2_SparseGrid  = L2_SparseBasis
      IF (L2_SparseBasis == huge(1) .AND. L2_SparseGrid  < huge(1)) L2_SparseBasis = L2_SparseGrid
      IF (L2_SparseBasis > L2_SparseGrid) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  L2_SparseBasis > L2_SparseGrid : it is impossible'
          write(out_unit,*) '  L2_SparseBasis,L2_SparseGrid',L2_SparseBasis,L2_SparseGrid
          STOP
      END IF
      basis_temp%para_SGType2%L2_SparseGrid  = L2_SparseGrid
      basis_temp%para_SGType2%L2_SparseBasis = L2_SparseBasis

      !write(out_unit,*) 'Li_SparseGrid',Li_SparseGrid
      IF (count(Li_SparseGrid < huge(1)) > 0 .AND. (L1_SparseGrid < huge(1) .OR. L2_SparseGrid < huge(1))) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  Li_SparseGrid and L1_SparseGrid or L2_SparseGrid are defined'
        write(out_unit,*) '  Li_SparseGrid',Li_SparseGrid
        write(out_unit,*) '  L1_SparseGrid,L2_SparseGrid',L1_SparseGrid,L2_SparseGrid
        write(out_unit,*) '  Only Li_SparseGrid OR (L1_SparseGrid or L2_SparseGrid) can be defined'
        write(out_unit,*) '  Check your data'
        STOP
      ELSE IF (count(Li_SparseGrid < huge(1)) == 0) THEN
        Li_SparseGrid(1:2) = [L1_SparseGrid,L2_SparseGrid]
      END IF

      !write(out_unit,*) 'Li_SparseBasis',Li_SparseBasis
      IF (count(Li_SparseBasis < huge(1)) > 0 .AND. (L1_SparseBasis < huge(1) .OR. L2_SparseBasis < huge(1))) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  Li_SparseBasis and L1_SparseBasis or L2_SparseBasis are defined'
        write(out_unit,*) '  Li_SparseBasis',Li_SparseBasis
        write(out_unit,*) '  L1_SparseBasis,L2_SparseBasis',L1_SparseBasis,L2_SparseBasis
        write(out_unit,*) '  Only Li_SparseBasis OR (L1_SparseBasis or L2_SparseBasis) can be defined'
        write(out_unit,*) '  Check your data'
        STOP
      ELSE IF (count(Li_SparseBasis < huge(1)) == 0) THEN
        Li_SparseBasis(1:2) = [L1_SparseBasis,L2_SparseBasis]
      END IF

      IF (count(Li_SparseBasis < huge(1)) /= count(Li_SparseGrid < huge(1))) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  The number of defined Li in Li_SparseBasis(:) and Li_SparseBasis(:) must be the same'
        nbLi = count(Li_SparseBasis < huge(1))
        write(out_unit,*) '  Li_SparseBasis',Li_SparseBasis(1:nbLi)
        nbLi = count(Li_SparseGrid < huge(1))
        write(out_unit,*) '  Li_SparseGrid',Li_SparseGrid(1:nbLi)
        write(out_unit,*) '  Check your data'
        STOP
      END IF

      nbLi = count(Li_SparseBasis < huge(1))
      basis_temp%para_SGType2%Li_SparseBasis = Li_SparseBasis(1:nbLi)
      nbLi = count(Li_SparseGrid < huge(1))
      basis_temp%para_SGType2%Li_SparseGrid  = Li_SparseGrid(1:nbLi)

      !write(out_unit,*) 'Li_SparseBasis',Li_SparseBasis
      !write(out_unit,*) 'Li_SparseGrid ',Li_SparseGrid
      write(out_unit,*) 'para_SGType2%Li_SparseBasis',basis_temp%para_SGType2%Li_SparseBasis
      write(out_unit,*) 'para_SGType2%Li_SparseGrid ',basis_temp%para_SGType2%Li_SparseGrid

      !IF (Num_OF_Lmax < 0 .OR. Num_OF_Lmax > 2) Num_OF_Lmax = 0 ! the condition (Num_OF_Lmax > 2) is not possible anymore
      IF (Num_OF_Lmax < 0) Num_OF_Lmax = 0
      basis_temp%para_SGType2%Num_OF_Lmax    = Num_OF_Lmax

      IF (max_nb > max_nq ) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  max_nb MUST be <= max_nq'
          write(out_unit,*) '  max_nb, max_nq',max_nb, max_nq
          STOP
      END IF


        IF (read_L_TO_n) THEN
          CALL alloc_NParray(Tab_L_TO_n,[10],"Tab_L_TO_n",name_sub,[0])

          read(in_unit,*,IOSTAT=err_io) dummy_name,Tab_L_TO_n
          IF (err_io /= 0) THEN
            write(out_unit,basis_nD)
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  while reading  "Tab_L_TO_n" for L_TO_nb'
            write(out_unit,*) ' Probably, you some integers are missing ...'
            write(out_unit,*) ' => The line has to be like that (with 11 integers):'
            write(out_unit,*) ' l_to_nb 1 2 3 4 5 6   6 6 6 6 6'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF
          write(out_unit,'(2a,11(1x,i0))') 'dummy_name ',dummy_name,Tab_L_TO_n
          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nb,max_n=max_nb,Tab_L_TO_n=Tab_L_TO_n)

          read(in_unit,*,IOSTAT=err_io) dummy_name,Tab_L_TO_n
          IF (err_io /= 0) THEN
            write(out_unit,basis_nD)
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  while reading "Tab_L_TO_n" for L_TO_nq'
            write(out_unit,*) ' Probably, you some integers are missing ...'
            write(out_unit,*) ' => The line has to like that (with 11 integers):'
            write(out_unit,*) ' l_to_nq 1 2 3 4 5 6   6 6 6 6 6'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF
          write(out_unit,'(2a,11(1x,i0))') 'dummy_name ',dummy_name,Tab_L_TO_n
          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nq,max_n=max_nq,Tab_L_TO_n=Tab_L_TO_n)

          CALL dealloc_NParray(Tab_L_TO_n,"Tab_L_TO_n",name_sub)

        ELSE
          IF (L_TO_nb_A   == -1) L_TO_nb_A   = L_TO_nq_A
          IF (L_TO_nb_B   == -1) L_TO_nb_B   = L_TO_nq_B
          IF (Lexpo_TO_nb == -1) Lexpo_TO_nb = Lexpo_TO_nq

          IF (L_TO_nb_B > L_TO_nq_B .OR. L_TO_nb_A > L_TO_nq_A .OR. Lexpo_TO_nb > Lexpo_TO_nq) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  L_TO_nb_B MUST be <= L_TO_nq_B'
            write(out_unit,*) '  L_TO_nb_B,L_TO_nq_B',L_TO_nb_B,L_TO_nq_B

            write(out_unit,*) '       and'

            write(out_unit,*) '  L_TO_nb_A MUST be <= L_TO_nq_A'
            write(out_unit,*) '  L_TO_nb_A,L_TO_nq_A',L_TO_nb_A,L_TO_nq_A

            write(out_unit,*) '       and'

            write(out_unit,*) '  Lexpo_TO_nb MUST be <= Lexpo_TO_nq'
            write(out_unit,*) '  Lexpo_TO_nb,Lexpo_TO_nq',Lexpo_TO_nb,Lexpo_TO_nq

            STOP
          END IF
          IF (L_TO_nq_C == 0) L_TO_nq_C = L_TO_nq_B

          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nb,L_TO_nb_A,L_TO_nb_B, &
                  expo=Lexpo_TO_nb,max_n=max_nb,L_TO_n_type=L_TO_n_type)
          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nq,L_TO_nq_A,L_TO_nq_B, &
             L_TO_nq_C,Lexpo_TO_nq,max_n=max_nq,L_TO_n_type=L_TO_n_type)


        END IF

        basis_temp%SparseGrid_With_Cuba    = SparseGrid_With_Cuba
        basis_temp%SparseGrid_With_Smolyak = SparseGrid_With_Smolyak
        basis_temp%SparseGrid_With_DP      = SparseGrid_With_DP

        basis_temp%With_L                  = With_L
        IF (SparseGrid_type == 2) basis_temp%With_L = .TRUE.
        IF (SparseGrid_type == 4) basis_temp%With_L = .TRUE.

        IF (SparseGrid_type == 1 .OR. SparseGrid_type == 2 .OR. SparseGrid_type == 4) SGtype = SparseGrid_type


        IF (print_level > 1) THEN
          write(out_unit,*) ' Parameters: nq(L) = A + B * L**expo'
          CALL Write_Basis_L_TO_n(basis_temp%L_TO_nq)
          write(out_unit,*) ' Parameters: nb(L) = A + B * L**expo'
          CALL Write_Basis_L_TO_n(basis_temp%L_TO_nb)
        END IF


      basis_temp%Type_OF_nDindB           = Type_OF_nDindB
      basis_temp%Norm_OF_nDindB           = Norm_OF_nDindB
      basis_temp%weight_OF_nDindB         = weight_OF_nDindB
      basis_temp%nDinit_OF_nDindB         = nDinit_OF_nDindB
      basis_temp%contrac_WITH_nDindB      = contrac_WITH_nDindB
      IF (nb_basis == 0) MaxCoupling_OF_nDindB = 1
      IF (nb_basis > 0 .AND. MaxCoupling_OF_nDindB < 1) MaxCoupling_OF_nDindB = nb_basis
      basis_temp%MaxCoupling_OF_nDindB    = MaxCoupling_OF_nDindB

      basis_temp%dnBBRep                  = dnBBRep
      basis_temp%packed                   = packed
      basis_temp%contrac                  = contrac .OR. auto_contrac .OR. contrac_analysis .OR. contrac_RVecOnly
      basis_temp%contrac_analysis         = contrac_analysis
      basis_temp%contrac_RVecOnly         = contrac_RVecOnly
      basis_temp%auto_contrac             = auto_contrac
      basis_temp%make_cubature            = make_cubature
      basis_temp%restart_make_cubature    = restart_make_cubature
      basis_temp%read_para_cubature       = read_para_cubature
      basis_temp%max_ene_contrac          = convRWU_TO_R_WITH_WorkingUnit(max_ene_contrac)
      basis_temp%max_nbc                  = max_nbc
      basis_temp%min_nbc                  = min_nbc
      basis_temp%auto_contrac_type1_TO    = auto_contrac_type1_TO
      basis_temp%auto_contrac_type21_TO   = auto_contrac_type21_TO
      basis_temp%POGridRep                = POGridRep
      basis_temp%POGridRep_polyortho      = POGridRep_polyortho
      basis_temp%read_contrac_file        = read_contrac_file
      basis_temp%file_contrac%name        = name_contrac_file
      basis_temp%nb_basis                 = nb_basis


      CALL Set_ReadsymabOFSymAbelian(basis_temp%P_SymAbelian,symab(1))

      IF (ndim == 0) THEN
         IF (.NOT. associated(basis_temp%nDindB)) THEN
           CALL alloc_array(basis_temp%nDindB,"basis_temp%nDindB",name_sub)
         END IF
         basis_temp%nDindB%packed          = basis_temp%packed .OR. (Type_OF_nDindB == 0)
      ELSE
        CALL alloc_init_basis(basis_temp)
        basis_temp%nDindB%packed          = basis_temp%packed .OR. (Type_OF_nDindB == 0)
        ! now basis_temp%iQdyn(:)
        DO i=1,ndim
          basis_temp%iQdyn(i) = mole%ActiveTransfo%list_QactTOQdyn(iQact(i))
          IF (basis_temp%iQdyn(i) > mole%nb_var) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' iQdyn(i) is larger than nb_var',               &
                                    i,basis_temp%iQdyn(i),mole%nb_var
            STOP
          END IF
        END DO

        IF (print_level > 1) THEN
          write(out_unit,*) 'auto_basis:    ',auto_basis
          write(out_unit,*) 'A,B,Q0,scaleQ: ',A(1:ndim),B(1:ndim),Q0(1:ndim),scaleQ(1:ndim)
          write(out_unit,*) 'opt of A,B,Q0,scaleQ: ',opt_A(1:ndim),opt_B(1:ndim),opt_Q0(1:ndim),opt_scaleQ(1:ndim)
        END IF
        basis_temp%A(:)      = ZERO
        basis_temp%B(:)      = ZERO
        basis_temp%Q0(:)     = ZERO
        basis_temp%scaleQ(:) = ONE

        DO i=1,ndim
          IF (scaleQ(i) == ZERO .AND. k_HO(i) > ZERO .AND. &
              (m_HO(i) > ZERO .OR. G_HO(i) > ZERO)) THEN
            IF (m_HO(i) > ZERO) scaleQ(i) = sqrt(sqrt(k_HO(i)*m_HO(i)))
            IF (G_HO(i) > ZERO) scaleQ(i) = sqrt(sqrt(k_HO(i)/G_HO(i)))
            write(out_unit,*) ' scaleQ(i)',i,scaleQ(i)
          END IF

          IF (A(i) /= B(i) .AND. scaleQ(i) > ZERO) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' You give the range ("A" and "B"):',A(i),B(i)
            write(out_unit,*) '   and also the scaling factors ("scaleQ" and "Q0"):',scaleQ(i),Q0(i)
            write(out_unit,*) ' You have to chose the range or the scaling factors'
            write(out_unit,*) ' CHECK your data'
            write(out_unit,basis_nD)
            STOP
          ELSE IF (A(i) > B(i) .AND. scaleQ(i) == ZERO) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' The range ("A" and "B") is :',A(i),B(i)
            write(out_unit,*) '   "A" MUST be < "B" '
            write(out_unit,*) ' CHECK your data'
            write(out_unit,basis_nD)
            STOP
          ELSE IF (A(i) /= B(i) .AND. scaleQ(i) == ZERO) THEN
            basis_temp%A(i)      = A(i)
            basis_temp%B(i)      = B(i)
            basis_temp%Q0(i)     = ZERO
            basis_temp%scaleQ(i) = ONE
            basis_temp%opt_A(i)  = opt_A(i)
            basis_temp%opt_B(i)  = opt_B(i)
          ELSE IF (A(i) == B(i) .AND. scaleQ(i) > ZERO) THEN
            basis_temp%A(i)          = ZERO
            basis_temp%B(i)          = ZERO
            basis_temp%Q0(i)         = Q0(i)
            basis_temp%scaleQ(i)     = scaleQ(i)
            basis_temp%opt_Q0(i)     = opt_Q0(i)
            basis_temp%opt_scaleQ(i) = opt_scaleQ(i)
          ELSE IF (A(i) == B(i) .AND. scaleQ(i) == ZERO .AND. Q0(i) /= ZERO) THEN
            basis_temp%A(i)          = ZERO
            basis_temp%B(i)          = ZERO
            basis_temp%Q0(i)         = Q0(i)
            basis_temp%scaleQ(i)     = ONE
            basis_temp%opt_Q0(i)     = opt_Q0(i)
            basis_temp%opt_scaleQ(i) = opt_scaleQ(i)
          ELSE
            basis_temp%A(i)      = ZERO
            basis_temp%B(i)      = ZERO
            basis_temp%Q0(i)     = ZERO
            basis_temp%scaleQ(i) = ONE
          END IF
          basis_temp%TD_Q0(i)     = TD_Q0(i)
          basis_temp%TD_scaleQ(i) = TD_scaleQ(i)
        END DO

        basis_temp%auto_basis = (auto_basis .AND. ndim == 1)
        IF (auto_basis .AND. ndim == 1 .AND. associated(mole%NMTransfo)) THEN
          IF (.NOT. mole%tab_Qtransfo(mole%itNM)%skip_transfo) THEN
          IF (associated(mole%NMTransfo%Q0_HObasis) .AND. associated(mole%NMTransfo%scaleQ_HObasis)) THEN

            basis_temp%Q0(1)     = mole%NMTransfo%Q0_HObasis(basis_temp%iQdyn(1))
            basis_temp%scaleQ(1) = mole%NMTransfo%scaleQ_HObasis(basis_temp%iQdyn(1))
            write(out_unit,*) 'Q0,scaleQ (from auto_basis): ',basis_temp%Q0,basis_temp%scaleQ
            basis_temp%auto_basis = .FALSE.
          ELSE
            write(out_unit,*) ' WARNING in ',name_sub
            write(out_unit,*) '  auto_basis=t and ...'
            write(out_unit,*) '   NMTransfo%Q0_HObasis or NMTransfo%scaleQ_HObasis are not associated'
            write(out_unit,*) '      The Q0 and scaleQ are not modified !!'
            write(out_unit,*) ' CHECK your data'
            !write(out_unit,basis_nD)
            !STOP
          END IF
          END IF
        ELSE IF (auto_basis .AND. ndim == 1 .AND. .NOT. associated(mole%NMTransfo)) THEN
          write(out_unit,*) ' WARNING in ',name_sub
          write(out_unit,*) '  auto_basis=t and mole%NMTransfo is not associated'
          write(out_unit,*) ' CHECK your data'
          !write(out_unit,basis_nD)
          !STOP
        END IF

        basis_temp%opt_param = count(basis_temp%opt_A /= 0) +           &
          count(basis_temp%opt_B /= 0) + count(basis_temp%opt_Q0 /= 0) +&
                                       count(basis_temp%opt_scaleQ /= 0)
        IF (print_level > 1)                                            &
                  write(out_unit,*) 'Parameter(s) to be optimized?: ', &
                                                   basis_temp%opt_param

      END IF


      basis_temp%nq_extra = nq_extra
      IF (basis_temp%nq_extra > 0 .AND. ndim > 0) THEN
        write(out_unit,*) 'Read nq_extra (',basis_temp%nq_extra,') grid points'
        CALL alloc_NParray(basis_temp%x_extra,[ndim,basis_temp%nq_extra],      &
                          'basis_temp%x_extra',name_sub)
        DO i=1,basis_temp%nq_extra
          read(in_unit,*)  dummy_name,basis_temp%x_extra(:,i)
        END DO
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)
!---------------------------------------------------------------------
      end subroutine read5_basis_nD
!================================================================
!       analysis of a string
!       output : nb_word,word (the i th word of name)
!================================================================
      SUBROUTINE analysis_name(name,word,i,nb_word)
      USE EVR_system_m
      IMPLICIT NONE


!     - analysis of the basis name --------------------------
      character (len=*)  :: word,name
      character          :: ch
      integer            :: i,nb_word
      integer            :: iw,ic,icw
      logical            :: blank


      iw = 0
      icw = 0
      blank = .TRUE.
!     write(out_unit,*) 'analysis_name: ',name,len(name)
      DO ic=1,len(name)
        ch = name(ic:ic)
        IF (ch .EQ. " ") THEN
          IF (.NOT. blank) THEN
            iw = iw + 1
            blank = .TRUE.
          END IF
        ELSE
          IF (iw .EQ. i-1) THEN
            icw = icw + 1
            word(icw:icw) = ch
          END IF
          blank = .FALSE.
        END IF
!       write(out_unit,*) 'analysis_name: ',ic,ch,blank,iw
      END DO

      nb_word = iw
!     write(out_unit,*) 'analysis_name: ',name,':',nb_word
!     write(out_unit,*) 'analysis_name: ',i,word


      end subroutine analysis_name
