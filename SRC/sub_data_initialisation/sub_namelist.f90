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
! ++    read inactive parameters (for cHAC or HADA methods)
!
!       nb_quadrature : grid quadrature points
!       max_excit     : (or max_herm) max degree for the hermite polynomias
!       max_ene_h     : energy limit (in cm-1 for the imput then in au)
!       num           : .TRUE. if numerical derivatives
!       step          : step for numerical derivatives
!       auTOcm_inv    : conversion factor hartree (au) to cm-1
!
!================================================================
!
      SUBROUTINE read_inactive(para_AllBasis,mole,QMLib_in)
      USE EVR_system_m
      USE mod_nDindex
      USE mod_Constant,  only: REAL_WU, convRWU_TO_R_WITH_WorkingUnit
      use mod_Coord_KEO, only: CoordType, alloc_array, dealloc_array,   &
                               set_rphtransfo, Tnum, Init_degenerate_freq
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis), intent(inout) :: para_AllBasis
      TYPE (CoordType),      intent(inout) :: mole
      logical,               intent(in)    :: QMLib_in

!-----------------------------------------------------------

!----- working variables -----------------------------------

      integer           :: nb_quadrature,max_excit,max_coupling,n_h
      real (kind=Rkind) :: step,max_ene_h_val
      logical           :: num,ADA
      TYPE (REAL_WU)    :: max_ene_h



      integer, parameter :: max_inact2n = 500
      integer       :: i,k,nb_inact21
      integer       :: tab_nq(max_inact2n)
      integer       :: tab_nb(max_inact2n)
      integer       :: nDinit(max_inact2n)
      logical       :: SparseGrid
      integer       :: isort,L_SparseGrid

      logical       :: H0_sym,gradTOpot0,diabatic_freq
      integer       :: Qinact2n_sym(max_inact2n)
      integer       :: Qinact2n_eq(max_inact2n,max_inact2n)

      logical       :: non_adia,contrac_ba_ON_HAC
      integer       :: max_nb_ba_ON_HAC
      logical       :: QMLib

      NAMELIST /inactives/nb_quadrature,max_excit,max_coupling,         &
                          isort,                                        &
                          max_ene_h,num,step,n_h,                       &
                          ADA,tab_nq,tab_nb,                            &
                          SparseGrid,L_SparseGrid,                      &
                          H0_sym,diabatic_freq,Qinact2n_sym,Qinact2n_eq,&
                          gradTOpot0,QMLib,                             &
                          non_adia,contrac_ba_ON_HAC,max_nb_ba_ON_HAC


!------- read the inactive namelist ----------------------------
      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
      character (len=*), parameter :: name_sub='read_inactive'
!      -----------------------------------------------------------------
      write(out_unit,*) 'INACTIVES PARAMETERS'
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'QMLib_in',QMLib_in
      END IF

      ! get nb_inact21 from mole%nb_inact2n or mole%RPHTransfo%nb_inact21
      IF (associated(mole%RPHTransfo)) THEN
        nb_inact21 = mole%RPHTransfo%nb_inact21
      ELSE
        nb_inact21 = mole%nb_inact2n
      END IF

      IF (mole%nb_inact2n <= 0) THEN
        read(in_unit,inactives) ! for nagfor we must read the inactive namelist
        RETURN
      END IF
      IF (nb_inact21 > max_inact2n) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'max_inact2n is too small',max_inact2n
        write(out_unit,*) 'It should be larger than "nb_inact21"',nb_inact21
        STOP
      END IF

      tab_nq(:) = 0
      tab_nb(:) = 0

      Qinact2n_sym(:)  = 0
      Qinact2n_eq(:,:) = 0
      H0_sym           = .FALSE.
      diabatic_freq    = .FALSE.
      gradTOpot0       = .FALSE.

      step             = ONETENTH**5
      num              = .TRUE. ! not used anymore
      ADA              = .FALSE.

      max_ene_h        = REAL_WU(huge(ONE),   'cm-1','E') ! huge (cm-1)
      nb_quadrature    = 10
      max_excit        = -1
      max_coupling     = mole%nb_inact2n
      n_h              = -1 ! all channels
      SparseGrid       = .FALSE.
      L_SparseGrid     = -1
      isort            = -1 ! 1: sort the HADA basis (energy)
                            ! 2: sort the HADA basis in term of excitation

      contrac_ba_ON_HAC = .FALSE.
      max_nb_ba_ON_HAC  = huge(1)

      QMLib             = QMLib_in

      read(in_unit,inactives)

      IF (step < epsilon(ONE)*TEN**3) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'step is too small',step
        write(out_unit,*) 'It should be larger than',epsilon(ONE)*TEN**3
        STOP
      END IF


      IF (.NOT. associated(mole%RPHTransfo)) THEN

        max_ene_h_val = convRWU_TO_R_WITH_WorkingUnit(max_ene_h)

        IF (max_ene_h_val == huge(ONE) .AND. max_excit  > -1) isort = 2
        IF (max_ene_h_val /= huge(ONE) .AND. max_excit == -1) isort = 1
        IF (isort == -1) isort = 1


        IF (product(tab_nq(1:nb_inact21)) == 0)  tab_nq(:) = nb_quadrature

        IF (product(tab_nb(1:nb_inact21)) == 0 ) tab_nb(:) = max_excit + 1
        ! here tab_nb can be 0 if max_excit=-1, therefore, we test again
        IF (product(tab_nb(1:nb_inact21)) == 0 ) tab_nb(:) = 1


        IF (max_excit == -1) THEN
          max_excit = sum(tab_nb(1:nb_inact21))-nb_inact21
        END IF

        write(out_unit,*) 'max_excit',max_excit
        write(out_unit,*) 'tab_nb(:)',tab_nb(1:nb_inact21)

        nDinit(:) = 0
        CALL alloc_array(para_AllBasis%Basis2n%nDindB,                   &
                        'para_AllBasis%Basis2n%nDindB',name_sub)
        IF (isort == 1) THEN ! sort with energy
          CALL init_nDindexPrim(para_AllBasis%Basis2n%nDindB,nb_inact21, &
                      nDsize=tab_nb(1:nb_inact21),nDinit=nDinit(1:nb_inact21), &
                      MaxNorm=convRWU_TO_R_WITH_WorkingUnit(max_ene_h),  &
                      type_OF_nDindex=0,With_nDindex=.FALSE.)
        ELSE IF (isort == 2) THEN ! sort with excitation
          CALL init_nDindexPrim(para_AllBasis%Basis2n%nDindB,nb_inact21,                 &
                                nDsize=tab_nb(1:nb_inact21),nDinit=nDinit(1:nb_inact21), &
                                Lmax=max_excit,type_OF_nDindex=0,With_nDindex=.FALSE.)
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '   the isort value ',isort,' cannot be used.'
          write(out_unit,*) '   check your data!!'
          STOP
        END IF
        para_AllBasis%Basis2n%nDindB%Max_nDI = n_h
        !CALL Write_nDindex(para_AllBasis%Basis2n%nDindB)

        CALL init_nDindexPrim(para_AllBasis%Basis2n%nDindG,nb_inact21,  &
                              nDsize=tab_nq(1:nb_inact21),type_OF_nDindex=0,With_nDindex=.FALSE.)
        !CALL Write_nDindex(para_AllBasis%Basis2n%nDindG)

        para_AllBasis%basis_ext2n%ADA               = ADA
        para_AllBasis%basis_ext2n%contrac_ba_ON_HAC = contrac_ba_ON_HAC
        para_AllBasis%basis_ext2n%max_nb_ba_ON_HAC  = max_nb_ba_ON_HAC
        para_AllBasis%basis_ext2n%max_ene_ON_HAC    = convRWU_TO_R_WITH_WorkingUnit(max_ene_h)


        para_AllBasis%Basis2n%nb_basis            = nb_inact21
        IF (SparseGrid) THEN
          para_AllBasis%Basis2n%SparseGrid_type   = 3
        ELSE
          para_AllBasis%Basis2n%SparseGrid_type   = 0
        END IF
        para_AllBasis%Basis2n%L_SparseGrid        = L_SparseGrid

        CALL alloc_tab_Pbasis_OF_basis(para_AllBasis%Basis2n)

        IF (SparseGrid .AND. isort /= 2) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' With SparseGrid, you MUST use isort=2'
          write(out_unit,*) '   SparseGrid,isort: ',SparseGrid,isort
          write(out_unit,*) 'Check your data!'
          STOP
        END IF

        IF (para_AllBasis%Basis2n%SparseGrid_type == 3 .AND.            &
            para_AllBasis%Basis2n%L_SparseGrid < 1) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' SparseGrid_type /= 3 and L_SparseGrid < 1'
          write(out_unit,*) ' You should increase L_SparseGrid!'
          STOP
        END IF
      END IF

      IF (nb_inact21 > 0 .AND. .NOT. associated(mole%RPHTransfo)) THEN
        ! here it is for the HADA/cHAC models (not RPH)
        IF (associated(mole%RPHTransfo_inact2n)) THEN
          CALL dealloc_array(mole%RPHTransfo_inact2n,                   &
                            'mole%RPHTransfo_inact2n',name_sub)
        END IF

        CALL alloc_array(mole%RPHTransfo_inact2n,                       &
                        'mole%RPHTransfo_inact2n',name_sub)

        CALL Set_RPHTransfo(mole%RPHTransfo_inact2n,                    &
                            mole%ActiveTransfo%list_act_OF_Qdyn,        &
                            gradTOpot0,diabatic_freq,step,              &
                            H0_sym,H0_sym,Qinact2n_sym(1:nb_inact21),   &
                            Qinact2n_eq(1:nb_inact21,1:nb_inact21),     &
                            QMLib=QMLib,list_QMLMapping=mole%ActiveTransfo%list_QMLMapping)

        CALL Init_degenerate_freq(mole%RPHTransfo_inact2n%degenerate_freq,nb_inact21)

      ELSE IF (nb_inact21 > 0 .AND. associated(mole%RPHTransfo)) THEN
        STOP 'In read_inactive: this option is not used anymore ?'

        IF (mole%RPHTransfo%option == 0) THEN
          ! here we do not set up the list_act_OF_Qdyn because it is already
          !     done in type_var_analysis

          CALL Set_RPHTransfo(mole%RPHTransfo,                          &
            gradTOpot0=gradTOpot0,diabatic_freq=diabatic_freq,step=step,&
                                    purify_hess=H0_sym,eq_hess=H0_sym,  &
                              Qinact2n_sym=Qinact2n_sym(1:nb_inact21),  &
                     Qinact2n_eq=Qinact2n_eq(1:nb_inact21,1:nb_inact21),&
                     QMLib=QMLib,list_QMLMapping=mole%ActiveTransfo%list_QMLMapping)
        END IF
      END IF


      IF (debug) THEN
        CALL Write_nDindex(para_AllBasis%Basis2n%nDindB,'para_AllBasis%Basis2n%nDindB')
        CALL Write_nDindex(para_AllBasis%Basis2n%nDindG,'para_AllBasis%Basis2n%nDindG')
        write(out_unit,*) 'END ',name_sub
      END IF

      END SUBROUTINE read_inactive
!
!================================================================
! ++    read the active parameters
!
!       max_ene1       : energy limit (in cm-1 for the imput then in au)
!       test           : .TRUE.    if a test on ONE active grid point
!       auTOcm_inv     : conversion factor hartree (au) to cm-1
!
!================================================================
!

