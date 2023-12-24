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
      SUBROUTINE vib(max_mem,test_mem,intensity_only)
      USE mod_system
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis

      USE mod_psi

      USE mod_propa
      USE mod_FullPropa
      USE mod_FullControl
      USE mod_Davidson
      USE mod_Filter
      USE mod_Arpack
      USE mod_CRP
      USE mod_Op
      USE mod_analysis
      USE mod_fullanalysis
      USE mod_Auto_Basis
      USE mod_MPI_aux
      IMPLICIT NONE

!---------------------------------------------------------------------------------------
!     variables
!---------------------------------------------------------------------------------------
!----- variables for the dynamical memory allocation -----------------------------------
      logical  :: intensity_only
      integer  :: nio_res_int

      integer (kind=ILkind)  :: max_mem
      logical   test_mem

!----- physical and mathematical constants ---------------------------------------------
      TYPE (constant)  :: const_phys

!----- for the CoordType and Tnum --------------------------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

!----- variables for the construction of H ---------------------------------------------
      TYPE (param_AllOp), target  :: para_AllOp
      TYPE (param_Op),    pointer :: para_H      => null()
      TYPE (param_Op),    pointer :: para_S      => null()
      TYPE (param_Op),    pointer :: para_Dip(:) => null()
      integer                     :: iOp
      real (kind=Rkind)           :: max_Sii,max_Sij

!----- variables pour la namelist analyse ----------------------------------------------
      TYPE (param_ana)            :: para_ana
      TYPE (param_intensity)      :: para_intensity
      TYPE (param_ana_psi)        :: ana_WP0

!----- variables for the WP propagation ------------------------------------------------
      TYPE (param_propa)          :: para_propa
      TYPE (param_psi)            :: WP0(1),WP0tmp,MuWP0
      complex (kind=Rkind)        :: s,c

!----- for Davidson diagonalization ----------------------------------------------------
      integer                        :: nb_diago
      integer                        :: max_diago
      TYPE (param_psi), pointer      :: Tab_Psi(:) => null()
      real (kind=Rkind),allocatable  :: Ene0(:)

!----- variables for the namelist actives ----------------------------------------------
!----- for the basis set ---------------------------------------------------------------
      TYPE (param_AllBasis) :: para_AllBasis

      integer :: Get_nbPERsym_FROM_SymAbelianOFAllBasis ! function

!----- variables divers ----------------------------------------------------------------
      integer           :: i,ip,i_baie,f_baie,id,nb_ScalOp
      integer           :: nb_baie_sym,nb_bi,nb_be,nb_ba,nb_bie
      real (kind=Rkind) :: T,DE,Ep,Em,Q,fac,zpe,pop
      logical           :: print_mat
      integer           :: err
      integer           :: err_mem,memory
      real (kind=Rkind) :: part_func ! function


!para_mem%mem_debug=.TRUE.


!---------------------------------------------------------------------------------------
!=====================================================================
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: initialization of variables'

      ! Disassociate for initialization
      nullify(para_H)
      nullify(para_S)
      nullify(para_Dip)

      write(out_unitp,*) ' END VIB: initialization of variables'
      write(out_unitp,*) '================================================='

      write(out_unitp,*)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: BEGINNING ini_data'
      CALL write_mem_file('ini_data ini')
      CALL time_perso('ini_data ini')
      write(out_unitp,*)

      !CALL system_mem_usage(memory_RSS,'before ini_data')

      CALL     ini_data(const_phys,para_Tnum,mole,                      &
                        para_AllBasis,para_AllOp,                       &
                        para_ana,para_intensity,intensity_only,         &
                        para_propa)

      !CALL system_mem_usage(memory_RSS,'after ini_data')
      para_H => para_AllOp%tab_Op(1)

      write(out_unitp,*)
      CALL time_perso('ini_data end')
      CALL write_mem_file('ini_data end')
      write(out_unitp,*)
      write(out_unitp,*) ' VIB: END ini_data'
      write(out_unitp,*) '================================================='
      write(out_unitp,*)
!---------------------------------------------------------------------------------------

      IF (intensity_only) GOTO 100

!---------------------------------------------------------------------------------------
!      Grids (V, T, Dip) calculations
!---------------------------------------------------------------------------------------
      !> turn off the allocation of grid, to be done in action
      !> sub_qa_bhe should be refined later
      !para_AllOp%tab_Op(1)%Grid_save_ac
!#if(run_MPI)
!      Grid_allco=.FALSE.
!#endif
      write(out_unitp,*)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: BEGINNING sub_qa_bhe'
      CALL write_mem_file('sub_qa_bhe ini')
      CALL time_perso('sub_qa_bhe ini')
      write(out_unitp,*)

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Read_FileGrid) THEN ! test only for H
        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*) ' The grid of operators (S V Veff T1 and T2) will be read'
          write(out_unitp,*)
        ENDIF
        IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid) THEN
          IF(MPI_id==0) write(out_unitp,*) 'Test_Grid=.TRUE. => STOP in vib'
          STOP
        END IF
      END IF

      IF (para_ana%VibRot) THEN
         para_Tnum%JJ = para_ana%JJmax
         para_Tnum%With_Cart_Transfo = (para_Tnum%JJ>0) .AND. mole%Cart_transfo
      END IF

      ! calculate potential in action
      CALL sub_qa_bhe(para_AllOp)

      IF (para_ana%VibRot) THEN
         para_Tnum%JJ = 0
         para_Tnum%With_Cart_Transfo = (para_Tnum%JJ>0) .AND. mole%Cart_transfo
      END IF

      write(out_unitp,*)
      CALL time_perso('sub_qa_bhe end')
      CALL write_mem_file('sub_qa_bhe end')
      write(out_unitp,*) ' VIB: END sub_qa_bhe'
      write(out_unitp,*) '================================================='
      write(out_unitp,*)

      IF (para_EVRT_calc%Grid_only) STOP 'EVRT stop: Grid_only=.TRUE.'
!#if(run_MPI)
!      Grid_allco=.True.
!#endif
!---------------------------------------------------------------------------------------
!      contraction of the active basis set with HADA basis
!---------------------------------------------------------------------------------------
      IF (para_H%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC .AND. para_H%nb_bi>1) THEN
        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING HADA contraction',          &
                                   para_AllBasis%BasisnD%nb
          CALL write_mem_file('HADA contraction ini')
          CALL time_perso('HADA contraction ini')
          write(out_unitp,*)
          write(out_unitp,*)
        ENDIF

        CALL sub_MatOp_HADA(para_H,para_ana,para_intensity,para_AllOp,  &
                            const_phys)

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('HADA contraction end')
          CALL write_mem_file('HADA contraction end')
          write(out_unitp,*) ' VIB: END HADA contraction'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF

        para_AllOp%tab_Op(:)%nb_tot = para_H%para_ReadOp%nb_elec *      &
                   sum(para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(:))
        para_AllOp%tab_Op(:)%nb_tot_ini = para_AllOp%tab_Op(:)%nb_tot
      END IF

!---------------------------------------------------------------------------------------
!     propagation:  Test if propagation will be done
!---------------------------------------------------------------------------------------

      IF (para_ana%propa) THEN
!---------------------------------------------------------------------------------------
!       => Time-dependent calculation
!---------------------------------------------------------------------------------------
        write(out_unitp,*) 'Propogation start'
        CALL write_mem_file('psi0 ini')

        CALL init_psi(WP0(1),para_H,cplx=.TRUE.)
        write(out_unitp,*) 'Propogation initialized'

        ! building of WP0 for WP propagation ------------------------------------------
        IF (.NOT. para_ana%control .AND. para_propa%type_WPpropa /=100) THEN
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: Generate Psi0',WP0%nb_tot
            CALL time_perso('psi0 ini')
            write(out_unitp,*)
          ENDIF

          IF ( .NOT. abs(para_propa%type_WPpropa) == 33 .AND.           &
               .NOT. abs(para_propa%type_WPpropa) == 34) THEN

            !CALL alloc_psi(WP0(1),BasisRep=.TRUE.,GridRep=.FALSE.)
            !WP0(1)%CvecB(:) = CZERO
            !WP0(1)%CvecB(WP0(1)%nb_ba+1) = CONE

            CALL init_psi0(WP0,para_propa%para_WP0,mole)

            ip = para_propa%para_WP0%WP0_dip
            IF (ip > 0 .AND. ip < 4) THEN
              iOp = 2
              para_Dip => para_AllOp%tab_Op(iOp+1:iOp+3)

              write(out_unitp,*) ' calculation Mu|WP0> with: ',para_Dip(ip)%name_Op
              write(out_unitp,*) ' th_WP0: ',para_propa%para_WP0%th_WP0
              flush(out_unitp)

              WP0tmp = WP0(1)
              CALL sub_OpPsi(WP0tmp,MuWP0,para_Dip(ip))
              MuWP0  = cos(para_propa%para_WP0%th_WP0)*MuWP0
              WP0tmp = sin(para_propa%para_WP0%th_WP0)*WP0tmp
              WP0(1) = MuWP0 + WP0tmp

              CALL dealloc_psi(MuWP0)
              CALL dealloc_psi(WP0tmp)

            END IF ! for ip > 0 .AND. ip <4

            IF (para_propa%para_WP0%WP0_nb_CleanChannel > 0) THEN
              write(out_unitp,*) ' clean channel of |WP0>'
              DO i=1,para_propa%para_WP0%WP0_nb_CleanChannel
                i_baie = 1 +                                            &
                (para_propa%para_WP0%WP0_CleanChannellist(i)-1)*WP0(1)%nb_ba
                f_baie = i_baie-1 + WP0(1)%nb_ba
                IF (WP0(1)%cplx) THEN
                  WP0(1)%CvecB(i_baie:f_baie) = ZERO
                ELSE
                  WP0(1)%RvecB(i_baie:f_baie) = ZERO
                END IF
              END DO
            END IF ! for para_propa%para_WP0%WP0_nb_CleanChannel > 0

            IF(keep_MPI) CALL norm2_psi(WP0(1))
            CALL MPI_Bcast_(WP0(1)%norm2,size1_MPI,root_MPI)

            IF(keep_MPI) THEN
              write(out_unitp,*) ' Norm^2 of |WP0>',WP0(1)%norm2
              CALL renorm_psi_With_norm2(WP0(1))
              write(out_unitp,*) ' Analysis of |WP0> or Mu|WP0>'
              write(out_unitp,*)
            ENDIF

            IF(keep_MPI) THEN
              ana_WP0               = para_propa%ana_psi
              CALL modif_ana_psi(ana_WP0,Ene=ZERO,T=ZERO,ZPE=ZERO)

              ana_WP0%file_Psi%name = trim(para_propa%file_WP%name) // '_WP0'

              CALL sub_analyze_tab_psi(WP0(:),ana_WP0,adia=.FALSE.,Write_psi=.FALSE.)
              CALL dealloc_ana_psi(ana_WP0)
            END IF
            ! spectral tranformation cannot be done here,
            !   because the matrix representation is not done yet

            write(out_unitp,*)

          ELSE ! for optimal control

            write(out_unitp,*) ' WP0 will be read later!'

          END IF ! for .NOT. abs(para_propa%type_WPpropa) == 33

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            CALL time_perso('psi0 end')
            CALL write_mem_file('psi0 end')
            write(out_unitp,*) ' VIB: END Generate Psi0'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
          flush(out_unitp)
        END IF ! .NOT. para_ana%control .AND. para_propa%type_WPpropa /=100

        !================================================================
        !================================================================
        !================================================================
        !===== Tune the number of threads (for SG4) =====================
        !================================================================
        ! for only one WP (complex)
        CALL write_mem_file('Tune_SG4threads_HPsi end')
        CALL Tune_SG4threads_HPsi(.TRUE.,1,para_H)
        CALL write_mem_file('Tune_SG4threads_HPsi end')

        !================================================================
        !================================================================
        !================================================================
        !===== build S and/or H if necessary ============================
        !================================================================
        IF (para_H%Make_Mat .AND. .NOT. para_H%spectral) THEN

          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING sub_matOp',para_H%nb_tot
          CALL write_mem_file('sub_matOp: H and S ini')
          CALL time_perso('sub_matOp: H and S ini')
          write(out_unitp,*) 'para_S...%comput_S',para_AllOp%tab_Op(2)%para_ReadOp%comput_S
          write(out_unitp,*)


          IF (para_AllOp%tab_Op(2)%para_ReadOp%comput_S) THEN
            para_S => para_AllOp%tab_Op(2)
            CALL sub_MatOp(para_S,para_ana%print)

            !> analysis of the overlap matrix
            CALL sub_ana_S(para_S%Rmat,para_S%nb_tot,max_Sii,max_Sij,.TRUE.)

            CALL dealloc_para_Op(para_S)
            nullify(para_S)
          END IF

          CALL sub_MatOp(para_H,para_ana%print)

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_matOp: H and S end')
            CALL write_mem_file('sub_matOp: H and S end')
            write(out_unitp,*) ' VIB: END sub_matOp'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        ELSE IF (para_H%Make_Mat .AND. para_H%spectral) THEN
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING Spectral Representation',para_H%nb_tot
            CALL write_mem_file('Spectral: All Op ini')
            CALL time_perso('Spectral: All Op ini')
            write(out_unitp,*)
          ENDIF

          DO i=1,size(para_AllOp%tab_Op)
            IF (i == 2) CYCLE
            CALL sub_MatOp(para_AllOp%tab_Op(i),para_ana%print)
          END DO

          IF (para_AllOp%tab_Op(2)%para_ReadOp%comput_S) THEN

            para_S => para_AllOp%tab_Op(2)
            CALL sub_MatOp(para_S,para_ana%print)

            !- analysis of the overlap matrix
            CALL sub_ana_S(para_S%Rmat,para_S%nb_tot,max_Sii,max_Sij,.TRUE.)
            CALL dealloc_para_Op(para_S)
            nullify(para_S)

          END IF

          IF (para_H%Spectral .AND. para_H%spectral_done .AND. .NOT. para_ana%control) THEN
            write(out_unitp,*) 'WP0 spectral representation'
            DO i=1,size(WP0)
              IF (para_H%cplx) THEN
                ! 1st: project psi on the spectral basis
                WP0(i)%CvecB(:) = matmul(transpose(para_H%Cvp),WP0(i)%CvecB)
              ELSE
                ! 1st: project psi on the spectral basis
                WP0(i)%CvecB(:) = matmul(transpose(para_H%Rvp),WP0(i)%CvecB)
              END IF
            END DO
          END IF

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            CALL time_perso('Spectral: All Op')
            CALL write_mem_file('Spectral: All Op end')
            write(out_unitp,*) ' VIB: END Spectral Representation'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        END IF ! IF (para_H%Make_Mat .AND. .NOT. para_H%spectral)
        !================================================================
        !================================================================

        !================================================================
        !===== Hmax calculation
        !================================================================
        write(out_unitp,*)
        write(out_unitp,*) '================================================'
        write(out_unitp,*) ' VIB: Hmax and Hmin calculation'
        CALL write_mem_file('sub_Hmax ini')
        CALL time_perso('sub_Hmax ini')

        write(out_unitp,*)
        write(out_unitp,*)

        IF ( .NOT. abs(para_propa%type_WPpropa) == 33 .AND.           &
             .NOT. abs(para_propa%type_WPpropa) == 34 .AND.           &
             .NOT. abs(para_propa%type_WPpropa) == 100 ) THEN

!         - Hmax and Hmin calculation ---------------------------------
          !IF(.NOT. openmpi) CALL sub_Hmax(para_propa,para_H)
          IF(.NOT. (para_H%para_ReadOp%para_FileGrid%Type_FileGrid==4))                &
            CALL sub_Hmax(para_propa,para_H)

        ELSE
          write(out_unitp,*) ' Calculation of Hmax is skiped'
        ENDIF

        write(out_unitp,*)
        CALL time_perso('sub_Hmax end')
        CALL write_mem_file('sub_Hmax end')
        write(out_unitp,*) ' VIB: END Hmax and Hmin calculation'
        write(out_unitp,*) '================================================'
        write(out_unitp,*)
        !================================================================
        !================================================================

        !================================================================
        !==== propagation or control
        !================================================================
        IF (para_ana%control) THEN
!         - control -----------------------------------------------
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================'
            write(out_unitp,*) ' VIB: control'
            write(out_unitp,*) ' VIB: propagation'
            CALL write_mem_file('sub_control ini')
            CALL time_perso('sub_control ini')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          CALL sub_Opt_control(para_AllOp,para_propa)
!         CALL sub_nonOpt_control(para_AllOp,para_propa)

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_control end')
            CALL write_mem_file('sub_control end')
            write(out_unitp,*) ' VIB: END control'
            write(out_unitp,*) '================================================'
            write(out_unitp,*)
          ENDIF
        ELSE
!         - propagation -----------------------------------------------
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================'
            write(out_unitp,*) ' VIB: propagation'
            CALL write_mem_file('sub_propagation ini')
            CALL time_perso('sub_propagation ini')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          CALL sub_propagation(WP0,para_AllOp,para_propa)

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_propagation end')
            CALL write_mem_file('sub_propagation end')
            write(out_unitp,*) ' VIB: END propagation'
            write(out_unitp,*) '================================================'
            write(out_unitp,*)
          ENDIF
        END IF
        !================================================================
        !================================================================

!=====================================================================
!       END Time-dependent calculation
!=====================================================================
      ELSE ! for not para_ana%propa
!=====================================================================
!       => Time-independent calculation
!=====================================================================

        !================================================================
        !================================================================
        !================================================================
        !===== Tune the number of threads (for SG4) =====================
        !================================================================
        CALL write_mem_file('Tune_SG4threads_HPsi ini')
        max_diago = max(10,para_propa%para_Davidson%nb_WP,para_H%nb_tot/10)
        max_diago = min(max_diago,10,para_H%nb_tot)
        CALL Tune_SG4threads_HPsi(para_H%cplx,max_diago,para_H)
        CALL write_mem_file('Tune_SG4threads_HPsi end')

        !================================================================
        !===== build S and/or H if necessary ============================
        !================================================================
        IF (para_H%Make_Mat) THEN
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING sub_matOp',para_H%nb_tot
            CALL write_mem_file('sub_matOp: H and S ini')
            CALL time_perso('sub_matOp: H and S')
            write(out_unitp,*)
            write(out_unitp,*) 'para_S...%comput_S',para_AllOp%tab_Op(2)%para_ReadOp%comput_S
            write(out_unitp,*)
          ENDIF

          IF (para_AllOp%tab_Op(2)%para_ReadOp%comput_S) THEN

            para_S => para_AllOp%tab_Op(2)
            CALL sub_MatOp(para_S,para_ana%print)

            !- analysis of the overlap matrix
            CALL sub_ana_S(para_S%Rmat,para_S%nb_tot,max_Sii,max_Sij,.TRUE.)

            CALL dealloc_para_Op(para_S)
            nullify(para_S)
          END IF ! for para_AllOp%tab_Op(2)%para_ReadOp%comput_S

          CALL sub_MatOp(para_H,para_ana%print)

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_matOp: H and S')
            CALL write_mem_file('sub_matOp: H and S end')
            write(out_unitp,*) ' VIB: END sub_matOp'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        END IF ! for para_H%Make_Mat
        !================================================================
        !================================================================

        !================================================================
        !===== Hmax calculation (for filter diagonalization)
        !================================================================
        IF (para_ana%filter .AND. para_propa%auto_Hmax) THEN
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================'
            write(out_unitp,*) ' VIB: Hmax and Hmin calculation'
            CALL write_mem_file('sub_Hmax ini2')
            CALL time_perso('sub_Hmax ini2')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          !IF(.NOT. openmpi) CALL sub_Hmax(para_propa,para_H)
          IF(.NOT. (para_H%para_ReadOp%para_FileGrid%Type_FileGrid==4))                &
            CALL sub_Hmax(para_propa,para_H)

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_Hmax end2')
            CALL write_mem_file('sub_Hmax end2')
            write(out_unitp,*) ' VIB: END Hmax and Hmin calculation'
            write(out_unitp,*) '================================================'
            write(out_unitp,*)
          ENDIF
        END IF ! for para_ana%filter .AND. para_propa%auto_Hmax

        !================================================================
        !================================================================

        IF (.NOT. para_H%cplx .AND. para_ana%CRP > 0) THEN
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================'
            write(out_unitp,*) ' VIB: BEGINNING sub_CRP',para_H%nb_tot
            CALL time_perso('sub_CRP ini')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          CALL sub_CRP(para_AllOp%tab_Op,size(para_AllOp%tab_Op),para_ana%print,para_ana%para_CRP)

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_CRP end')
            write(out_unitp,*) ' VIB: END sub_CRP'
            write(out_unitp,*) '================================================'
            write(out_unitp,*)
          ENDIF
        END IF ! for .NOT. para_H%cplx .AND. para_ana%CRP > 0
        flush(out_unitp)

        !================================================================
        !===== Diagonalisation ==========================================
        !================================================================
        IF (para_ana%davidson .OR. para_ana%arpack .OR. para_ana%filter) THEN
          IF (para_H%Partial_MatOp) STOP 'STOP the Matrix is incomplete'
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING ITERATIVE DIAGONALIZATION'
            CALL time_perso('sub_Iterative_Diago')
            write(out_unitp,*)
          ENDIF


          IF (para_propa%para_Davidson%max_WP == 0) THEN
            max_diago = max(1000,para_propa%para_Davidson%nb_WP,          &
                            para_H%nb_tot/10)
          ELSE
            max_diago = para_propa%para_Davidson%max_WP
          END IF

          IF (Get_nbPERsym_FROM_SymAbelianOFAllBasis(para_AllBasis,       &
                               para_propa%para_Davidson%symab) == 0) THEN ! (test on -1 ???)
            max_diago = min(max_diago,para_H%nb_tot)
          ELSE
            nb_baie_sym = Get_nbPERsym_FROM_SymAbelianOFAllBasis(       &
                           para_AllBasis,para_propa%para_Davidson%symab)
            nb_be = para_H%para_ReadOp%nb_elec
            nb_bi = get_nb_FROM_basis(para_AllBasis%Basis2n)
            nb_baie_sym = nb_baie_sym * nb_bi*nb_be
            max_diago = min(max_diago,para_H%nb_tot,nb_baie_sym)
          END IF
          para_propa%para_Davidson%max_WP = max_diago

          nb_diago = min(para_propa%para_Davidson%nb_WP,para_H%nb_tot,max_diago)
!#if(run_MPI)
!          CALL MPI_Bcast(nb_diago,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
          IF(openmpi .AND. MPI_scheme/=1) THEN
            CALL MPI_Bcast_(nb_diago,size1_MPI,root_MPI)
            CALL MPI_Bcast_(max_diago,size1_MPI,root_MPI)
          ENDIF

          nullify(Tab_Psi)
          CALL alloc_array(Tab_Psi,[max_diago],"Tab_Psi","vib")
          CALL alloc_NParray(Ene0,[max_diago],"Ene0","vib")

          IF (para_ana%davidson) THEN

            CALL sub_propagation_Davidson(Tab_Psi,Ene0,nb_diago,max_diago,             &
                                          para_H,para_propa%para_Davidson,para_propa)

          ELSE IF (para_ana%arpack) THEN ! arpack=t
            !CALL sub_propagation_Arpack(Tab_Psi,Ene0,nb_diago,max_diago,  &
            !                            para_H,para_propa)
            CALL sub_propagation_Arpack_Sym(Tab_Psi,Ene0,nb_diago,max_diago,&
                                           para_H,para_propa)

          ELSE ! filter diagonalization
            CALL sub_GaussianFilterDiagonalization(Tab_Psi,Ene0,nb_diago,max_diago,&
                                           para_H,para_propa)

            !CALL sub_GaussianFilterDiagonalization_v0(Tab_Psi,Ene0,nb_diago,max_diago,&
            !                               para_H,para_propa)

            !CALL sub_FilterDiagonalization(Tab_Psi,Ene0,nb_diago,max_diago,&
            !                               para_H,para_propa)

          END IF

          CALL dealloc_NParray(Ene0,"Ene0","vib")
          para_H%spectral = .TRUE.

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_Iterative_Diago')
            write(out_unitp,*) ' VIB: END ITERATIVE DIAGONALIZATION'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF

        ELSE IF (para_ana%CRP == 0) THEN ! for para_ana%davidson .OR. para_ana%arpack .OR. para_ana%filter
          IF (para_H%Partial_MatOp) STOP 'STOP the Matrix is incomplete'
          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING DIAGONALIZATION'
            CALL time_perso('sub_diago_H')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          nb_diago  = para_H%nb_tot
          max_diago = para_H%nb_tot

          nullify(Tab_Psi)
          IF(keep_MPI) CALL alloc_array(Tab_Psi,[max_diago],"Tab_Psi","vib")

          IF(keep_MPI) THEN
            DO i=1,max_diago
              CALL init_psi(Tab_psi(i),para_H,para_H%cplx)
            END DO
          ENDIF

          para_H%diago = .TRUE.
          IF(keep_MPI) CALL alloc_para_Op(para_H,Grid=.FALSE.,Mat=.TRUE.)
          IF (para_H%cplx) THEN
            IF(keep_MPI) THEN
              CALL sub_diago_CH(para_H%Cmat,para_H%Cdiag,para_H%Cvp,      &
                                para_H%nb_tot)

              nb_diago = count(abs(para_H%Cdiag(:)-para_H%Cdiag(1))< para_ana%max_ene)
              nb_diago = min(nb_diago,para_ana%max_ana)
              IF (nb_diago < 3) nb_diago = 10
              IF (nb_diago > para_H%nb_tot) nb_diago = para_H%nb_tot
              DO i=1,nb_diago
                CALL alloc_psi(Tab_Psi(i))
                Tab_Psi(i)%CvecB(:)    = para_H%Cvp(:,i)
                Tab_psi(i)%CAvOp       = para_H%Cdiag(i)
                Tab_psi(i)%IndAvOp     = para_H%n_Op  ! it should be 0
                Tab_psi(i)%convAvOp    = .TRUE.
              END DO
              para_H%para_AllBasis%basis_ext%Cvp_spec    => para_H%Cvp
            ENDIF
          ELSE
            IF(keep_MPI) THEN
              CALL sub_diago_H(para_H%Rmat,para_H%Rdiag,para_H%Rvp,       &
                               para_H%nb_tot,para_H%sym_Hamil)

              write(out_unitp,*) 'HMin,HMax (ua)  : ',[minval(para_H%Rdiag),maxval(para_H%Rdiag)]
              write(out_unitp,*) 'HMin,HMax (cm-1): ',   &
                  [minval(para_H%Rdiag),maxval(para_H%Rdiag)]*get_Conv_au_TO_unit('E','cm-1')

              nb_diago = count((para_H%Rdiag(:)-para_H%Rdiag(1))< para_ana%max_ene)
              IF (para_ana%max_ana > 0) nb_diago = min(nb_diago,para_ana%max_ana)
              IF (nb_diago < 3) nb_diago = 10
              IF (nb_diago > para_H%nb_tot) nb_diago = para_H%nb_tot
              DO i=1,nb_diago
                CALL alloc_psi(Tab_Psi(i))
                Tab_Psi(i)%RvecB(:)    = para_H%Rvp(:,i)
                Tab_psi(i)%CAvOp       = para_H%Rdiag(i)
                Tab_psi(i)%IndAvOp     = para_H%n_Op  ! it should be 0
                Tab_psi(i)%convAvOp    = .TRUE.
                CALL Set_symab_OF_psiBasisRep(Tab_psi(i))
              END DO
              para_H%para_AllBasis%basis_ext%Rvp_spec    => para_H%Rvp
            ENDIF ! for MPI_id=0
          END IF ! for para_H%cplx
          para_H%para_AllBasis%basis_ext%nb_vp_spec  = nb_diago

          CALL alloc_NParray(para_AllBasis%basis_ext%liste_spec,        &
                 [nb_diago],"para_AllBasis%basis_ext%liste_spec","vib")
          para_AllBasis%basis_ext%liste_spec(:) = [ (i,i=1,nb_diago) ]

          !CALL sub_analyse(Tab_Psi,nb_diago,para_H,para_ana,             &
          !                 para_intensity,para_AllOp,const_phys)

          !  IF (.NOT. para_H%cplx .AND. para_ana%VibRot) THEN
          !    CALL sub_VibRot(Tab_Psi,para_ana%max_ana,para_H,para_ana)
          !  END IF

          IF (allocated(para_H%Rmat)) THEN
            CALL dealloc_NParray(para_H%Rmat,"para_H%Rmat","vib")
          END IF
          IF (allocated(para_H%Cmat)) THEN
            CALL dealloc_NParray(para_H%Cmat,"para_H%Cmat","vib")
          END IF
          para_ana%max_ana = nb_diago

          IF(MPI_id==0) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_diago_H')
            write(out_unitp,*) ' VIB: END DIAGONALIZATION'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        END IF ! for para_ana%davidson .OR. para_ana%arpack .OR. para_ana%filter
        flush(out_unitp)
        !===============================================================
        !===============================================================


        !===============================================================
        !===============================================================
        IF (para_ana%CRP == 0) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING WAVE FUNCTION ANALYSIS'
          CALL time_perso('sub_analyse ini')
          write(out_unitp,*)

          IF(keep_MPI) THEN
             CALL sub_analyse(Tab_Psi,nb_diago,para_H,                                 &
                              para_ana,para_intensity,para_AllOp,const_phys)
          ENDIF

          flush(out_unitp)

          IF (.NOT. para_H%cplx .AND. para_ana%VibRot) THEN
            CALL sub_VibRot(Tab_Psi,para_ana%max_ana,para_H,para_ana)
          END IF
          flush(out_unitp)

          !===============================================================
          ! Spectral representation of operator
          !===============================================================
          print_mat = (para_ana%MaxWP_TO_Write_MatOp >= nb_diago)
          IF (para_ana%Spectral_ScalOp) THEN
            DO iOp=1,size(para_AllOp%tab_Op)

              IF (para_AllOp%tab_Op(iOp)%n_Op == -1) CYCLE ! S
              IF (para_AllOp%tab_Op(iOp)%spectral) THEN
                write(out_unitp,*) '==========================================='
                write(out_unitp,*) ' Spectral representation of: ',        &
                                   trim(para_AllOp%tab_Op(iOp)%name_Op)

                CALL sub_build_MatOp(Tab_Psi,nb_diago,                     &
                                     para_AllOp%tab_Op(iOp),.TRUE.,print_mat)

                write(out_unitp,*) '==========================================='
                flush(out_unitp)
              END IF
            END DO
          END IF

          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_analyse end')
          write(out_unitp,*) ' VIB: END WAVE FUNCTION ANALYSIS'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        END IF
        !===============================================================
        !===============================================================

        !===============================================================
        ! Deallocation of Tab_Psi
        !===============================================================
        IF(keep_MPI) THEN
        IF (associated(Tab_Psi)) THEN
          DO i=1,size(Tab_Psi)
            CALL dealloc_psi(Tab_Psi(i))
          END DO
          CALL dealloc_array(Tab_Psi,"Tab_Psi","vib")
        END IF
        ENDIF
        !===============================================================
        !===============================================================

!=====================================================================
!       => Time-independent calculation
!=====================================================================
      END IF
!=====================================================================
!=====================================================================


!=====================================================================
!=====================================================================
100   CONTINUE
      IF (.NOT. para_H%cplx .AND. para_ana%intensity) THEN
        write(out_unitp,*)
        write(out_unitp,*) '================================================'
        write(out_unitp,*) ' VIB: BEGINNING sub_intensity',para_H%nb_tot,para_ana%max_ana
        CALL time_perso('sub_intensity')
        write(out_unitp,*)
        write(out_unitp,*)

        !         ----- file for restart or for changed the temprature --------------
        IF(MPI_id==0) THEN
        CALL file_open(para_intensity%file_resart_int,nio_res_int)

        IF (intensity_only) THEN
          read(nio_res_int,*,IOSTAT=err) para_H%nb_tot,para_ana%max_ana
          read(nio_res_int,*,IOSTAT=err)

          para_H%diago = .TRUE.
          CALL alloc_para_Op(para_H)

          CALL Read_Vec(para_H%Rdiag,nio_res_int,5,err)
          IF (err /= 0) THEN
            write(out_unitp,*) 'ERROR in vib'
            write(out_unitp,*) ' reading the vector "para_H%Rdiag"'
            STOP
          END IF

          IF (intensity_only) THEN
            write(out_unitp,*)
            Q =  part_func(para_H%Rdiag,size(para_H%Rdiag),para_ana%Temp,const_phys)
            fac = const_phys%Eh / (const_phys%k * para_ana%Temp)
            zpe = minval(para_H%Rdiag(:))
            write(out_unitp,*) 'population at T, Q',para_ana%Temp,Q
            write(out_unitp,*) 'Energy level (',const_phys%ene_unit,') pop:'
            DO i=1,size(para_H%Rdiag)
              pop = exp(-(para_H%Rdiag(i)-zpe)*fac)
              write(out_unitp,10) i,para_H%Rdiag(i) * const_phys%auTOenergy,&
                    (para_H%Rdiag(i) - zpe) * const_phys%auTOenergy,pop/Q
10             format(i4,x,3f20.5)
            END DO
          END IF

          read(nio_res_int,*,IOSTAT=err)

          CALL Read_Mat(para_H%Rvp,nio_res_int,5,err)
          IF (err /= 0) THEN
            write(out_unitp,*) 'ERROR in vib'
            write(out_unitp,*) ' reading the matrix "para_H%Rvp"'
            STOP
          END IF
        ELSE ! for intensity_only
          write(out_unitp,*) 'write restart file for intensity: ',                     &
                 para_intensity%file_resart_int%name
          flush(out_unitp)
          write(nio_res_int,*) para_H%nb_tot,para_ana%max_ana
          write(nio_res_int,*) 'ene'
          CALL Write_Vec(para_H%Rdiag,nio_res_int,5,Rformat='e30.23')
          write(nio_res_int,*) 'psi'
          flush(out_unitp)
          CALL Write_Mat(para_H%Rvp,nio_res_int,5,Rformat='e30.23')
          flush(nio_res_int)
        END IF ! for intensity_only
        flush(out_unitp)
!         -------------------------------------------------------------------

        iOp = 2
        para_Dip => para_AllOp%tab_Op(iOp+1:iOp+3)
        CALL sub_intensity(para_Dip,                                  &
                           para_ana%print,para_H,para_ana%max_ana,    &
                           para_intensity,const_phys,intensity_only,nio_res_int)
        ENDIF ! MPI_id==0

        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('sub_intensity')
        write(out_unitp,*) ' VIB: END sub_intensity'
        write(out_unitp,*) '================================================'
        write(out_unitp,*)

        IF(MPI_id==0) CALL file_close(para_intensity%file_resart_int)
        nullify(para_Dip)
      END IF !for .NOT. para_H%cplx .AND. para_ana%intensity
      flush(out_unitp)

      IF (.NOT. para_H%cplx .AND. para_ana%Psi_ScalOp) THEN
        write(out_unitp,*)
        write(out_unitp,*) '================================================'
        write(out_unitp,*) ' VIB: BEGINNING sub_AnalysePsy_ScalOp',para_H%nb_tot,para_ana%max_ana
        CALL time_perso('sub_AnalysePsy_ScalOp')
        write(out_unitp,*)
        write(out_unitp,*) 'nb_scalar_Op',para_H%para_ReadOp%nb_scalar_Op

        IF(MPI_id==0) THEN
          iOp = 2
          nb_ScalOp = para_H%para_ReadOp%nb_scalar_Op
          para_Dip => para_AllOp%tab_Op(iOp+1:iOp+nb_ScalOp)
          CALL sub_AnalysePsy_ScalOp(para_Dip,nb_ScalOp,para_H,para_ana%max_ana)
        ENDIF

        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('sub_AnalysePsy_ScalOp')
        write(out_unitp,*) ' VIB: END sub_AnalysePsy_ScalOp'
        write(out_unitp,*) '================================================'
        write(out_unitp,*)

        nullify(para_Dip)
      END IF
      flush(out_unitp)

      IF (.NOT. para_H%cplx .AND. para_ana%NLO) THEN
        write(out_unitp,*)
        write(out_unitp,*) '================================================'
        write(out_unitp,*) ' VIB: BEGINNING sub_NLO',para_H%nb_tot,para_ana%max_ana
        CALL time_perso('sub_NLO')
        write(out_unitp,*)
        write(out_unitp,*)

        IF(MPI_id==0) THEN
          iOp = 2
          para_Dip => para_AllOp%tab_Op(iOp+1:iOp+3)
          CALL sub_NLO(para_Dip,para_ana%print,para_H,para_ana%max_ana, &
                      para_intensity)
        ENDIF

        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('sub_NLO')
        write(out_unitp,*) ' VIB: END sub_NLO'
        write(out_unitp,*) '================================================'
        write(out_unitp,*)

        nullify(para_Dip)
      END IF
      flush(out_unitp)

!=====================================================================
!=====================================================================
!=====================================================================
!       deallocated memories
!=====================================================================
      CALL dealloc_table_at(const_phys%mendeleev)

      CALL dealloc_CoordType(mole)
      IF (associated(para_Tnum%Gref)) THEN
        CALL dealloc_array(para_Tnum%Gref,"para_Tnum%Gref","vib")
      END IF
      CALL dealloc_Tnum(para_Tnum)

      CALL dealloc_para_AllOp(para_AllOp)
      CALL dealloc_para_ana(para_ana)
      CALL dealloc_param_propa(para_propa)
      CALL dealloc_psi(WP0(1))

      IF ( associated(Tab_Psi) ) THEN
        DO i=1,size(Tab_Psi)
          CALL dealloc_psi(Tab_Psi(i))
        END DO
        CALL dealloc_array(Tab_Psi,"Tab_Psi","vib")
      END IF

      CALL dealloc_AllBasis(para_AllBasis)

      write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
      write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
      write(out_unitp,*) '================================================'
      IF(openmpi) THEN
        write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!', ' from ', MPI_id
      ELSE
        write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
      ENDIF
      write(out_unitp,*) '================================================'

      END SUBROUTINE vib
!=======================================================================================
!=======================================================================================
!=======================================================================================
  SUBROUTINE sub_Opt_CAP_basis()
    USE mod_system
    USE mod_Constant
    USE mod_Coord_KEO
    USE mod_PrimOp
    USE mod_basis

    USE mod_psi

    USE mod_propa

    USE mod_CRP
    USE mod_Op
    USE mod_analysis
    USE mod_fullanalysis
    USE mod_Auto_Basis
    USE mod_MPI_aux
    USE mod_Optimization

    IMPLICIT NONE

  !---------------------------------------------------------------------------------------
  !     variables
  !---------------------------------------------------------------------------------------
  !----- variables for the dynamical memory allocation -----------------------------------
        logical  :: intensity_only
        integer  :: nio_res_int
  
        integer (kind=ILkind)  :: max_mem
        logical   test_mem
  
  !----- physical and mathematical constants ---------------------------------------------
        TYPE (constant)  :: const_phys
  
  !----- for the CoordType and Tnum --------------------------------------------------------
        TYPE (CoordType) :: mole,mole_100
        TYPE (Tnum)      :: para_Tnum
  
  !----- variables for the construction of H ---------------------------------------------
        TYPE (param_AllOp), target  :: para_AllOp
        TYPE (param_Op),    pointer :: para_H      => null()
        TYPE (param_Op),    pointer :: para_S      => null()
        TYPE (param_Op),    pointer :: para_Dip(:) => null()
        integer                     :: iOp
        real (kind=Rkind)           :: max_Sii,max_Sij
  
  !----- variables pour la namelist analyse ----------------------------------------------
        TYPE (param_ana)            :: para_ana
        TYPE (param_intensity)      :: para_intensity
  
  !----- variables for the WP propagation ------------------------------------------------
        TYPE (param_propa)          :: para_propa
  
  !----- for Davidson diagonalization ----------------------------------------------------
        integer                        :: nb_diago
        integer                        :: max_diago
        real (kind=Rkind),allocatable  :: Ene0(:)
  
  !----- variables for the namelist actives ----------------------------------------------
  !----- for the basis set ---------------------------------------------------------------
        TYPE (param_AllBasis) :: para_AllBasis
  
        integer :: Get_nbPERsym_FROM_SymAbelianOFAllBasis ! function
  
  !----- variables divers ----------------------------------------------------------------
        integer           :: i,ip,i_baie,f_baie,id,nb_ScalOp
        integer           :: nb_baie_sym,nb_bi,nb_be,nb_ba,nb_bie
        real (kind=Rkind) :: T,DE,Ep,Em,Q,fac,zpe,pop
        logical           :: print_mat
        integer           :: err
        integer           :: err_mem,memory
        real (kind=Rkind) :: part_func ! function
  
  
        integer :: iOp_CAPR,iOp_CAPP,iQs,iQopt,iQdyns,iact
        integer :: nb_Opt
        real (kind=Rkind), allocatable :: Qopt(:)

        real (kind=Rkind), allocatable :: Qact(:)
        TYPE (param_Optimization)      :: para_Optimization

        real (kind=Rkind), allocatable :: QTS(:),QR(:),QP(:)
        real (kind=Rkind), allocatable :: Qdyn_TS(:),Qdyn_P(:),Qdyn_R(:)

        real (kind=Rkind), allocatable :: hess_TS(:,:),hess_R(:,:),hess_P(:,:)
        real (kind=Rkind), allocatable :: d0G_TS(:,:),d0G_R(:,:),d0G_P(:,:)
        real (kind=Rkind) :: Ene_TS,Ene_P,Ene_R,Ene_Asymp,Ene_Col,E,LP,LR,mass
        real (kind=Rkind) :: A,B
        real (kind=Rkind), parameter :: c = 2.62206_Rkind

        real (kind=Rkind) :: Ene
        TYPE (param_dnMatOp) :: dnMatOp(1)

        !para_mem%mem_debug=.TRUE.
     !---------------------------------------------------------------------
        logical,parameter :: debug= .FALSE.
        !logical,parameter :: debug= .TRUE.
        character (len=*), parameter :: name_sub = 'sub_Opt_CAP_basis'
     !---------------------------------------------------------------------

  
  !---------------------------------------------------------------------------------------
  !=====================================================================
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: initialization of variables'
  
        nullify(para_H)
        nullify(para_S)
        nullify(para_Dip)
  
        write(out_unitp,*) ' END VIB: initialization of variables'
        write(out_unitp,*) '================================================='
  
        write(out_unitp,*)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: BEGINNING ini_data'
        CALL write_mem_file('ini_data ini')
        CALL time_perso('ini_data ini')
        write(out_unitp,*)

        CALL ini_data(const_phys,para_Tnum,mole,                      &
                      para_AllBasis,para_AllOp,                       &
                      para_ana,para_intensity,intensity_only,         &
                      para_propa)
  
        !CALL system_mem_usage(memory_RSS,'after ini_data')
        para_H => para_AllOp%tab_Op(1)
  
        write(out_unitp,*)
        CALL time_perso('ini_data end')
        CALL write_mem_file('ini_data end')
        write(out_unitp,*)
        write(out_unitp,*) ' VIB: END ini_data'
        write(out_unitp,*) '================================================='
        write(out_unitp,*)
    !---------------------------------------------------------------------------------------
    CALL alloc_NParray(Qact,   [mole%nb_var],'Qact',   name_sub)
    CALL alloc_NParray(Qdyn_TS,[mole%nb_var],'Qdyn_TS',name_sub)
    CALL alloc_NParray(Qdyn_R ,[mole%nb_var],'Qdyn_R', name_sub)
    CALL alloc_NParray(Qdyn_P ,[mole%nb_var],'Qdyn_P', name_sub)

    !get s index (iQs)), QO+, Q0-, QTS (ref geom)
    iOp_CAPR = para_ana%para_CRP%iOp_CAP_Reactif
    iOp_CAPP = para_ana%para_CRP%iOp_CAP_Product
    IF (iOp_CAPR < 1 .OR. iOp_CAPP < 1) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' iOp_CAP_Reactif or iOp_CAP_Product are not defined',iOp_CAPR,iOp_CAPP
      write(out_unitp,*) ' CRP calculation must be defined and CAPs must but read'
      STOP 'ERROR in sub_Opt_CAP_basis: CAP and CRP are not defined'
    END IF
    IF (para_AllOp%tab_Op(iOp_CAPR)%para_ReadOp%tab_CAP(1)%ind_Q /=  para_AllOp%tab_Op(iOp_CAPR)%para_ReadOp%tab_CAP(2)%ind_Q) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Two different s coordinates:'
      write(out_unitp,*) 'iQs',para_AllOp%tab_Op(iOp_CAPR)%para_ReadOp%tab_CAP(:)%ind_Q
      write(out_unitp,*) '  The feature is not implemented yet!'
      STOP 'ERROR in sub_Opt_CAP_basis: Two different s coordinates are not possible'
    END IF
    iqS = para_AllOp%tab_Op(iOp_CAPR)%para_ReadOp%tab_CAP(1)%ind_Q
    write(out_unitp,*) 'iQs',iqS
    write(out_unitp,*)
    write(out_unitp,*)
    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '======= TS geometry optimisation ================================='
    write(out_unitp,*) '=================================================================='

    CALL Read_param_Optimization(para_Optimization,mole,para_AllBasis%BasisnD,read_nml=.TRUE.)
    IF (debug) CALL Write_param_Optimization(para_Optimization)

    allocate(QTS(para_Optimization%nb_Opt))

    CALL Sub_Optimization(para_AllBasis%BasisnD,para_Tnum,mole,para_H%para_ReadOp%PrimOp_t, &
                          QTS,para_Optimization)

    write(out_unitp,*) 'Q_TS',QTS
    CALL Init_Tab_OF_dnMatOp(dnMatOp,para_Optimization%nb_Opt,para_H%para_ReadOp%PrimOp_t%nb_elec,nderiv=2)
    CALL alloc_NParray(hess_TS,[para_Optimization%nb_Opt,para_Optimization%nb_Opt],'hess',name_sub)
    CALL alloc_NParray(d0G_TS,[para_Optimization%nb_Opt,para_Optimization%nb_Opt],'d0G_TS',name_sub)

    CALL get_Qact0(Qact,mole%ActiveTransfo)
    Qact(1:para_Optimization%nb_Opt) = QTS(:)
    write(out_unitp,*) 'Qact_TS',Qact
    CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn_TS,mole%ActiveTransfo) ! here with mole_100
    write(out_unitp,*) 'Qdyn_TS',Qdyn_TS

    CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_H%para_ReadOp%PrimOp_t)
    CALL Get_Hess_FROM_Tab_OF_dnMatOp(hess_TS,dnMatOp) ! for the ground state
    CALL Write_Mat(hess_TS, out_unitp, 5, info='hess_TS')
    Ene_TS = Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)
    write(out_unitp,*) 'Ene_TS',Ene_TS
    CALL get_d0GG(Qact,para_Tnum,mole,d0G_TS,def=.TRUE.)
    CALL Write_Mat(d0G_TS, out_unitp, 5, info='d0G_TS')
    CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

    write(out_unitp,*)
    write(out_unitp,*)
    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '======= Reactif geometry optimisation ============================'
    write(out_unitp,*) '=================================================================='
    ! optimisation for CAP (iOp_CAPR) (without s)
    CALL dealloc_param_Optimization(para_Optimization)
    iQdyns = mole%liste_QactTOQdyn(iQs)
    mole_100 = mole
    mole_100%ActiveTransfo%list_act_OF_Qdyn(iqS) = 100
    mole_100%opt_Qdyn(iQdyns) = 0
    CALL type_var_analysis_OF_CoordType(mole_100,print_lev=.FALSE.)
    CALL get_Qact0(Qact,mole_100%ActiveTransfo)
    iQs =  mole_100%liste_QactTOQdyn(iQdyns)

    para_Optimization%Optimization_method = 'bfgs'
    para_Optimization%Optimization_param = 'coordbasis'
    para_Optimization%para_BFGS%nb_neg = 0
    para_Optimization%para_BFGS%calc_hessian_always = .TRUE.

    CALL Read_param_Optimization(para_Optimization,mole_100,para_AllBasis%BasisnD,read_nml=.FALSE.)
    IF (debug) CALL Write_param_Optimization(para_Optimization)

    Qact(iqS) = para_AllOp%tab_Op(iOp_CAPR)%para_ReadOp%tab_CAP(1)%Q0
    mole_100%ActiveTransfo%Qact0(iqS)    = para_AllOp%tab_Op(iOp_CAPR)%para_ReadOp%tab_CAP(1)%Q0
    mole_100%ActiveTransfo%Qdyn0(iQdyns) = para_AllOp%tab_Op(iOp_CAPR)%para_ReadOp%tab_CAP(1)%Q0
    QR = Qact(1:mole_100%nb_act)

    CALL Sub_Optimization(para_AllBasis%BasisnD,para_Tnum,mole_100,para_H%para_ReadOp%PrimOp_t, &
                          QR,para_Optimization)
    write(out_unitp,*) 'QR',QR
    Qact(1:mole_100%nb_act) = QR(:)
    CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn_R,mole_100%ActiveTransfo) ! here with mole_100
    CALL Qdyn_TO_Qact_FROM_ActiveTransfo(Qdyn_R,Qact,mole%ActiveTransfo)  !    then with mole
    write(out_unitp,*) 'Qact_R',Qact
    write(out_unitp,*) 'Qdyn_R',Qdyn_R

    CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,para_H%para_ReadOp%PrimOp_t%nb_elec,nderiv=2)
    CALL alloc_NParray(hess_R,[mole%nb_act,mole%nb_act],'hess_R',name_sub)
    CALL alloc_NParray(d0G_R, [mole%nb_act,mole%nb_act],'d0G_R',name_sub)

    CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_H%para_ReadOp%PrimOp_t)
    CALL Get_Hess_FROM_Tab_OF_dnMatOp(hess_R,dnMatOp) ! for the ground state
    CALL Write_Mat(hess_R, out_unitp, 5, info='hess_R')
    Ene_R = Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)
    write(out_unitp,*) 'Ene_R',Ene_R
    CALL get_d0GG(Qact,para_Tnum,mole,d0G_R,def=.TRUE.)
    CALL Write_Mat(d0G_R, out_unitp, 5, info='d0G_R')
    CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

    write(out_unitp,*)
    write(out_unitp,*)
    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '======= Product geometry optimisation ============================'
    write(out_unitp,*) '=================================================================='

    Qact(iqS) = para_AllOp%tab_Op(iOp_CAPP)%para_ReadOp%tab_CAP(2)%Q0
    mole_100%ActiveTransfo%Qact0(iqS)    = para_AllOp%tab_Op(iOp_CAPP)%para_ReadOp%tab_CAP(2)%Q0
    mole_100%ActiveTransfo%Qdyn0(iQdyns) = para_AllOp%tab_Op(iOp_CAPP)%para_ReadOp%tab_CAP(2)%Q0
    QP = Qact(1:mole_100%nb_act)

    CALL Sub_Optimization(para_AllBasis%BasisnD,para_Tnum,mole_100,para_H%para_ReadOp%PrimOp_t, &
                          QP,para_Optimization)
    write(out_unitp,*) 'QP',QP
    Qact(1:mole_100%nb_act) = QP(:)
    CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn_P,mole_100%ActiveTransfo) ! here with mole_100
    CALL Qdyn_TO_Qact_FROM_ActiveTransfo(Qdyn_P,Qact,mole%ActiveTransfo)  !    then with mole
    write(out_unitp,*) 'Qact_P',Qact
    write(out_unitp,*) 'Qdyn_P',Qdyn_P

    CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,para_H%para_ReadOp%PrimOp_t%nb_elec,nderiv=2)
    CALL alloc_NParray(hess_P,[mole%nb_act,mole%nb_act],'hess_P',name_sub)
    CALL alloc_NParray(d0G_P, [mole%nb_act,mole%nb_act],'d0G_P',name_sub)

    CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_H%para_ReadOp%PrimOp_t)
    CALL Get_Hess_FROM_Tab_OF_dnMatOp(hess_P,dnMatOp) ! for the ground state
    CALL Write_Mat(hess_P, out_unitp, 5, info='hess_P')
    Ene_P = Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)
    write(out_unitp,*) 'Ene_P',Ene_P
    CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
    CALL get_d0GG(Qact,para_Tnum,mole,d0G_P,def=.TRUE.)
    CALL Write_Mat(d0G_P, out_unitp, 5, info='d0G_P')


    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '=================================================================='
    write(out_unitp,*) '======= CAP and basis parameters ================================='
    write(out_unitp,*) '=================================================================='
    Ene_Col    = 0.01_Rkind / get_Conv_au_TO_unit('E','eV')

    IF (Ene_P < Ene_R) THEN
      Ene_Asymp  = Ene_R
      write(out_unitp,*) " Higher asymptote (R): ",Qdyn_R(iQdyns),Ene_R
    ELSE
      Ene_Asymp  = Ene_P
      write(out_unitp,*) " Higher asymptote (P): ",Qdyn_P(iQdyns),Ene_P
    END IF
    write(out_unitp,*) " Higher asymptote: ",Ene_Asymp

    ! for the product
    mass = ONE/d0G_P(iQs,iQs)
    E = Ene_Asymp-Ene_P + Ene_Col
    LP = TWO*PI/sqrt(TWO*mass*E)
    ! for the reactif
    mass = ONE/d0G_R(iQs,iQs)
    E = Ene_Asymp-Ene_R + Ene_Col
    LR = TWO*PI/sqrt(TWO*mass*E)

    ! CAP parameters
    DO iact=1,mole%nb_act
      IF (iact == iQs) THEN
        IF (Qdyn_R(iQdyns) < 0) THEN
          A = Qdyn_R(iQdyns) - c*LR
          B = Qdyn_P(iQdyns) + c*LP
        ELSE
          A = Qdyn_P(iQdyns) - c*LP
          B = Qdyn_R(iQdyns) + c*LR
        END IF
        write(out_unitp,'(a,i0,a,f7.2,a,f7.2,a)') "&basis_nD iQdyn=",iQdyns," name='boxAB'  A=",A," B=",B," nb=$nb nq=$nq /"
      ELSE
        write(out_unitp,'(a,i0,a,3(f9.2,a))') "&basis_nD iQact=",iact," name='HO' Q0=",Qdyn_TS(iact), &
                      " k_HO=",hess_TS(iact,iact)," m_HO=",ONE/d0G_TS(iact,iact)," nb=$nb2 nq=$nq2 /"
      END IF
    END DO

    ! CAP parameters
    ! for the product
    E = Ene_Asymp-Ene_P + Ene_Col
    IF (Qdyn_P(iQdyns) < 0) THEN
      !write(out_unitp,*) 'mass (-)',ONE/d0G_P(iQs,iQs)
      write(out_unitp,'(a,e9.3,a,f7.2,a,f7.2,a,i0,a)') &
        "&CAP Name_Cap='-inv' A=",E,' Q0=',Qdyn_P(iQdyns),' LQ=',LP,' ind_Q=',iQs,' /'
    ELSE
      !write(out_unitp,*) 'mass (+)',ONE/d0G_P(iQs,iQs)
      write(out_unitp,'(a,e9.3,a,f7.2,a,f7.2,a,i0,a)') &
        "&CAP Name_Cap='+inv' A=",E,' Q0=',Qdyn_P(iQdyns),' LQ=',LP,' ind_Q=',iQs,' /'
    END IF
    ! for the reactif
    E = Ene_Asymp-Ene_R + Ene_Col

    IF (Qdyn_R(iQdyns) < 0) THEN
      !write(out_unitp,*) 'mass (-)',ONE/d0G_R(iQs,iQs)
      write(out_unitp,'(a,e9.3,a,f7.2,a,f7.2,a,i0,a)') &
        "&CAP Name_Cap='-inv' A=",E,' Q0=',Qdyn_R(iQdyns),' LQ=',LR,' ind_Q=',iQs,' /'
    ELSE
      !write(out_unitp,*) 'mass (+)',ONE/d0G_R(iQs,iQs)
      write(out_unitp,'(a,e9.3,a,f7.2,a,f7.2,a,i0,a)') &
        "&CAP Name_Cap='+inv' A=",E,' Q0=',Qdyn_R(iQdyns),' LQ=',LR,' ind_Q=',iQs,' /'
    END IF
  
  !=====================================================================
  !=====================================================================
  !=====================================================================
  !       deallocated memories
  !=====================================================================
        CALL dealloc_table_at(const_phys%mendeleev)
  
        CALL dealloc_CoordType(mole)
        CALL dealloc_Tnum(para_Tnum)
  
        CALL dealloc_para_AllOp(para_AllOp)
        CALL dealloc_para_ana(para_ana)
        CALL dealloc_param_propa(para_propa)
  
        CALL dealloc_AllBasis(para_AllBasis)
  
        write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
        write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
        write(out_unitp,*) '================================================'
        IF(openmpi) THEN
          write(out_unitp,*) ' sub_Opt_CAP_basis AU REVOIR!!!', ' from ', MPI_id
        ELSE
          write(out_unitp,*) ' sub_Opt_CAP_basis AU REVOIR!!!'
        ENDIF
        write(out_unitp,*) '================================================'
  
  END SUBROUTINE sub_Opt_CAP_basis
      SUBROUTINE sub_GridTOBasis_test(max_mem)
      USE mod_system
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_psi
      USE mod_Op
      USE mod_analysis
      USE mod_propa
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================
!

!----- variables for the dynamical memory allocation -----------------
      logical  :: intensity_only
      integer  :: nio_res_int

      integer (kind=ILkind)  :: max_mem
      logical   test_mem



!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

!----- variables for the construction of H ----------------------------
      TYPE (param_AllOp), target  :: para_AllOp


!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)           :: para_ana
      TYPE (param_intensity)     :: para_intensity

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: WP0,WP1,WP3

!----- variables for the namelist actives ----------------------------
!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis), target :: para_AllBasis

!----- variables divers ----------------------------------------------
       integer           :: i,ibb,nb_it = 1
!      real (kind=Rkind) :: T,DE,Ep,Em
!      logical           :: print_mat
!      integer           :: err
!      integer :: err_mem,memory


!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================
      write(out_unitp,*)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: BEGINNING ini_data'
      CALL time_perso('ini_data')
      write(out_unitp,*)
      write(out_unitp,*)
      CALL     ini_data(const_phys,para_Tnum,mole,                      &
                        para_AllBasis,para_AllOp,                       &
                        para_ana,para_intensity,intensity_only,         &
                        para_propa)

      write(out_unitp,*)
      write(out_unitp,*)
      CALL time_perso('ini_data')
      write(out_unitp,*) ' VIB: END ini_data'
      write(out_unitp,*) '================================================='
      write(out_unitp,*)
!=====================================================================




      CALL init_psi(WP0,para_AllOp%tab_Op(1),.FALSE.)
      CALL alloc_psi(WP0,BasisRep=.TRUE.,GridRep=.TRUE.)
      nb_mult_BTOG = 0
      nb_mult_GTOB = 0

para_mem%mem_debug = .FALSE.

      write(out_unitp,*)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: BEGINNING GridTOBasis_test'
      write(out_unitp,*)

!DO ibb=1,WP0%nb_tot
  ibb=2
  WP0 = ONETENTH
  WP0%RvecB(ibb) = ONE
  WP1 = WP0 !for the comparison
  !CALL ecri_psi(psi=WP0)


      CALL time_perso('B=>G')

      DO i=1,nb_it
        IF (mod(i,100) == 0) write(out_unitp,*) i ; flush(out_unitp)
        CALL sub_PsiBasisRep_TO_GridRep(WP0)
      END DO

      CALL time_perso('B=>G')

      !CALL ecri_psi(psi=WP0)
      WP0%RvecB(:) = ZERO

      write(out_unitp,*)
      write(out_unitp,*)
      CALL time_perso('G=>B')
      DO i=1,nb_it
        IF (mod(i,100) == 0) write(out_unitp,*) i ; flush(out_unitp)
        CALL sub_PsiGridRep_TO_BasisRep(WP0)
      END DO
      CALL time_perso('G=>B')
      !CALL ecri_psi(psi=WP0)

      write(out_unitp,*)
      write(out_unitp,*) 'nb_mult_BTOG,nb_mult_GTOB',nb_mult_BTOG,nb_mult_GTOB
      write(out_unitp,*)
      write(out_unitp,*) ' VIB: END GridTOBasis_test'
      write(out_unitp,*) '================================================='
      write(out_unitp,*)

      WP3 = WP1-WP0
      write(out_unitp,*) 'max diff',maxval(abs(WP3%RvecB))
!END DO
      write(out_unitp,*) '================================================'
      write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
      write(out_unitp,*) '================================================'


para_mem%mem_debug = .FALSE.


      END SUBROUTINE sub_GridTOBasis_test

      SUBROUTINE Sub_OpPsi_test(max_mem)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_psi

      USE mod_propa
      USE mod_FullPropa
      USE mod_FullControl
      USE mod_Davidson
      USE mod_Filter
      USE mod_Arpack
      USE mod_Op
      USE mod_analysis
      USE mod_fullanalysis
      USE mod_Auto_Basis
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the dynamical memory allocation -----------------
      logical  :: intensity_only
      integer  :: nio_res_int

      integer (kind=ILkind)  :: max_mem
      logical   test_mem



!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum


!----- variables for the construction of H ----------------------------
      TYPE (param_AllOp), target  :: para_AllOp
      TYPE (param_Op), pointer    :: para_H      => null()
      TYPE (param_Op), pointer    :: para_S      => null()
      TYPE (param_Op), pointer    :: para_Dip(:) => null()
      integer                     :: iOp
      real (kind=Rkind)           :: max_Sii,max_Sij

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)           :: para_ana
      TYPE (param_intensity)     :: para_intensity

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: WP0,OpWP0

!----- for Davidson diagonalization ----------------------------
      integer                    :: nb_diago,max_diago
      TYPE (param_psi),pointer   :: Tab_Psi(:) => null()
      TYPE (param_psi),pointer   :: Tab_OpPsi(:) => null()
      real (kind=Rkind),allocatable  :: Ene0(:)

!----- variables for the namelist actives ----------------------------
!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis) :: para_AllBasis

      integer :: Get_nbPERsym_FROM_SymAbelianOFAllBasis ! function
!----- variables divers ----------------------------------------------
      integer           :: i,ib,ip,i_baie,f_baie,id,nb_ScalOp
      real (kind=Rkind) :: T,DE,Ep,Em,Q,fac,zpe,pop,a,b
      logical           :: print_mat
      integer           :: err
      integer :: err_mem,memory
      real (kind=Rkind) :: part_func ! function

      integer           :: nb_it = 10
      integer           :: PSG4_maxth_save
      logical           :: cplx  = .FALSE.

!para_mem%mem_debug=.TRUE.


!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================
      IF(MPI_id==0) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: initialization of variables'
      ENDIF

      nullify(para_H)
      nullify(para_S)
      nullify(para_Dip)

      IF(MPI_id==0) THEN
        write(out_unitp,*) ' END VIB: initialization of variables'
        write(out_unitp,*) '================================================='

        write(out_unitp,*)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: BEGINNING ini_data'
        CALL time_perso('ini_data ini')
        write(out_unitp,*)
        write(out_unitp,*)
      ENDIF
      CALL     ini_data(const_phys,para_Tnum,mole,                      &
                        para_AllBasis,para_AllOp,                       &
                        para_ana,para_intensity,intensity_only,         &
                        para_propa)

      para_H => para_AllOp%tab_Op(1)

      IF(MPI_id==0) THEN
        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('ini_data end')
        write(out_unitp,*) ' VIB: END ini_data'
        write(out_unitp,*) '================================================='
        write(out_unitp,*)
!=====================================================================

!=====================================================================
!      Grids (V, T, Dip) calculations
!=====================================================================


        write(out_unitp,*)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: BEGINNING sub_qa_bhe'
        CALL time_perso('sub_qa_bhe')
        write(out_unitp,*)
        write(out_unitp,*)
      ENDIF  ! for MPI_id=0

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Read_FileGrid) THEN ! test only for H
        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*) ' The grid of operators (S V Veff T1 and T2) will be read'
          write(out_unitp,*)
        ENDIF
        IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid) THEN
          IF(MPI_id==0) write(out_unitp,*) 'Test_Grid=.TRUE. => STOP in vib'
          STOP
        END IF
      ELSE
        IF (para_ana%VibRot) THEN
           para_Tnum%JJ = para_ana%JJmax
           para_Tnum%With_Cart_Transfo = (para_Tnum%JJ>0) .AND. mole%Cart_transfo
        END IF

        CALL sub_qa_bhe(para_AllOp)

        IF (para_ana%VibRot) THEN
           para_Tnum%JJ = 0
           para_Tnum%With_Cart_Transfo = (para_Tnum%JJ>0) .AND. mole%Cart_transfo
        END IF
      END IF

      IF(MPI_id==0) THEN
        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('sub_qa_bhe')
        write(out_unitp,*) ' VIB: END sub_qa_bhe'
        write(out_unitp,*) '================================================='
        write(out_unitp,*)
      ENDIF

!=====================================================================
!===== contraction of the active basis set with HADA basis ===========
!=====================================================================

      IF (para_H%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC .AND. para_H%nb_bi>1) THEN

        write(out_unitp,*)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: BEGINNING HADA contraction',          &
                                   para_AllBasis%BasisnD%nb
        CALL time_perso('HADA contraction')
        write(out_unitp,*)
        write(out_unitp,*)

        CALL sub_MatOp_HADA(para_H,para_ana,para_intensity,para_AllOp,  &
                            const_phys)

        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('HADA contraction')
        write(out_unitp,*) ' VIB: END HADA contraction'
        write(out_unitp,*) '================================================='
        write(out_unitp,*)

        para_AllOp%tab_Op(:)%nb_tot     =                               &
                sum(para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(:)) * &
                para_H%para_ReadOp%nb_elec
        para_AllOp%tab_Op(:)%nb_tot_ini = para_AllOp%tab_Op(:)%nb_tot
      END IF


      !================================================================
      !===== build S and/or H if necessary ============================
      !================================================================
      IF (para_H%Make_Mat) THEN

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING sub_matOp',para_H%nb_tot
          CALL time_perso('sub_matOp: H and S')
          write(out_unitp,*)
          write(out_unitp,*) 'para_S...%comput_S',para_AllOp%tab_Op(2)%para_ReadOp%comput_S
          write(out_unitp,*)
        ENDIF

        IF (para_AllOp%tab_Op(2)%para_ReadOp%comput_S) THEN

          para_S => para_AllOp%tab_Op(2)
          CALL sub_MatOp(para_S,para_ana%print)

          !- analysis of the overlap matrix
          CALL sub_ana_S(para_S%Rmat,para_S%nb_tot,max_Sii,max_Sij,.TRUE.)

          CALL dealloc_para_Op(para_S)
          nullify(para_S)
          CALL time_perso('sub_matOp: S')
        END IF

        CALL sub_MatOp(para_H,para_ana%print)

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_matOp: H and S')
          write(out_unitp,*) ' VIB: END sub_matOp'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF
      END IF

      cplx = (para_ana%propa .OR. para_H%cplx)

      !@chen  different in my version
      IF (cplx) nb_it = 1
      allocate(Tab_Psi(nb_it))
      allocate(Tab_OpPsi(nb_it))
      DO i=1,nb_it
        CALL init_psi(Tab_Psi(i),para_H,cplx)
        CALL alloc_psi(Tab_Psi(i),BasisRep=.TRUE.,GridRep=.FALSE.)
        Tab_Psi(i)   = ZERO
        Tab_OpPsi(i) = Tab_Psi(i)

        CALL Set_psi_With_index(Tab_Psi(i),ONE,ind_aie=i)

      END DO

para_mem%mem_debug = .FALSE.

write(out_unitp,*)
write(out_unitp,*) '================================================='
write(out_unitp,*) ' VIB: BEGINNING sub_OpPsi_test'

      write(out_unitp,*)
      write(out_unitp,*) ' Number of psi:',nb_it
      write(out_unitp,*) ' Number of threads (SG4):',SG4_maxth

      write(out_unitp,*)
      CALL time_perso('HPsi')

      IF (cplx) THEN
        CALL sub_OpPsi(Tab_Psi(1),Tab_OpPsi(1),para_H)
      ELSE
        CALL sub_TabOpPsi(Tab_Psi,Tab_OpPsi,para_H)
      END IF

      CALL time_perso('HPsi')
      write(out_unitp,*)
      write(out_unitp,*) 'nb_mult_BTOG,nb_mult_GTOB',nb_mult_BTOG,nb_mult_GTOB
      write(out_unitp,*)

      write(out_unitp,*) 'OpPsi states:'
      para_propa%file_WP%name = 'file_WP'
      CALL sub_save_psi(Tab_OpPsi,size(Tab_Psi),para_propa%file_WP)

write(out_unitp,*) ' VIB: END sub_OpPsi_test'
write(out_unitp,*) '================================================='
write(out_unitp,*)


CALL Tune_SG4threads_HPsi(cplx,nb_it,para_H)


      !CALL ecri_psi(psi=OpWP0)

      write(out_unitp,*) '================================================'
      write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
      write(out_unitp,*) '================================================'

para_mem%mem_debug = .FALSE.


      CALL dealloc_table_at(const_phys%mendeleev)

      CALL dealloc_CoordType(mole)
      IF (associated(para_Tnum%Gref)) THEN
        CALL dealloc_array(para_Tnum%Gref,"para_Tnum%Gref","vib")
      END IF
      CALL dealloc_Tnum(para_Tnum)

      CALL dealloc_para_AllOp(para_AllOp)

      CALL dealloc_para_ana(para_ana)


      CALL dealloc_param_propa(para_propa)

      CALL dealloc_psi(WP0)

      IF ( associated(Tab_Psi) ) THEN
        DO i=1,size(Tab_Psi)
          CALL dealloc_psi(Tab_Psi(i))
        END DO
        CALL dealloc_array(Tab_Psi,"Tab_Psi","vib")
      END IF
      IF ( associated(Tab_OpPsi) ) THEN
        DO i=1,size(Tab_OpPsi)
          CALL dealloc_psi(Tab_OpPsi(i))
        END DO
        CALL dealloc_array(Tab_OpPsi,"Tab_OpPsi","vib")
      END IF

      CALL dealloc_AllBasis(para_AllBasis)

      write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
      write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
      write(out_unitp,*) '================================================'
      write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
      write(out_unitp,*) '================================================'


END SUBROUTINE Sub_OpPsi_test


SUBROUTINE Tune_SG4threads_HPsi(cplx,nb_psi,para_H)
USE mod_system
USE mod_psi,     ONLY : param_psi,alloc_psi,dealloc_NParray,Set_Random_psi
USE mod_Op
IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================
!

!----- variables for the construction of H ----------------------------
 TYPE (param_Op)    :: para_H


!----- for Davidson diagonalization ----------------------------
 integer            :: nb_psi
 logical            :: cplx

 TYPE (param_psi), allocatable   :: Tab_Psi(:)
 TYPE (param_psi), allocatable   :: Tab_OpPsi(:)

 TYPE (Time_t) :: HPsiTime

 integer           :: nb_psi_loc,i,ib,PSG4_maxth_save,opt_PSG4_maxth
 logical           :: Make_Mat_save
 real(kind=Rkind)  :: a,b,Opt_RealTime,RealTime(SG4_maxth)

!----- for debuging --------------------------------------------------
 character (len=*), parameter :: name_sub='Tune_SG4threads_HPsi'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

IF(openmpi) RETURN

IF (.NOT. Tune_SG4_omp) RETURN
IF (para_H%BasisnD%SparseGrid_type /= 4) RETURN

para_mem%mem_debug = .FALSE.

Make_Mat_save = para_H%Make_Mat
para_H%Make_Mat = .FALSE.

nb_psi_loc = nb_psi
IF (cplx) nb_psi_loc = 1

allocate(Tab_Psi(nb_psi_loc))
allocate(Tab_OpPsi(nb_psi_loc))
DO i=1,nb_psi_loc
  CALL init_psi(Tab_Psi(i),para_H,cplx)
  CALL alloc_psi(Tab_Psi(i),BasisRep=.TRUE.,GridRep=.FALSE.)
  Tab_Psi(i)   = ZERO
  Tab_OpPsi(i) = Tab_Psi(i)
  CALL Set_Random_psi(Tab_Psi(i))
END DO

write(out_unitp,*) '============================================'
write(out_unitp,*) '== Tuning the number of OMP threads ========'

write(out_unitp,*)
write(out_unitp,*) ' Number of psi:',size(Tab_Psi)
write(out_unitp,*) ' cplx?:',cplx

write(out_unitp,*) ' Number of threads (SG4):',SG4_maxth

RealTime(1) = Delta_RealTime(HPsiTime)

PSG4_maxth_save = SG4_maxth
Opt_RealTime    = huge(ONE)
opt_PSG4_maxth  = SG4_maxth
DO SG4_maxth=1,PSG4_maxth_save

  IF (Tab_Psi(1)%cplx) THEN
    CALL sub_OpPsi(Tab_Psi(1),Tab_OpPsi(1),para_H)
  ELSE
    CALL sub_TabOpPsi(Tab_Psi,Tab_OpPsi,para_H)
  END IF

  RealTime(SG4_maxth) = Delta_RealTime(HPsiTime)
  write(out_unitp,*) 'With ',SG4_maxth,'threads, Delta Real Time',RealTime(SG4_maxth)
  IF (RealTime(SG4_maxth) < Opt_RealTime) THEN
    IF (RealTime(SG4_maxth) > 0) THEN
      Opt_RealTime   = RealTime(SG4_maxth)
      opt_PSG4_maxth = SG4_maxth
    END IF
  ELSE
    IF (SG4_maxth > 1) EXIT
  END IF


END DO

SG4_maxth = opt_PSG4_maxth
write(out_unitp,*) 'Optimal threads: ',SG4_maxth,' Delta Real Time',RealTime(SG4_maxth)
write(out_unitp,*) '    => Speed-up: ',RealTime(1)/RealTime(SG4_maxth)

write(out_unitp,*) '============================================'

para_H%Make_Mat = Make_Mat_save

CALL dealloc_NParray(Tab_OpPsi,'Tab_OpPsi',name_sub)
CALL dealloc_NParray(Tab_Psi,  'Tab_Psi',  name_sub)

END SUBROUTINE Tune_SG4threads_HPsi

      SUBROUTINE sub_Analysis_Only(max_mem)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_psi

      USE mod_Op
      USE mod_analysis
      USE mod_fullanalysis
      USE mod_propa
      USE mod_Auto_Basis

      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================
!

!----- variables for the dynamical memory allocation -----------------
      logical  :: intensity_only
      integer  :: nio_res_int

      integer (kind=ILkind)  :: max_mem
      logical   test_mem



!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

!----- variables for the construction of H ----------------------------
      TYPE (param_AllOp), target  :: para_AllOp
      TYPE (param_Op),    pointer :: para_H      => null()



!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)           :: para_ana
      TYPE (param_intensity)     :: para_intensity

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa)            :: para_propa
      TYPE (param_WP0)              :: para_WP0
      TYPE (param_psi), allocatable :: Tab_WP(:)
      TYPE (param_psi)              :: WP0
      integer                       :: nioWP

!----- variables for the namelist actives ----------------------------
!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis) :: para_AllBasis
      integer :: Get_nbPERsym_FROM_SymAbelianOFAllBasis ! function

!----- variables divers ----------------------------------------------
       integer           :: i,max_diago
!      real (kind=Rkind) :: T,DE,Ep,Em
!      logical           :: print_mat
!      integer           :: err


!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================
      write(out_unitp,*)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: BEGINNING ini_data'
      CALL time_perso('ini_data')
      write(out_unitp,*)
      write(out_unitp,*)
      CALL     ini_data(const_phys,para_Tnum,mole,                      &
                        para_AllBasis,para_AllOp,                       &
                        para_ana,para_intensity,intensity_only,         &
                        para_propa)

      para_H => para_AllOp%tab_Op(1)
      CALL init_psi(WP0,para_H,para_propa%para_WP0%WP0cplx)

      para_ana%intensity    = .FALSE.
      para_H%para_ReadOp%nb_scalar_Op   = 0
      write(out_unitp,*)
      write(out_unitp,*)
      CALL time_perso('ini_data')
      write(out_unitp,*) ' VIB: END ini_data'
      write(out_unitp,*) '================================================='
      write(out_unitp,*)
!=====================================================================


      write(out_unitp,*)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' VIB: BEGINNING WP analysis only'
      CALL time_perso('WP analysis only')
      write(out_unitp,*)
      write(out_unitp,*)


      ! first the number of WP
      IF (para_ana%davidson) THEN


        IF (para_propa%para_Davidson%max_WP == 0) THEN
          max_diago = max(1000,para_propa%para_Davidson%nb_WP,          &
                          para_H%nb_tot/10)
        ELSE
          max_diago = para_propa%para_Davidson%max_WP
        END IF
        IF (Get_nbPERsym_FROM_SymAbelianOFAllBasis(para_AllBasis,       &
                             para_propa%para_Davidson%symab) == 0) THEN
          max_diago = min(max_diago,para_H%nb_tot)
        ELSE
          max_diago = min(max_diago,para_H%nb_tot,                      &
                  Get_nbPERsym_FROM_SymAbelianOFAllBasis(para_AllBasis, &
                                        para_propa%para_Davidson%symab))
        END IF

        para_WP0%nb_WP0              = 0
        para_WP0%read_file           = .TRUE.
        para_propa%file_WP%formatted = para_propa%para_Davidson%formatted_file_readWP
        para_WP0%file_WP0            = para_propa%file_WP
        para_WP0%read_listWP0        = para_propa%para_Davidson%read_listWP
        para_WP0%WP0cplx             = para_H%cplx
        para_propa%file_WP%name      = para_propa%para_Davidson%name_file_saveWP
        para_WP0%file_WP0%name       = para_propa%para_Davidson%name_file_readWP
      ELSE
        write(out_unitp,*) ' ERROR in sub_Analysis_Only'
        write(out_unitp,*) ' Not yet without Davidson'
        STOP
      END IF


      ! allocation
      CALL alloc_NParray(Tab_WP,[max_diago],"Tab_WP","sub_Analysis_Only")
      DO i=1,size(Tab_WP)
        CALL init_psi(Tab_WP(i),para_H,para_WP0%WP0cplx)
      END DO


      ! read the WPs
      CALL sub_read_psi0(tab_WP,para_WP0,max_diago)

      ! Analysis
      IF(keep_MPI) CALL sub_analyse(Tab_WP,para_WP0%nb_WP0,para_H,                   &
                       para_ana,para_intensity,para_AllOp,const_phys)

      write(out_unitp,*)
      write(out_unitp,*)
      CALL time_perso('WP analysis only')
      write(out_unitp,*) ' VIB: BEGINNING WP analysis only'
      write(out_unitp,*) '================================================='
      write(out_unitp,*)

      IF ( allocated(tab_WP) ) THEN
        DO i=1,size(tab_WP)
          CALL dealloc_psi(tab_WP(i))
        END DO
        CALL dealloc_NParray(tab_WP,"tab_WP","sub_Analysis_Only")
      END IF

      write(out_unitp,*) '================================================'
      write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
      write(out_unitp,*) '================================================'

      END SUBROUTINE sub_Analysis_Only
