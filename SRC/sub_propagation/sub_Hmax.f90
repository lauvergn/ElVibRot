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
      SUBROUTINE sub_Hmax(para_propa,para_H)
      USE EVR_system_m
      USE mod_Op
      USE mod_psi,      ONLY : param_psi,Set_psi_With_index,renorm_psi, &
                       alloc_psi,dealloc_psi,alloc_array,dealloc_array, &
                       Write_ana_psi
      USE mod_ana_psi_MPI
      USE mod_propa
      USE mod_FullPropa
      USE mod_Davidson
      USE mod_Hmax_MPI
      USE mod_MPI_aux
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)    :: para_H

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)   :: psi,Hpsi
      TYPE (param_propa) :: para_propa

!------ active Matrix H ------------------------------------------
      real (kind=Rkind)    :: Hinter
      complex (kind=Rkind) :: Emax

!----- quadrature points and weight -----------------------------
      real (kind=Rkind)    :: WnD

!----- working parameters --------------------------------------------
      integer           :: i_qa
      integer           :: i_h,i1_h,i,k_term
      integer           :: nio
      real (kind=Rkind) :: a
      !logical           :: relax = .TRUE.
      logical           :: relax = .FALSE.

      real (kind=Rkind) :: Qdyn(para_H%mole%nb_var)
      real (kind=Rkind) :: Qact(para_H%mole%nb_act1)
      TYPE (param_d0MatOp) :: d0MatOp
      integer              :: type_Op


      real (kind=Rkind), pointer :: Grid(:,:,:)

      integer                    :: nb_diago,max_diago
      TYPE (param_psi),pointer   :: Tab_Psi(:) => null()
      real (kind=Rkind),allocatable  :: Ene0(:)
      integer                    :: ii,temp_Hmin

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'sub_Hmax'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub,' ',para_H%nb_tot
        write(out_unit,*) 'Hmin,Hmax',para_H%Hmin,para_H%Hmax
        flush(out_unit)
        CALL Write_ana_psi(para_propa%ana_psi)
      END IF
!-----------------------------------------------------------

!     --------------------------------------------------------
!     -  Hmin ------------------------------------------------
!     --------------------------------------------------------
      !IF (para_H%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN
      IF(para_H%OpGrid(1)%para_FileGrid%Save_MemGrid_done) THEN
!       IF(openmpi) THEN
!         CALL get_Hmin_MPI(para_H)
        IF(para_H%para_ReadOp%para_FileGrid%Type_FileGrid==4) THEN
          IF(openmpi) THEN
            CALL get_Hmin_MPI(para_H)
          ELSE
            para_H%Hmin = HUGE(ONE)
            DO ii=1,para_H%BasisnD%para_SGType2%nb_SG
              temp_Hmin=minval(para_H%OpGrid(1)%SRep%SmolyakRep(ii)%V)
              IF(temp_Hmin<para_H%Hmin) para_H%Hmin=temp_Hmin
            ENDDO
          ENDIF
        ELSE
          ! Minimal value of Veff
          para_H%Hmin = para_H%OpGrid(1)%Op_min

        ENDIF ! para_H%BasisnD%SparseGrid_type==4

      ELSE ! the grids are on a file

!       -  Hmin ------------------------------------------------
        para_H%Hmin = huge(ONE)

        IF (para_H%para_ReadOp%para_FileGrid%Type_FileGrid == 0) THEN
          type_Op = para_H%para_ReadOp%Type_HamilOp ! H
          IF (type_Op /= 1) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '    Type_HamilOp MUST be equal to 1 for the usual SHADA file'
            write(out_unit,*) '    CHECK your data!!'
            STOP
          END IF

          CALL Init_d0MatOp(d0MatOp,para_H%param_TypeOp,para_H%nb_bie)

          !- loop on i_qa -----------------------------------------
          DO i_qa=1,para_H%nb_qa

            CALL sub_reading_Op(i_qa,para_H%nb_qa,d0MatOp,para_H%n_Op,  &
                           Qdyn,para_H%mole%nb_var,para_H%mole%nb_act1, &
                                Qact,WnD,para_H%file_grid)

            k_term = d0MatOp%derive_term_TO_iterm(0,0)
            DO i1_h=1,para_H%nb_bie
              para_H%Hmin = min(para_H%Hmin,d0MatOp%ReVal(i1_h,i1_h,k_term))
            END DO

          END DO
          CALL dealloc_d0MatOp(d0MatOp)

        ELSE
          nullify(Grid)
          CALL alloc_array(Grid,                                        &
                          [para_H%nb_qa,para_H%nb_bie,para_H%nb_bie], &
                          'Grid',name_sub)

          DO k_term=1,para_H%nb_term
            IF (para_H%OpGrid(k_term)%derive_termQact(1) /= 0 .OR.      &
                para_H%OpGrid(k_term)%derive_termQact(2) /= 0) CYCLE

            IF (para_H%OpGrid(k_term)%grid_cte) THEN
               para_H%Hmin = minval(para_H%OpGrid(k_term)%Mat_cte(:,:))
            ELSE

              !$OMP critical(CRIT_sub_Hmax)
              IF (para_H%OpGrid(k_term)%file_Grid%seq) THEN   ! sequential acces file
                CALL sub_ReadSeq_Grid_iterm(Grid,para_H%OpGrid(k_term))
              ELSE  ! direct acces file
                CALL sub_ReadDir_Grid_iterm(Grid,para_H%OpGrid(k_term))
              END IF
              !$OMP end critical(CRIT_sub_Hmax)
            END IF

            EXIT ! we exit the loop because we find the potential (k_term)
          END DO
          para_H%Hmin = minval(Grid)

          CALL dealloc_array(Grid,'Grid',name_sub)

        END IF

!       --------------------------------------------------------
      END IF ! for para_H%para_ReadOp%para_FileGrid%Save_MemGrid_done
      IF (debug) write(out_unit,*) 'Hmin: ',para_H%Hmin
      flush(out_unit)

!     --------------------------------------------------------
!     -  Hmax ------------------------------------------------
!     --------------------------------------------------------
      para_H%Hmax = -huge(ONE)

      CALL init_psi(psi,para_H,para_H%cplx)
      IF (debug) write(out_unit,*) 'nb_tot',psi%nb_tot
      CALL alloc_psi(psi,BasisRep =.TRUE.)

      IF(keep_MPI) psi = ZERO
      IF(keep_MPI) CALL Set_psi_With_index(psi,R=ONE,ind_aie=psi%nb_tot)
      psi%symab = -1

      IF(keep_MPI) CALL renorm_psi(psi,BasisRep=.TRUE.)

      IF (debug) write(out_unit,*) 'psi%symab',psi%symab

      IF(keep_MPI) Hpsi = psi

      CALL sub_PsiOpPsi(Emax,psi,Hpsi,para_H)
      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(Emax,size1_MPI,root_MPI)
      para_H%Hmax = Real(Emax,kind=Rkind)

      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(para_H%Hmax,size1_MPI,root_MPI)

      IF (debug) write(out_unit,*) 'Hmax: ',para_H%Hmax
      flush(out_unit)

      CALL dealloc_psi(psi)
      CALL dealloc_psi(Hpsi)
!       --------------------------------------------------------

      write(out_unit,*) 'non-auto : Hmin,Hmax',para_H%Hmin,para_H%Hmax
      write(out_unit,*) 'nb_tot',para_H%nb_tot
      flush(out_unit)

      IF (para_H%read_Op) THEN
        IF (para_H%cplx) THEN
          para_H%Hmax = para_H%Cmat(para_H%nb_tot,para_H%nb_tot)
          para_H%Hmin = para_H%Cmat(1,1)
        ELSE
          para_H%Hmax = para_H%Rmat(para_H%nb_tot,para_H%nb_tot)
          para_H%Hmin = para_H%Rmat(1,1)
        END IF
        para_propa%Hmax = para_H%Hmax
        para_propa%Hmin = para_H%Hmin
        write(out_unit,*) 'read_Op : Hmin,Hmax',para_H%Hmin,para_H%Hmax
        RETURN
      END IF

      IF (para_H%spectral) THEN
        IF (para_H%cplx) THEN
          IF (.NOT. associated(para_H%Cdiag)) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  spectral=t and para_H%Cdiag is not allocated'
            STOP 'ERROR in sub_Hmax: pectral=t and para_H%Cdiag is not allocated'
          END IF

          para_H%Hmax = para_H%Cdiag(para_H%para_AllBasis%basis_ext%nb_vp_spec)
          para_H%Hmin = para_H%Cdiag(1)
        ELSE
           IF (.NOT. associated(para_H%Rdiag)) THEN
                 write(out_unit,*) ' ERROR in ',name_sub
                 write(out_unit,*) '  spectral=t and para_H%Rdiag is not allocated'
                 STOP 'ERROR in sub_Hmax: pectral=t and para_H%Rdiag is not allocated'
               END IF
           para_H%Hmax = para_H%Rdiag(para_H%para_AllBasis%basis_ext%nb_vp_spec)
          para_H%Hmin = para_H%Rdiag(1)
        END IF
        para_propa%Hmax = para_H%Hmax
        para_propa%Hmin = para_H%Hmin
        write(out_unit,*) 'spectral : Hmin,Hmax',para_H%Hmin,para_H%Hmax
        RETURN
      END IF

      !CALL sub_Auto_Hmax_cheby(para_propa,para_H)

      IF (para_propa%auto_Hmax) THEN

        IF (para_H%cplx) relax = .TRUE. ! Because Davidson doesn't work (yet) with complex H
relax = .TRUE.
        IF (relax) THEN
          CALL sub_Auto_HmaxHmin_relax(para_propa,para_H)
        ELSE
          write(out_unit,*) 'Davidson Hmin Hmax' ; flush(out_unit)
          para_propa%para_Davidson%nb_WP            = 0
          para_propa%para_Davidson%lower_states     = .TRUE.
          para_propa%para_Davidson%project_WP0      = .FALSE.
          para_propa%para_Davidson%all_lower_states = .FALSE.
          para_propa%para_Davidson%max_it           = 200
          para_propa%para_Davidson%num_resetH       = para_propa%para_Davidson%max_it
          para_propa%para_Davidson%num_checkS       = para_propa%para_Davidson%max_it
          para_propa%para_Davidson%NewVec_type      = 4
          para_propa%para_Davidson%read_WP          = .FALSE.
          para_propa%para_Davidson%conv_ene         = ONETENTH**6
          para_propa%para_Davidson%conv_resi        = ONETENTH**4
          para_propa%para_Davidson%symab            = -1


          nb_diago = 0
          max_diago = para_propa%para_Davidson%max_it + 1
          nullify(Tab_Psi)
          CALL alloc_array(Tab_Psi,[max_diago],"Tab_Psi",name_sub)
          CALL alloc_NParray(Ene0,[max_diago],"Ene0",name_sub)


          para_propa%para_Davidson%Hmin_propa       = .TRUE.
          para_propa%para_Davidson%Hmax_propa       = .FALSE.
          para_propa%para_Davidson%name_file_saveWP = 'file_WP_Hmin'

          CALL sub_propagation_Davidson(Tab_Psi,Ene0,nb_diago,max_diago,&
                              para_H,para_propa%para_Davidson,para_propa)

          nb_diago = 0
          para_propa%para_Davidson%Hmin_propa       = .FALSE.
          para_propa%para_Davidson%Hmax_propa       = .TRUE.
          para_propa%para_Davidson%name_file_saveWP = 'file_WP_Hmax'

          CALL sub_propagation_Davidson(Tab_Psi,Ene0,nb_diago,max_diago,&
                              para_H,para_propa%para_Davidson,para_propa)

          CALL dealloc_array(Tab_Psi,"Tab_Psi",name_sub)
          CALL dealloc_NParray(Ene0,"Ene0",name_sub)
          nullify(Tab_Psi)

          write(out_unit,*) 'END Davidson Hmin Hmax' ; flush(out_unit)


        END IF
      END IF

      para_propa%Hmax = para_H%Hmax
      para_propa%Hmin = para_H%Hmin
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
      END IF
      write(out_unit,*) 'END ',name_sub
!-----------------------------------------------------------

      END SUBROUTINE sub_Hmax
!=======================================================================================

!=======================================================================================
      SUBROUTINE sub_Auto_HmaxHmin_relax(para_propa,para_H)
      USE EVR_system_m
      USE mod_psi,      ONLY : param_psi,renorm_psi,alloc_psi,dealloc_psi,&
                               Set_psi_With_index,Write_ana_psi
      USE mod_Op
      USE mod_propa
      USE mod_FullPropa
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_H

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)   :: WP0(1),WP(1)
      TYPE (param_propa) :: para_propa

      TYPE (param_propa) :: para_propa_loc

!----- working parameters --------------------------------------------
      complex (kind=Rkind) :: Emax
      integer           :: i
      real (kind=Rkind) :: a

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'sub_Auto_HmaxHmin_relax'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      write(out_unit,*) 'BEGINNING ',name_sub,' ',para_H%nb_tot
      write(out_unit,*) 'Hmin,Hmax',para_H%Hmin,para_H%Hmax
      IF (debug) THEN
        CALL Write_ana_psi(para_propa%ana_psi)
      END IF
!-----------------------------------------------------------

        para_propa_loc%Hmin           = para_H%Hmin
        para_propa_loc%Hmax           = para_H%Hmax
        para_propa_loc%para_poly%Hmin = para_H%Hmin
        para_propa_loc%para_poly%Hmax = para_H%Hmax

        para_propa_loc%WPTmax         = TEN**6
        para_propa_loc%WPdeltaT       = ONE
        para_propa_loc%nb_micro       = 1

        para_propa_loc%para_poly%max_poly     = 20
        para_propa_loc%para_poly%npoly        = 20

        CALL alloc_array(para_propa_loc%para_poly%coef_poly,[2],      &
                        "para_propa_loc%para_poly%coef_poly",name_sub)

        para_propa_loc%para_poly%poly_tol     = ONETENTH**6
        para_propa_loc%para_poly%DHmax        = ZERO

        para_propa_loc%write_iter        = (print_level > 1)
        para_propa_loc%WriteWP_nDeltaT   = 1
        para_propa_loc%WPpsi2            = .FALSE.
        para_propa_loc%WPpsi             = .FALSE.
        para_propa_loc%file_autocorr%name= 'WP_auto'
        para_propa_loc%ana_psi           = para_propa%ana_psi

        CALL init_psi(WP(1),para_H,para_H%cplx)
        WP(1)%GridRep = .TRUE.
        CALL alloc_psi(WP(1))

        CALL init_psi(WP0(1),para_H,para_H%cplx)
        WP0(1)%BasisRep = .TRUE.
        CALL alloc_psi(WP0(1))

        !---- for Hmax -----------------------------------------
        para_propa_loc%name_WPpropa      = 'Emax'
        para_propa_loc%file_WP%name      = make_EVRTFileName('file_WP_Hmax')
        para_propa_loc%ana_psi%file_Psi  = para_propa_loc%file_WP
        para_propa_loc%ana_psi%propa     = .TRUE.
        para_propa_loc%ana_psi%Write_psi = .FALSE.

        WP0(1) = ZERO
        CALL Set_psi_With_index(WP0(1),R=ONE,ind_aie=WP0(1)%nb_tot)
        WP0(1)%symab = -1
        CALL renorm_psi(WP0(1),BasisRep=.TRUE.)

        para_propa_loc%type_WPpropa = -3
        CALL sub_propagation3(Emax,WP0,WP,[para_H],para_propa_loc)
        para_H%Hmax = Real(Emax,kind=Rkind)

        !---- END for Hmax -------------------------------------


        !---- for Hmin -----------------------------------------
        para_propa_loc%name_WPpropa       = 'Emin'
        para_propa_loc%file_WP%name       = make_EVRTFileName('file_WP_Hmin')
        para_propa_loc%ana_psi%file_Psi   = para_propa_loc%file_WP
        para_propa_loc%para_poly%poly_tol = ONETENTH**8
        para_propa_loc%WPdeltaT           = ONE

        WP0(1) = ZERO
        CALL Set_psi_With_index(WP0(1),R=ONE,ind_aie=1)
        WP0(1)%symab = -1
        CALL renorm_psi(WP0(1),BasisRep=.TRUE.)

        para_propa_loc%type_WPpropa = 3
        CALL sub_propagation3(Emax,WP0,WP,[para_H],para_propa_loc)
        para_H%Hmin = Real(Emax,kind=Rkind)
        !---- END for Hmin -------------------------------------

        CALL dealloc_psi(WP0(1))
        CALL dealloc_psi(WP(1))

        CALL dealloc_param_propa(para_propa_loc)

        para_propa%para_poly%DHmax        = ZERO

        write(out_unit,*) 'auto : Hmin,Hmax',para_H%Hmin,para_H%Hmax

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
      END IF
      write(out_unit,*) 'END ',name_sub
!-----------------------------------------------------------

      END SUBROUTINE sub_Auto_HmaxHmin_relax

!=======================================================================================
      ! we are using the fact that chebychev propagation is very sensitive to the spectral range
      SUBROUTINE sub_Auto_Hmax_cheby(para_propa,para_H)
      USE EVR_system_m
      USE mod_psi,      ONLY : param_psi,renorm_psi,alloc_psi,dealloc_psi,&
                               Set_Random_psi,Write_ana_psi
      USE mod_Op
      USE mod_march
      USE mod_propa
      USE mod_FullPropa
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_H

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)   :: WP0(1),WP(1)
      TYPE (param_propa) :: para_propa

!----- working parameters --------------------------------------------
      integer           :: i,nb_HPsi,nb_HPsi_opt
      real (kind=Rkind) :: a,WPdeltaT,T,DeltaT_opt
      integer           :: type_WPpropa
      logical           :: With_field
      character (len=Name_len) :: name_WPpropa


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'sub_Auto_Hmax_cheby'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      write(out_unit,*) 'BEGINNING ',name_sub,' ',para_H%nb_tot
      write(out_unit,*) 'Hmin,Hmax',para_H%Hmin,para_H%Hmax
      IF (debug) THEN

      END IF
!-----------------------------------------------------------

     T = ZERO

     type_WPpropa            = para_propa%type_WPpropa
     para_propa%type_WPpropa = 1
     name_WPpropa            = para_propa%name_WPpropa
     para_propa%name_WPpropa = 'cheby'
     With_field              = para_propa%With_field
     para_propa%With_field   = .FALSE.

     para_propa%para_poly%Hmax = para_H%Hmax
     para_propa%para_poly%Hmin = para_H%Hmin

     WPdeltaT                = para_propa%WPdeltaT
     para_propa%WPdeltaT     = min(ONE,WPdeltaT)

     CALL init_psi(WP(1),para_H,cplx=.TRUE.)
     WP(1)%GridRep = .TRUE.
     CALL alloc_psi(WP(1))
     WP0(1) = WP(1)

     nb_HPsi_opt = huge(1)

     DO

        CALL Set_Random_psi(WP0(1))
        WP0(1)%symab = -1
        CALL renorm_psi(WP0(1),BasisRep=.TRUE.)
        WP(1) = WP0(1)

        CALL march_gene(T,WP,WP0,1,.FALSE.,para_H,para_propa)

        write(out_unit,*) 'March_cheby with WPdeltaT',para_propa%WPdeltaT
        write(out_unit,*) 'Opt npoly',para_propa%para_poly%npoly_Opt

        nb_HPsi = para_propa%para_poly%npoly_Opt*(ONE+para_propa%WPTmax/para_propa%WPdeltaT)
        IF (nb_HPsi < nb_HPsi_opt) THEN
          nb_HPsi_opt = nb_HPsi
          DeltaT_opt  = para_propa%WPdeltaT
        END IF
        write(out_unit,*) '# HPsi',nb_HPsi
        write(out_unit,*) 'opt nb_HPsi,DeltaT_opt',nb_HPsi_opt,DeltaT_opt


        IF (para_propa%WPdeltaT == WPdeltaT) EXIT

        para_propa%WPdeltaT = para_propa%WPdeltaT*TEN
        IF (para_propa%WPdeltaT > WPdeltaT) para_propa%WPdeltaT = WPdeltaT

     END DO
     para_propa%WPdeltaT     = WPdeltaT

      write(out_unit,*) 'Optimal WPdeltaT',DeltaT_opt
      write(out_unit,*) '# HPsi (for the whole propagation)',nb_HPsi_opt
      para_propa%WPdeltaT     = DeltaT_opt



      ! back to the old values
      !para_propa%WPdeltaT     = WPdeltaT
      para_propa%type_WPpropa = type_WPpropa
      para_propa%name_WPpropa = name_WPpropa
      para_propa%With_field   = With_field


      CALL dealloc_psi(WP0(1))
      CALL dealloc_psi(WP(1))


      para_propa%para_poly%DHmax        = ZERO


     para_H%Hmax = para_propa%para_poly%Hmax
     para_H%Hmin = para_propa%para_poly%Hmin
      write(out_unit,*) 'auto : Hmin,Hmax',para_H%Hmin,para_H%Hmax

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
      END IF
      write(out_unit,*) 'END ',name_sub
!-----------------------------------------------------------

      END SUBROUTINE sub_Auto_Hmax_cheby
