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
MODULE mod_psi
  USE mod_param_WP0
  USE mod_type_ana_psi
  USE mod_psi_set_alloc
  USE mod_ana_psi
  USE mod_psi_io
  USE mod_psi_B_TO_G
  USE mod_psi_Op
  USE mod_ana_psi_MPI, ONLY:share_psi_nodes_MPI
  IMPLICIT NONE

  PRIVATE :: psi0_gaussGridRep
CONTAINS
      SUBROUTINE init_psi0(WP0,para_WP0,mole)
      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_psi_Op
      USE mod_psi_B_TO_G
      USE mod_param_WP0
      USE mod_ana_psi
      USE mod_psi_io
      USE mod_WP0
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_WP0), intent(inout) :: para_WP0
      TYPE (param_psi), intent(inout) :: WP0(1)

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(in)    :: mole

!------ working variables ---------------------------------
      integer  :: nio,nb_WPdum
      integer  :: ecri_numi,ecri_nume

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING init_psi0'
        write(out_unit,*) 'WP0n_h,WP0nb_elec',para_WP0%WP0n_h,para_WP0%WP0nb_elec
        write(out_unit,*)
        write(out_unit,*) 'WP0(:)%symab',WP0(:)%symab

        CALL ecri_init_psi(WP0(1))

        write(out_unit,*)
        write(out_unit,*) 'WP0BasisRep',para_WP0%WP0BasisRep
        write(out_unit,*) 'lect_WP0GridRep,lect_WP0BasisRep',         &
                  para_WP0%lect_WP0GridRep,para_WP0%lect_WP0BasisRep
        write(out_unit,*)
        write(out_unit,*) 'nb_basis_act1',WP0(1)%BasisnD%nb_basis
        CALL RecWrite_basis(WP0(1)%BasisnD)
        write(out_unit,*)
        flush(out_unit)
      END IF

!-----------------------------------------------------------

      IF (para_WP0%New_Read_WP0) THEN
        IF(keep_MPI) CALL sub_read_psi0(WP0,para_WP0,max_WP=1)
      ELSE
        ecri_numi = para_WP0%WP0n_h
        ecri_nume = para_WP0%WP0nb_elec


        ! alloc the WP0 ----------------------------------------------------------------
        ! only allocated for MPI_id=0 inside the subroutine
        ! for psi%CvecB,psi%RvecB,psi%CvecG,psi%RvecG
        CALL alloc_psi(WP0(1))

        IF (.NOT. para_WP0%lect_WP0BasisRep .AND.                         &
            .NOT. para_WP0%WP0restart) THEN

          CALL alloc_psi(WP0(1),GridRep=.TRUE.)

          IF (para_WP0%lect_WP0GridRep) THEN
            !- read WP0 on the grid ----------------------------
            CALL lect_psiBasisRep(WP0(1)%CvecG,para_WP0%WP0cplx,          &
                             WP0(1)%nb_qa,WP0(1)%nb_bi,WP0(1)%nb_be,      &
                             para_WP0%WP0n_h,para_WP0%WP0nb_elec)
          ELSE IF (para_WP0%new_GWP0) THEN
            CALL New_WP0_GridRep(WP0(1),para_WP0,mole)
          ELSE
            !- initialisation with a gaussian WP ---------------
            CALL psi0_gaussGridRep(WP0(1),para_WP0,mole)
          END IF



          IF (debug) THEN
            write(out_unit,*) 'psiGridRep ini'
            CALL ecri_psi(ZERO,WP0(1),out_unit,.TRUE.,.FALSE.)
          END IF

          CALL norm2_psi(WP0(1),GridRep=.TRUE.)
          IF(openmpi) CALL MPI_Bcast_(WP0(1)%norm2,size1_MPI,root_MPI)
          write(out_unit,*) 'norm2WP GridRep',WP0(1)%norm2
          flush(out_unit)


          IF (WP0(1)%norm2 < ZERO) THEN
            write(out_unit,*) 'WARNING norm2WP < 0',WP0(1)%norm2
            write(out_unit,*) ' ... Bug in the Smolyak implementation!'
            write(out_unit,*) ' The renormalization is not performed'
          ELSE
            CALL renorm_psi_WITH_norm2(WP0(1),GridRep=.TRUE.)
            IF(openmpi) CALL MPI_Bcast_(WP0(1)%norm2,size1_MPI,root_MPI)
            write(out_unit,*) 'norm2WP GridRep',WP0(1)%norm2
            flush(out_unit)

            IF (debug) THEN
              write(out_unit,*) 'psiGridRep normalized'
              CALL ecri_psi(ZERO,WP0(1),out_unit,.TRUE.,.FALSE.)
            END IF
          END IF

          !- GridRep=>BasisRep -------------------------------------------------
          IF (para_WP0%WP0BasisRep) THEN

            CALL sub_PsiGridRep_TO_BasisRep(WP0(1))

            IF(keep_MPI) CALL norm2_psi(WP0(1),BasisRep=.TRUE.)
            IF(openmpi) CALL MPI_Bcast_(WP0(1)%norm2,size1_MPI,root_MPI)
            write(out_unit,*) 'norm2WP BasisRep',WP0(1)%norm2

            IF (debug) THEN
              write(out_unit,*) 'psiBasisRep unnormalized'
              CALL ecri_psi(ZERO,WP0(1),out_unit,.FALSE.,.TRUE.)
            END IF

            IF (abs(ONE-WP0(1)%norm2) >= ONETENTH**5) THEN
              write(out_unit,*) ' WARNNIG in psi0'
              write(out_unit,*) ' the transformation GridRep to BasisRep is NOT exact'
              write(out_unit,*) ' => used more basis functions'
              !CALL ecri_psi(psi=WP0(1),                                      &
              !              ecri_BasisRep=.TRUE.,ecri_GridRep=.TRUE.)
            END IF
            IF(keep_MPI) CALL renorm_psi(WP0(1),BasisRep=.TRUE.)
          END IF

        ELSE IF (para_WP0%WP0restart ) THEN
          !write(out_unit,*) 'WP0%cplx',WP0%cplx
          CALL file_open(para_WP0%file_WP0,nio)
          read(nio,*) nb_WPdum
          CALL lect_psiBasisRepnotall_nD(WP0(1),nio,WP0(1)%cplx,para_WP0%file_WP0%formatted)
          CALL file_close(para_WP0%file_WP0)

        ELSE IF (para_WP0%lect_WP0BasisRep .AND. para_WP0%lect_WP0BasisRepall) THEN

          IF (WP0(1)%nb_baie > WP0(1)%nb_tot) THEN
          CALL lect_psiBasisRep(WP0(1)%CvecB,para_WP0%WP0cplx,                  &
                           WP0(1)%nb_tot,1,1,                                   &
                           para_WP0%WP0n_h,para_WP0%WP0nb_elec)
          ELSE
          CALL lect_psiBasisRep(WP0(1)%CvecB,para_WP0%WP0cplx,                  &
                           WP0(1)%nb_ba,WP0(1)%nb_bi,WP0(1)%nb_be,                 &
                           para_WP0%WP0n_h,para_WP0%WP0nb_elec)
          END IF

        ELSE IF (para_WP0%lect_WP0BasisRep                                     &
                                .AND. .NOT. para_WP0%lect_WP0BasisRepall) THEN

          CALL lect_psiBasisRepnotall(WP0(1),para_WP0%WP0cplx)
        ELSE
          write(out_unit,*) ' ERROR in psi0'
          write(out_unit,*) ' I do not know what to do!!!'
          STOP
        END IF ! for .NOT. para_WP0%lect_WP0BasisRep

!     IF (para_WP0%WP0_DIP .GT. 0) THEN
!       IF (para_WP0%WP0BasisRep) THEN
!       ELSE
!       ENDIF
!     END IF

      END IF ! for para_WP0%New_Read_WP0

      IF(keep_MPI) CALL renorm_psi(WP0(1),BasisRep=.TRUE.)
      ! ENDIF

      IF(keep_MPI) write(out_unit,*) 'norm2WP BasisRep',WP0(1)%norm2

      !- clear WP0%...GridRep, if not need ------------------
      IF (para_WP0%WP0BasisRep) THEN
        WP0(1)%GridRep = .FALSE.
        CALL alloc_psi(WP0(1)) ! deallocate the grid representation
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'WP0(:)%symab',WP0(:)%symab
        write(out_unit,*) 'WP0BasisRep',WP0(1)%norm2
        IF (para_WP0%WP0BasisRep) THEN
          CALL ecri_psi(ZERO,WP0(1),out_unit,.FALSE.,.TRUE.)
        ELSE
          CALL ecri_psi(ZERO,WP0(1),out_unit,.TRUE.,.FALSE.)
        END IF
        write(out_unit,*) 'END psi0'
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE init_psi0
!==============================================================
!     Wave packet initialisation
!
!     for each active dimension :
!          psi0(qi) =  exp[-((Q-Qeq)/sigma)^2]*exp[i*imp_k*(Q-Qeq)+i*phase]
!
!     If WP0_DIP=1,2,3     => WP0 = WP0*dip(i) [i=1,2,3 => x,y,z]
!             =0 (default) => WP0 = WP0
!
!==============================================================
      SUBROUTINE psi0_gaussGridRep(WP0,para_WP0,mole)
      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_param_WP0
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_WP0), intent(in)      :: para_WP0
      TYPE (param_psi), intent(inout)   :: WP0


!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(in)      :: mole

      ! local variable
      integer            :: WP0nb_elec,WP0n_h

      real (kind=Rkind)      :: Qact(WP0%nb_act1)



      real (kind=Rkind)      :: z,ze,zk
      real (kind=Rkind)      :: czk,szk
      integer                :: i_act
      integer                :: i_qa,i_qaie
      TYPE (OldParam)        :: OldPara
!------ dipole moment ----------------------------------------------------
      real (kind=Rkind)      :: rhonD

!----- for debuging --------------------------------------------------
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING psi0_gaussGridRep'
        write(out_unit,*) 'nb_act1',WP0%nb_act1
        write(out_unit,*) 'nq_a,n_h,nb_elec',WP0%nb_qa,WP0%nb_bi,WP0%nb_be
        write(out_unit,*) 'WP0n_h,WP0nb_elec',para_WP0%WP0n_h,para_WP0%WP0nb_elec
        write(out_unit,*) 'sigma,Qeq,imp_k',                                   &
           para_WP0%WP0sigma,para_WP0%WP0Qeq,para_WP0%WP0imp_k
        write(out_unit,*) 'WP0_DIP',para_WP0%WP0_DIP
        write(out_unit,*)
        flush(out_unit)
      END IF
!-----------------------------------------------------------


      WP0n_h     = para_WP0%WP0n_h
      WP0nb_elec = para_WP0%WP0nb_elec
      IF (WP0n_h <1)     WP0n_h = 1
      IF (WP0nb_elec <1) WP0nb_elec = 1



!     - initisalisation ----------------------------------
      IF (WP0%cplx) THEN
        WP0%CvecG(:) = cmplx(ZERO,ZERO,kind=Rkind)
      ELSE
        WP0%RvecG(:) = ZERO
      END IF

      DO i_qa=1,WP0%nb_qa

!       - calculation of Qact -------------------------------
        CALL Rec_Qact(Qact,WP0%BasisnD,i_qa,mole,OldPara)
!       - calculation of WrhonD ------------------------------
        rhonD = Rec_rhonD(WP0%BasisnD,i_qa,OldPara)

!       ------------------------------------------------
        ze = ZERO
        zk = ZERO
        DO i_act=1,WP0%nb_act1
           z = (Qact(i_act)-para_WP0%WP0Qeq(i_act))/                    &
                                     para_WP0%WP0sigma(i_act)
           ze = ze + z*z
           zk = zk + para_WP0%WP0phase(i_act) +                         &
             (Qact(i_act)-para_WP0%WP0Qeq(i_act))*para_WP0%WP0imp_k(i_act)
        END DO
        zk = mod(zk,TWO*pi)
        ze = exp(-ze)/sqrt(sqrt(pi/TWO))**WP0%nb_act1
        ze = ze /sqrt(product(para_WP0%WP0sigma))
        IF (para_WP0%WP0nrho == 1) ze = ze /sqrt(rhonD)
        czk = cos(zk)
        szk = sin(zk)

        i_qaie = i_qa + ( (WP0n_h-1) +                                  &
               (WP0nb_elec-1)*WP0%nb_bi ) * WP0%nb_qa

        IF (WP0%cplx) THEN
          WP0%CvecG(i_qaie) = cmplx(czk*ze,szk*ze,kind=Rkind)
          IF (debug) write(out_unit,11) Qact,WP0%CvecG(i_qaie)
 11       format(10f15.5)
        ELSE
          WP0%RvecG(i_qaie) = ze
          IF (debug) write(out_unit,11) Qact,WP0%RvecG(i_qaie)
        END IF

      END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'END psi0_gaussGridRep'
         flush(out_unit)
       END IF
!-----------------------------------------------------------


      END SUBROUTINE psi0_gaussGridRep
      SUBROUTINE New_WP0_GridRep(WP0,para_WP0,mole)
      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_param_WP0
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_WP0), intent(in)      :: para_WP0
      TYPE (param_psi), intent(inout)   :: WP0


!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(in)      :: mole

      ! local variable
      real (kind=Rkind)      :: Qact(WP0%nb_act1)

      complex (kind=Rkind)   :: WP0_at_Q
      integer                :: iGWP,iQdyn,i
      integer                :: i_qa,i_qaie
      TYPE (OldParam)        :: OldPara
!------ dipole moment ----------------------------------------------------
      real (kind=Rkind)      :: rhonD

!----- for debuging --------------------------------------------------
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING New_WP0_GridRep'
        write(out_unit,*)
        CALL Write_Tab_GWP(para_WP0%tab_GWP0)
        write(out_unit,*) 'WP0nrho',para_WP0%WP0nrho
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      IF (WP0%cplx) THEN
        IF(keep_MPI) WP0%CvecG(:) = CZERO
      ELSE
        IF(keep_MPI) WP0%RvecG(:) = ZERO
      END IF

      DO i_qa=1,WP0%nb_qa

        !- calculation of Qact -------------------------------
        CALL Rec_Qact(Qact,WP0%BasisnD,i_qa,mole,OldPara)
        !- calculation of WrhonD ------------------------------
        IF (para_WP0%WP0nrho == 1) rhonD = Rec_rhonD(WP0%BasisnD,i_qa,OldPara)

        DO iGWP=1,size(para_WP0%tab_GWP0)
          WP0_at_Q = calc_GWP(para_WP0%tab_GWP0(iGWP),Qact)
          IF (para_WP0%WP0nrho == 1) WP0_at_Q = WP0_at_Q /sqrt(rhonD)

          i_qaie = i_qa + ( (para_WP0%tab_GWP0(iGWP)%I_HAChannel -1) +          &
               (para_WP0%tab_GWP0(iGWP)%I_ElecChannel-1)*WP0%nb_bi ) * WP0%nb_qa

          IF (WP0%cplx) THEN
            IF(keep_MPI) WP0%CvecG(i_qaie) = WP0%CvecG(i_qaie) + WP0_at_Q
          ELSE
            IF(keep_MPI) WP0%RvecG(i_qaie) = WP0%RvecG(i_qaie) + real(WP0_at_Q,kind=Rkind)
          END IF
        END DO

        IF (debug) THEN
          IF (WP0%cplx) THEN
            write(out_unit,*) 'i_qa,i_qaie,Qact',i_qa,i_qaie,Qact,WP0%CvecG(i_qaie)
          ELSE
            write(out_unit,*) 'i_qa,i_qaie,Qact',i_qa,i_qaie,Qact,WP0%RvecG(i_qaie)
          END IF
        END If

      END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'END New_WP0_GridRep'
         flush(out_unit)
       END IF
!-----------------------------------------------------------


END SUBROUTINE New_WP0_GridRep
END MODULE mod_psi
