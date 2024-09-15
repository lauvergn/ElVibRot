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
!     RoVibrational levels calculation
!================================================================
      SUBROUTINE sub_VibRot(Tab_Psi,nb_psi,para_H,para_ana)
      USE EVR_system_m
      USE mod_Constant
      use mod_Coord_KEO, only: CoordType, tnum
      use mod_PrimOp,    only: param_d0matop, init_d0matop, write_d0matop, dealloc_d0matop
      USE mod_basis
      USE mod_psi,       ONLY : param_psi,Overlap_psi1_psi2
      USE mod_Op
      USE mod_analysis
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      integer            :: nb_psi
      TYPE (param_psi)   :: Tab_Psi(nb_psi)

!----- Operator variables ----------------------------------------------
      TYPE (param_Op)  :: para_H
      logical          :: print_Op

!----- variables pour la namelist analyse ------------------------------
      TYPE (param_ana)           :: para_ana



!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum

!----- working variables -----------------------------
      TYPE (param_d0MatOp)              :: MatRV
      TYPE (param_psi)                  :: OpPsi
      integer                           :: iv,jv

      integer           :: nb_bVR
      real (kind=Rkind), allocatable ::    H_VR(:,:)
      real (kind=Rkind), allocatable ::    Ene_VR(:)
      real (kind=Rkind), allocatable ::    Vec_VR(:,:)
      real (kind=Rkind) ::    rho_V(nb_psi)

      integer              :: nb_bRot,JRot,iterm_Op,iterm_BasisRot
      integer              :: J1,J2,i1,i2,f1,f2,ibRot,jbRot,nb_ana
      integer              :: i,jR,nb_shift,type_Op

      real (kind=Rkind)    :: Val_BasisRot,non_hermitic,auTOcm_inv
      complex (kind=Rkind) :: C_over



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_VibRot'
!-----------------------------------------------------------
      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum

      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nb_psi',nb_psi

      IF (debug) THEN
        write(out_unit,*)
      END IF
      flush(out_unit)

      !-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
      JRot = para_ana%JJmax
      para_H%Mat_done = .FALSE.

      type_Op = para_H%para_ReadOp%Type_HamilOp ! H
      IF (type_Op /= 1) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '    Type_HamilOp MUST be equal to 1 here!!'
        write(out_unit,*) '    CHECK your data!!'
        STOP
      END IF

      CALL Init_d0MatOp(MatRV,type_Op,0,nb_psi,iQact=0,JRot=JRot,cplx=para_H%cplx)
      para_H%Make_Mat = .FALSE.
      DO iterm_Op=1,MatRV%nb_term
        DO iv=1,nb_psi
          IF (debug) THEN
            write(out_unit,*) '======================================='
            write(out_unit,*) '======================================='
            write(out_unit,*) '======================================='
            write(out_unit,*) "iterm_Op,iv",iterm_Op,iv
            flush(out_unit)
          END IF

          CALL sub_OpPsi(Tab_Psi(iv),OpPsi,para_H,MatRV%derive_termQact(:,iterm_Op))
          DO jv=1,nb_psi
            CALL Overlap_psi1_psi2(C_over,Tab_Psi(jv),OpPsi)
            !write(out_unit,*) 'jv,iv,C_over',jv,iv,C_over
            MatRV%ReVal(jv,iv,iterm_Op) = real(C_over,kind=Rkind)
            IF (MatRV%cplx) MatRV%ImVal(jv,iv) = aimag(C_over)
          END DO
        END DO
      END DO
      IF (debug) CALL Write_d0MatOp(MatRV)
      flush(out_unit)

      !SET Rotational basis
      CALL init_RotBasis_Param(para_H%BasisnD%RotBasis,Jrot)
      IF (debug) CALL Write_RotBasis_Param(para_H%BasisnD%RotBasis)
      flush(out_unit)

      nb_bRot = para_H%BasisnD%RotBasis%nb_Rot
      nb_bVR  = nb_psi*nb_bRot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL alloc_NParray(H_VR,  [nb_bVR,nb_bVR],'H_VR',  name_sub)
      CALL alloc_NParray(Vec_VR,[nb_bVR,nb_bVR],'Vec_VR',name_sub)
      CALL alloc_NParray(Ene_VR,[nb_bVR],       'Ene_VR',name_sub)

      H_VR(:,:) = ZERO

      DO iterm_Op=1,MatRV%nb_term

          ! Rotational contribution
          J1       = MatRV%derive_termQact(1,iterm_Op)
          J2       = MatRV%derive_termQact(2,iterm_Op)
          iterm_BasisRot = para_H%BasisnD%RotBasis%tab_der_TO_iterm(J1,J2)
          !write(out_unit,*) 'J1,J2',J1,J2,'iterm_Op,iterm_BasisRot',iterm_Op,iterm_BasisRot

          DO ibRot=1,nb_bRot
          DO jbRot=1,nb_bRot
            Val_BasisRot = para_H%BasisnD%RotBasis%tab_RotOp(ibRot,jbRot,iterm_BasisRot)
            IF (abs(Val_BasisRot) < ONETENTH**9) CYCLE

            i1 = (ibRot-1)*nb_psi
            i2 = (jbRot-1)*nb_psi
            f1 = nb_psi
            f2 = nb_psi

            !IF (debug) THEN
            !  write(out_unit,*) 'J1,J2',J1,J2
            !  write(out_unit,*) 'ibRot,i1+1:i1+f1',ibRot,i1+1,i1+f1
            !  write(out_unit,*) 'jbRot,i2+1:i2+f2',jbRot,i2+1,i2+f2
            !END IF

            H_VR(i1+1:i1+f1 , i2+1:i2+f2) =                             &
                                        H_VR(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
          END DO
          END DO
          !---- END LOOP on the rotational basis function


      END DO

      IF (debug) THEN
        write(out_unit,*) 'H_VR'
        CALL Write_Mat(H_VR,out_unit,5)
      END IF

      CALL sub_hermitic_H(H_VR,nb_bVR,non_hermitic,para_H%sym_Hamil)

      IF (non_hermitic >= FOUR/TEN**4) THEN
        If(MPI_id==0) write(out_unit,*) 'WARNING: non_hermitic is BIG'
        If(MPI_id==0) write(out_unit,31) non_hermitic
 31     format(' Hamiltonien: ',f16.12,' au')
      ELSE
        IF (print_level>-1 .AND. MPI_id==0) write(out_unit,21) non_hermitic*auTOcm_inv
 21     format(' Hamiltonien: ',f16.12,' cm-1')
      END IF
      flush(out_unit)


      IF (para_H%sym_Hamil) THEN
        CALL diagonalization(H_VR,Ene_VR,Vec_VR,nb_bVR,3,1,.TRUE.)
      ELSE
        CALL diagonalization(H_VR,Ene_VR,Vec_VR,nb_bVR,4,1,.TRUE.)
      END IF
      nb_shift = count(Ene_VR(:) <= para_H%ZPE)
      IF (nb_shift > 0) THEN
        If(MPI_id==0) write(out_unit,*) 'WARNING the vectors 1 to ',nb_shift,'have negative energies',Ene_VR(1:nb_shift)
        If(MPI_id==0)  write(out_unit,*) '=> They will be shifted'
        Vec_VR = cshift(Vec_VR,shift=nb_shift,dim=2)
        Ene_VR = cshift(Ene_VR,shift=nb_shift)
      END IF

      IF (debug) THEN
        write(out_unit,*) 'Vec_VR (in column)'
        CALL Write_Mat(Vec_VR,out_unit,5)
      END IF


      write(out_unit,*) 'ZPE',para_H%ZPE*auTOcm_inv
      nb_ana = min(nb_bVR,nb_bRot*2)
      write(out_unit,'(A,i4,30f15.6)') 'Ene RV',JRot,                  &
                         (Ene_VR(1:nb_bRot)-para_H%ZPE)*auTOcm_inv
      write(out_unit,'(A,i4,30f15.6)') 'Ene RV',JRot,                  &
          (Ene_VR(nb_bRot+1:nb_ana)-real(Tab_psi(2)%CAvOp,kind=Rkind))* &
                                                             auTOcm_inv

      write(out_unit,*) 'Ene RV (all), J:',JRot
      DO i=1,min(nb_bVR,10*nb_bRot)
        write(out_unit,'(A,i5,f15.6)') ' levR:',i,                     &
                                 (Ene_VR(i)-para_H%ZPE)*auTOcm_inv
        rho_V(:) = ZERO
        DO jR=1,nb_bRot
          i1 = (jR-1)*nb_psi + 1
          f1 = (jR-1)*nb_psi + nb_psi
          rho_V(:) = rho_V(:) + abs(Vec_VR( i1:f1 ,i))**2
        END DO
        write(out_unit,'(A,10f5.2)') ' densVib:',rho_V(1:min(10,nb_psi))

      END DO
      flush(out_unit)

      CALL dealloc_NParray(H_VR,  'H_VR',  name_sub)
      CALL dealloc_NParray(Vec_VR,'Vec_VR',name_sub)
      CALL dealloc_NParray(Ene_VR,'Ene_VR',name_sub)

      CALL dealloc_d0MatOp(MatRV)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)

!----------------------------------------------------------

      end subroutine sub_VibRot
