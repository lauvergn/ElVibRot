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
!=====================================================================
!
!  Op(Q).d0b(Q,i).W(Q) and d0b(Q,i) calculations for the nD quadrature point k
!
!=====================================================================
      SUBROUTINE calc_td0b_OpRVd0bW(iq,k,td0b,d0MatOpd0bWrho,WnD,kmem,  &
                                    d0MatOp,para_Op,BasisnD)
      USE mod_system
      USE mod_PrimOp, only: Write_d0MatOp
      USE mod_basis
      USE mod_SetOp
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (Basis) :: BasisnD

!----- Operator variables --------------------------------------------
      TYPE (param_Op)      :: para_Op
      integer              :: k,kmem,iq
      real (kind=Rkind)    :: td0b(para_Op%nb_ba,kmem)
      TYPE (param_d0MatOp) :: d0MatOpd0bWrho(kmem,para_Op%nb_ba)


      TYPE (param_d0MatOp) :: d0MatOp

!------ quadrature weights ----------------------------------------
      real (kind=Rkind) :: WnD


!------ working variables ---------------------------------
      integer  :: ib
      integer  :: j_act,i_act,i,j

      real (kind=Rkind) :: d0bnD
      real (kind=Rkind) :: d1bnD(para_Op%nb_act1)
      real (kind=Rkind) :: d2bnD(para_Op%nb_act1,para_Op%nb_act1)

      real (kind=Rkind) :: dnbnD

      integer ::  i_term,i_Op

!----- for debuging --------------------------------------------------
       character (len=*), parameter :: name_sub = 'calc_td0b_OpRVd0bW'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq,k,kmem',iq,k,kmem
         write(out_unitp,*) 'nb_ba',para_Op%nb_ba
         write(out_unitp,*) 'nb_act1',para_Op%nb_act1
         write(out_unitp,*)
         write(out_unitp,*) 'WnD',WnD
         ! CALL write_param_Op(para_Op)
          write(out_unitp,*) 'd0MatOp:'
          CALL Write_d0MatOp(d0MatOp)
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      IF (d0MatOp%cplx) THEN
        d0MatOp%ImVal(:,:)    = d0MatOp%ImVal(:,:) * WnD
      END IF
      d0MatOp%ReVal(:,:,:)    = d0MatOp%ReVal(:,:,:) * WnD

      DO ib=1,para_Op%nb_ba

        CALL d0d1d2bnDQact(d0bnD,d1bnD,d2bnD,BasisnD,iq,ib,para_Op%mole)
        td0b(ib,k)            = d0bnD

        !- initialisation ----
        d0MatOpd0bWrho(k,ib)%ReVal(:,:,:) = ZERO

        DO i_term=1,d0MatOp%nb_term
          i_act = d0MatOp%derive_termQact(1,i_term)
          j_act = d0MatOp%derive_termQact(2,i_term)

          !- 2d order derivatives ------------------
          IF ( j_act > 0 .AND. i_act > 0 ) THEN
            dnbnD = d2bnD(j_act,i_act)

          !- first order derivatives ---------------
          ELSE IF ( j_act > 0 .AND. i_act <= 0 ) THEN
            dnbnD = d1bnD(j_act)

          ELSE IF ( j_act <= 0 .AND. i_act > 0 ) THEN
            dnbnD = d1bnD(i_act)

          !- no derivative  of deformation part ----------------------
          ELSE
            dnbnD = d0bnD

          END IF

          i    = min(0,i_act)
          j    = min(0,j_act)
          i_Op = d0MatOpd0bWrho(k,ib)%derive_term_TO_iterm(j,i)

          d0MatOpd0bWrho(k,ib)%ReVal(:,:,i_Op) =                        &
                                d0MatOpd0bWrho(k,ib)%ReVal(:,:,i_Op) +  &
                                       dnbnD * d0MatOp%ReVal(:,:,i_term)
        END DO

        IF (d0MatOp%cplx) THEN
          d0MatOpd0bWrho(k,ib)%ImVal(:,:) = d0bnD * d0MatOp%ImVal(:,:)
        END IF

      END DO
      !-----------------------------------------------------------------


!-----------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' td0b(:,k)',k,para_Op%nb_ba
        CALL Write_Vec(td0b(:,k),out_unitp,8)
        write(out_unitp,*)
        write(out_unitp,*) ' d0MatOpd0bWrho(:,:)'
        DO i=1,ubound(d0MatOpd0bWrho,dim=2)
          write(out_unitp,*) 'k,i',k,i
          CALL Write_d0MatOp(d0MatOpd0bWrho(k,i))
        END DO
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF

      END SUBROUTINE calc_td0b_OpRVd0bW
