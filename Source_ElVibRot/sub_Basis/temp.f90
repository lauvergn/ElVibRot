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
SUBROUTINE OpRV_TO_RV_basis(Vin,Vout,lBin,lBout,base,der)
  USE EVR_system_m
  IMPLICIT NONE
    
  real (kind=Rkind), intent(in)           :: Vin(:)
  real (kind=Rkind), intent(inout)        :: Vout(:)
  logical,           intent(in)           :: lBin,lBout

  TYPE (basis),      intent(in), target   :: base
  integer,           intent(in), optional :: der(2)

  integer :: nb,nq,n_Vin,n_Vout,der_basis(2)
  real (kind=Rkind), pointer       :: Mat(:,:)

  !----- for debuging --------------------------------------------------
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------
  nq     = get_nq_FROM_basis(base)
  nb     = get_nb_FROM_basis(base)
  n_Vin  = size(Vin)
  n_Vout = size(Vout)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING OpRV_TO_RV_basis'
    write(out_unit,*) 'nb,nq',nb,nq
    write(out_unit,*) 'n_Vin',n_Vin
    write(out_unit,*) 'n_Vout',n_Vout
    write(out_unit,*) 'Vin',Vin(:)
  END IF
  !-----------------------------------------------------------

  IF (NewBasisEl .AND. base%ndim == 0) THEN
    G(:) = B
  ELSE
    IF (nb_basis < nb_B .OR. nq /= size(G)) THEN
      write(out_unit,*) ' ERROR in RB_TO_RG_basis'
      write(out_unit,*) ' nb_basis is inconsistent with nb_B',nb_basis,nb_B
      write(out_unit,*) ' nq is inconsistent with size(G)',nq,size(G)
      STOP ' ERROR in RB_TO_RG_basis: inconsistent sizes'
    END IF
    nb_mult_BTOG = nb_mult_BTOG + int(nb_B*nq,kind=ILkind)

    IF (present(der)) THEN
      der_basis = base%Tabder_Qdyn_TO_Qbasis(der(:))
    ELSE
      der_basis = 0
    END IF

    IF (base%dnBBRep_done) THEN
      IF (der_basis(1) == 0 .AND. der_basis(2) == 0) THEN
        Binter = B
      ELSE IF (der_basis(1) > 0 .AND. der_basis(2) == 0) THEN
        Mat => base%dnRBB%d1(:,:,der_basis(1))
        Binter = matmul(Mat,B)
      ELSE IF (der_basis(1) == 0 .AND. der_basis(2) > 0) THEN
        Mat => base%dnRBB%d1(:,:,der_basis(2))
        Binter = matmul(Mat,B)
      ELSE ! der_basis(1) > 0 and der_basis(2) > 0
        Mat => base%dnRBB%d2(:,:,der_basis(1),der_basis(2))
        Binter = matmul(Mat,B)
      END IF

      Mat => base%dnRGB%d0(:,1:nb_B)
      G = matmul(Mat,Binter)

      deallocate(Binter)

    ELSE
      IF (der_basis(1) == 0 .AND. der_basis(2) == 0) THEN
        Mat => base%dnRGB%d0(:,1:nb_B)
      ELSE IF (der_basis(1) > 0 .AND. der_basis(2) == 0) THEN
        Mat => base%dnRGB%d1(:,1:nb_B,der_basis(1))
      ELSE IF (der_basis(1) == 0 .AND. der_basis(2) > 0) THEN
        Mat => base%dnRGB%d1(:,1:nb_B,der_basis(2))
      ELSE ! der_basis(1) > 0 and der_basis(2) > 0
        Mat => base%dnRGB%d2(:,1:nb_B,der_basis(1),der_basis(2))
      END IF

      G = matmul(Mat,B)
    END IF

  END IF

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'Vout',Vout(:)
    write(out_unit,*) 'END OpRV_TO_RV_basis'
  END IF
  !-----------------------------------------------------------
END SUBROUTINE OpRV_TO_RV_basis