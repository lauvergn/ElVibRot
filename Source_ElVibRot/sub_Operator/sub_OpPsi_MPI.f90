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

MODULE mod_OpPsi_MPI
  PUBLIC :: sub_scaledOpPsi_SR_MPI

  CONTAINS
!=======================================================================================
!> OpPsi = (OpPsi - E0*Psi) / Esc, for working on full Smolyak rep.
!=======================================================================================
    SUBROUTINE sub_scaledOpPsi_SR_MPI(Psi,OpPsi,E0,Esc)
      USE EVR_system_m
      USE mod_psi,ONLY:param_psi,ecri_psi
      IMPLICIT NONE

      TYPE(param_psi)                          :: OpPsi
      TYPE(param_psi)                          :: Psi
      Real(kind=Rkind)                         :: E0
      Real(kind=Rkind)                         :: Esc

      IF(Psi%SRG_MPI .AND. OpPsi%SRG_MPI) THEN
        OpPsi%SR_G(:,:)=(OpPsi%SR_G(:,:)-E0*Psi%SR_G(:,:))/Esc
      ELSE IF(Psi%SRB_MPI .AND. OpPsi%SRB_MPI) THEN
        OpPsi%SR_B(:,:)=(OpPsi%SR_B(:,:)-E0*Psi%SR_B(:,:))/Esc
      ELSE
        STOP 'error in sub_scaledOpPsi_SR_MPI'
      ENDIF
      
      IF(OpPsi%symab/=Psi%symab) THEN
        IF(OpPsi%symab==-2) THEN
          OpPsi%symab=Psi%symab
        ELSE
          OpPsi%symab=-1
        ENDIF
      ENDIF

    ENDSUBROUTINE sub_scaledOpPsi_SR_MPI
!=======================================================================================

END MODULE mod_OpPsi_MPI
