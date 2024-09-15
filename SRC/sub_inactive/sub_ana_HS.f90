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
!=====================================================================
!
!      hermiticity analysis
!      and symetrisation
!
!=====================================================================
      SUBROUTINE sub_hermitic_H(H,nb_bases,non_hermitic,sym)
      USE EVR_system_m
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer   nb_bases
      real (kind=Rkind) ::    H(nb_bases,nb_bases)
      real (kind=Rkind) ::    non_hermitic
      logical           ::    sym


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING sub_hermitic_H'
      END IF
!-----------------------------------------------------------

      non_hermitic = maxval(abs(H-transpose(H)))*HALF
      IF (sym) H = (H+transpose(H))*HALF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'non_hermitique H',non_hermitic
        write(out_unit,*) 'END sub_hermitic_H'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_hermitic_H
!=====================================================================
!
!      hermiticity analysis
!      and symetrisation
!
!=====================================================================
      SUBROUTINE sub_hermitic_cplxH(H,nb_bases,non_hermitic,sym)
      USE EVR_system_m
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer       nb_bases
      complex (kind=Rkind) ::    H(nb_bases,nb_bases)
      real (kind=Rkind)    ::    val,non_hermitic
      logical              ::    sym


      integer   i,j

!----- for debuging --------------------------------------------------
      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING sub_hermitic_cplxH'
      END IF
!-----------------------------------------------------------

      non_hermitic = ZERO
      DO i=1,nb_bases
        DO j=i+1,nb_bases

          val = abs( H(i,j) - H(j,i) )
          IF (sym) THEN
            H(i,j) = (H(i,j) + H(j,i) )*cmplx(HALF,ZERO,kind=Rkind)
            H(j,i) = H(i,j)
          END IF

          IF ( val > non_hermitic) non_hermitic = val
!         write(out_unit,*) 'sub_hermitic_cplxH',val,non_hermitic
        END DO
      END DO

      non_hermitic = non_hermitic * HALF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'non_hermitique cplxH',non_hermitic
        write(out_unit,*) 'END sub_hermitic_cplxH'
      END IF
!-----------------------------------------------------------

      RETURN
      end subroutine sub_hermitic_cplxH
!=====================================================================
!
!       analysis of the overlap matrix (upper part)
!
!=====================================================================
      SUBROUTINE sub_ana_S(S,nb_bases,max_Sii,max_Sij,write_maxS)
      USE EVR_system_m
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer   nb_bases
      real (kind=Rkind) :: S(nb_bases,nb_bases)
      real (kind=Rkind) :: max_Sii,max_Sij
      logical           :: write_maxS


      integer   i,j

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING sub_ana_S'
        CALL Write_Mat(S,out_unit,5)
      END IF
!-----------------------------------------------------------

!     - analysis of the overlap matrix
      max_Sii = ZERO
      max_Sij = ZERO
      DO i=1,nb_bases
        max_Sii=max(max_Sii,abs(S(i,i)-ONE))
        DO j=i+1,nb_bases
          max_Sij = max(max_Sij,abs(S(i,j)))
        END DO
      END DO

      IF (write_maxS) THEN
         write(out_unit,"(' Max Overlap:',2e11.3)") max_Sii,max_Sij
         flush(out_unit)
      END IF
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) ' Max Overlap:',max_Sii,max_Sij
        write(out_unit,*) 'END sub_ana_S'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_ana_S
!=====================================================================
!
!       analysis of the overlap matrix (upper part)
!
!=====================================================================
      SUBROUTINE sub_ana_cplxS(S,nb_bases,max_Sii,max_Sij,write_maxS)
      USE EVR_system_m
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer              :: nb_bases
      complex (kind=Rkind) :: S(nb_bases,nb_bases)
      real (kind=Rkind)    :: max_Sii,max_Sij
      logical              :: write_maxS


      integer   :: i,j

!----- for debuging --------------------------------------------------
      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING sub_ana_cplxS'
      END IF
!-----------------------------------------------------------

!     - analysis of the overlap matrix
      max_Sii = ZERO
      max_Sij = ZERO
      DO i=1,nb_bases
        max_Sii=max(max_Sii,abs(S(i,i)-CONE))
        DO j=i+1,nb_bases
          max_Sij = max(max_Sij,abs(S(i,j)))
        END DO
      END DO

      IF (write_maxS) THEN
         write(out_unit,"(' Max Overlap:',2e11.3)") max_Sii,max_Sij
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) max_Sii,max_Sij
        write(out_unit,*) 'sub_ana_cplxS'
      END IF
!-----------------------------------------------------------

      RETURN
      end subroutine sub_ana_cplxS

