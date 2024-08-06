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
!     NLO
!================================================================
      SUBROUTINE sub_NLO(para_Dip,print_Op,para_H,nb_ana,para_intensity)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      USE mod_analysis
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      TYPE (param_Op)  :: para_Dip(3),para_H
      logical          :: print_Op
      integer          :: nb_ana

!----- variables pour la namelist analyse ------------------------------
      TYPE (param_intensity) :: para_intensity



!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      real (kind=Rkind), allocatable ::    Mat_Aif(:,:)
      integer       ::    i,k



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_NLO'
!-----------------------------------------------------------
      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum

      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unitp,*) 'nb_ana',nb_ana
        write(out_unitp,*) 'Rvp',shape(para_H%Rvp)
!       CALL Write_Mat(para_H%Rvp,out_unitp,5)
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------


      IF (nb_ana > 0) THEN
        CALL alloc_NParray(Mat_Aif,[nb_ana,nb_ana],'Mat_Aif',name_sub)
        Mat_Aif(:,:) = ZERO
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_ana <=0',nb_ana
        STOP
      END IF
      flush(out_unitp)



!     -------------------------------------------
      write(out_unitp,*) ' ene (ua): ',nb_ana
      DO i=1,nb_ana
        write(out_unitp,*) i,para_H%Rdiag(i)
      END DO
      write(out_unitp,*) ' END ene',nb_ana
      write(out_unitp,*)

      write(out_unitp,*) '==================================================='
      write(out_unitp,*) '==================================================='
      write(out_unitp,*) ' Calculation of "Mat_Aif(:,:)": '

      Mat_Aif(:,:) = ZERO
      DO k=1,3
        Mat_Aif(:,:) = Mat_Aif(:,:) + para_Dip(k)%Rmat(1:nb_ana,1:nb_ana)**2
      END DO
      CALL Write_Mat(Mat_Aif,out_unitp,5,Rformat='e30.23')
      write(out_unitp,*) '==================================================='
      write(out_unitp,*) '==================================================='
      flush(out_unitp)


      CALL dealloc_NParray(Mat_Aif,'Mat_Aif',name_sub)


!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
!----------------------------------------------------------

      end subroutine sub_NLO
