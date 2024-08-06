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
MODULE mod_EVR
 USE mod_system
!$ USE omp_lib, only : omp_get_max_threads
 USE mod_Constant
 USE mod_Coord_KEO
 USE mod_PrimOp
 USE mod_basis

 USE mod_psi

 USE mod_propa
 USE mod_Op
 USE mod_analysis
 IMPLICIT NONE

 TYPE param_EVRT

   !----- physical and mathematical constants ---------------------------
   TYPE (constant) :: const_phys

   !----- for the CoordType and Tnum --------------------------------------
   TYPE (CoordType) :: mole
   TYPE (Tnum)      :: para_Tnum

   !----- for the basis set ----------------------------------------------
   TYPE (param_AllBasis) :: para_AllBasis

   !----- variables for the construction of H ----------------------------
   TYPE (param_AllOp)  :: para_AllOp


   !----- variables pour la namelist analyse ----------------------------
   TYPE (param_ana)           :: para_ana
   TYPE (param_intensity)     :: para_intensity
   logical                    :: intensity_only
   integer                    :: nio_res_int

   !----- variables for the WP propagation ----------------------------
   TYPE (param_propa) :: para_propa
   TYPE (param_psi)   :: WP0
 END TYPE param_EVRT

 TYPE (param_EVRT)              :: para_EVRT
 TYPE (param_EVRT), allocatable :: tab_EVRT(:) ! for the openmp


END MODULE mod_EVR

