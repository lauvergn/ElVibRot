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
      SUBROUTINE ElVibRot_version(write_version)
      USE mod_system
      IMPLICIT NONE

      logical :: write_version

      character (len=*), parameter :: EVR_name='ElVibRot'



      IF (write_version .AND. MPI_id==0) THEN
        write(out_unitp,*) '==============================================='
        write(out_unitp,*) '==============================================='
        write(out_unitp,*) 'Working with ',EVR_name,trim(adjustl(EVR_version))

        write(out_unitp,*) 'Compiled on "',trim(compile_host), '" the ',trim(compile_date)
        write(out_unitp,*) 'Compiler version: ',trim(compiler_ver)
        write(out_unitp,*) 'Compiler options: ',trim(compiler_opt)
        write(out_unitp,*) 'Compiler libs: ',trim(compiler_libs)

        write(out_unitp,*) 'EVRT_path: ',trim(EVRT_path)
        write(out_unitp,*) 'git ',trim(git_branch)

        write(out_unitp,*) '-----------------------------------------------'

        write(out_unitp,*) EVR_name,' is written by David Lauvergnat [1] '
        write(out_unitp,*) '  with contributions of'
        write(out_unitp,*) '     Josep Maria Luis (optimization) [2]'
        write(out_unitp,*) '     Ahai Chen (MPI) [1,4]'
        write(out_unitp,*) '     Lucien Dupuy (CRP) [5]'

        write(out_unitp,*) EVR_name,' is under MIT license.'
        write(out_unitp,*)

        write(out_unitp,*)
        write(out_unitp,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
        write(out_unitp,*) '[2]: Institut de Química Computacional and Departament de Química',&
                                   ' Universitat de Girona, Catalonia, Spain'
        write(out_unitp,*) '[4]: Maison de la Simulation USR 3441, CEA Saclay, France'
        write(out_unitp,*) '[5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,', &
                                   ' Université de Montpellier, France'
        write(out_unitp,*) '==============================================='
        write(out_unitp,*) '==============================================='
      END IF
      END SUBROUTINE ElVibRot_version
