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
!
!================================================================
!    read the parameters of the fits
!================================================================
      SUBROUTINE read_para0d(F,nn,max_points,nom1,exist)
      USE EVR_system_m
      IMPLICIT NONE

       integer :: nn
       character (len=*) :: nom1
       integer :: max_points
       real (kind=Rkind) :: F(max_points)
       logical :: exist

       integer :: no,ios,kl

       write(out_unit,*) nom1
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)

       IF (ios .EQ. 0) THEN
         read(no,*) nn
         IF (nn .GT. max_points) THEN
           write(out_unit,*) ' ERROR : nb de points du fit (',nn,') >'
           write(out_unit,*) '         a max_points (',max_points,')'
           write(out_unit,*) '         STOP in read_para0d'
           STOP
         END IF
         DO kl=1,nn
           read(no,*) F(kl)
!          write(out_unit,*) F(kl)
         END DO
         CLOSE(no) ! CALL file_close cannot be used
         exist = .TRUE.
       ELSE
         write(out_unit,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


!      write(out_unit,*) 'nom1,exist ',nom1,exist
!      write(out_unit,*) 'nn,max_points',nn,max_points

       RETURN
       end subroutine read_para0d
!
!================================================================
!    read the parameters of the fits
!================================================================
      SUBROUTINE read_para1d(F,nn,max_points,nb_fit,nom1,exist)
      USE EVR_system_m
      IMPLICIT NONE

       integer :: nn,nb_fit
       character (len=*) :: nom1
       integer :: max_points
       real (kind=Rkind) :: F(max_points)
       logical :: exist

       integer :: no,ios,kl,i

!      write(out_unit,*) nom1
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)

       IF (ios .EQ. 0) THEN
         DO i=1,nb_fit
           read(no,*) nn
           IF (nn .GT. max_points) THEN
             write(out_unit,*) ' ERROR : nb de points du fit (',nn,') >'
             write(out_unit,*) '         a max_points (',max_points,')'
             write(out_unit,*) '         STOP in read_para0d'
             STOP
           END IF
           DO kl=1,nn
            read(no,*) F(kl)
!           write(out_unit,*) F(kl)
           END DO
         END DO
         CLOSE(no)  ! CALL file_close cannot be used
         exist = .TRUE.
       ELSE
         write(out_unit,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


!      write(out_unit,*) 'nom1,exist ',nom1,exist
!      write(out_unit,*) 'nn,max_points',nn,max_points

       RETURN
       end subroutine read_para1d
!
!================================================================
!    read parameters for the fit
!    (several fits)
!================================================================
      SUBROUTINE read_para2d(F,nn,nb_fit,max_fit,max_points,nom1,exist)
      USE EVR_system_m
      IMPLICIT NONE

       integer :: max_points,max_fit,nb_fit
       integer :: nn(max_fit)
       character (len=*) :: nom1
       real (kind=Rkind) :: F(max_points,max_fit)
       logical :: exist

       integer :: no,ios,kl,i

       write(out_unit,*) nom1,nb_fit,max_fit,max_points
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)
       IF (ios .EQ. 0) THEN

         IF (nb_fit == 0) read(no,*) nb_fit

         IF (nb_fit > max_fit) THEN
           write(out_unit,*) ' ERROR in read_para2d'
           write(out_unit,*) ' The index nfit (',nb_fit,                &
                     ') is greater than max_fit',max_fit
           STOP
         END IF
         write(out_unit,*) nom1,nb_fit
         IF (nb_fit <= 0) THEN
           write(out_unit,*) ' ERROR in read_para2d'
           write(out_unit,*) ' nb_fit is < 1 !!',nb_fit,' in the file :',nom1
           STOP
         END IF
         DO i=1,nb_fit
           read(no,*) nn(i)
           write(out_unit,*) 'nom1,nb_fit,i,nn ',nom1,nb_fit,i,nn(i)
           IF (nn(i) .GT. max_points) THEN
             write(out_unit,*) ' ERROR : nb de points du fit (',nn(i),') >'
             write(out_unit,*) '         a max_points (',max_points,')'
             write(out_unit,*) '         STOP in read_para2d'
             STOP
           END IF
           DO kl=1,nn(i)
            read(no,*) F(kl,i)
!           write(out_unit,*) F(kl,i)
           END DO
         END DO
         CLOSE(no) ! CALL file_close cannot be used
         exist = .TRUE.
       ELSE
         write(out_unit,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


!      write(out_unit,*) 'nom1,exist ',nom1,exist
!      write(out_unit,*) 'nn,max_points',nn,max_points

       end subroutine read_para2d
!
!================================================================
!    read the parameters of the fits
!================================================================
      SUBROUTINE read_para3d(F,n,ndim,nb_fit,max_fit,max_points,        &
                              nom1,exist)
      USE EVR_system_m
      IMPLICIT NONE

       integer :: max_points,max_fit,ndim
       integer :: n(0:ndim,max_fit),nb_fit
       character (len=*) :: nom1
       real (kind=Rkind) :: F(max_points,max_fit)
       logical :: exist

       integer :: no,ios,kl,i

       write(out_unit,*) 'read_para3d: nom1,nb_fit,max_fit,max_points: ',&
                                      nom1,nb_fit,max_fit,max_points
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)

       IF (ios == 0) THEN

         IF (nb_fit == 0) read(no,*) nb_fit

         IF (nb_fit > max_fit) THEN
           write(out_unit,*) ' ERROR in read_para3d'
           write(out_unit,*) ' The index nfit (',nb_fit,               &
                     ') is greater than max_fit',max_fit
           STOP
         END IF
         write(out_unit,*) 'nom1,nb_fit,ndim: ',nom1,nb_fit,ndim
         IF (nb_fit <= 0) THEN
           write(out_unit,*) ' ERROR in read_para3d'
           write(out_unit,*) ' nb_fit is < 1 !!',nb_fit,' in the file :',nom1
           STOP
         END IF
         DO i=1,nb_fit
           read(no,*) n(0:ndim,i)
           write(out_unit,*) 'nom1,nb_fit,i,n ',nom1,nb_fit,i,n(0:ndim,i)
           IF (n(0,i) > max_points) THEN
             write(out_unit,*) ' ERROR : nb de points du fit (',n(0,i),') >'
             write(out_unit,*) '         a max_points (',max_points,')'
             write(out_unit,*) '         STOP in read_para3d'
             STOP
           END IF
           DO kl=1,n(0,i)
            read(no,*) F(kl,i)
!           write(out_unit,*) F(kl,i)
           END DO
         END DO
         CLOSE(no) ! CALL file_close cannot be used
         exist = .TRUE.
       ELSE
         write(out_unit,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF

       end subroutine read_para3d
!================================================================
!    read the parameters of the fits
!    + the type of the functions
!================================================================
      SUBROUTINE read_para4d(F,n,ndim,nt,max_points,nom1,exist)
      USE EVR_system_m
      IMPLICIT NONE

       integer :: max_points,ndim,nt
       integer :: n(0:ndim)
       character (len=*) :: nom1
       real (kind=Rkind) :: F(max_points)
       logical :: exist

       integer :: no,ios,kl,i

       write(out_unit,*) 'read_para4d: nom1,max_points: ',nom1,max_points

       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)
       IF (ios == 0) THEN

         read(no,*) nt
         read(no,*) i ! for nb_fit (not used)

         write(out_unit,*) 'nom1,nt,ndim: ',nom1,nt,ndim
         read(no,*) n(0:ndim)
         write(out_unit,*) 'nom1,n ',nom1,n(0:ndim)
         IF (n(0) > max_points) THEN
             write(out_unit,*) ' ERROR : nb de points du fit (',n(0),') >'
             write(out_unit,*) '         a max_points (',max_points,')'
             write(out_unit,*) '         STOP in read_para4d'
             STOP
           END IF
           DO kl=1,n(0)
            read(no,*) F(kl)
!           write(out_unit,*) F(kl)
           END DO
         CLOSE(no) ! CALL file_close cannot be used
         exist = .TRUE.
       ELSE
         write(out_unit,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


       end subroutine read_para4d

