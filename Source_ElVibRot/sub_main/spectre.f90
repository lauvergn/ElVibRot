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
  PROGRAM spectre
      implicit none

      integer :: npts,nptE
      real (kind=8) :: t0,tmax,dt
      real (kind=8) :: E,Emin,Emax,dE
      real (kind=8) :: conv
      real (kind=8), pointer :: t(:)
      complex (kind=8), pointer :: auto(:)
      complex (kind=8), pointer :: Expiwt(:)
      complex (kind=8), pointer :: funcE(:)
      character (len=50) :: file_auto

      real (kind=8) :: a,b
      integer :: i,nio

      complex (kind=8), parameter :: EYE = (0,1)
      real (kind=8), parameter ::                                       &
       pi = 3.14159265358979323846264338327950288419716939937511d0


      namelist / param / Emin,Emax,conv,file_auto


      file_auto = 'file_auto'
      conv = 1.d0
      Emin = 0.d0
      Emax = -1.d0
      read(in_unitp,param)
      write(out_unitp,param)

!     read the time function (autocorrelation...)

       open(newunit=nio,file=file_auto)
       read(nio,*) npts

          allocate(auto(npts))
          allocate(t(npts))

          DO i=1,npts
            read(nio,*) t(i),a,b
            auto(i) = cmplx(a,b,kind=8)
          END DO
          dt = t(2)-t(1)
          tmax = t(npts)
          close(nio) ! CALL file_close cannot be used
          write(out_unitp,*) 'npts t0,tmax,dt: ',npts,t0,tmax,dt
!     END read the time function

!     Check the funcErgy grid...
          dE = 2.d0*Pi/tmax
          IF (Emax <0) Emax=2.d0*Pi/dt
          IF (Emax> 2.d0*Pi/dt) THEN
            write(out_unitp,*) 'Emax> 2.d0*Pi/dt',Emax,2.d0*Pi/tmax
            STOP
          END IF
          nptE = int((Emax-Emin)/dE)

          allocate(funcE(nptE))
          allocate(Expiwt(npts))
          DO i=1,nptE
            E=Emin+real(i-1,kind=8)*dE
            Expiwt(:) = exp(EYE*E*t(:))

            funcE(i) = 2.d0*dt*sum(Expiwt(:)*auto(:))
            write(out_unitp,*) E,funcE(i)
          END DO



          deallocate(funcE)
          deallocate(Expiwt)
          deallocate(t)
          deallocate(auto)
  END PROGRAM spectre

