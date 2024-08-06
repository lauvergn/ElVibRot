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
 SUBROUTINE TF_autocorr(para_propa)
   USE mod_system
   USE mod_Constant, ONLY: get_Conv_au_TO_unit
   USE mod_propa
   IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      real (kind=Rkind)  :: microDeltaT


      complex (kind=Rkind),pointer :: C(:),W(:)
      real (kind=Rkind) :: E,DE,T,Tmax
      real (kind=Rkind) :: A,B,cmc


      integer :: no,ni

      integer :: NPT,NPT2,I,max_iE
      integer :: necr

      character (len=Name_len) :: name
      real (kind=Rkind)        :: ca,cb,cc,auTOcm_inv
      integer :: err_mem,memory

      CALL file_open(para_propa%file_autocorr,ni)
      IF(MPI_id==0) THEN
        CALL file_open(para_propa%file_spectrum,no)
        write(out_unitp,*) 'ni,no',ni,no
      ENDIF

      microDeltaT = para_propa%WPdeltaT/para_propa%nb_micro
      auTOcm_inv  = get_Conv_au_TO_unit('E','cm-1')

      NPT  = para_propa%WPTmax/microDeltaT
      NPT2 = 2**para_propa%TFnexp2
      IF (NPT2 < NPT) THEN
        write(out_unitp,*) ' ERROR in TF_autocorr'
        write(out_unitp,*) ' NPT2(=2^TFnexp2) is inferior to npt',NPT2,NPT
        write(out_unitp,*) ' Increase "TFnexp2" in the "&propa" namlist'
        STOP
      END IF

      nullify(C)
      nullify(W)
      CALL alloc_array(C,[NPT2],"C","TF_autocorr")
      CALL alloc_array(W,[NPT2],"C","TF_autocorr")

      write(out_unitp,*)  'NPT,NPT2,microWPdeltaT',NPT,NPT2,microDeltaT
      flush(out_unitp)

      DO I=1,NPT
        CALL Read_AutoCorr(ni,t,C(I))
        !write(out_unitp,*) t,C(I) ; flush(out_unitp)
      END DO
      C(1) = C(1) * HALF
      CALL file_close(para_propa%file_autocorr)
      C(NPT+1:NPT2) = cmplx(ZERO,ZERO,kind=Rkind)


      TMAX = NPT2*microDeltaT
      DE   = TWO*pi/TMAX * auTOcm_inv

      write(out_unitp,*) 'DeltaE(TF)  (cm-1)',DE
      write(out_unitp,*) 'TFmaxE  (cm-1)', para_propa%TFmaxE*auTOcm_inv
      write(out_unitp,*) 'TFminE  (cm-1)', para_propa%TFminE*auTOcm_inv

      CALL PREFFT(NPT2,1,para_propa%TFnexp2,W)
      CALL FFT(NPT2,1,microDeltaT,para_propa%TFnexp2,W,C)

      E = para_propa%TFminE*auTOcm_inv

      max_iE = min(NPT2,int(para_propa%TFmaxE * auTOcm_inv/DE))

      IF (max_iE == 0) max_iE = NPT2

      DO i=1,max_iE
        IF(MPI_id==0) write(no,*) E,real(C(i),kind=Rkind),aimag(C(i)),abs(C(i))
        E = E + DE
      END DO

      IF(MPI_id==0) CALL file_close(para_propa%file_autocorr)


  33  FORMAT(5(E13.6,' '))

      CALL dealloc_array(C,"C","TF_autocorr")
      CALL dealloc_array(W,"C","TF_autocorr")

  END SUBROUTINE TF_autocorr

