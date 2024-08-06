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
MODULE mod_basis_BtoG_GtoB_MPI
  USE mod_system
  USE mod_basis_BtoG_GtoB_SGType4
  USE mod_MPI_aux
  IMPLICIT NONE

  CONTAINS
!=======================================================================================
  SUBROUTINE CVecB_TO_CVecG_R_MPI(SRep,CVecG)
    IMPLICIT NONE

    TYPE(Type_SmolyakRep),            Intent(in) :: SRep
    Complex(kind=Rkind),           Intent(inout) :: CVecG(:)

    Complex(kind=Rkind),allocatable              :: Cvec_temp(:)
    Integer                                      :: Cvec_length(0:MPI_np-1)
    Integer                                      :: itabR
    Integer                                      :: iG
    Integer                                      :: nR
    Integer                                      :: d1,d2

#if(run_MPI)

    Cvec_length=0

    IF(MPI_scheme==1 .OR. MPI_scheme==3) THEN
      itabR=0
      IF(MPI_scheme==1) THEN
        d1=iGs_MPI(1,MPI_id)
        d2=iGs_MPI(2,MPI_id)
        DO iG=d1,d2
          nR=size(SRep%SmolyakRep(iG)%V)
          itabR=itabR+nR
        ENDDO
        Cvec_length(MPI_id)=itabR

      ELSEIF(MPI_scheme==3 .AND. MPI_nodes_p0) THEN
        itabR=0
        d1=iGs_MPI(1,MPI_id)
        d2=iGs_MPI(2,MPI_sub_id(2))
        DO iG=d1,d2
          nR=size(SRep%SmolyakRep(iG)%V)
          itabR=itabR+nR
        ENDDO
        Cvec_length(MPI_id)=itabR
      ENDIF

      CALL MPI_collect_info(Cvec_length,bcast=.True.,MS=MPI_scheme)

      ! IF(MPI_scheme==1 .OR. (MPI_scheme==3 .AND. MPI_nodes_p0)) THEN
      IF(keep_MPI) THEN
        CALL allocate_array(Cvec_temp,1,Cvec_length(MPI_id))
        itabR=0
        DO iG=d1,d2
          nR=size(SRep%SmolyakRep(iG)%V)
          Cvec_temp(itabR+1:itabR+nR)=cmplx(SRep%SmolyakRep(iG)%V,kind=Rkind)
          itabR=itabR+nR
        ENDDO

        CALL MPI_combine_array(CVecG,Cvec_temp,lengths=Cvec_length,MS=MPI_scheme)
        IF(allocated(Cvec_temp)) deallocate(Cvec_temp)
      ENDIF

      d1=Sum(Cvec_length(0:MPI_np-1))
      IF(MPI_scheme==1) THEN
        CALL MPI_Bcast_(CVecG,d1,root_MPI)
      ELSEIF(MPI_scheme==3) THEN
        CALL MPI_Bcast(CVecG,d1,Cplx_MPI,root_MPI,MPI_NODE_0_COMM,MPI_err)
      ENDIF
    ENDIF

#endif
  ENDSUBROUTINE CVecB_TO_CVecG_R_MPI
!=======================================================================================

!=======================================================================================
  SUBROUTINE CVecB_TO_CVecG_C_MPI(SRep,CVecG)
    IMPLICIT NONE

    TYPE(Type_SmolyakRep),            Intent(in) :: SRep
    Complex(kind=Rkind),           Intent(inout) :: CVecG(:)

    Complex(kind=Rkind),allocatable              :: Cvec_temp1(:)
    Complex(kind=Rkind),allocatable              :: Cvec_temp2(:)

    Integer                                      :: Cvec_length(0:MPI_np-1)
    Integer                                      :: itabR
    Integer                                      :: iG
    Integer                                      :: nR
    Integer                                      :: d1,d2

#if(run_MPI)

    Cvec_length=0

    IF(MPI_scheme==1 .OR. MPI_scheme==3) THEN
      IF(MPI_scheme==1) THEN
        itabR=0
        d1=iGs_MPI(1,MPI_id)
        d2=iGs_MPI(2,MPI_id)

        DO iG=d1,d2
          nR=size(SRep%SmolyakRep(iG)%V)
          itabR=itabR+nR
        ENDDO
        Cvec_length(MPI_id)=itabR

      ELSEIF(MPI_scheme==3 .AND. MPI_nodes_p0) THEN
        itabR=0
        d1=iGs_MPI(1,MPI_id)
        d2=iGs_MPI(2,MPI_sub_id(2))

        DO iG=d1,d2
          nR=size(SRep%SmolyakRep(iG)%V)
          itabR=itabR+nR
        ENDDO
        Cvec_length(MPI_id)=itabR

      ENDIF

      CALL MPI_collect_info(Cvec_length,bcast=.True.,MS=MPI_scheme)

      ! IF(MPI_scheme==1 .OR. (MPI_scheme==3 .AND. MPI_nodes_p0)) THEN
      IF(keep_MPI) THEN
        CALL allocate_array(Cvec_temp1,1,Cvec_length(MPI_id))
        itabR=0
        DO iG=d1,d2
          nR=size(SRep%SmolyakRep(iG)%V)
          Cvec_temp1(itabR+1:itabR+nR)=EYE*cmplx(SRep%SmolyakRep(iG)%V,kind=Rkind)
          itabR=itabR+nR
        ENDDO
      ENDIF

      d1=Sum(Cvec_length(0:MPI_np-1))
      IF(keep_MPI) THEN
        CALL allocate_array(Cvec_temp2,1,d1)
        CALL MPI_combine_array(Cvec_temp2,Cvec_temp1,lengths=Cvec_length,MS=MPI_scheme)
        IF(allocated(Cvec_temp1)) deallocate(Cvec_temp1)
      ENDIF

      IF(MPI_scheme==1) THEN
        CALL MPI_Bcast_(Cvec_temp2,d1,root_MPI)
      ELSEIF(MPI_scheme==3 .AND. MPI_nodes_p0) THEN
        CALL MPI_Bcast(Cvec_temp2,d1,Cplx_MPI,root_MPI,MPI_NODE_0_COMM,MPI_err)
      ENDIF

      IF(keep_MPI) THEN
        CVecG(1:d1)=CVecG(1:d1)+Cvec_temp2(1:d1)
        IF(allocated(Cvec_temp2)) deallocate(Cvec_temp2)
      ENDIF
    ENDIF

#endif
  ENDSUBROUTINE CVecB_TO_CVecG_C_MPI
!=======================================================================================

END MODULE mod_basis_BtoG_GtoB_MPI
