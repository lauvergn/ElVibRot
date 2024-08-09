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
!
!      Determination of the nD-harmonic functions: nb_bi
!      and the corresponding table: ind_harm(nb_var+1,nb_bi)
!      as a function of max_excit (the degree maximal of the hermite polynomia)
!      and the harmonic frequencies (at the Qdyn0).
!
!      The nD-functions are sorted (or not) relatively to the harmonic energy.
!
!=====================================================================
      SUBROUTINE sub2_ind_harm(Basis2n,PrimOp,para_Tnum,mole)
      use EVR_system_m
      USE mod_nDindex
      use mod_Constant,  only: get_conv_au_to_unit
      USE mod_Coord_KEO, only: CoordType, Tnum, gaussian_width, get_Qact0
      use mod_PrimOp,    only: PrimOp_t, sub_dnfreq
      USE mod_basis
      IMPLICIT NONE

!----- for the Basis2n -------------------------------------------------
      TYPE (basis)   :: Basis2n

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

!----- for the primitive operators: PES --------------------------------
      TYPE (PrimOp_t)   :: PrimOp

!------ working variables ----------------------------------------------
      TYPE (Type_RPHpara_AT_Qact1) :: RPHpara_AT_Qact1
      real (kind=Rkind)            :: pot0_corgrad
      real (kind=Rkind)            :: Qact(mole%nb_var)

      real (kind=Rkind)            :: A(mole%nb_inact2n,mole%nb_inact2n)
      integer                      :: i,k,n_h, min_i(mole%nb_inact2n)

      real (kind=Rkind)            :: ene_freq,ZPE,ene,auTOcm_inv

!---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='sub2_ind_harm'
      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) ' BEGINNING ',name_sub
        write(out_unit,*) '  Qdyn0',mole%ActiveTransfo%Qdyn0
        write(out_unit,*) '  nb_inact2n',mole%nb_inact2n
        write(out_unit,*) '  nb_var',mole%nb_var
      END IF
!---------------------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      IF (mole%nb_inact2n == 0 .AND. PrimOp%nb_elec == 1) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_inact2n (nb_inact21+nb_inact22)= 0'
        write(out_unit,*) ' no harmonic inactive variables !!'
        STOP
      END IF
!------------------------------------------------------------

      CALL get_Qact0(Qact,mole%ActiveTransfo)

      CALL alloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1,mole%nb_act1,        &
                                               mole%nb_inact2n,nderiv=0)
      RPHpara_AT_Qact1%RPHQact1(:) = Qact(1:mole%nb_act1)
      ! frequencies at RPHpara_AT_Qact1%Qact1

      CALL sub_dnfreq(RPHpara_AT_Qact1,pot0_corgrad,                    &
                      para_Tnum,mole,mole%RPHTransfo_inact2n,nderiv=0,  &
                      test=.FALSE.,cHAC=.TRUE.)

      write(out_unit,*) '------------------------------------------'
      write(out_unit,*) '------------------------------------------'
      write(out_unit,*) ' Parameters for the Basis2n (HADA or cHAC)'
      write(out_unit,*) 'freq',RPHpara_AT_Qact1%dneHess%d0(:)*auTOcm_inv

      write(out_unit,*) 'd0Qeq',RPHpara_AT_Qact1%dnQopt%d0
      write(out_unit,*) 'd0c'
      CALL Write_Mat(RPHpara_AT_Qact1%dnC%d0,out_unit,5)
      write(out_unit,*)
!-----------------------------------------------------------------
      n_h = Basis2n%nDindB%Max_nDI

      CALL gaussian_width(mole%nb_inact2n,A,RPHpara_AT_Qact1%dnC%d0)

      IF (.NOT. Basis2n%nDindB%With_L) THEN
        Basis2n%nDindB%nDweight = RPHpara_AT_Qact1%dneHess%d0 ! change the weight with the frequnecies
      END IF
      n_h = Basis2n%nDindB%Max_nDI
      CALL init_nDindexPrim(Basis2n%nDindB,mole%nb_inact2n,min_i,With_init=.FALSE.)
      CALL sort_nDindex(Basis2n%nDindB)
      IF (n_h > 0) Basis2n%nDindB%Max_nDI = min(n_h,Basis2n%nDindB%Max_nDI)
      !CALL write_nDindex(Basis2n%nDindB)
!-----------------------------------------------------------------


!     -----------------------------------------------------
!     write label of the anharmonic basis functions
      CALL alloc_NParray(Basis2n%EneH0,[Basis2n%nDindB%max_nDI],        &
                        'Basis2n%nDindB%max_nDI',name_sub)
      ZPE = HALF*sum(RPHpara_AT_Qact1%dneHess%d0)
      DO i=1,Basis2n%nDindB%max_nDI
        Ene = ZPE + sum(real(Basis2n%nDindB%Tab_nDval(:,i),kind=Rkind) *&
                                          RPHpara_AT_Qact1%dneHess%d0(:))
        Basis2n%EneH0(i) = Ene ! this enables the NewVec_type=4 in Davidson
        write(out_unit,'(i4,2(1x,f16.4))',advance='no') i,             &
                                     Ene*auTOcm_inv,(Ene-ZPE)*auTOcm_inv
        DO k=1,mole%nb_inact2n
          write(out_unit,'(i4)',advance='no') Basis2n%nDindB%Tab_nDval(k,i)
        END DO
        write(out_unit,*)
      END DO
!     -----------------------------------------------------
      write(out_unit,*) ' number of nD inactive harmonic functions, nb_bi: ', &
                    Basis2n%nDindB%max_nDI
      IF (Basis2n%nDindB%max_nDI < 1) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_bi < 1 !!',Basis2n%nDindB%max_nDI
        STOP
      END IF

      write(out_unit,*) '------------------------------------------'
      write(out_unit,*) '------------------------------------------'

      CALL dealloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1)


!---------------------------------------------------------------------
      IF (debug) THEN
        CALL write_nDindex(Basis2n%nDindB,'Basis2n%nDindB')
        write(out_unit,*) ' END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE sub2_ind_harm
