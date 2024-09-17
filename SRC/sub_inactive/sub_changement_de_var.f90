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
!=============================================================
!
!      on calcule H*d0f
!      decompose en
!      une partie active et inactive
!
!      le vecteur i : c(.,i)
!
!=============================================================
      SUBROUTINE calc_Tinact_new(T0inact,Veff,T1,T2,                    &
                             d0c,d1c,d1xa,d2xaa,                        &
                             d0cd0c,d1xad0c,d1xad1xa,                   &
                             d1lnN,d2lnN,                               &
                             f2Qaa,f2Qii,f2Qai,                         &
                             f1Qa,f1Qi,                                 &
                             d0herm_ij,d1herm_ij,d2herm_ij,             &
                             nb_inact2,nb_act)
      USE EVR_system_m
      IMPLICIT NONE

      integer       :: nb_inact2,nb_act

      real (kind=Rkind) :: T0inact,Veff
      real (kind=Rkind) :: T1(nb_act),T2(nb_act,nb_act)

!------ f2Q and f1Q in programm order ------------------
      real (kind=Rkind) :: f1Qa(nb_act),f1Qi(nb_inact2)
      real (kind=Rkind) :: f2Qaa(nb_act,nb_act)
      real (kind=Rkind) :: f2Qii(nb_inact2,nb_inact2)
      real (kind=Rkind) :: f2Qai(nb_act,nb_inact2)


      real (kind=Rkind) :: d0herm_ij(nb_inact2)
      real (kind=Rkind) :: d1herm_ij(nb_inact2)
      real (kind=Rkind) :: d2herm_ij(nb_inact2)

      real (kind=Rkind) :: d1xa(nb_inact2,nb_act)
      real (kind=Rkind) :: d2xaa(nb_inact2,nb_act,nb_act)
      real (kind=Rkind) :: d0c(nb_inact2,nb_inact2)
      real (kind=Rkind) :: d1c(nb_inact2,nb_inact2,nb_act)
      real (kind=Rkind) :: d1lnN(nb_act),d2lnN(nb_act,nb_act)

      real (kind=Rkind) :: d0cd0c(nb_inact2,nb_inact2,nb_inact2)
      real (kind=Rkind) :: d1xad1xa(nb_inact2,nb_act,nb_act)
      real (kind=Rkind) :: d1xad0c(nb_inact2,nb_inact2,nb_act)


      real (kind=Rkind) :: d0f,d1f(nb_inact2),d2f(nb_inact2,nb_inact2)
      real (kind=Rkind) :: d1fQa(nb_act),d1fQi(nb_inact2)
      real (kind=Rkind) :: d2fQaa(nb_act,nb_act)
      real (kind=Rkind) :: d2fQii(nb_inact2,nb_inact2)
      real (kind=Rkind) :: d2fQai(nb_act,nb_inact2)

!     - local variables -------------
      integer       :: i,j,k,l
      real (kind=Rkind) :: A
!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*) 'BEGINNING calc_Tinact_new_old'
      write(out_unit,*) 'd1d2herm',d1herm_ij,d2herm_ij
      write(out_unit,*) 'f1Qi f2Qii ',f1Qi,f2Qii
      write(out_unit,*) 'f1Qa f2Qaa ',f1Qa,f2Qaa
      write(out_unit,*) 'f2Qai ',f2Qai
      write(out_unit,*) 'd1d2lnN ',d1lnN,d2lnN
      write(out_unit,*) 'nb_inact2',nb_inact2
      write(out_unit,*) 'd0c :',d0c
      write(out_unit,*) 'd1c :',d1c
      END IF
!---------------------------------------------------------------------


!-------------------------------------------------------
!-------------------------------------------------------
!     derivatives relatively to x
!     d0f,d1f(:),d2f(:,:)
      CALL d0d1d2f_harmo(d0f,d1f,d2f,                                   &
                         d0herm_ij,d1herm_ij,d2herm_ij,                 &
                         nb_inact2)
!-------------------------------------------------------
!-------------------------------------------------------


!-------------------------------------------------------
!-------------------------------------------------------
!     derivatives relatively to Qi inactives and actives
!     d0f
!     d1fQi(:),d1fQa(:),d2fQii(:,:)
!     d2fQii(:,:),d2fQaa(:,:),d2fQai(:,:)

!     -------------------------------------------------
!     - d./dQi inactif terms --------------------------
!     - Rq : d1x/dQi   = d0c(i,k)
!     - Rq : d1N/dQi   = ZERO
!     -------------------------------------------------
      DO i=1,nb_inact2
        d1fQi(i) = dot_product( d0c(i,:) , d1f(:) )
!       write(out_unit,*) 'd1psi_Qi',i,d1fQi(i)
      END DO

!     -------------------------------------------------
!     - d./dQi actif terms ----------------------------
!     - Rq : d1x/dQi   = d1xa(:,i)
!     - Rq : d1N/dQi   = d1lnN(i)
!     -------------------------------------------------
      DO i=1,nb_act
        d1fQa(i) = d1lnN(i) * d0f +                                     &
                   dot_product( d1xa(:,i) , d1f(:) )
!       write(out_unit,*) 'd1psi_Qa',i,d1fQa(i)
      END DO

!     -------------------------------------------------
!     - d2./dQiQj inactif inactif terms ---------------
!     - Rq : d1x/dQi     = d0c(i,k)
!     - Rq : d2x/dQidQj  = ZERO
!     - Rq : d1lnN/dQi   = ZERO
!     - Rq : d2lnN/dQiQj = ZERO
!     -------------------------------------------------
      DO i=1,nb_inact2
      DO j=i,nb_inact2
        d2fQii(i,j) = ZERO
        DO k=1,nb_inact2
        DO l=1,nb_inact2
          d2fQii(i,j) = d2fQii(i,j) + d0c(i,k) * d0c(j,l) * d2f(k,l)
        END DO
        END DO
        d2fQii(j,i) = d2fQii(i,j)
!       write(out_unit,*) 'd2psi_Qii',i,j,d2fQii(i,j)
      END DO
      END DO
!     -------------------------------------------------

!     -------------------------------------------------
!     - d2./dQiQj actif actif terms -------------------
!     - Rq : d1x/dQi     = d1xa(:,i)
!     - Rq : d2x/dQidQj  = d2xaa(:,i,j)
!     - Rq : d1lnN/dQi   = d1lnN(i)
!     - Rq : d2lnN/dQiQj = d2lnN(i,j)
!     -------------------------------------------------
      DO i=1,nb_act
      DO j=i,nb_act
        d2fQaa(i,j) = (d2lnN(i,j)+d1lnN(i)*d1lnN(j)) * d0f +            &
                      d1lnN(i)*dot_product( d1xa(:,j) , d1f(:) ) +      &
                      d1lnN(j)*dot_product( d1xa(:,i) , d1f(:) ) +      &
                      dot_product( d2xaa(:,i,j) , d1f(:) )
        DO k=1,nb_inact2
        DO l=1,nb_inact2
          d2fQaa(i,j) = d2fQaa(i,j) + d1xa(k,i) * d1xa(l,j) * d2f(k,l)
        END DO
        END DO
!       write(out_unit,*) 'd2psi_Qaa',i,j,d2fQaa(i,j)
        d2fQaa(j,i) = d2fQaa(i,j)
      END DO
      END DO
!     -------------------------------------------------

!     -------------------------------------------------
!     - d2./dQiQj inactif actif terms -----------------
!     - Rq : d1x/dQi     = d1xa(:,i)
!     - Rq : d1x/dQj     = d0c(j,k)
!     - Rq : d2x/dQidQj  = d1c(j,k,i)
!     - Rq : d1lnN/dQi   = d1lnN(i)
!     - Rq : d1lnN/dQj   = ZERO
!     - Rq : d2lnN/dQiQj = ZERO
!     -------------------------------------------------
      DO i=1,nb_act
      DO j=1,nb_inact2
        d2fQai(i,j) = d1lnN(i) * dot_product( d0c(j,:) , d1f(:) ) +     &
                      dot_product( d1c(j,:,i) , d1f(:) )
        DO k=1,nb_inact2
        DO l=1,nb_inact2
          d2fQai(i,j) = d2fQai(i,j) + d1xa(k,i) * d0c(j,l) * d2f(k,l)
        END DO
        END DO
!       write(out_unit,*) 'd2psi_Qai',i,j,d2fQai(i,j)
      END DO
      END DO
!     --------------------------------------------------

!-------------------------------------------------------
!-------------------------------------------------------



!--------------------------------------------------------
!--------------------------------------------------------

!     -------------------------------------------------
!     - d./dQi inactif terms --------------------------
!     - Veff and T0inact ------------------------------
!     -------------------------------------------------
      Veff = dot_product( f1Qi(:) , d1fQi(:) )
      T0inact = Veff
!     -------------------------------------------------

!     -------------------------------------------------
!     - d./dQi actif terms ----------------------------
!     - Veff and T1(:) kinetic energy part ------------
!     -------------------------------------------------
      DO i=1,nb_act
        Veff = Veff + f1Qa(i)*d1fQa(i)
        T1(i) = f1Qa(i) * d0f
      END DO
!     -------------------------------------------------

!     -------------------------------------------------
!     - d2./dQiQj inactif inactif terms ---------------
!     - Veff and T0inact ------------------------------
!     -------------------------------------------------
      A = ZERO
      DO i=1,nb_inact2
      DO j=i,nb_inact2
        A = A + f2Qii(i,j) * d2fQii(i,j)
      END DO
      END DO
      T0inact = T0inact + A
      Veff    = Veff    + A
!     -------------------------------------------------

!     -------------------------------------------------
!     - d2./dQiQj actif actif terms -------------------
!     - Veff, T2(:,:) and T1(:) kinetic energy part ---
!     -------------------------------------------------

      DO i=1,nb_act
      DO j=i,nb_act
        T2(i,j) = f2Qaa(i,j) * d0f
        T2(j,i) = f2Qaa(j,i) * d0f

        Veff  = Veff  + f2Qaa(i,j) * d2fQaa(i,j)
        T1(i) = T1(i) + f2Qaa(i,j) * d1fQa(j)
        T1(j) = T1(j) + f2Qaa(i,j) * d1fQa(i)
      END DO
      END DO
!     -------------------------------------------------

!     -------------------------------------------------
!     - d2./dQiQj inactif actif terms -----------------
!     - Veff and T1(:) kinetic energy part ------------
!     -------------------------------------------------
      DO i=1,nb_act
      DO j=1,nb_inact2
        Veff  = Veff  + f2Qai(i,j) * d2fQai(i,j)
        T1(i) = T1(i) + f2Qai(i,j) * d1fQi(j)
      END DO
      END DO
!     -------------------------------------------------

!     END IF

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*) 'T0inact ',T0inact
      write(out_unit,*) 'Veff ',Veff
      write(out_unit,*) 'T1 ',T1
      write(out_unit,*) 'T2 ',T2
      write(out_unit,*) 'END calc_Tinact'
      END IF
!---------------------------------------------------------------------
!     STOP
      RETURN
      end subroutine calc_Tinact_new
!
!=============================================================
!
!      Calculation of [(Tcor2 + Tcor1 + Trot) Psi_h]/Psi_h
!
!=============================================================
      SUBROUTINE calc_Tcorrot(Tcor2,Tcor1,Trot,                         &
                              d0c,d1c,d1xa,d2xaa,                       &
                              d0cd0c,d1xad0c,d1xad1xa,                  &
                              d1lnN,d2lnN,                              &
                              tcor2a,tcor2i,tcor1a,trota,               &
                              d0herm_ij,d1herm_ij,d2herm_ij,            &
                              gi,ga,                                    &
                              nb_inact2,nb_act                          &
                              )
      USE EVR_system_m
      IMPLICIT NONE

      integer       :: nb_inact2,nb_act

      real (kind=Rkind) :: Tcor2(nb_act,3),Tcor1(3),Trot(3,3)

!------ Tcor2, Tcor1 and Trot in programm order --------------
      real (kind=Rkind) :: tcor2a(nb_act,3),tcor2i(nb_inact2,3)
      real (kind=Rkind) :: tcor1a(3),trota(3,3)


      real (kind=Rkind) :: d0herm_ij(nb_inact2)
      real (kind=Rkind) :: d1herm_ij(nb_inact2)
      real (kind=Rkind) :: d2herm_ij(nb_inact2)

      real (kind=Rkind) :: gi(nb_inact2),ga(nb_act)
      real (kind=Rkind) :: Gii,Gaa,Gai
      real (kind=Rkind) :: d1xa(nb_inact2,nb_act)
      real (kind=Rkind) :: d2xaa(nb_inact2,nb_act,nb_act)
      real (kind=Rkind) :: d0c(nb_inact2,nb_inact2)
      real (kind=Rkind) :: d1c(nb_inact2,nb_inact2,nb_act)
      real (kind=Rkind) :: d1lnN(nb_act),d2lnN(nb_act,nb_act)

      integer       :: i,j,k

      real (kind=Rkind) :: d0cd0c(nb_inact2,nb_inact2,nb_inact2)
      real (kind=Rkind) :: d1xad1xa(nb_inact2,nb_act,nb_act)
      real (kind=Rkind) :: d1xad0c(nb_inact2,nb_inact2,nb_act)

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*) 'BEGINNING calc_Tcorrot'
      write(out_unit,*) 'd1d2herm',d1herm_ij,d2herm_ij
      write(out_unit,*) 'tcor2a,tcor2i',tcor2a,tcor2i
      write(out_unit,*) 'tcor1a',tcor1a
      write(out_unit,*) 'trota',trota
      write(out_unit,*) 'd1d2lnN ',d1lnN,d2lnN
      write(out_unit,*) 'nb_inact2',nb_inact2
      write(out_unit,*) 'd0c :',d0c
      write(out_unit,*) 'd1c :',d1c
      END IF
!---------------------------------------------------------------------


!     -------------------------------------------------
!     - initialisation --------------------------------
!     -------------------------------------------------
      Tcor1(:)   = tcor1a(:)
      Tcor2(:,:) = Tcor2a(:,:)

      Trot(:,:)  = trota(:,:)



!     -------------------------------------------------
!     - d./dQi inactif terms --------------------------
!     - Tcor1(:) part ---------------------------------
!     - Rq : d1x/dQi   = d0c(i,k)
!     - Rq : d1lnN/dQi = ZERO
!     -------------------------------------------------
      DO i=1,nb_inact2
        gi(i) = ZERO
        DO k=1,nb_inact2
          gi(i) = gi(i) + d0c(i,k) * d1herm_ij(k)
        END DO
        Tcor1(1) = Tcor1(1) + gi(i) * tcor2i(i,1)
        Tcor1(2) = Tcor1(2) + gi(i) * tcor2i(i,2)
        Tcor1(3) = Tcor1(3) + gi(i) * tcor2i(i,3)
      END DO
!     -------------------------------------------------


!     -------------------------------------------------
!     - d./dQi actif terms ----------------------------
!     - Tcor1 and Tcor2 (done in the initialisation) --
!     -------------------------------------------------
      DO i=1,nb_act
        ga(i) = d1lnN(i)
        DO k=1,nb_inact2
          ga(i) = ga(i) + d1xa(k,i) * d1herm_ij(k)
        END DO
        Tcor1(1) = Tcor1(1) + ga(i) * tcor2a(i,1)
        Tcor1(2) = Tcor1(2) + ga(i) * tcor2a(i,2)
        Tcor1(3) = Tcor1(3) + ga(i) * tcor2a(i,3)
      END DO
!     -------------------------------------------------
!     write(out_unit,*) 'd1f d1lnN',ga,d1lnN


!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*) 'Trot',Trot
      write(out_unit,*) 'Tcor1',Tcor1
      write(out_unit,*) 'Tcor2',Tcor2
      write(out_unit,*) 'END calc_Tcorrot'
      END IF
!---------------------------------------------------------------------
!     STOP
      RETURN
      end subroutine calc_Tcorrot
!
!=============================================================
!
!      d0x determination in function of the quadrature points
!      and their indexes (ind_quadra)
!
!=============================================================
      SUBROUTINE calc2_d0x(d0x,nb_inact2n,ind_quadra,tab_Pbasis)
      USE EVR_system_m
      USE mod_basis
      IMPLICIT NONE

       integer       :: nb_inact2n
       real (kind=Rkind) :: d0x(nb_inact2n)

       integer       :: ind_quadra(nb_inact2n)

!----- for the inactive basis sets ----------------------------------
      TYPE (P_basis)  :: tab_Pbasis(nb_inact2n)



       integer      :: k
!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING calc_d0x'
        write(out_unit,*) 'nb_inact2n',nb_inact2n
        write(out_unit,*) 'ind_quadra',ind_quadra

      END IF
!---------------------------------------------------------------------

       DO k=1,nb_inact2n
        d0x(k) = tab_Pbasis(k)%Pbasis%x(1,ind_quadra(k))
       END DO

!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'd0x',d0x
       write(out_unit,*) 'END calc_d0x'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc2_d0x
!
!=============================================================
!
!      calculs of Q and deltaQ in fuction of x
!
!=============================================================
      SUBROUTINE calc_d0xTOQ(Q,Qeq,deltaQ,x,c_inv,nb_inact2)
      USE EVR_system_m
      IMPLICIT NONE


       integer           :: nb_inact2
       real (kind=Rkind) :: Qeq(nb_inact2)
       real (kind=Rkind) :: Q(nb_inact2)
       real (kind=Rkind) :: deltaQ(nb_inact2)
       real (kind=Rkind) :: x(nb_inact2)
       real (kind=Rkind) :: c_inv(nb_inact2,nb_inact2)


       integer       :: i

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'BEGINNING calc_d0xTOQ'
       write(out_unit,*) 'Qeq',Qeq
       write(out_unit,*) 'x',x
       write(out_unit,*) 'c_inv'
       CALL Write_Mat_MPI(c_inv,out_unit,5)
      END IF
!---------------------------------------------------------------------

       DO i=1,nb_inact2
         deltaQ(i) = dot_product(c_inv(:,i),x(:))
         Q(i)      = Qeq(i) + deltaQ(i)
       END DO

!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'Q',Q
       write(out_unit,*) 'deltaQ',deltaQ
       write(out_unit,*) 'END calc_d0xTOQ'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc_d0xTOQ
!
!=============================================================
!
!      calculs of d1x and d2x in function deltaQ d0c
!
!=============================================================
      SUBROUTINE calc_d1d2x(                                            &
                            d1xa,d2xaa,nb_inact2,                       &
                            deltaQ,d0c,                                 &
                            d1Qeq,d1c,                                  &
                            d2Qeq,d2c,                                  &
                            nb_act)
      USE EVR_system_m
      IMPLICIT NONE

       integer       :: nb_act,nb_inact2
       real (kind=Rkind) :: deltaQ(nb_inact2)
       real (kind=Rkind) :: d1Qeq(nb_inact2,nb_act)
       real (kind=Rkind) :: d2Qeq(nb_inact2,nb_act,nb_act)
       real (kind=Rkind) :: d1xa(nb_inact2,nb_act)
       real (kind=Rkind) :: d2xaa(nb_inact2,nb_act,nb_act)
       real (kind=Rkind) :: d0c(nb_inact2,nb_inact2)
       real (kind=Rkind) :: d1c(nb_inact2,nb_inact2,nb_act)
       real (kind=Rkind) :: d2c(nb_inact2,nb_inact2,nb_act,nb_act)


       integer       :: i,j,k,l

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'BEGINNING calc_d1d2x'
       write(out_unit,*) 'deltaQ',deltaQ
       write(out_unit,*) 'd1Qeq',d1Qeq
       write(out_unit,*) 'd2Qeq',d2Qeq
       write(out_unit,*) 'd0c'
       CALL Write_Mat_MPI(d0c,out_unit,5)
       DO i=1,nb_inact2
         write(out_unit,*) 'd1c',i
         CALL Write_Mat_MPI(d1c(:,:,i),out_unit,5)
       END DO
       DO i=1,nb_inact2
       DO j=1,nb_inact2
       write(out_unit,*) 'd2c',i,j
       CALL Write_Mat_MPI(d2c(:,:,i,j),out_unit,5)
       END DO
       END DO
      END IF
!---------------------------------------------------------------------

       DO j=1,nb_act
       DO k=1,nb_inact2
         d1xa(k,j) = dot_product(d1c(:,k,j),deltaQ) -                   &
                     dot_product(d0c(:,k),d1Qeq(:,j))
       END DO
       END DO

       DO i=1,nb_act
       DO j=1,nb_act
       DO k=1,nb_inact2

         d2xaa(k,i,j) = dot_product(d2c(:,k,i,j),deltaQ(:)) -           &
                        dot_product(d1c(:,k,j)  ,d1Qeq(:,i)) -          &
                        dot_product(d1c(:,k,i)  ,d1Qeq(:,j)) -          &
                        dot_product(d0c(:,k)    ,d2Qeq(:,i,j) )

       END DO
       END DO
       END DO

!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'd1xa',d1xa
       write(out_unit,*) 'd2xaa',d2xaa
       write(out_unit,*) 'END calc_d1d2x'
      END IF
!---------------------------------------------------------------------


      RETURN
      end subroutine calc_d1d2x
!
!=====================================================================
!
! ++   the weight in nD
!
!=====================================================================
!
      SUBROUTINE weight2_nD(wnDh,tab_wnDh,i_modif,ind_quadra,           &
                            tab_Pbasis,nb_inact2n)
      USE EVR_system_m
      USE mod_basis
      IMPLICIT NONE

      integer       :: nb_inact2n
      real (kind=Rkind) :: wnDh,tab_wnDh(nb_inact2n)
      integer       :: i_modif
      integer       :: ind_quadra(nb_inact2n)


!----- for the inactive basis sets ----------------------------------
      TYPE (P_basis)  :: tab_Pbasis(nb_inact2n)

      integer       :: i

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'BEGINNING weight_nD'
       write(out_unit,*) 'nb_inact2n',nb_inact2n
       write(out_unit,*) 'i_modif',i_modif
       write(out_unit,*) 'ind_quadra',ind_quadra
      END IF
!---------------------------------------------------------------------

!     GOTO 100
      IF ( i_modif == 0 ) THEN

        tab_wnDh(1) = tab_Pbasis(1)%Pbasis%w(1)

        DO i=2,nb_inact2n-1
          tab_wnDh(i) = tab_wnDh(i-1) * tab_Pbasis(i)%Pbasis%w(1)
!         multi = multi + 1
        END DO

      ELSE IF ( i_modif == 1 ) THEN

        tab_wnDh(1) = tab_Pbasis(1)%Pbasis%w(ind_quadra(1))
        DO i=2,nb_inact2n
          tab_wnDh(i) = tab_wnDh(i-1) * tab_Pbasis(i)%Pbasis%w(ind_quadra(i))
!         multi = multi + 1
        END DO

      ELSE
        DO i=i_modif,nb_inact2n
          tab_wnDh(i) = tab_wnDh(i-1) * tab_Pbasis(i)%Pbasis%w(ind_quadra(i))
!         multi = multi + 1
        END DO

      END IF

      wnDh = tab_wnDh(nb_inact2n)


!100  CONTINUE
!     -----------------------------------------------------
!     the weight in nD (old one)
!     i=1
!     write(out_unit,*) i,tab_Pbasis(i)%Pbasis%w(ind_quadra(i))
!     wnDh = tab_Pbasis(1)%Pbasis%w(ind_quadra(1))
!     DO i=2,nb_inact2n
!       write(out_unit,*) i,tab_Pbasis(i)%Pbasis%w(ind_quadra(i))
!       wnDh = wnDh * tab_Pbasis(i)%Pbasis%w(ind_quadra(i))
!     END DO
!     -----------------------------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'wnDh',wnDh
       write(out_unit,*) 'END weight_nD'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine weight2_nD
!
!=====================================================================
!
! ++   d0f = Produit de hermite
!
!=====================================================================
!
      SUBROUTINE d0f_harmo(d0f,d0herm_ij,nb_inact2)
      USE EVR_system_m
      IMPLICIT NONE


      integer       :: nb_inact2
      real (kind=Rkind) :: d0herm_ij(nb_inact2)
      real (kind=Rkind) :: d0f

      integer       :: i

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'BEGINNING d0f_harmo'
       write(out_unit,*) 'd0f_herm_ij',d0herm_ij
      END IF
!---------------------------------------------------------------------

      d0f = product(d0herm_ij)

!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'END d0f_harmo'
       write(out_unit,*) 'd0f',d0f
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine d0f_harmo

!=====================================================================
!
! ++   d0f = Produit de hermite
!
!=====================================================================
!
      SUBROUTINE d0d1d2f_harmo(d0f,d1f,d2f,                             &
                               d0herm_ij,d1herm_ij,d2herm_ij,           &
                               nb_inact2)
      USE EVR_system_m
      IMPLICIT NONE


      integer           :: nb_inact2
      real (kind=Rkind) :: d0herm_ij(nb_inact2)
      real (kind=Rkind) :: d1herm_ij(nb_inact2)
      real (kind=Rkind) :: d2herm_ij(nb_inact2)
      real (kind=Rkind) :: d0f
      real (kind=Rkind) :: d1f(nb_inact2)
      real (kind=Rkind) :: d2f(nb_inact2,nb_inact2)

      integer :: i,j,k

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING d0d1d2f_harmo'
        write(out_unit,*) 'd0f_herm_ij',d0herm_ij
        write(out_unit,*) 'd1f_herm_ij',d1herm_ij
        write(out_unit,*) 'd2f_herm_ij',d2herm_ij
      END IF
!---------------------------------------------------------------------

      DO i=1,nb_inact2
      DO j=i,nb_inact2
        d2f(i,j) = ONE
        DO k=1,nb_inact2
          IF (k /= i .AND. k /= j)                                      &
                    d2f(i,j) = d2f(i,j) * d0herm_ij(k)
        END DO
        d1f(i) = d2f(i,i)
      END DO
      END DO
      d0f = d1f(1)*d0herm_ij(1)

      DO i=1,nb_inact2
        d1f(i)   = d1f(i)   * d1herm_ij(i)
        d2f(i,i) = d2f(i,i) * d2herm_ij(i)
        DO j=i+1,nb_inact2
          d2f(i,j) = d2f(i,j) * d1herm_ij(i)*d1herm_ij(j)
          d2f(j,i) = d2f(i,j)
        END DO
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
       write(out_unit,*) 'd0f',d0f
       write(out_unit,*) 'd1f',d1f
       write(out_unit,*) 'd2f',d2f
       write(out_unit,*) 'END d0d1d2f_harmo'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine d0d1d2f_harmo
!
!=====================================================================
!
! ++   calculation of d0herm_ij,d1herm_ij,d2herm_ij
!
!=====================================================================
!
      SUBROUTINE calc_d0herm_ij(ind_quadra,herm_ij,ind_herm_ij,         &
                                tab_Pbasis,nb_inact2)

      USE EVR_system_m
      USE mod_basis
      IMPLICIT NONE



      integer  ::  nb_inact2
      integer  ::  ind_quadra(nb_inact2),ind_herm_ij(nb_inact2)
      real (kind=Rkind) ::     herm_ij(nb_inact2)

!----- for the inactive basis sets ----------------------------------
      TYPE (P_basis)  :: tab_Pbasis(nb_inact2)




      integer  ::  i

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!     ----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING calc_d0herm_ij'
        write(out_unit,*) 'ind_quadra',ind_quadra
        write(out_unit,*) 'ind_herm_ij',ind_herm_ij
        write(out_unit,*) 'nb_inact2',nb_inact2
      END IF
!---------------------------------------------------------------------

      DO i=1,nb_inact2
        herm_ij(i) = tab_Pbasis(i)%Pbasis%dnRGB%d0(ind_quadra(i),ind_herm_ij(i)+1)
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*) 'herm_ij ',herm_ij
      write(out_unit,*) 'END calc_d0herm_ij'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc_d0herm_ij
      SUBROUTINE calc_d1herm_ij(ind_quadra,herm_ij,ind_herm_ij,         &
                                tab_Pbasis,nb_inact2)

      USE EVR_system_m
      USE mod_basis
      IMPLICIT NONE



      integer :: nb_inact2
      integer :: ind_quadra(nb_inact2),ind_herm_ij(nb_inact2)
      real(kind=Rkind) :: herm_ij(nb_inact2)

!----- for the inactive basis sets ----------------------------------
      TYPE (P_basis)  :: tab_Pbasis(nb_inact2)




      integer :: i

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!     ----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING calc_d1herm_ij'
        write(out_unit,*) 'ind_quadra',ind_quadra
        write(out_unit,*) 'ind_herm_ij',ind_herm_ij
        write(out_unit,*) 'nb_inact2',nb_inact2
      END IF
!---------------------------------------------------------------------

      DO i=1,nb_inact2
        herm_ij(i) = tab_Pbasis(i)%Pbasis%dnRGB%d1(ind_quadra(i),ind_herm_ij(i)+1,1)
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*) 'herm_ij ',herm_ij
      write(out_unit,*) 'END calc_d1herm_ij'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc_d1herm_ij
      SUBROUTINE calc_d2herm_ij(ind_quadra,herm_ij,ind_herm_ij,         &
                                tab_Pbasis,nb_inact2)

      USE EVR_system_m
      USE mod_basis
      IMPLICIT NONE



      integer       :: nb_inact2
      integer       :: ind_quadra(nb_inact2),ind_herm_ij(nb_inact2)
      real (kind=Rkind) :: herm_ij(nb_inact2)

!----- for the inactive basis sets ----------------------------------
      TYPE (P_basis)  :: tab_Pbasis(nb_inact2)




      integer    :: i

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!     ----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING calc_d2herm_ij'
        write(out_unit,*) 'ind_quadra',ind_quadra
        write(out_unit,*) 'ind_herm_ij',ind_herm_ij
        write(out_unit,*) 'nb_inact2',nb_inact2
      END IF
!---------------------------------------------------------------------

      DO i=1,nb_inact2
      herm_ij(i) = tab_Pbasis(i)%Pbasis%dnRGB%d2(ind_quadra(i),ind_herm_ij(i)+1,1,1)
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unit,*) 'herm_ij ',herm_ij
      write(out_unit,*) 'END calc_d2herm_ij'
      END IF
!---------------------------------------------------------------------

      end subroutine calc_d2herm_ij

!=====================================================================
!
! ++   calc d0cd0c
!
!
!=====================================================================
!
      SUBROUTINE calc_d0cd0c(d0cd0c,d0c,n)
      USE EVR_system_m
      IMPLICIT NONE

      integer           :: n
      real (kind=Rkind) :: d0c(n,n)
      real (kind=Rkind) :: d0cd0c(n,n,n)

      integer           :: i,j,k

!     -------------------------------------------------
!     - d0c * d0c
!     -------------------------------------------------
      DO i=1,n
      DO j=i,n
      DO k=1,n
        d0cd0c(k,j,i) = d0c(i,k) * d0c(j,k)
      END DO
      END DO
      END DO
!     -------------------------------------------------
!     -------------------------------------------------


      end subroutine calc_d0cd0c
