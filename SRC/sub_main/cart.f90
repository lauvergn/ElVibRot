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
      SUBROUTINE sub_cart(max_mem)
      USE EVR_system_m
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_cart
      IMPLICIT NONE

      TYPE (Type_cart) :: para_cart


!
!=====================================================================
!
!     variables
!
!=====================================================================
!

!----- variables for the dynamical memory allocation -----------------
      integer   max_mem



!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

      integer :: i,j,jat,ivA,ivB
      real(kind=Rkind) :: x,y,z
      character (len=Line_len) :: line

        character (len=Name_len) :: name_transfo
        integer :: nat,nb_vect,nb_G,iv_tot
        logical :: cos_th
        real(kind=Rkind) :: conv
        namelist /Coord_transfo/ name_transfo,nat,nb_G,nb_vect,cos_th,conv

      integer, parameter :: max_atG = 100
      integer :: tab_At_TO_G(max_atG)
      NAMELIST /recenterG / tab_At_TO_G
      real (kind=Rkind) :: Mtot


      logical :: Frame,cart
      character (len=Name_len) :: name_d,name_th,name_dih,              &
                                  name_x,name_y,name_z,                 &
                                  name_alpha,name_beta,name_gamma

      TYPE (Type_dnVec) :: UnitVect_F0(3) ! unit vectors (ndim=3)

      real (kind=Rkind) :: Ri,ui,sthi,thi,phi
      NAMELIST /vector/ nb_vect,Frame,                                  &
                        cos_th,name_d,name_th,name_dih,                 &
                        cart,name_x,name_y,name_z,                      &
                        name_alpha,name_beta,name_gamma
!=====================================================================
      write(out_unit,*) '================================================='
      write(out_unit,*) ' CART: initialization of variables'

      CALL sub_constantes(const_phys,.FALSE.)

      write(out_unit,*) ' END CART: initialization of variables'
      write(out_unit,*) '================================================='

      DO
        nat     = 0
        nb_vect = 0
        nb_G    = 0
        conv    = ONE
        read(in_unit,Coord_transfo,err=999,end=999)
        write(out_unit,Coord_transfo)
        IF (nb_vect < 1) nb_vect = nat-1

        write(out_unit,*) '================================================='
        write(out_unit,*) ' CART: ',trim(name_transfo)

        SELECT CASE (name_transfo)
        CASE ('cart')
          write(out_unit,*) nat
          CALL alloc_Type_cart(para_cart,nb_at=nat+nb_G)
          DO i=1,nat
            read(in_unit,41) line
            !write(out_unit,*) line
            flush(out_unit)
 41         format(72a)
            read(line,*) name_Z,x,y,z
            para_cart%masses(i) = get_mass_Tnum(const_phys%mendeleev,name=name_Z)
            para_cart%dnAt(i)%d0(1:3) = [x,y,z]
            para_cart%dnAt(i)%d0(:) = para_cart%dnAt(i)%d0(:) * conv
            !write(out_unit,*) para_cart%masses(i),para_cart%dnAt(i)%d0
          END DO
          para_cart%Mtot = sum(para_cart%masses)

          ! for the center-of-mass
          DO i=1,nb_G
            tab_At_TO_G(:) = 0
            read(in_unit,recenterG)
            write(out_unit,recenterG)
            para_cart%dnAt(nat+i)%d0(1:3) = ZERO
            Mtot = ZERO
            DO j=1,count(tab_At_TO_G > 0)
              jat = tab_At_TO_G(j)
              IF (jat < 1 .OR. jat > nat) THEN
                write(out_unit,*) ' ERROR in sub_cart'
                write(out_unit,*) ' The index of the center-of-mass is out of range: [1,',nat,']'
                write(out_unit,*) ' tab_At_TO_G(:) ',tab_At_TO_G(1:count(tab_At_TO_G > 0))
                STOP
              END IF
              Mtot = Mtot + para_cart%masses(jat)
              para_cart%dnAt(nat+i)%d0(1:3) = para_cart%dnAt(nat+i)%d0(1:3) + &
                                    para_cart%masses(jat)*para_cart%dnAt(jat)%d0(1:3)
            END DO
            para_cart%dnAt(nat+i)%d0(1:3) = para_cart%dnAt(nat+i)%d0(1:3) / Mtot
          END DO

          ! for the vectors
          CALL alloc_Type_cart(para_cart,nb_vect=nb_vect+1)  ! one vector is added for the center-of-mass
          DO i=1,nb_vect
            read(in_unit,*) ivA,ivB
            IF (ivA < 1 .OR. ivA > para_cart%nb_at .OR.                 &
                ivB < 1 .OR. ivB > para_cart%nb_at) THEN
                write(out_unit,*) ' ERROR in sub_cart'
                write(out_unit,*) ' The index of the vectors are out of range: [1,',nat,']'
                write(out_unit,*) ' ivA,ivB ',ivA,ivb
                STOP
            END IF
            para_cart%dnVect(i)%d0(1:3) = para_cart%dnAt(ivB)%d0(1:3) - &
                                          para_cart%dnAt(ivA)%d0(1:3)
          END DO
          ! the last vector: the COM
          i = nb_vect+1
          para_cart%dnVect(i)%d0(1:3) = ZERO
          DO j=1,nat
            para_cart%dnVect(i)%d0(1:3) = para_cart%dnVect(i)%d0(1:3) + &
                           para_cart%masses(j)*para_cart%dnAt(j)%d0(1:3)
          END DO
          para_cart%dnVect(i)%d0(1:3) = para_cart%dnVect(i)%d0(1:3) / para_cart%Mtot

          CALL write_Type_cart(para_cart)
        CASE ('vector')
          DO i=1,nb_vect
            write(out_unit,*) 'vect',i,para_cart%dnVect(i)%d0
          END DO

          CALL alloc_dnSVM(UnitVect_F0(1),nb_var_vec=3,nderiv=0)
          UnitVect_F0(1)%d0(1) = ONE

          CALL alloc_dnSVM(UnitVect_F0(2),nb_var_vec=3,nderiv=0)
          UnitVect_F0(2)%d0(3) = ONE

          CALL alloc_dnSVM(UnitVect_F0(3),nb_var_vec=3,nderiv=0)
          UnitVect_F0(3)%d0(3) = ONE

          iv_tot = 0
          CALL RecGet_Vec_Fi(para_cart%dnVect,nb_vect,iv_tot,1,UnitVect_F0)
          write(out_unit,*) 'iv_tot',iv_tot

          CALL dealloc_dnSVM(UnitVect_F0(1))
          CALL dealloc_dnSVM(UnitVect_F0(2))
          CALL dealloc_dnSVM(UnitVect_F0(3))
        CASE default ! ERROR: wrong transformation !
          write(out_unit,*) ' ERROR in sub_cart'
          write(out_unit,*) ' The transformation is UNKNOWN: ',trim(name_transfo)
          STOP
        END SELECT


      write(out_unit,*) ' CART: END ',trim(name_transfo)
      write(out_unit,*) '================================================='
      write(out_unit,*)
      flush(out_unit)
     END DO
 999 CONTINUE ! no more namelist (end-of-file)



      write(out_unit,*) '================================================'
      write(out_unit,*) ' CART AU REVOIR!!!'
      write(out_unit,*) '================================================'


      END SUBROUTINE sub_cart

      RECURSIVE SUBROUTINE RecGet_Vec_Fi(tab_Vect_Fi,nb_vect_tot,iv_tot,iv_Fi,UnitVect_Fi)
      USE EVR_system_m
      USE mod_dnSVM
      IMPLICIT NONE

      integer, intent(in)    :: nb_vect_tot,iv_Fi
      integer, intent(inout) :: iv_tot
      TYPE (Type_dnVec), intent(in) :: tab_Vect_Fi(nb_vect_tot) ! table of vectors (ndim=3)
      TYPE (Type_dnVec), intent(in) :: UnitVect_Fi(3) ! unit vectors (ndim=3)

      TYPE (Type_dnVec) :: UnitVect_Fij(3) ! unit vectors if frame=t (ndim=3)
      integer :: iv_Fij

      integer :: nb_vect
      logical :: Frame,cart,cos_th,zmat_order
      character (len=Name_len) :: name_d,name_th,name_dih,              &
                                  name_x,name_y,name_z,                 &
                                  name_alpha,name_beta,name_gamma

      real (kind=Rkind) :: Riv,px,py,pz
      real (kind=Rkind) :: alphaiv,betaiv,ubetaiv,gammaiv,sgamma,cgamma
      real (kind=Rkind) :: uiv,thiv,phiv

      NAMELIST /vector/ nb_vect,Frame,                                  &
                        cos_th,name_d,name_th,name_dih,                 &
                        cart,name_x,name_y,name_z,                      &
                        name_alpha,name_beta,name_gamma,                &
                        zmat_order
      real (kind=Rkind), parameter :: radTOdeg = 180._Rkind / pi

      iv_tot = iv_tot + 1
      !write(out_unit,*) 'RecGet_Vec_Fi: nb_vect_tot,iv_Fi',nb_vect_tot,iv_Fi,iv_tot
      !write(out_unit,*) 'vect:',tab_Vect_Fi(iv_tot)%d0
      nb_vect = 0
      Frame = .FALSE.
      zmat_order = .TRUE.
      read(in_unit,vector)
      !write(out_unit,vector)

      IF (Frame) THEN
        write(out_unit,*) '============================================='
        write(out_unit,*) ' Frame = T, iv_tot',iv_tot
        ! norm of the vector (distance)
        Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0))

        ! The 3 unit vectors in the new frame Fij

        CALL alloc_dnSVM(UnitVect_Fij(1),nb_var_vec=3,nderiv=0)
        CALL alloc_dnSVM(UnitVect_Fij(2),nb_var_vec=3,nderiv=0)
        CALL alloc_dnSVM(UnitVect_Fij(3),nb_var_vec=3,nderiv=0)

        UnitVect_Fij(3)%d0 = tab_Vect_Fi(iv_tot)%d0 / Riv ! ez_Fij
        IF (nb_vect > 0) THEN
          pz = dot_product(tab_Vect_Fi(iv_tot+1)%d0,UnitVect_Fij(3)%d0)
          UnitVect_Fij(1)%d0 = tab_Vect_Fi(iv_tot+1)%d0 - pz * UnitVect_Fij(3)%d0
          UnitVect_Fij(1)%d0 = UnitVect_Fij(1)%d0 /                       &
               sqrt(dot_product(UnitVect_Fij(1)%d0,UnitVect_Fij(1)%d0))

          UnitVect_Fij(2)%d0(1) =                                         &
                       UnitVect_Fij(3)%d0(2) * UnitVect_Fij(1)%d0(3) -  &
                       UnitVect_Fij(3)%d0(3) * UnitVect_Fij(1)%d0(2)
          UnitVect_Fij(2)%d0(2) =                                         &
                       UnitVect_Fij(3)%d0(3) * UnitVect_Fij(1)%d0(1) -  &
                       UnitVect_Fij(3)%d0(1) * UnitVect_Fij(1)%d0(3)
          UnitVect_Fij(2)%d0(3) =                                         &
                       UnitVect_Fij(3)%d0(1) * UnitVect_Fij(1)%d0(2) -  &
                       UnitVect_Fij(3)%d0(2) * UnitVect_Fij(1)%d0(1)
           !write(out_unit,*) 'ex_Fij',UnitVect_Fij(1)%d0
           !write(out_unit,*) 'ey_Fij',UnitVect_Fij(2)%d0
           !write(out_unit,*) 'ez_Fij',UnitVect_Fij(3)%d0
        END IF

        ! Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0)) !already calculated

        pz = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(3)%d0)
        px = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(1)%d0)
        py = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(2)%d0)

        ubetaiv = pz / Riv  ! cos(beta)
        betaiv = acos(ubetaiv)
        alphaiv = atan2(py,px)

        write(out_unit,*) 'R       : ',iv_Fi,':',Riv,Riv

        IF (zmat_order) THEN
          IF (iv_Fi == 2) THEN
            write(out_unit,*) 'beta (u): ',iv_Fi,':',betaiv*radTOdeg,ubetaiv
          ELSE ! iv_Fi /= 2
            write(out_unit,*) 'beta (u): ',iv_Fi,':',betaiv*radTOdeg,ubetaiv
            write(out_unit,*) 'alpha   : ',iv_Fi,':',alphaiv*radTOdeg,alphaiv
          END IF

        END IF

        ! loop on the vectors in frame Fij
        DO iv_Fij=2,nb_vect+1
           CALL RecGet_Vec_Fi(tab_Vect_Fi,nb_vect_tot,iv_tot,iv_Fij,UnitVect_Fij)

          IF (zmat_order .AND. iv_Fij==2) THEN
            ! for gamma we use ex_BF projected on the SF
            px = dot_product(UnitVect_Fij(1)%d0,UnitVect_Fi(3)%d0)
            py = dot_product(UnitVect_Fij(2)%d0,UnitVect_Fi(3)%d0)
            !write(out_unit,*) 'for gamma, px,py',px,py
            gammaiv = atan2(py,-px)
            write(out_unit,*) 'gamma   : ',iv_Fi,':',gammaiv*radTOdeg,gammaiv
          END IF
        END DO

        IF (.NOT. zmat_order) THEN
          IF (iv_Fi == 2) THEN
            write(out_unit,*) 'beta (u): ',iv_Fi,':',betaiv*radTOdeg,ubetaiv
          ELSE ! iv_Fi /= 2
            write(out_unit,*) 'alpha   : ',iv_Fi,':',alphaiv*radTOdeg,alphaiv
            write(out_unit,*) 'beta (u): ',iv_Fi,':',betaiv*radTOdeg,ubetaiv
          END IF

          IF (nb_vect > 0) THEN
            ! for gamma we use ex_BF projected on the SF
            px = dot_product(UnitVect_Fij(1)%d0,UnitVect_Fi(3)%d0)
            py = dot_product(UnitVect_Fij(2)%d0,UnitVect_Fi(3)%d0)
            !write(out_unit,*) 'for gamma, px,py',px,py
            gammaiv = atan2(py,-px)
            write(out_unit,*) 'gamma   : ',iv_Fi,':',gammaiv*radTOdeg,gammaiv
          END IF
        END IF

        CALL dealloc_dnSVM(UnitVect_Fij(1))
        CALL dealloc_dnSVM(UnitVect_Fij(2))
        CALL dealloc_dnSVM(UnitVect_Fij(3))

        write(out_unit,*) '============================================='
      ELSE

!        write(out_unit,*) 'ex_Fi',UnitVect_Fi(1)%d0
!        write(out_unit,*) 'ey_Fi',UnitVect_Fi(2)%d0
!        write(out_unit,*) 'ez_Fi',UnitVect_Fi(3)%d0

        IF (iv_Fi == 1) THEN
          write(out_unit,*) ' ERROR in RecGet_Vec_Fi'
          write(out_unit,*) ' iv_Fi=1 and frame=F is NOT possible'
          write(out_unit,*) ' check the fortran!'
          STOP
        ELSE IF (iv_Fi == 2) THEN
          ! norm of the vector (distance)
          Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0))

          pz = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(3)%d0)
          uiv = pz / Riv  ! cos(thi)
          thiv = acos(uiv)

          write(out_unit,*) 'R       : ',iv_Fi,':',Riv,Riv
          write(out_unit,*) 'th (u)  : ',iv_Fi,':',thiv*radTOdeg,uiv
        ELSE
          ! norm of the vector (distance)
          Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0))

          pz = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(3)%d0)
          uiv = pz / Riv  ! cos(thi)
          thiv = acos(uiv)

          px = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(1)%d0)
          py = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(2)%d0)
          phiv = atan2(py,px)

          write(out_unit,*) 'R       : ',iv_Fi,':',Riv,Riv
          write(out_unit,*) 'th (u)  : ',iv_Fi,':',thiv*radTOdeg,uiv
          write(out_unit,*) 'phi     : ',iv_Fi,':',phiv*radTOdeg,phiv

        END IF

      END IF

      END SUBROUTINE RecGet_Vec_Fi
