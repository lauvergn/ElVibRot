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
MODULE mod_fullanalysis
  IMPLICIT NONE
CONTAINS
!================================================================
!
!     write the energy ene(i) and the vectors i psi(.,i)
!
!================================================================
      SUBROUTINE sub_analyse(Tab_Psi,nb_psi_in,para_H,para_ana,         &
                             para_intensity,para_AllOp,const_phys)
      USE EVR_system_m
      USE mod_Constant
      USE mod_Coord_KEO

      USE mod_basis
      USE mod_param_RD

      USE mod_psi,    ONLY : param_psi,Write_Psi_nDBasis,               &
                             sub_PsiGridRep_TO_BasisRep,sub_analyze_psi,&
                             Write_header_saveFile_psi

      USE mod_Op
      USE mod_analysis
      USE mod_ana_psi
      IMPLICIT NONE

!----- variables for the construction of H ----------------------------
      TYPE (param_AllOp)  :: para_AllOp
      TYPE (param_Op)     :: para_H

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)                  :: para_ana
      TYPE (param_intensity)            :: para_intensity
      TYPE (param_ana_psi), allocatable :: ana_psi(:)


!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),pointer     :: mole
      TYPE (Tnum),pointer          :: para_Tnum

!----- physical and mathematical constants ---------------------------
      TYPE (constant)            :: const_phys

!----- for the basis set ----------------------------------------------
      TYPE (Basis), pointer      :: BasisnD

!----- variables for the WP propagation ----------------------------
      integer            :: nb_psi_in
      TYPE (param_psi)   :: Tab_Psi(nb_psi_in)

!------ working variables ---------------------------------

      real (kind=Rkind), allocatable :: ene(:)
      real (kind=Rkind), allocatable :: Mat_psi(:,:)
      character(len=:),  allocatable :: info
      TYPE (string_t),   allocatable :: tab_PsiAna(:)

      logical                        :: cube = .FALSE.

      integer                        :: i,nb_col,ib,maxth,iana
      real (kind=Rkind)              :: Q,E,DE
      TYPE (File_t)                  :: file_WPspectral
      integer                        :: nioWP
      character (len=Name_longlen)   :: lformat
      TYPE(REAL_WU)                  :: RWU_ZPE,RWU_E,RWU_DE

      real (kind=Rkind), allocatable :: AllPsi_max_RedDensity(:)

      integer  :: Version_File,nb_psi,nb_tot
      namelist / headerFile / Version_File,nb_psi,nb_tot


!----- FUNCTION --------------------------------------------------
      real (kind=Rkind) :: part_func
!---- FUNCTION ---------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_analyse'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (para_ana%ana_level == 0) RETURN

      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum
      BasisnD    => para_H%para_AllBasis%BasisnD

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_ba,nb_qa',para_H%nb_ba,para_H%nb_qa
        write(out_unit,*) 'nb_bi',para_H%nb_bi
        write(out_unit,*) 'nb_bie',para_H%nb_bie
        write(out_unit,*) 'nb_act1',mole%nb_act1
        write(out_unit,*) 'max_ana,max_ene',para_ana%max_ana,para_ana%max_ene
        write(out_unit,*) 'nb_psi_in',nb_psi_in
        flush(out_unit)
      END IF
!-----------------------------------------------------------
      CALL alloc_NParray(ene,shape(Tab_psi),'ene',name_sub)

      ene(:) = real(Tab_psi(:)%CAvOp,kind=Rkind)
      CALL Set_ZPE_OF_Op(para_H,Ene=ene)

      IF (count(ene(:)-para_H%ZPE <= para_ana%max_ene) == 0) RETURN

      RWU_ZPE = REAL_WU(para_H%ZPE,'au','E')
      RWU_E   = REAL_WU(sum(ene) / real(nb_psi_in,kind=Rkind),'au','E')
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*)
      write(out_unit,*) 'ZPE        : ',RWU_Write(RWU_ZPE,WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
      write(out_unit,*) 'Average_ene: ',RWU_Write(RWU_E,  WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
      RWU_E   = REAL_WU(sum(ene),'au','E') ! trace
      write(out_unit,*) 'trace_ene  : ',RWU_Write(RWU_E,  WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
      flush(out_unit)
      IF (para_ana%max_ana > nb_psi_in) para_ana%max_ana = nb_psi_in

      IF (para_ana%intensity .AND. para_intensity%l_IntVR) THEN
        CALL alloc_NParray(para_intensity%ABC,[3,nb_psi_in],            &
                          "para_intensity%ABC",name_sub)
      END IF

       write(out_unit,*)
       Q =  part_func(ene,nb_psi_in,para_ana%Temp)

      file_WPspectral%name = make_EVRTFileName(para_ana%name_file_spectralWP)
      IF(MPI_id==0) CALL file_open(file_WPspectral,nioWP,lformatted=para_ana%formatted_file_WP)

      ! For the header of the file
      nb_psi        = count((ene(:)-para_H%ZPE) <= para_ana%max_ene)
      !write(out_unit,*) 'nb_psi',nb_psi
      !write(out_unit,*) 'ZPE,max_ene',para_H%ZPE,para_ana%max_ene


      IF(MPI_id==0) CALL Write_header_saveFile_psi(tab_Psi,nb_psi,file_WPspectral)

      ! write the energy level + save the psi
      allocate(ana_psi(nb_psi_in))
      allocate(tab_PsiAna(nb_psi_in))

      write(out_unit,*) 'population at T, Q',para_ana%Temp,Q
      write(out_unit,*) 'Energy level (',const_phys%ene_unit,') pop and means :'
      flush(out_unit)
      DO i=1,nb_psi_in
        ana_psi(i) = para_ana%ana_psi

        ana_psi(i)%Ene     = ene(i)
        ana_psi(i)%num_psi = i

        ana_psi(i)%ZPE        = para_H%ZPE
        ana_psi(i)%Part_Func  = Q
        ana_psi(i)%Temp       = para_ana%Temp

        IF (ene(i)-para_H%ZPE > para_ana%max_ene) CYCLE

        RWU_E  = REAL_WU(ene(i),'au','E')
        RWU_DE = REAL_WU(ene(i)-para_H%ZPE,'au','E')
        E  = convRWU_TO_R_WITH_WritingUnit(RWU_E)
        DE = convRWU_TO_R_WITH_WritingUnit(RWU_DE)


        IF (i < 10000) THEN
          lformat = '("lev0: ",i4,i4,l3,3(1x,' // trim(adjustl(EneIO_format)) // '))'
        ELSE
          lformat = '("lev0: ",i0,i0,l3,3(1x,' // trim(adjustl(EneIO_format)) // '))'
        END IF

        write(out_unit,lformat) i,0,tab_Psi(i)%convAvOp,E,DE

        IF(MPI_id==0) CALL Write_Psi_nDBasis(tab_Psi(i),nioWP,i,ZERO,file_WPspectral%formatted,FilePsiVersion)

      END DO
      IF(MPI_id==0) CALL file_close(file_WPspectral)


      IF (.NOT. allocated(AllPsi_max_RedDensity)) THEN
        CALL alloc_NParray(AllPsi_max_RedDensity,[BasisnD%nDindB%ndim],"AllPsi_max_RedDensity",name_sub)
        AllPsi_max_RedDensity(:) = ZERO
      END IF

      maxth              = 1
      !$ maxth           = omp_get_max_threads()

      IF (.NOT. openmp) maxth = 1

      para_ana%ana_psi%ZPE        = para_H%ZPE
      para_ana%ana_psi%Part_Func  = Q
      para_ana%ana_psi%Temp       = para_ana%Temp
      iana                        = 0
      IF (Ana_maxth > 1)  THEN
        write(out_unit,'(a)')              'Psi(:) analysis: (%): [--0-10-20-30-40-50-60-70-80-90-100]'
        write(out_unit,'(a)',ADVANCE='no') 'Psi(:) analysis: (%): ['
        flush(out_unit)
      END IF
!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(nb_psi_in,ene,para_H,para_AllOp,para_ana,tab_Psi,ana_psi,AllPsi_max_RedDensity,const_phys) &
!$OMP   SHARED(out_unit,para_intensity,tab_PsiAna,iana,MPI_id,Ana_maxth) &
!$OMP   PRIVATE(i,info) &
!$OMP   NUM_THREADS(Ana_maxth)

!$OMP   DO SCHEDULE(STATIC)
      DO i=1,nb_psi_in
        !OMP ATOMIC
        iana = iana + 1

        IF (ana_psi(i)%Ene-ana_psi(i)%ZPE > para_ana%max_ene) CYCLE

        IF (.NOT. tab_Psi(i)%BasisRep) THEN
          CALL sub_PsiGridRep_TO_BasisRep(tab_Psi(i))
        END IF

        CALL SET_string(info," ",TO_string(ene(i)*const_phys%auTOenergy,"f12.6")," : ")

        CALL sub_analyze_psi(tab_Psi(i),ana_psi(i),adia=.FALSE.,PsiAna=tab_PsiAna(i)%str)
        IF (Ana_maxth == 1) THEN
          write(out_unit,'(a)') tab_PsiAna(i)%str
          flush(out_unit)
        END IF

        !$OMP CRITICAL (sub_analyse_CRIT)
        IF (allocated(para_ana%ana_psi%max_RedDensity)) THEN
          DO ib=1,size(AllPsi_max_RedDensity)
            AllPsi_max_RedDensity(ib) = max(AllPsi_max_RedDensity(ib),para_ana%ana_psi%max_RedDensity(ib))
          END DO
        END IF
        !$OMP END CRITICAL (sub_analyse_CRIT)

        IF (para_ana%intensity .AND. para_intensity%l_IntVR) THEN
          CALL sub_moyABC(tab_Psi(i),i,info,para_intensity%ABC(:,i),para_H,PsiAna=tab_PsiAna(i)%str)
        ELSE IF (para_AllOp%tab_Op(1)%para_ReadOp%nb_scalar_Op > 0 .AND. &
                 para_ana%ana_psi%AvScalOp) THEN
          CALL sub_moyScalOp(tab_Psi(i),i,info,para_AllOp%tab_Op,PsiAna=tab_PsiAna(i)%str)
        END IF

        IF (para_ana%ana_psi%AvHiterm) THEN
          CALL sub_psiHitermPsi(tab_Psi(i),i,info,para_H,PsiAna=tab_PsiAna(i)%str)
        END IF

        deallocate(info)

        IF (mod(iana,nb_psi_in) == 0 .AND. MPI_id == 0 .AND. Ana_maxth > 1) write(out_unit,'(a)',ADVANCE='no') '---'
      END DO
!$OMP   END DO
!$OMP   END PARALLEL
      IF (Ana_maxth > 1) write(out_unit,'(a)',ADVANCE='yes') '----]'

      IF (Ana_maxth > 1) THEN
        write(out_unit,*) '=============================================================='
        write(out_unit,*) '=============================================================='
        write(out_unit,*) '=============================================================='
        write(out_unit,*) 'population at T, Q',para_ana%Temp,Q
        write(out_unit,*) 'Energy level (',const_phys%ene_unit,') pop and averages :'
        flush(out_unit)
        DO i=1,nb_psi_in
          write(out_unit,'(a)') tab_PsiAna(i)%str
          flush(out_unit)
        END DO
        write(out_unit,*) '=============================================================='
        write(out_unit,*) '=============================================================='
        write(out_unit,*) '=============================================================='
        flush(out_unit)
      END IF


      IF (allocated(AllPsi_max_RedDensity)) THEN
        CALL Write_Vec(AllPsi_max_RedDensity,out_unit,6,Rformat='e10.3',info='For all psi max_RedDensity ')
        !write(out_unit,*) 'For all psi max_RedDensity ',AllPsi_max_RedDensity(:)
        CALL dealloc_NParray(AllPsi_max_RedDensity,"AllPsi_max_RedDensity",name_sub)
      END IF

      flush(out_unit)

!----------------------------------------------------------


!----------------------------------------------------------
!     writing the eigenvectors
      para_ana%print_psi = min(para_ana%print_psi,nb_psi_in)
      IF (debug) para_ana%print_psi = nb_psi_in
      write(out_unit,*) 'para_ana%print_psi',para_ana%print_psi
      flush(out_unit)

      IF (cube) CALL write_cube(Tab_Psi)


      IF (para_ana%print_psi > 0 .OR. debug) THEN


        CALL alloc_NParray(Mat_psi,[tab_Psi(1)%nb_tot,para_ana%print_psi], &
                        "Mat_psi",name_sub)
        DO i=1,para_ana%print_psi
          Mat_psi(:,i) = tab_Psi(i)%RvecB(:)
        END DO
        IF (.NOT. para_H%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC    &
                                            .AND. mole%nb_act1 < 3) THEN
          CALL write_psi(Mat_psi,para_ana%psi2,para_ana%print_psi,      &
                          tab_Psi(1)%nb_tot,                            &
                          para_H%nb_ba,para_H%nb_qa,para_H%nb_bie,      &
                          para_Tnum,mole,BasisnD,para_H)
        !ELSE IF (para_H%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC    &
        !                                   .AND. mole%nb_act1 < 3) THEN
        !  CALL write_psi2_new(tab_Psi)
        END IF
        flush(out_unit)
        nb_col = 5
        write(out_unit,*) 'eigenvectors in column'
        write(out_unit,*) nb_col,para_ana%print_psi,tab_Psi(1)%nb_tot
        CALL Write_Mat(Mat_psi(:,1:para_ana%print_psi),out_unit,nb_col)
        CALL dealloc_NParray(Mat_psi,"Mat_psi",name_sub)

      END IF
!----------------------------------------------------------

      ! write basis
      !if(para_H%para_ReadOp%para_FileGrid%Type_FileGrid==4) CALL write_basis_biqi(tab_Psi(1))

      IF (allocated(ene))     CALL dealloc_NParray(ene,'ene',name_sub)
      IF (allocated(Mat_psi)) CALL dealloc_NParray(Mat_psi,'Mat_psi',name_sub)
      IF (allocated(AllPsi_max_RedDensity))                              &
          CALL dealloc_NParray(AllPsi_max_RedDensity,'AllPsi_max_RedDensity',name_sub)
      IF (allocated(info)) deallocate(info)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unit,*) 'END ',name_sub
!----------------------------------------------------------


  END SUBROUTINE sub_analyse

!================================================================
!
!     write psi or psi^2 on the grid point
!
!================================================================
      SUBROUTINE write_psi(Mat_psi,psi2,nb_psi,nb_tot,                  &
                           nb_ba,nb_qa,n_h,                             &
                           para_Tnum,mole,                              &
                           BasisnD,para_Op)

      USE EVR_system_m
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- variables for the construction of H ----------------------------
      TYPE (param_Op)   :: para_Op


      integer           :: nb_ba,nb_qa,n_h
      integer           :: nb_psi,nb_tot
      real (kind=Rkind) :: Mat_psi(nb_tot,nb_psi)
      logical           :: psi2

!----- for the basis set ----------------------------------------------
      TYPE (Basis) :: BasisnD


!------ working variables ---------------------------------
      integer       :: i,k,l,i_q
      integer       :: ih,ihk
      real (kind=Rkind), allocatable :: Qact1(:)
      real (kind=Rkind), allocatable :: psi_q(:,:)
      integer       :: ih_print
      integer, parameter :: max_print = 200

      real (kind=Rkind), allocatable :: d0b(:)
      real (kind=Rkind), allocatable :: psid0b_k(:)

      real (kind=Rkind)              :: WrhonD
      real (kind=Rkind), save        :: Q1


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING write_psi'
        write(out_unit,*) 'nb_psi,psi2',nb_psi,psi2
        write(out_unit,*) 'nb_ba',nb_ba
        write(out_unit,*) 'nb_qa',nb_qa
        write(out_unit,*) 'n_h',n_h
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
      write(out_unit,*)
      write(out_unit,*) 'eigenvectors on a grid',nb_psi

!     - initisalisation ----------------------------------
      CALL alloc_NParray(d0b,     [nb_ba],       'd0b','write_psi')
      CALL alloc_NParray(psid0b_k,[nb_ba],       'psid0b_k','write_psi')
      CALL alloc_NParray(Qact1,   [mole%nb_act1],'Qact1','write_psi')
      CALL alloc_NParray(psi_q,   [nb_psi,n_h],  'psi_q','write_psi')


!      - check the phase of psi(:,i) ----------
       DO i=1,nb_psi
         IF (Mat_psi(1,i) < ZERO) Mat_psi(:,i) = -Mat_psi(:,i)
       END DO


       DO i_q=1,nb_qa

        CALL Rec_Qact(Qact1,BasisnD,i_q,mole)

        IF (i_q == 1) Q1 = Qact1(1)

!       - calculation of WrhonD ------------------------------
         WrhonD = Rec_WrhonD(BasisnD,i_q)
         CALL calc_d0b(d0b,BasisnD,i_q)

         DO i=1,nb_psi

           psi_q(i,:) = ZERO
           IF (psi2) THEN


             DO ih=1,n_h
               ihk = (ih-1)*nb_ba
               DO k=1,nb_ba
                 ihk = ihk + 1
                 psid0b_k(k) = Mat_psi(ihk,i) * d0b(k)
               END DO

               DO k=1,nb_ba
               DO l=1,nb_ba
                 psi_q(i,1) = psi_q(i,1) + psid0b_k(k) * psid0b_k(l)
               END DO
               END DO

             END DO


           ELSE


             DO ih=1,n_h
               ihk = (ih-1)*nb_ba
               psi_q(i,ih)=dot_product(Mat_psi(ihk+1:ihk+nb_ba,i),d0b(:))
             END DO

           END IF

         END DO

         ih_print = n_h
         IF (psi2) ih_print = 1

          IF (mole%nb_act1 == 2 .AND. abs(Q1-Qact1(1)) > ONETENTH**6) THEN
             write(out_unit,*)
             Q1 = Qact1(1)
          END IF

         write(out_unit,31) Qact1(:),WrhonD,psi_q(1:min(max_print,nb_psi),1:ih_print)
 31      format(3f20.10,200f20.10)

       END DO

      CALL dealloc_NParray(d0b,     'd0b',     'write_psi')
      CALL dealloc_NParray(psid0b_k,'psid0b_k','write_psi')
      CALL dealloc_NParray(Qact1,   'Qact1',   'write_psi')
      CALL dealloc_NParray(psi_q,   'psi_q',   'write_psi')
!----------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*) 'END write_psi'
        END IF
!----------------------------------------------------------


      end subroutine write_psi

      SUBROUTINE write_psi2_new(Tab_Psi)

      USE EVR_system_m
      USE mod_basis
      USE mod_psi,      ONLY : param_psi,sub_PsiBasisRep_TO_GridRep
      IMPLICIT NONE


      TYPE (param_psi)   :: Tab_Psi(:)


      integer           :: nb_psi


!------ working variables ---------------------------------
      integer       :: i,i_q,ie,ieq,nb_bi,nb_qa
      real (kind=Rkind), allocatable :: Qact1(:)
      real (kind=Rkind), allocatable :: psi_q(:,:)
      integer, parameter :: max_print = 200

      real (kind=Rkind)              :: WrhonD

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='write_psi2_new'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nb_psi = size(Tab_Psi)
      IF (nb_psi < 1) RETURN

      nb_bi  = Tab_Psi(1)%nb_bi
      nb_qa  = Tab_Psi(1)%nb_qa
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_psi',nb_psi
      END IF
!-----------------------------------------------------------
      write(out_unit,*)
      write(out_unit,*) 'eigenvectors on a grid',nb_psi


      CALL alloc_NParray(psi_q,[nb_qa,nb_psi],'psi_q',name_sub)
      CALL alloc_NParray(Qact1,[Tab_Psi(1)%nb_act1],'Qact1',name_sub)


!-----------------------------------------------------------


       DO i=1,nb_psi
         CALL sub_PsiBasisRep_TO_GridRep(Tab_Psi(i))

         DO i_q=1,nb_qa
           psi_q(i_q,i)  = ZERO
           DO ie=1,Tab_Psi(i)%nb_be*nb_bi
             ieq = (ie-1)*nb_qa + i_q
             IF (Tab_Psi(i)%cplx) THEN
               psi_q(i_q,i) = psi_q(i_q,i) + abs(Tab_Psi(i)%CvecG(ieq))**2
             ELSE
               psi_q(i_q,i) = psi_q(i_q,i) + Tab_Psi(i)%RvecG(ieq)**2
             END IF
           END DO

         END DO
       END DO


       DO i_q=1,nb_qa

         CALL Rec_x(Qact1,Tab_Psi(1)%BasisnD,i_q)

         !- calculation of WrhonD ------------------------------
         WrhonD = Rec_WrhonD(Tab_Psi(1)%BasisnD,i_q)

         write(out_unit,31) Qact1(:),WrhonD,psi_q(i_q,1:nb_psi)
 31      format(3f20.10,200f20.10)

       END DO

      CALL dealloc_NParray(psi_q,'psi_q',name_sub)
      CALL dealloc_NParray(Qact1,'Qact1',name_sub)


!----------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*) 'END ',name_sub
        END IF
!----------------------------------------------------------


      end subroutine write_psi2_new

      SUBROUTINE write_cube(Tab_Psi)

      USE EVR_system_m
      USE mod_psi,      ONLY : param_psi,sub_PsiBasisRep_TO_GridRep
      IMPLICIT NONE


      TYPE (param_psi)   :: Tab_Psi(:)


      integer           :: nb_psi


!------ working variables ---------------------------------
      integer       :: i,i_q,ie,ieq,nb_bi,nb_qa
      real (kind=Rkind), allocatable :: psi_q(:,:)
      integer, parameter :: max_print = 200


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='write_cube'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nb_psi = size(Tab_Psi)
      IF (nb_psi < 1) RETURN

      nb_bi  = Tab_Psi(1)%nb_bi
      nb_qa  = Tab_Psi(1)%nb_qa
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_psi',nb_psi
      END IF
!-----------------------------------------------------------
      write(out_unit,*)
      write(out_unit,*) 'eigenvectors on a cube (test)',nb_psi

!-----------------------------------------------------------

      CALL alloc_NParray(psi_q,[nb_qa,nb_psi],'psi_q',name_sub)


       DO i=1,nb_psi
         CALL sub_PsiBasisRep_TO_GridRep(Tab_Psi(i))

         DO i_q=1,nb_qa
           psi_q(i_q,i)  = ZERO
           DO ie=1,Tab_Psi(i)%nb_be*nb_bi
             ieq = (ie-1)*nb_qa + i_q
             IF (Tab_Psi(i)%cplx) THEN
               psi_q(i_q,i) = psi_q(i_q,i) + abs(Tab_Psi(i)%CvecG(ieq))**2
             ELSE
               psi_q(i_q,i) = psi_q(i_q,i) + Tab_Psi(i)%RvecG(ieq)**2
             END IF
           END DO

         END DO
         !write cube file one for each Tab_Psi(i)
         write(100+i,*)  psi_q(:,i)
       END DO

      CALL dealloc_NParray(psi_q,'psi_q',name_sub)

!----------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*) 'END ',name_sub
        END IF
!----------------------------------------------------------


      end subroutine write_cube

      ! write basis b^n(Q)
      SUBROUTINE write_basis_biqi(psi)
        USE EVR_system_m
        USE mod_nDindex
        USE mod_psi,                         ONLY: param_psi
        IMPLICIT NONE

        TYPE(param_psi),              intent(in) :: psi
        Integer                                  :: nDval(psi%BasisnD%nDindB%ndim)
        Integer                                  :: ndim_AT_ib(psi%BasisnD%nDindB%ndim)

        Integer                                  :: ib,nD
        Integer                                  :: iQ,nQ
        Integer                                  :: ibiq
        Integer                                  :: LG


        ! write first 3 basis b_ib(Q)
        LG=psi%BasisnD%L_SparseGrid

        ndim_AT_ib=0
        DO ib=1,psi%nb_ba
          CALL calc_nDindex(psi%BasisnD%nDindB,ib,nDval)
          DO nD=1,psi%BasisnD%nDindB%ndim
            ibiq=nDval(nD)
            ndim_AT_ib(nD)=max(ibiq,ndim_AT_ib(nD))
          ENDDO
        ENDDO

        DO nD=1,3 ! psi%BasisnD%nDindB%ndim in total
          nQ=size(Psi%BasisnD%tab_basisPrimSG(LG,nD)%dnRGB%d0,1)
          write(out_unit,*) "b(Q): ",nD,ndim_AT_ib(nD),nQ,LG ! dimension: b*Q
          ! tab_basisPrimSG(0 ~ L_SparseGrid,nD)
          DO iQ=1,nQ
            write(out_unit,51) psi%BasisnD%tab_basisPrimSG(LG,nD)&
                                   %dnRGB%d0(iQ,1:ndim_AT_ib(nD))
          ENDDO

          DO iQ=1,nQ
            write(out_unit,51) Psi%BasisnD%tab_basisPrimSG(LG,nD)%x(:,iQ)
          ENDDO
        ENDDO

51      format(100(e16.6))
      END SUBROUTINE

END MODULE mod_fullanalysis
