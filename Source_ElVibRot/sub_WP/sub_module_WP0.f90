!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
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
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!
!===========================================================================
!===========================================================================
      MODULE mod_WP0
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: sub_read_psi0
      CONTAINS

      SUBROUTINE sub_read_psi0(psi0,para_WP0,max_WP,symab,ortho)
        USE mod_system
        USE mod_psi_set_alloc
        USE mod_ana_psi
        USE mod_psi_io
        USE mod_psi_Op
        USE mod_param_WP0
        IMPLICIT NONE
      
      !----- variables for the WP propagation ----------------------------
        integer,             intent(in)               :: max_WP
        TYPE (param_WP0),    intent(inout)            :: para_WP0
        TYPE (param_psi),    intent(inout)            :: psi0(max_WP)
        integer,             intent(in),   optional   :: symab
        logical,             intent(in),   optional   :: ortho
      
      !------ working parameters --------------------------------
        logical                  :: cplx,ortho_loc
        integer                  :: i,symab_loc
      
        integer                  :: nioWP
        character (len=Line_len) :: line
      
        integer  :: nb_readWP_file,nb_tot_file,n1,n2,nb_WP0
        integer  :: ilist
        integer, allocatable :: list_nDindBasis1_TO_nDindBasis2(:)
        integer, allocatable :: list_readWP(:)
      
        integer   :: Version_File,nb_psi,nb_tot
        character :: ch
        namelist / headerFile / Version_File,nb_psi,nb_tot
      
      !----- for debuging --------------------------------------------------
        integer :: err_mem,memory,ioerr
        character (len=*), parameter ::name_sub='sub_read_psi0'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) ' para_WP0%nb_WP0            ',para_WP0%nb_WP0
          write(out_unitp,*) ' para_WP0%read_listWP0      ',para_WP0%read_listWP0
          write(out_unitp,*) ' para_WP0%read_file         ',para_WP0%read_file
          write(out_unitp,*) ' para_WP0%file_WP0          ',trim(adjustl(para_WP0%file_WP0%name))
          write(out_unitp,*) ' para_WP0%WP0cplx           ',para_WP0%WP0cplx
          write(out_unitp,*) ' para_WP0%file_WP0%formatted',para_WP0%file_WP0%formatted
          write(out_unitp,*) ' max_WP                     ',max_WP
          write(out_unitp,*)
        END IF
      
      !------ initialization -------------------------------------
        IF (para_WP0%nb_WP0 < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' para_WP0%nb_WP0 < 0!!',para_WP0%nb_WP0
          STOP
        END IF
      
        IF (present(symab)) THEN
          symab_loc = symab
        ELSE
          symab_loc = -1
        END IF
      
        IF (present(ortho)) THEN
          ortho_loc = ortho
        ELSE
          ortho_loc = .FALSE.
        END IF
      
      !------ read guess vectors ---------------------------------
        cplx = para_WP0%WP0cplx
      
        IF (para_WP0%read_file) THEN
          CALL file_open(para_WP0%file_WP0,nioWP,                         &
                         lformatted=para_WP0%file_WP0%formatted,old=.TRUE.)
      
          ! first we read the header of the file. There are 3 options:
          ! - one integer:  nb_readWP_file,                 hence Version_File=0, formatted only
          ! - two integers: nb_readWP_file and nb_tot_file, hence Version_File=0, formatted only
          ! - the &headerFile namelist,                     hence Version_File=1, formatted and unformatted (not a namelist)
          CALL Read_header_saveFile_psi(psi0,nb_readWP_file,              &
                                        list_nDindBasis1_TO_nDindBasis2,  &
                                        para_WP0%file_WP0,Version_File)
      
          write(out_unitp,*) ' file name for WP0: ',trim(para_WP0%file_WP0%name)
          write(out_unitp,*) ' nb_readWP_file: ',nb_readWP_file
          write(out_unitp,*) ' nb_tot_file   : ',size(list_nDindBasis1_TO_nDindBasis2)
          write(out_unitp,*) ' Version_File  : ',Version_File
          flush(out_unitp)
      
          IF (nb_readWP_file < para_WP0%nb_WP0 .OR. para_WP0%nb_WP0 == 0) THEN
             write(out_unitp,*) ' WARNING in ',name_sub
             write(out_unitp,*) ' The number of WP0 in the file (',               &
                trim(para_WP0%file_WP0%name),') is smaller than the requested number'
             write(out_unitp,*) ' OR the request number is zero'
             write(out_unitp,*) 'nb_WP0,nb_readWP_file',                          &
                          para_WP0%nb_WP0,nb_readWP_file
             para_WP0%nb_WP0 = nb_readWP_file
          END IF
      
          CALL alloc_NParray(list_readWP,[para_WP0%nb_WP0],"list_readWP",name_sub)
          write(out_unitp,*) 'nb_WP0,nb_readWP_file',para_WP0%nb_WP0,nb_readWP_file
          IF (para_WP0%read_listWP0) THEN
            read(in_unitp,*) list_readWP(1:para_WP0%nb_WP0)
          ELSE
            list_readWP(:) = [(i,i=1,para_WP0%nb_WP0)]
          END IF
      
          IF (debug) write(out_unitp,*) 'list_readWP:',list_readWP(1:para_WP0%nb_WP0)
          flush(out_unitp)
      
      
          IF (.NOT. allocated(list_nDindBasis1_TO_nDindBasis2) ) THEN ! Version_File=0, option=1
            ilist = 1
            DO i=1,nb_readWP_file
              IF (debug .OR. print_level > 1) write(out_unitp,*) i,ilist
              flush(out_unitp)
              CALL lect_psiBasisRepnotall_nD(psi0(ilist),nioWP,cplx,      &
                                      para_WP0%file_WP0%formatted)
              IF (list_readWP(ilist) == i) ilist = ilist + 1
              IF (ilist > para_WP0%nb_WP0 .OR. ilist > max_WP) EXIT
            END DO
            para_WP0%nb_WP0 = ilist - 1
            IF (debug) write(out_unitp,*) ' read with lect_psiBasisRepnotall_nD ' // &
                                       '(Version_File=0, option=1): done'
          ELSE ! Version_File=0, option=2 or Version_File=1 (with the namelist)
            nb_tot_file = size(list_nDindBasis1_TO_nDindBasis2)
      
            ilist = 1
            DO i=1,nb_readWP_file
              IF (debug .OR. print_level > 1) write(out_unitp,*) 'i,ilist',i,ilist
              flush(out_unitp)
      
              CALL Read_psi_nDBasis(psi0(ilist),nioWP,                     &
                                 para_WP0%file_WP0%formatted,Version_File, &
                               list_nDindBasis1_TO_nDindBasis2,nb_tot_file)
      
              IF (list_readWP(ilist) == i) ilist = ilist + 1
              IF (ilist > para_WP0%nb_WP0 .OR. ilist > max_WP) EXIT
            END DO
            para_WP0%nb_WP0 = ilist - 1
      
            CALL dealloc_NParray(list_nDindBasis1_TO_nDindBasis2,         &
                                "list_nDindBasis1_TO_nDindBasis2",name_sub)
      
            IF (debug) write(out_unitp,*) ' read with Read_psi_nDBasis ' //&
                      '(Version_File=0, option=2 or Version_File=1): done'
          END IF
      
          CALL file_close(para_WP0%file_WP0)
          CALL dealloc_NParray(list_readWP,"list_readWP",name_sub)
        ELSE
          DO i=1,para_WP0%nb_WP0
            write(out_unitp,*) ' Read in',in_unitp,i
            flush(out_unitp)
            CALL alloc_psi(psi0(i),BasisRep=.TRUE.)
            CALL lect_psiBasisRepnotall_nD(psi0(i),in_unitp,cplx,.TRUE.)
            write(out_unitp,*) ' write in',out_unitp,i
            flush(out_unitp)
            IF(keep_MPI) CALL ecri_psiBasisRepnotall_nD(psi0(i),out_unitp,ONETENTH**4,.TRUE.,i)
            IF (debug) write(out_unitp,*) ' read with lect_psiBasisRepnotall_nD ' // &
                                '(Version_File=0, option=1 ???): done'
          END DO
        END IF
      
      !-----Check the norm and Check and set the symmetry--------------------
        DO i=1,para_WP0%nb_WP0
          CALL Set_symab_OF_psiBasisRep(psi0(i))
          CALL renorm_psi(psi0(i))
          IF (debug .OR. print_level > 1) write(out_unitp,*) ' norm psi0(i),symab',i,psi0(i)%norm2,psi0(i)%symab
          flush(out_unitp)
        END DO
      
      !-----Check the norm and Check and set the symmetry--------------------
        ! keep only the right symmetry if symab > -1
        IF (symab_loc > -1) THEN
          nb_WP0 = 0
      
          DO i=1,para_WP0%nb_WP0
            IF (psi0(i)%symab == symab_loc) THEN
              nb_WP0 = nb_WP0 + 1
              ! copy the psi0 only if nb_WP0 /= i
              IF (nb_WP0 /= i) psi0(nb_WP0) = psi0(i)
            END IF
          END DO
      
          DO i=nb_WP0+1,para_WP0%nb_WP0
            CALL dealloc_psi(psi0(i))
          END DO
      
          para_WP0%nb_WP0 = nb_WP0
      
        END IF
        write(out_unitp,*) 'Number of read vector(s) after the symmetry analysis:',para_WP0%nb_WP0
        ! orthogonalisation
        IF (ortho_loc) CALL sub_Schmidt(psi0,para_WP0%nb_WP0)
      
        write(out_unitp,*) 'Number of read vector(s):',para_WP0%nb_WP0
      
        IF (debug) THEN
          write(out_unitp,*) ' ',para_WP0%nb_WP0,' read WP0s'
          DO i=1,para_WP0%nb_WP0
            write(out_unitp,*) 'psiO(i)',i
            write(out_unitp,*) 'symab, bits(symab)',WriteTOstring_symab(psi0(i)%symab)
            CALL ecri_psi(psi=psi0(i))
          END DO
          write(out_unitp,*)
          write(out_unitp,*) 'END ',name_sub
        END IF
        END SUBROUTINE sub_read_psi0
      !=======================================================================================

      END MODULE mod_WP0
