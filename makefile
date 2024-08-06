#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx)
# F90 = mpifort
 FC = gfortran
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 1
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK = 1
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT = 4
## extension for the "sub_system." file. Possible values: f; f90
extf = f
#
## how to get external libraries;  "loc" (default): from local zip file, Empty or something else (v0.5): from github
EXTLIB_TYPE = loc
#=================================================================================
#=================================================================================
ifeq ($(FC),)
  FFC      := gfortran
else
  FFC      := $(FC)
endif
ifeq ($(OPT),)
  OOPT      := 1
else
  OOPT      := $(OPT)
endif
ifeq ($(OMP),)
  OOMP      := 1
else
  OOMP      := $(OMP)
endif
ifeq ($(LAPACK),)
  LLAPACK      := 1
else
  LLAPACK      := $(LAPACK)
endif
#===============================================================================
# setup for mpifort
ifeq ($(FFC),mpifort)
  ## MPI compiled with: gfortran or ifort
  MPICORE := $(shell ompi_info | grep 'Fort compiler:' | awk '{print $3}')
  OOMP = 0
endif
#===============================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)

# about EVRT, path, versions ...:
LOC_path:= $(shell pwd)

# Extension for the object directory and the library
ifeq ($(FFC),mpifort)
  extlibwi_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
endif
extlib_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)

OBJ_DIR = obj/obj$(extlibwi_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA=libEVR$(extlibwi_obj).a
#=================================================================================
# cpp preprocessing
CPPSHELL_LAPACK  = -D__LAPACK="$(LLAPACK)"

#===============================================================================
#
#===============================================================================
# external lib (QML, AD_dnSVM ...)
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(LOC_path)/Ext_Lib
endif

QD_DIR            = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR         = $(QD_DIR)/OBJ/obj$(extlib_obj)
QDLIBA            = $(QD_DIR)/libQD$(extlib_obj).a

AD_DIR            = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR         = $(AD_DIR)/OBJ/obj$(extlib_obj)
ADLIBA            = $(AD_DIR)/libAD_dnSVM$(extlib_obj).a

QML_DIR           = $(ExtLibDIR)/QuantumModelLib
QMLMOD_DIR        = $(QML_DIR)/OBJ/obj$(extlib_obj)
QMLLIBA           = $(QML_DIR)/libQMLib$(extlib_obj).a

nDindex_DIR       = $(ExtLibDIR)/nDindex
nDindexMOD_DIR    = $(nDindex_DIR)/obj/obj$(extlib_obj)
nDindexLIBA       = $(nDindex_DIR)/libnDindex$(extlib_obj).a

EVRTdnSVM_DIR     = $(ExtLibDIR)/EVRT_dnSVM
EVRTdnSVMMOD_DIR  = $(EVRTdnSVM_DIR)/obj/obj$(extlib_obj)
EVRTdnSVMLIBA     = $(EVRTdnSVM_DIR)/libEVRT_dnSVM$(extlib_obj).a

FOREVRT_DIR       = $(ExtLibDIR)/FOR_EVRT
FOREVRTMOD_DIR    = $(FOREVRT_DIR)/obj/obj$(extlibwi_obj)
FOREVRTLIBA       = $(FOREVRT_DIR)/libFOR_EVRT$(extlibwi_obj).a

CONSTPHYS_DIR     = $(ExtLibDIR)/ConstPhys
CONSTPHYSMOD_DIR  = $(CONSTPHYS_DIR)/obj/obj$(extlibwi_obj)
CONSTPHYSLIBA     = $(CONSTPHYS_DIR)/libPhysConst$(extlibwi_obj).a

TNUMTANA_DIR      = $(ExtLibDIR)/Tnum-Tana
TNUMTANAMOD_DIR   = $(TNUMTANA_DIR)/obj/obj$(extlibwi_obj)
TNUMTANALIBA      = $(TNUMTANA_DIR)/libTnum-Tana$(extlibwi_obj).a

#EXTLib_pot        = /Users/lauvergn/trav/ElVibRot-work/exa_work/exa_C2H3p/TEST_EVRT/Lauvergnat_PA/tests/libpotFull.a
EXTLib_pot        = 

EXTLib     = $(EXTLib_pot) $(TNUMTANALIBA) $(CONSTPHYSLIBA) $(FOREVRTLIBA) $(QDLIBA)  $(ADLIBA) $(EVRTdnSVMLIBA) $(nDindexLIBA) $(QMLLIBA)
EXTMod     = -I$(TNUMTANAMOD_DIR) -I$(CONSTPHYSMOD_DIR) -I$(FOREVRTMOD_DIR) -I$(nDindexMOD_DIR) \
             -I$(EVRTdnSVMMOD_DIR) -I$(QMLMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR)
#===============================================================================
#
#===============================================================================
# gfortran (osx and linux)
#ifeq ($(F90),gfortran)
#===============================================================================
ifeq ($(F90),$(filter $(F90),gfortran gfortran-8))

  # opt management
  ifeq ($(OOPT),1)
    FFLAGS = -O5 -g -fbacktrace -funroll-loops -ftree-vectorize -falign-loops=16
  else
    FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
  endif

  # integer kind management
  ifeq ($(INT),8)
    FFLAGS += -fdefault-integer-8 -Dint8=1
  endif

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -fopenmp
  endif
  FFLAGS0 := $(FFLAGS)


  # where to store the .mod files
  FFLAGS +=-J$(MOD_DIR)

  # where to look the .mod files
  FFLAGS += $(EXTMod)

  # some cpreprocessing
  FFLAGS += -cpp $(CPPSHELL)
  ifeq ($(OMP),1)
    FFLAGS += -Drun_openMP=1
  endif

  FLIB   = $(EXTLib)
  # OS management
  ifeq ($(LLAPACK),1)
    ifeq ($(OS),Darwin)    # OSX
      # OSX libs (included lapack+blas)
      FLIB += -framework Accelerate
    else                   # Linux
      # linux libs
      FLIB += -llapack -lblas
      #
      # linux libs with mkl and with openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
      # linux libs with mkl and without openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
    endif
  endif

  FC_VER = $(shell $(FFC) --version | head -1 )

endif
#=================================================================================
#=================================================================================
#=================================================================================
# ifort compillation v17 v18 with mkl
#=================================================================================
ifeq ($(FFC),$(filter $(FFC),ifort ifx))

  # opt management
  ifeq ($(OOPT),1)
      FFLAGS = -O  -g -traceback -heap-arrays
  else
      FFLAGS = -O0 -check all -g -traceback
  endif

  # integer kind management
  ifeq ($(INT),8)
    FFLAGS += -i8 -Dint8=1
  endif

  # where to store the modules
  FFLAGS +=-module $(MOD_DIR)

 # omp management
  ifeq ($(OOMP),1)
    ifeq ($(FFC),ifort)
      FFLAGS += -qopenmp -parallel
    else # ifx
      FFLAGS += -qopenmp
    endif
  endif

  # where to look the .mod files
  FFLAGS += $(EXTMod)

  # some cpreprocessing
  FFLAGS += -cpp $(CPPSHELL)
  ifeq ($(OMP),1)
    FFLAGS += -Drun_openMP=1
  endif

  FLIB    = $(EXTLib)
  ifneq ($(LLAPACK),1)
    ifeq ($(FFC),ifort)
      FLIB += -mkl -lpthread
    else # ifx
      FLIB += -qmkl -lpthread
    endif
  else
    FLIB += -lpthread
  endif

  FC_VER = $(shell $(F90) --version | head -1 )

endif
#=================================================================================
#=================================================================================
#===============================================================================
#===============================================================================
$(info ************************************************************************)
$(info ***********OS:               $(OS))
$(info ***********COMPILER:         $(FFC))
$(info ***********OPTIMIZATION:     $(OOPT))
$(info ***********COMPILER VERSION: $(FC_VER))
ifeq ($(FFC),mpifort)
$(info ***********COMPILED with:    $(MPICORE))
endif
$(info ***********OpenMP:           $(OOMP))
$(info ***********Lapack:           $(LLAPACK))
$(info ***********FFLAGS0:          $(FFLAGS0))
$(info ***********FLIB:             $(FLIB))
$(info ***********ExtLibDIR:        $(ExtLibDIR))
$(info ***********extf [sub_system]:$(extf))
$(info ************************************************************************)
$(info ************************************************************************)

#==========================================
VPATH = Source_ElVibRot/sub_Basis Source_ElVibRot/sub_Basis/sub_Basis_SG4 \
  Source_ElVibRot/sub_Basis/sub_ReducedDensity Source_ElVibRot/sub_Basis/sub_SymAbelian \
  Source_ElVibRot/sub_CRP Source_ElVibRot/sub_GWP Source_ElVibRot/sub_Operator \
  Source_ElVibRot/sub_Optimization Source_ElVibRot/sub_Smolyak_test Source_ElVibRot/sub_WP \
  Source_ElVibRot/sub_active Source_ElVibRot/sub_analysis Source_ElVibRot/sub_data_initialisation Source_ElVibRot/sub_inactive \
  Source_ElVibRot/sub_main Source_ElVibRot/sub_propagation Source_ElVibRot/sub_rotation sub_pot sub_operator_T



#SRCFILES=  $(basis_SRCFILES) $(main_SRCFILES) $(EVR-Mod_SRCFILES) 
include ./fortranlist.mk

OBJ0=${SRCFILES:.f90=.o}
OBJ0 += QMRPACK_lib.o

OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
OBJ0_EXT= sub_system.o calc_f2_f1Q.o Sub_X_TO_Q_ana.o Calc_Tab_dnQflex.o
OBJ_EXT=$(addprefix $(OBJ_DIR)/, $(OBJ0_EXT))
$(info ************ OBJ_EXT: $(OBJ_EXT))

#
#===============================================
#============= Several mains ===================
#===============================================
#===============================================
#==============================================
#ElVibRot:
VIBEXE  = vib.exe
VIBMAIN = EVR-T


#make all : EVR
.PHONY: all evr EVR libEVR libevr lib
evr EVR all :obj vib $(VIBEXE)
	@echo "EVR OK"
lib libEVR libevr: $(LIBA)
	@echo $(LIBA) " OK"
#
# vib script
.PHONY: vib
vib:
	@echo "make vib script"
	./scripts/make_vib.sh $(LOC_path) $(FFC) $(extf)
	chmod a+x vib
#
$(VIBEXE): $(OBJ_DIR)/$(VIBMAIN).o $(OBJ_EXT) $(LIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o $(VIBEXE) $(OBJ_DIR)/$(VIBMAIN).o $(LIBA) $(OBJ_EXT) $(FLIB)
	@echo EVR-T
#===============================================
#============= TESTS ===========================
#===============================================
.PHONY: ut UT
ut Ut:
	@echo "Unitary test"
	@cd UnitTests/HCN-WP_UT ; ./run_tests
	@cd UnitTests/HCN_UT    ; ./run_tests
	@cd UnitTests/HNO3_UT   ; ./run_tests
#===============================================
#============= Library: lib_FOR_EVRT.a  ========
#===============================================
$(LIBA): $(OBJ) | $(EXTLib)
	@echo "  LIBA from OBJ files"
	ar -cr $(LIBA) $(OBJ)
	@echo "  done Library: "$(LIBA)
#
#===============================================
#============= make sub_system =================
#=============  with the .f or .f90 extention ==
#===============================================
sub_pot/sub_system.$(extf): sub_pot/sub_system_save.$(extf)
	@echo "cp sub_system.... from sub_pot"
	cp sub_pot/sub_system_save.$(extf) sub_pot/sub_system.$(extf)
#
#===============================================
#============= compilation =====================
#===============================================
$(OBJ_DIR)/%.o: %.f90
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
$(OBJ_DIR)/%.o: %.f
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
$(OBJ_DIR)/sub_system.o: sub_pot/sub_system.$(extf)
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
#===============================================
#================ cleaning =====================
.PHONY: clean cleanall
clean:
	rm -f  $(OBJ_DIR)/*.o
	rm -f *.log 
	rm -f vib.exe
	@echo "  done cleaning"

cleanall : clean clean_extlib
	rm -fr obj/* build
	rm -f lib*.a
	rm -f *.exe
	rm -f TESTS/res* TESTS/*log
	cd UnitTests ; ./clean
	cd Examples ; ./clean
	cd Working_tests ; ./clean
	@echo "  done all cleaning"
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := EVR
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	$(ExtLibSAVEDIR)/makezip.sh $(BaseName)
	@echo "  done zip"
#===============================================
#=== external libraries ========================
# AD_dnSVM + QML Libs ...
#===============================================
#	@test -d $(TNUMTANA_DIR) || (cd $(ExtLibDIR) ; ./get_Tnum-Tana.sh $(EXTLIB_TYPE))
#	test -e $(TNUMTANALIBA) && ln -s $(TNUMTANALIBA_old) $(TNUMTANALIBA)

TNUMTANALIBA_old=$(TNUMTANA_DIR)/libCoord_KEO_PrimOp_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT).a
$(TNUMTANALIBA):
	@test -d $(ExtLibDIR)    || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(TNUMTANA_DIR) || (cd $(ExtLibDIR) ; ./get_Tnum-Tana.sh $(EXTLIB_TYPE))
	@test -d $(TNUMTANA_DIR) || (echo $(TNUMTANA_DIR) "does not exist" ; exit 1)
	cd $(TNUMTANA_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(TNUMTANA_DIR) " in "$(BaseName)
#
$(CONSTPHYSLIBA):
	@test -d $(ExtLibDIR)     || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(CONSTPHYS_DIR) || (cd $(ExtLibDIR) ; ./get_ConstPhys.sh $(EXTLIB_TYPE))
	@test -d $(CONSTPHYS_DIR) || (echo $(CONSTPHYS_DIR) "does not exist" ; exit 1)
	cd $(CONSTPHYS_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(CONSTPHYS_DIR) " in "$(BaseName)
#
$(nDindexLIBA):
	@test -d $(ExtLibDIR)   || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(nDindex_DIR) || (cd $(ExtLibDIR) ; ./get_nDindex.sh  $(EXTLIB_TYPE))
	@test -d $(nDindex_DIR) || (echo $(nDindex_DIR) "does not exist" ; exit 1)
	cd $(nDindex_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(nDindex_DIR) " in "$(BaseName)
#
$(EVRTdnSVMLIBA):
	@test -d $(ExtLibDIR)     || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(EVRTdnSVM_DIR) || (cd $(ExtLibDIR) ; ./get_EVRT_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(EVRTdnSVM_DIR) || (echo $(EVRTdnSVM_DIR) "does not exist" ; exit 1)
	cd $(EVRTdnSVM_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(EVRTdnSVM_DIR) " in "$(BaseName)
#
$(FOREVRTLIBA):
	@test -d $(ExtLibDIR)   || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(FOREVRT_DIR) || (cd $(ExtLibDIR) ; ./get_FOR_EVRT.sh $(EXTLIB_TYPE))
	@test -d $(FOREVRT_DIR) || (echo $(FOREVRT_DIR) "does not exist" ; exit 1)
	cd $(FOREVRT_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(FOREVRTLIBA) " in "$(BaseName)
#
$(QMLLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QML_DIR)   || (cd $(ExtLibDIR) ; ./get_QML.sh $(EXTLIB_TYPE))
	@test -d $(QML_DIR)   || (echo $(QML_DIR) "does not exist" ; exit 1)
	cd $(QML_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(AD_DIR)    || (cd $(ExtLibDIR) ; ./get_AD_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(AD_DIR)    || (echo $(AD_DIR) "does not exist" ; exit 1)
	cd $(AD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(AD_DIR) " in "$(BaseName)
#
$(QDLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QD_DIR)    || (cd $(ExtLibDIR) ; ./get_QDUtilLib.sh $(EXTLIB_TYPE))
	@test -d $(QD_DIR)    || (echo $(QD_DIR) "does not exist" ; exit 1)
	cd $(QD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
##
.PHONY: clean_extlib
clean_extlib:
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	cd $(ExtLibDIR) ; ./cleanlib
##
#=======================================================================================
#=======================================================================================
#add dependence for parallelization
#	@echo "OBJ with EXTLib"
$(OBJ) : | $(EXTLib)

#===============================================
#===============================================
#============= make dependencies =============
#===============================================
.PHONY: dep
dependencies.mk fortranlist.mk dep:
	./scripts/dependency.sh
#===============================================
include ./dependencies.mk
