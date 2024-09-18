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
## Arpack? Empty: default No Arpack; 0: without Arpack; 1 with Arpack
ARPACK = 0
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
ifeq ($(ARPACK),)
  AARPACK      := 0
else
  AARPACK      := $(ARPACK)
endif
#===============================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)
#
# about EVRT, path, versions ...:
MAIN_path:= $(shell pwd)
#
#===============================================================================
# We cannot use ARPACK without lapack
ifeq ($(LLAPACK),0)
  AARPACK = 0
endif
#===============================================================================
# turn off ARPACK when using pgf90
ifeq ($(F90),pgf90)
  AARPACK = 0
endif
#===============================================================================
# setup for mpifort
ifeq ($(FFC),mpifort)
  ## MPI compiled with: gfortran or ifort
  MPICORE := $(shell ompi_info | grep 'Fort compiler:' | awk '{print $3}')
  OOMP = 0
  ifeq ($(INT),8)
    AARPACK = 0 ## temp here, disable ARPACK for 64-bit case
  endif
endif

# Extension for the object directory and the library
ext_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
ifeq ($(FFC),mpifort)
  extlibwi_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj:=$(ext_obj)
endif

OBJ_DIR = obj/obj$(extlibwi_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA=libEVR$(extlibwi_obj).a
#=================================================================================
#
#===============================================================================
# external lib : QDUtilLib AD_dnSVM ConstPhys QuantumModelLib nDindex EVRT_dnSVM FOR_EVRT
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(MAIN_path)/Ext_Lib
endif
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA    = $(QD_DIR)/libQD$(ext_obj).a

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(ext_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(ext_obj).a

CONSTPHYS_DIR    = $(ExtLibDIR)/ConstPhys
CONSTPHYSMOD_DIR = $(CONSTPHYS_DIR)/obj/obj$(extlibwi_obj)
CONSTPHYSLIBA    = $(CONSTPHYS_DIR)/libPhysConst$(extlibwi_obj).a

QML_DIR    = $(ExtLibDIR)/QuantumModelLib
QMLMOD_DIR = $(QML_DIR)/OBJ/obj$(ext_obj)
QMLLIBA    = $(QML_DIR)/libQMLib$(ext_obj).a

nDindex_DIR    = $(ExtLibDIR)/nDindex
nDindexMOD_DIR = $(nDindex_DIR)/obj/obj$(ext_obj)
nDindexLIBA    = $(nDindex_DIR)/libnDindex$(ext_obj).a

EVRTdnSVM_DIR    = $(ExtLibDIR)/EVRT_dnSVM
EVRTdnSVMMOD_DIR = $(EVRTdnSVM_DIR)/obj/obj$(ext_obj)
EVRTdnSVMLIBA    = $(EVRTdnSVM_DIR)/libEVRT_dnSVM$(ext_obj).a

FOREVRT_DIR    = $(ExtLibDIR)/FOR_EVRT
FOREVRTMOD_DIR = $(FOREVRT_DIR)/obj/obj$(extlibwi_obj)
FOREVRTLIBA    = $(FOREVRT_DIR)/libFOR_EVRT$(extlibwi_obj).a


TNUMTANA_DIR      = $(ExtLibDIR)/Tnum-Tana
TNUMTANAMOD_DIR   = $(TNUMTANA_DIR)/obj/obj$(extlibwi_obj)
TNUMTANALIBA      = $(TNUMTANA_DIR)/libTnum-Tana$(extlibwi_obj).a
# AARPACK library
ARPACK_DIR        = $(ExtLibDIR)/ARPACK_EVR
ifeq ($(AARPACK),1)
  ARPACKLIBA = $(ARPACK_DIR)/libarpack$(extlibwi_obj).a
else
  ARPACKLIBA =
endif

#EXTLib_pot        = /Users/lauvergn/trav/ElVibRot-work/exa_work/exa_C2H3p/TEST_EVRT/Lauvergnat_PA/tests/libpotFull.a
EXTLib_pot        = 
#
EXTLib     = $(EXTLib_pot) $(ARPACKLIBA) \
             $(TNUMTANALIBA) $(CONSTPHYSLIBA) \
             $(FOREVRTLIBA) $(EVRTdnSVMLIBA) $(nDindexLIBA) \
             $(QMLLIBA) $(ADLIBA) $(QDLIBA)
EXTMod     = -I$(TNUMTANAMOD_DIR) -I$(CONSTPHYSMOD_DIR) -I$(FOREVRTMOD_DIR) -I$(nDindexMOD_DIR) \
             -I$(EVRTdnSVMMOD_DIR) -I$(QMLMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR)
#=================================================================================
#
#=================================================================================
# To deal with external compilers.mk file
CompilersDIR = $(MAIN_path)
ifeq ($(CompilersDIR),)
  include compilers.mk
else
  include $(CompilersDIR)/compilers.mk
endif
#=================================================================================
# cpp preprocessing
EVR_ver:=$(shell awk '/EVR/ {print $$3}' $(MAIN_path)/version-EVR)
ArpackFLAGS := $(FFLAGS)
ArpackFLAGS += -fallow-argument-mismatch
FFLAGS += -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
          -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
          -D__EVRTPATH="'$(MAIN_path)'" \
          -D__EVR_VER="'$(EVR_ver)'" \
          -D__ARPACK="$(AARPACK)"
#
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
$(info ***********Arpack:           $(AARPACK))
$(info ***********Arpack lib:       $(ARPACKLIBA))
#$(info ***********Arpack flags:     $(ArpackFLAGS))
$(info ***********ExtLibDIR:        $(ExtLibDIR))
$(info ***********extf [sub_system]:$(extf))
$(info ************************************************************************)
$(info ************************************************************************)

#==========================================
VPATH = APP SRC SRC/sub_Basis SRC/sub_Basis/sub_Basis_SG4 \
  SRC/sub_Basis/sub_ReducedDensity SRC/sub_Basis/sub_SymAbelian \
  SRC/sub_CRP SRC/sub_GWP SRC/sub_Operator \
  SRC/sub_Optimization SRC/sub_Smolyak_test SRC/sub_WP \
  SRC/sub_active SRC/sub_analysis SRC/sub_data_initialisation SRC/sub_inactive \
  SRC/sub_main SRC/sub_propagation SRC/sub_rotation sub_pot sub_operator_T



#SRCFILES=  $(basis_SRCFILES) $(main_SRCFILES) $(EVR-Mod_SRCFILES) 
include ./fortranlist.mk

OBJ0=${SRCFILES:.f90=.o}
OBJ0 += QMRPACK_lib.o

OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
#$(info ************ OBJ: $(OBJ))
#
OBJ0_EXT= read_para.o sub_system.o calc_f2_f1Q.o Sub_X_TO_Q_ana.o Calc_Tab_dnQflex.o
OBJ_EXT=$(addprefix $(OBJ_DIR)/, $(OBJ0_EXT))
#$(info ************ OBJ_EXT: $(OBJ_EXT))

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
	./scripts/make_vib.sh $(MAIN_path) $(FFC) $(extf)
	chmod a+x vib
#
$(VIBEXE): $(OBJ_DIR)/$(VIBMAIN).o $(OBJ_EXT) $(LIBA) | $(EXTLib)
	$(FFC) $(FFLAGS) -o $(VIBEXE) $(OBJ_DIR)/$(VIBMAIN).o $(LIBA) $(OBJ_EXT) $(EXTLib) $(FLIB)
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
	cd $(ARPACK_DIR) ; make clean
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
# QDUtilLib AD_dnSVM ConstPhys QuantumModelLib nDindex EVRT_dnSVM FOR_EVRT
#===============================================
#
DEV=dev
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib 		$(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM  		$(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh ConstPhys 		$(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh QuantumModelLib  $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh nDindex 			$(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh EVRT_dnSVM 		$(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh FOR_EVRT 		$(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh Tnum-Tana 		$(DEV)
#
$(QDLIBA): 
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib $(DEV)
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM $(DEV)
	cd $(ExtLibDIR)/AD_dnSVM ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(ADLIBA) || (echo $(ADLIBA) "does not exist" ; exit 1)
	@echo "  done " $(ADLIBA) " in "$(BaseName)
#
$(CONSTPHYSLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh ConstPhys $(DEV)
	cd $(ExtLibDIR)/ConstPhys ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(CONSTPHYSLIBA) || (echo $(CONSTPHYSLIBA) "does not exist" ; exit 1)
	@echo "  done " $(CONSTPHYSLIBA) " in "$(BaseName)
#
$(QMLLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh QuantumModelLib $(DEV)
	cd $(ExtLibDIR)/QuantumModelLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QMLLIBA) || (echo $(QMLLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QMLLIBA) " in "$(BaseName)
#
$(nDindexLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh nDindex $(DEV)
	cd $(ExtLibDIR)/nDindex ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(nDindexLIBA) || (echo $(nDindexLIBA) "does not exist" ; exit 1)
	@echo "  done " $(nDindexLIBA) " in "$(BaseName)
#
$(EVRTdnSVMLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh EVRT_dnSVM $(DEV)
	cd $(ExtLibDIR)/EVRT_dnSVM ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(EVRTdnSVMLIBA) || (echo $(EVRTdnSVMLIBA) "does not exist" ; exit 1)
	@echo "  done " $(EVRTdnSVMLIBA) " in "$(BaseName)
#
$(FOREVRTLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh FOR_EVRT $(DEV)
	cd $(ExtLibDIR)/FOR_EVRT ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(FOREVRTLIBA) || (echo $(FOREVRTLIBA) "does not exist" ; exit 1)
	@echo "  done " $(FOREVRTLIBA) " in "$(BaseName)
#
$(TNUMTANALIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh Tnum-Tana $(DEV)
	cd $(ExtLibDIR)/Tnum-Tana ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR)  CompilersDIR=$(CompilersDIR)
	@test -f $(TNUMTANALIBA) || (echo $(TNUMTANALIBA) "does not exist" ; exit 1)
	@echo "  done " $(TNUMTANALIBA) " in "$(BaseName)#
# ARPACK is needed
$(ARPACKLIBA):
	@test -d $(ARPACK_DIR)    || (echo $(ARPACK_DIR) "does not exist", exit 1)
	cd $(ARPACK_DIR) ; make lib home=$(ARPACK_DIR) FC=$(FFC) LIBEXT=$(ext_obj) FFLAGS="$(ArpackFLAGS)"
#
.PHONY: clean_extlib
clean_extlib:
	echo cleanlib, DIR=$(ExtLibDIR)
	cd $(ExtLibDIR) ; ./cleanlib
#
#=======================================================================================
.PHONY: fpm
fpm: getlib
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
