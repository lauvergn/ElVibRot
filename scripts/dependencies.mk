#===============================================
evr_system_m := $(OBJ_DIR)/EVR_system_m.o
mod_system := $(OBJ_DIR)/mod_system.o
mod_auto_basis := $(OBJ_DIR)/sub_Auto_Basis.o
mod_basis_btog_gtob_sgtype4 := $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4.o
mod_basis_btog_gtob_sgtype4_mpi := $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4_MPI.o
mod_basis_rcvec_sgtype4 := $(OBJ_DIR)/sub_module_basis_RCVec_SG4.o
mod_param_rd := $(OBJ_DIR)/sub_module_param_RD.o
mod_symabelian := $(OBJ_DIR)/sub_SymAbelian.o
basismakegrid := $(OBJ_DIR)/sub_module_BasisMakeGrid.o
mod_basis_l_to_n := $(OBJ_DIR)/sub_module_Basis_LTO_n.o
mod_rotbasis_param := $(OBJ_DIR)/sub_module_RotBasis.o
mod_basis := $(OBJ_DIR)/sub_module_basis.o
mod_basis_btog_gtob := $(OBJ_DIR)/sub_module_basis_BtoG_GtoB.o
mod_basis_btog_gtob_mpi := $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_MPI.o
mod_basis_btog_gtob_sgtype2 := $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SGType2.o
mod_basis_grid_param := $(OBJ_DIR)/sub_module_basis_Grid_Param.o
mod_basis_set_alloc := $(OBJ_DIR)/sub_module_basis_set_alloc.o
mod_param_sgtype2 := $(OBJ_DIR)/sub_module_param_SGType2.o
readbasis_m := $(OBJ_DIR)/sub_read_data.o
mod_crp := $(OBJ_DIR)/sub_CRP.o
mod_matop := $(OBJ_DIR)/sub_MatOp.o
mod_oppsi := $(OBJ_DIR)/sub_OpPsi.o
mod_oppsi_mpi := $(OBJ_DIR)/sub_OpPsi_MPI.o
mod_oppsi_sg4 := $(OBJ_DIR)/sub_OpPsi_SG4.o
mod_oppsi_sg4_mpi := $(OBJ_DIR)/sub_OpPsi_SG4_MPI.o
mod_op := $(OBJ_DIR)/sub_module_Op.o
mod_opgrid := $(OBJ_DIR)/sub_module_OpGrid.o
readop_m := $(OBJ_DIR)/sub_module_ReadOp.o
mod_setop := $(OBJ_DIR)/sub_module_SetOp.o
mod_bfgs := $(OBJ_DIR)/sub_module_BFGS.o
mod_optimization := $(OBJ_DIR)/sub_module_Optimization.o
mod_simulatedannealing := $(OBJ_DIR)/sub_module_SimulatedAnnealing.o
mod_smolyak_dind := $(OBJ_DIR)/sub_Smolyak_DInd.o
mod_smolyak_rdp := $(OBJ_DIR)/sub_Smolyak_RDP.o
mod_smolyak_ba := $(OBJ_DIR)/sub_Smolyak_ba.o
mod_smolyak_test := $(OBJ_DIR)/sub_Smolyak_module.o
mod_psi := $(OBJ_DIR)/mod_psi.o
mod_wp0 := $(OBJ_DIR)/sub_module_WP0.o
mod_ana_psi := $(OBJ_DIR)/sub_module_ana_psi.o
mod_ana_psi_mpi := $(OBJ_DIR)/sub_module_ana_psi_MPI.o
mod_param_wp0 := $(OBJ_DIR)/sub_module_param_WP0.o
mod_psi_b_to_g := $(OBJ_DIR)/sub_module_psi_B_TO_G.o
mod_psi_op := $(OBJ_DIR)/sub_module_psi_Op.o
mod_psi_op_mpi := $(OBJ_DIR)/sub_module_psi_Op_MPI.o
mod_psi_io := $(OBJ_DIR)/sub_module_psi_io.o
mod_psi_set_alloc := $(OBJ_DIR)/sub_module_psi_set_alloc.o
mod_type_ana_psi := $(OBJ_DIR)/sub_module_type_ana_psi.o
mod_set_pararph := $(OBJ_DIR)/sub_paraRPH.o
mod_fullanalysis := $(OBJ_DIR)/sub_analyse.o
mod_analysis := $(OBJ_DIR)/sub_module_analysis.o
mod_evr := $(OBJ_DIR)/EVR_Module.o
mod_ndgridfit := $(OBJ_DIR)/sub_main_nDfit.o
mod_hmax_mpi := $(OBJ_DIR)/sub_Hmax_MPI.o
mod_fullcontrol := $(OBJ_DIR)/sub_control.o
mod_arpack := $(OBJ_DIR)/sub_module_Arpack.o
mod_davidson := $(OBJ_DIR)/sub_module_Davidson.o
mod_davidson_mpi := $(OBJ_DIR)/sub_module_Davidson_MPI.o
mod_exactfact := $(OBJ_DIR)/sub_module_ExactFact.o
mod_filter := $(OBJ_DIR)/sub_module_Filter.o
mod_field := $(OBJ_DIR)/sub_module_field.o
mod_march := $(OBJ_DIR)/sub_module_propa_march.o
mod_march_mpi := $(OBJ_DIR)/sub_module_propa_march_MPI.o
mod_march_sg4 := $(OBJ_DIR)/sub_module_propa_march_SG4.o
mod_propa := $(OBJ_DIR)/sub_module_propagation.o
mod_propa_mpi := $(OBJ_DIR)/sub_module_propagation_MPI.o
mod_fullpropa := $(OBJ_DIR)/sub_propagation.o
#===============================================
#file+mod_name: SRC/EVR_system_m.f90 evr_system_m
$(OBJ_DIR)/EVR_system_m.o : \
          $(qdutil_m) \
          $(mod_mpi) \
          $(for_evrt_system_m) \
          $(iso_fortran_env) \
          $(tnumtana_system_m)
#file+mod_name: SRC/mod_system.f90 mod_system
$(OBJ_DIR)/mod_system.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_Basis/sub_Auto_Basis.f90 mod_auto_basis
$(OBJ_DIR)/sub_Auto_Basis.o : \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_mpi) \
          $(evr_system_m) \
          $(basismakegrid) \
          $(mod_constant) \
          $(mod_analysis) \
          $(mod_propa) \
          $(mod_psi) \
          $(mod_fullanalysis)
#file+mod_name: SRC/sub_Basis/sub_Basis_SG4/sub_module_basis_BtoG_GtoB_SG4.f90 mod_basis_btog_gtob_sgtype4
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux) \
          $(mod_basis_btog_gtob_sgtype4_mpi)
#file+mod_name: SRC/sub_Basis/sub_Basis_SG4/sub_module_basis_BtoG_GtoB_SG4_MPI.f90 mod_basis_btog_gtob_sgtype4_mpi
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4_MPI.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux) \
          $(mod_mpi)
#file+mod_name: SRC/sub_Basis/sub_Basis_SG4/sub_module_basis_RCVec_SG4.f90 mod_basis_rcvec_sgtype4
$(OBJ_DIR)/sub_module_basis_RCVec_SG4.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_Basis/sub_ReducedDensity/sub_module_param_RD.f90 mod_param_rd
$(OBJ_DIR)/sub_module_param_RD.o : \
          $(evr_system_m) \
          $(mod_ndindex)
#file+mod_name: SRC/sub_Basis/sub_SymAbelian/sub_SymAbelian.f90 mod_symabelian
$(OBJ_DIR)/sub_SymAbelian.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_Basis/sub_SymAbelian_OF_Basis.f90 
$(OBJ_DIR)/sub_SymAbelian_OF_Basis.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_basis_El.f90 
$(OBJ_DIR)/sub_basis_El.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_basis_NoGrid.f90 
$(OBJ_DIR)/sub_basis_NoGrid.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_module_BasisMakeGrid.f90 basismakegrid
$(OBJ_DIR)/sub_module_BasisMakeGrid.o : \
          $(evr_system_m) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_ndindex) \
          $(mod_dnsvm)
#file+mod_name: SRC/sub_Basis/sub_module_Basis_LTO_n.f90 mod_basis_l_to_n
$(OBJ_DIR)/sub_module_Basis_LTO_n.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_Basis/sub_module_RotBasis.f90 mod_rotbasis_param
$(OBJ_DIR)/sub_module_RotBasis.o : \
          $(evr_system_m) \
          $(mod_ndindex)
#file+mod_name: SRC/sub_Basis/sub_module_basis.f90 mod_basis
$(OBJ_DIR)/sub_module_basis.o : \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_rotbasis_param) \
          $(mod_basis_grid_param) \
          $(mod_basis_l_to_n) \
          $(mod_symabelian) \
          $(mod_param_sgtype2) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob) \
          $(evr_system_m) \
          $(mod_dnsvm) \
          $(mod_basis_btog_gtob_sgtype4)
#file+mod_name: SRC/sub_Basis/sub_module_basis_BtoG_GtoB.f90 mod_basis_btog_gtob
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB.o : \
          $(evr_system_m) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype2) \
          $(mod_param_sgtype2) \
          $(mod_basis_btog_gtob_mpi) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_Basis/sub_module_basis_BtoG_GtoB_MPI.f90 mod_basis_btog_gtob_mpi
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_MPI.o : \
          $(evr_system_m) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_Basis/sub_module_basis_BtoG_GtoB_SGType2.f90 mod_basis_btog_gtob_sgtype2
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SGType2.o : \
          $(evr_system_m) \
          $(mod_basis_set_alloc) \
          $(mod_module_dind)
#file+mod_name: SRC/sub_Basis/sub_module_basis_Grid_Param.f90 mod_basis_grid_param
$(OBJ_DIR)/sub_module_basis_Grid_Param.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_Basis/sub_module_basis_set_alloc.f90 mod_basis_set_alloc
$(OBJ_DIR)/sub_module_basis_set_alloc.o : \
          $(evr_system_m) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(qdutil_intvec_m) \
          $(mod_rotbasis_param) \
          $(mod_basis_grid_param) \
          $(mod_symabelian) \
          $(mod_param_sgtype2) \
          $(mod_basis_l_to_n) \
          $(mod_param_rd) \
          $(mod_mpi)
#file+mod_name: SRC/sub_Basis/sub_module_param_SGType2.f90 mod_param_sgtype2
$(OBJ_DIR)/sub_module_param_SGType2.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_Basis/sub_quadra_DirProd.f90 
$(OBJ_DIR)/sub_quadra_DirProd.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_dnsvm)
#file+mod_name: SRC/sub_Basis/sub_quadra_SincDVR.f90 
$(OBJ_DIR)/sub_quadra_SincDVR.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_SparseBasis.f90 
$(OBJ_DIR)/sub_quadra_SparseBasis.o : \
          $(evr_system_m) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_module_dind)
#file+mod_name: SRC/sub_Basis/sub_quadra_SparseBasis2n.f90 
$(OBJ_DIR)/sub_quadra_SparseBasis2n.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_Wigner.f90 
$(OBJ_DIR)/sub_quadra_Wigner.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_Ylm.f90 
$(OBJ_DIR)/sub_quadra_Ylm.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_box.f90 
$(OBJ_DIR)/sub_quadra_box.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_dfst.f90 
$(OBJ_DIR)/sub_quadra_dfst.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_fourier.f90 
$(OBJ_DIR)/sub_quadra_fourier.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_herm.f90 
$(OBJ_DIR)/sub_quadra_herm.o : \
          $(evr_system_m) \
          $(mod_basis) \
          $(addnsvm_m) \
          $(basismakegrid) \
          $(mod_ndindex) \
          $(mod_dnsvm)
#file+mod_name: SRC/sub_Basis/sub_quadra_inact.f90 
$(OBJ_DIR)/sub_quadra_inact.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_laguerre.f90 
$(OBJ_DIR)/sub_quadra_laguerre.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_quadra_legendre.f90 
$(OBJ_DIR)/sub_quadra_legendre.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_Basis/sub_read_data.f90 readbasis_m
$(OBJ_DIR)/sub_read_data.o : \
          $(evr_system_m) \
          $(mod_basis) \
          $(mod_coord_keo) \
          $(mod_constant)
#file+mod_name: SRC/sub_CRP/sub_CRP.f90 mod_crp
$(OBJ_DIR)/sub_CRP.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_psi) \
          $(mod_realwithunit) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(mod_mpi)
#file+mod_name: SRC/sub_Operator/sub_MatOp.f90 mod_matop
$(OBJ_DIR)/sub_MatOp.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_setop) \
          $(mod_psi) \
          $(mod_oppsi) \
          $(mod_ndindex) \
          $(mod_ana_psi) \
          $(mod_basis_btog_gtob_sgtype4)
#file+mod_name: SRC/sub_Operator/sub_OpPsi.f90 mod_oppsi
$(OBJ_DIR)/sub_OpPsi.o : \
          $(mod_oppsi_sg4) \
          $(mod_oppsi_sg4_mpi) \
          $(evr_system_m) \
          $(mod_psi) \
          $(mod_setop) \
          $(mod_mpi) \
          $(mod_mpi_aux) \
          $(mod_symabelian) \
          $(mod_basis_btog_gtob) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_opgrid) \
          $(mod_primop) \
          $(mod_param_sgtype2) \
          $(mod_constant)
#file+mod_name: SRC/sub_Operator/sub_OpPsi_MPI.f90 mod_oppsi_mpi
$(OBJ_DIR)/sub_OpPsi_MPI.o : \
          $(evr_system_m) \
          $(mod_psi)
#file+mod_name: SRC/sub_Operator/sub_OpPsi_SG4.f90 mod_oppsi_sg4
$(OBJ_DIR)/sub_OpPsi_SG4.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis_set_alloc) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_psi) \
          $(mod_setop) \
          $(mod_basis) \
          $(mod_symabelian) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_primop)
#file+mod_name: SRC/sub_Operator/sub_OpPsi_SG4_MPI.f90 mod_oppsi_sg4_mpi
$(OBJ_DIR)/sub_OpPsi_SG4_MPI.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_symabelian) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_psi) \
          $(mod_setop) \
          $(mod_mpi_aux) \
          $(mod_oppsi_sg4) \
          $(mod_basis_btog_gtob_sgtype4_mpi) \
          $(mod_psi_set_alloc) \
          $(mod_mpi)
#file+mod_name: SRC/sub_Operator/sub_lib_Op.f90 
$(OBJ_DIR)/sub_lib_Op.o : \
          $(evr_system_m) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_setop)
#file+mod_name: SRC/sub_Operator/sub_module_Op.f90 mod_op
$(OBJ_DIR)/sub_module_Op.o : \
          $(mod_setop) \
          $(readop_m) \
          $(mod_matop) \
          $(mod_oppsi)
#file+mod_name: SRC/sub_Operator/sub_module_OpGrid.f90 mod_opgrid
$(OBJ_DIR)/sub_module_OpGrid.o : \
          $(evr_system_m) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi)
#file+mod_name: SRC/sub_Operator/sub_module_ReadOp.f90 readop_m
$(OBJ_DIR)/sub_module_ReadOp.o : \
          $(evr_system_m) \
          $(mod_opgrid) \
          $(mod_primop) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_cap) \
          $(mod_hstep)
#file+mod_name: SRC/sub_Operator/sub_module_SetOp.f90 mod_setop
$(OBJ_DIR)/sub_module_SetOp.o : \
          $(evr_system_m) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_opgrid) \
          $(readop_m) \
          $(mod_mpi) \
          $(mod_param_sgtype2) \
          $(mod_psi)
#file+mod_name: SRC/sub_Optimization/sub_main_Optimization.f90 
$(OBJ_DIR)/sub_main_Optimization.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(basismakegrid) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_auto_basis) \
          $(mod_optimization)
#file+mod_name: SRC/sub_Optimization/sub_module_BFGS.f90 mod_bfgs
$(OBJ_DIR)/sub_module_BFGS.o : \
          $(evr_system_m) \
          $(mod_dnsvm) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis)
#file+mod_name: SRC/sub_Optimization/sub_module_Optimization.f90 mod_optimization
$(OBJ_DIR)/sub_module_Optimization.o : \
          $(evr_system_m) \
          $(mod_simulatedannealing) \
          $(mod_bfgs) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(basismakegrid) \
          $(mod_dnsvm) \
          $(mod_primop) \
          $(mod_op) \
          $(mod_auto_basis)
#file+mod_name: SRC/sub_Optimization/sub_module_SimulatedAnnealing.f90 mod_simulatedannealing
$(OBJ_DIR)/sub_module_SimulatedAnnealing.o : \
          $(evr_system_m) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis)
#file+mod_name: SRC/sub_Smolyak_test/sub_Smolyak_DInd.f90 mod_smolyak_dind
$(OBJ_DIR)/sub_Smolyak_DInd.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_Smolyak_test/sub_Smolyak_RDP.f90 mod_smolyak_rdp
$(OBJ_DIR)/sub_Smolyak_RDP.o : \
          $(evr_system_m) \
          $(mod_smolyak_ba) \
          $(mod_smolyak_dind)
#file+mod_name: SRC/sub_Smolyak_test/sub_Smolyak_ba.f90 mod_smolyak_ba
$(OBJ_DIR)/sub_Smolyak_ba.o : \
          $(mod_smolyak_dind) \
          $(evr_system_m)
#file+mod_name: SRC/sub_Smolyak_test/sub_Smolyak_module.f90 mod_smolyak_test
$(OBJ_DIR)/sub_Smolyak_module.o : \
          $(mod_smolyak_dind) \
          $(mod_smolyak_rdp) \
          $(mod_smolyak_ba) \
          $(evr_system_m)
#file+mod_name: SRC/sub_Smolyak_test/sub_Smolyak_test.f90 
$(OBJ_DIR)/sub_Smolyak_test.o : \
          $(evr_system_m) \
          $(mod_smolyak_dind) \
          $(mod_smolyak_rdp) \
          $(mod_smolyak_ba) \
          $(mod_smolyak_test)
#file+mod_name: SRC/sub_WP/mod_psi.f90 mod_psi
$(OBJ_DIR)/mod_psi.o : \
          $(mod_param_wp0) \
          $(mod_type_ana_psi) \
          $(mod_psi_set_alloc) \
          $(mod_ana_psi) \
          $(mod_psi_io) \
          $(mod_psi_b_to_g) \
          $(mod_psi_op) \
          $(mod_ana_psi_mpi) \
          $(evr_system_m) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_wp0) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_WP/sub_module_WP0.f90 mod_wp0
$(OBJ_DIR)/sub_module_WP0.o : \
          $(evr_system_m) \
          $(mod_psi_set_alloc) \
          $(mod_ana_psi) \
          $(mod_psi_io) \
          $(mod_psi_op) \
          $(mod_param_wp0)
#file+mod_name: SRC/sub_WP/sub_module_ana_psi.f90 mod_ana_psi
$(OBJ_DIR)/sub_module_ana_psi.o : \
          $(evr_system_m) \
          $(mod_psi_set_alloc) \
          $(mod_type_ana_psi) \
          $(mod_constant) \
          $(mod_psi_io) \
          $(mod_psi_b_to_g) \
          $(mod_dnsvm) \
          $(mod_basis) \
          $(mod_ndindex) \
          $(mod_param_rd) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_WP/sub_module_ana_psi_MPI.f90 mod_ana_psi_mpi
$(OBJ_DIR)/sub_module_ana_psi_MPI.o : \
          $(evr_system_m) \
          $(mod_psi_set_alloc) \
          $(mod_basis) \
          $(mod_param_sgtype2) \
          $(mod_ana_psi) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_WP/sub_module_param_WP0.f90 mod_param_wp0
$(OBJ_DIR)/sub_module_param_WP0.o : \
          $(evr_system_m) \
          $(mod_coord_keo)
#file+mod_name: SRC/sub_WP/sub_module_psi_B_TO_G.f90 mod_psi_b_to_g
$(OBJ_DIR)/sub_module_psi_B_TO_G.o : \
          $(mod_basis) \
          $(evr_system_m) \
          $(mod_psi_set_alloc)
#file+mod_name: SRC/sub_WP/sub_module_psi_Op.f90 mod_psi_op
$(OBJ_DIR)/sub_module_psi_Op.o : \
          $(mod_basis) \
          $(evr_system_m) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_ana_psi)
#file+mod_name: SRC/sub_WP/sub_module_psi_Op_MPI.f90 mod_psi_op_mpi
$(OBJ_DIR)/sub_module_psi_Op_MPI.o : \
          $(mod_basis) \
          $(evr_system_m) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_param_sgtype2) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_psi_op)
#file+mod_name: SRC/sub_WP/sub_module_psi_io.f90 mod_psi_io
$(OBJ_DIR)/sub_module_psi_io.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_psi_set_alloc)
#file+mod_name: SRC/sub_WP/sub_module_psi_set_alloc.f90 mod_psi_set_alloc
$(OBJ_DIR)/sub_module_psi_set_alloc.o : \
          $(evr_system_m) \
          $(mod_basis) \
          $(mod_type_ana_psi) \
          $(mod_mpi_aux) \
          $(mod_mpi)
#file+mod_name: SRC/sub_WP/sub_module_type_ana_psi.f90 mod_type_ana_psi
$(OBJ_DIR)/sub_module_type_ana_psi.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_active/sub_Grid_SG4.f90 
$(OBJ_DIR)/sub_Grid_SG4.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_op)
#file+mod_name: SRC/sub_active/sub_diago_H.f90 
$(OBJ_DIR)/sub_diago_H.o : \
          $(evr_system_m) \
          $(mod_constant)
#file+mod_name: SRC/sub_active/sub_ini_act_harm.f90 
$(OBJ_DIR)/sub_ini_act_harm.o : \
          $(evr_system_m) \
          $(mod_op) \
          $(mod_primop)
#file+mod_name: SRC/sub_active/sub_lib_act.f90 
$(OBJ_DIR)/sub_lib_act.o : \
          $(evr_system_m) \
          $(mod_primop) \
          $(mod_op) \
          $(mod_coord_keo) \
          $(mod_basis)
#file+mod_name: SRC/sub_active/sub_paraRPH.f90 mod_set_pararph
$(OBJ_DIR)/sub_paraRPH.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis)
#file+mod_name: SRC/sub_analysis/sub_NLO.f90 
$(OBJ_DIR)/sub_NLO.o : \
          $(evr_system_m) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_analysis)
#file+mod_name: SRC/sub_analysis/sub_VibRot.f90 
$(OBJ_DIR)/sub_VibRot.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_analysis)
#file+mod_name: SRC/sub_analysis/sub_analyse.f90 mod_fullanalysis
$(OBJ_DIR)/sub_analyse.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_param_rd) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_ana_psi) \
          $(mod_ndindex)
#file+mod_name: SRC/sub_analysis/sub_intensity.f90 
$(OBJ_DIR)/sub_intensity.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_analysis)
#file+mod_name: SRC/sub_analysis/sub_module_analysis.f90 mod_analysis
$(OBJ_DIR)/sub_module_analysis.o : \
          $(evr_system_m) \
          $(mod_psi) \
          $(mod_crp) \
          $(mod_constant) \
          $(mod_coord_keo)
#file+mod_name: SRC/sub_data_initialisation/ini_data.f90 
$(OBJ_DIR)/ini_data.o : \
          $(evr_system_m) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_cap) \
          $(mod_basis) \
          $(readbasis_m) \
          $(mod_set_pararph) \
          $(readop_m) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_propa) \
          $(mod_auto_basis) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_data_initialisation/nb_harm.f90 
$(OBJ_DIR)/nb_harm.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis)
#file+mod_name: SRC/sub_data_initialisation/sub_namelist.f90 
$(OBJ_DIR)/sub_namelist.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op)
#file+mod_name: SRC/sub_inactive/sub_HST_harm.f90 
$(OBJ_DIR)/sub_HST_harm.o : \
          $(evr_system_m) \
          $(mod_dnsvm) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_primop) \
          $(mod_constant) \
          $(mod_ndindex)
#file+mod_name: SRC/sub_inactive/sub_ana_HS.f90 
$(OBJ_DIR)/sub_ana_HS.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_inactive/sub_changement_de_var.f90 
$(OBJ_DIR)/sub_changement_de_var.o : \
          $(evr_system_m) \
          $(mod_basis)
#file+mod_name: SRC/sub_inactive/sub_inactive_harmo.f90 
$(OBJ_DIR)/sub_inactive_harmo.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_primop) \
          $(mod_basis)
#file+mod_name: SRC/sub_main/EVR_Module.f90 mod_evr
$(OBJ_DIR)/EVR_Module.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_op) \
          $(mod_analysis)
#file+mod_name: SRC/sub_main/EVR_driver.f90 
$(OBJ_DIR)/EVR_driver.o : \
          $(mod_evr) \
          $(mod_fullpropa) \
          $(mod_fullcontrol) \
          $(mod_davidson) \
          $(mod_filter) \
          $(mod_arpack) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_fullanalysis) \
          $(mod_auto_basis) \
          $(mod_psi) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_main/cart.f90 
$(OBJ_DIR)/cart.o : \
          $(evr_system_m) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_cart)
#file+mod_name: SRC/sub_main/sub_main_nDfit.f90 mod_ndgridfit
$(OBJ_DIR)/sub_main_nDfit.o : \
          $(evr_system_m) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(qdutil_intvec_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop)
#file+mod_name: SRC/sub_main/vib.f90 
$(OBJ_DIR)/vib.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_fullcontrol) \
          $(mod_davidson) \
          $(mod_filter) \
          $(mod_arpack) \
          $(mod_crp) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_fullanalysis) \
          $(mod_auto_basis) \
          $(mod_mpi_aux) \
          $(mod_optimization) \
          $(mod_wp0)
#file+mod_name: SRC/sub_propagation/sub_Hmax.f90 
$(OBJ_DIR)/sub_Hmax.o : \
          $(evr_system_m) \
          $(mod_op) \
          $(mod_psi) \
          $(mod_ana_psi_mpi) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_davidson) \
          $(mod_hmax_mpi) \
          $(mod_mpi_aux) \
          $(mod_march)
#file+mod_name: SRC/sub_propagation/sub_Hmax_MPI.f90 mod_hmax_mpi
$(OBJ_DIR)/sub_Hmax_MPI.o : \
          $(evr_system_m) \
          $(mod_setop) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_propagation/sub_TF_autocorr.f90 
$(OBJ_DIR)/sub_TF_autocorr.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_propa)
#file+mod_name: SRC/sub_propagation/sub_control.f90 mod_fullcontrol
$(OBJ_DIR)/sub_control.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_field) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_wp0) \
          $(mod_march)
#file+mod_name: SRC/sub_propagation/sub_module_Arpack.f90 mod_arpack
$(OBJ_DIR)/sub_module_Arpack.o : \
          $(mod_constant) \
          $(evr_system_m) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_mpi_aux) \
          $(mod_wp0)
#file+mod_name: SRC/sub_propagation/sub_module_Davidson.f90 mod_davidson
$(OBJ_DIR)/sub_module_Davidson.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(evr_system_m) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_propa_mpi) \
          $(mod_davidson_mpi) \
          $(mod_ana_psi_mpi) \
          $(mod_psi_op_mpi) \
          $(mod_mpi_aux) \
          $(mod_wp0) \
          $(mod_basis)
#file+mod_name: SRC/sub_propagation/sub_module_Davidson_MPI.f90 mod_davidson_mpi
$(OBJ_DIR)/sub_module_Davidson_MPI.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(evr_system_m) \
          $(mod_ana_psi_mpi) \
          $(mod_psi) \
          $(mod_psi_op_mpi) \
          $(mod_propa) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_propagation/sub_module_ExactFact.f90 mod_exactfact
$(OBJ_DIR)/sub_module_ExactFact.o : \
          $(evr_system_m) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_field)
#file+mod_name: SRC/sub_propagation/sub_module_Filter.f90 mod_filter
$(OBJ_DIR)/sub_module_Filter.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa)
#file+mod_name: SRC/sub_propagation/sub_module_field.f90 mod_field
$(OBJ_DIR)/sub_module_field.o : \
          $(evr_system_m)
#file+mod_name: SRC/sub_propagation/sub_module_propa_march.f90 mod_march
$(OBJ_DIR)/sub_module_propa_march.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_field) \
          $(mod_march_mpi) \
          $(mod_march_sg4) \
          $(mod_op) \
          $(mod_mpi_aux) \
          $(mod_basis) \
          $(mod_oppsi_sg4)
#file+mod_name: SRC/sub_propagation/sub_module_propa_march_MPI.f90 mod_march_mpi
$(OBJ_DIR)/sub_module_propa_march_MPI.o : \
          $(evr_system_m) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_march_sg4) \
          $(mod_op) \
          $(mod_psi_op_mpi) \
          $(mod_psi_set_alloc) \
          $(mod_oppsi_mpi) \
          $(mod_propa_mpi) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_basis_btog_gtob_sgtype4_mpi) \
          $(mod_mpi_aux) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_ndindex) \
          $(mod_symabelian) \
          $(mod_psi_op)
#file+mod_name: SRC/sub_propagation/sub_module_propa_march_SG4.f90 mod_march_sg4
$(OBJ_DIR)/sub_module_propa_march_SG4.o : \
          $(evr_system_m) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_op) \
          $(mod_oppsi_sg4) \
          $(mod_mpi_aux)
#file+mod_name: SRC/sub_propagation/sub_module_propagation.f90 mod_propa
$(OBJ_DIR)/sub_module_propagation.o : \
          $(evr_system_m) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_field) \
          $(mod_op) \
          $(mod_exactfact) \
          $(mod_mpi_aux) \
          $(mod_realwithunit) \
          $(mod_type_ana_psi) \
          $(mod_coord_keo)
#file+mod_name: SRC/sub_propagation/sub_module_propagation_MPI.f90 mod_propa_mpi
$(OBJ_DIR)/sub_module_propagation_MPI.o : \
          $(mod_propa) \
          $(mod_mpi_aux) \
          $(evr_system_m) \
          $(mod_op) \
          $(mod_psi_set_alloc) \
          $(mod_psi_op_mpi)
#file+mod_name: SRC/sub_propagation/sub_propagation.f90 mod_fullpropa
$(OBJ_DIR)/sub_propagation.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(evr_system_m) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_psi) \
          $(mod_field) \
          $(mod_march)
