#===============================================
mod_auto_basis = $(OBJ_DIR)/sub_Auto_Basis.o
mod_basis_btog_gtob_sgtype4 = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4.o
mod_basis_btog_gtob_sgtype4_mpi = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4_MPI.o
mod_basis_rcvec_sgtype4 = $(OBJ_DIR)/sub_module_basis_RCVec_SG4.o
mod_param_rd = $(OBJ_DIR)/sub_module_param_RD.o
mod_symabelian = $(OBJ_DIR)/sub_SymAbelian.o
basismakegrid = $(OBJ_DIR)/sub_module_BasisMakeGrid.o
mod_basis_l_to_n = $(OBJ_DIR)/sub_module_Basis_LTO_n.o
mod_rotbasis_param = $(OBJ_DIR)/sub_module_RotBasis.o
mod_basis = $(OBJ_DIR)/sub_module_basis.o
mod_basis_btog_gtob = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB.o
mod_basis_btog_gtob_mpi = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_MPI.o
mod_basis_btog_gtob_sgtype2 = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SGType2.o
mod_basis_grid_param = $(OBJ_DIR)/sub_module_basis_Grid_Param.o
mod_basis_set_alloc = $(OBJ_DIR)/sub_module_basis_set_alloc.o
mod_param_sgtype2 = $(OBJ_DIR)/sub_module_param_SGType2.o
mod_crp = $(OBJ_DIR)/sub_CRP.o
mod_matop = $(OBJ_DIR)/sub_MatOp.o
mod_oppsi = $(OBJ_DIR)/sub_OpPsi.o
mod_oppsi_mpi = $(OBJ_DIR)/sub_OpPsi_MPI.o
mod_oppsi_sg4 = $(OBJ_DIR)/sub_OpPsi_SG4.o
mod_oppsi_sg4_mpi = $(OBJ_DIR)/sub_OpPsi_SG4_MPI.o
mod_op = $(OBJ_DIR)/sub_module_Op.o
mod_opgrid = $(OBJ_DIR)/sub_module_OpGrid.o
mod_readop = $(OBJ_DIR)/sub_module_ReadOp.o
mod_setop = $(OBJ_DIR)/sub_module_SetOp.o
mod_bfgs = $(OBJ_DIR)/sub_module_BFGS.o
mod_optimization = $(OBJ_DIR)/sub_module_Optimization.o
mod_simulatedannealing = $(OBJ_DIR)/sub_module_SimulatedAnnealing.o
mod_smolyak_dind = $(OBJ_DIR)/sub_Smolyak_DInd.o
mod_smolyak_rdp = $(OBJ_DIR)/sub_Smolyak_RDP.o
mod_smolyak_ba = $(OBJ_DIR)/sub_Smolyak_ba.o
mod_smolyak_test = $(OBJ_DIR)/sub_Smolyak_module.o
mod_psi = $(OBJ_DIR)/mod_psi.o
mod_wp0 = $(OBJ_DIR)/sub_module_WP0.o
mod_ana_psi = $(OBJ_DIR)/sub_module_ana_psi.o
mod_ana_psi_mpi = $(OBJ_DIR)/sub_module_ana_psi_MPI.o
mod_param_wp0 = $(OBJ_DIR)/sub_module_param_WP0.o
mod_psi_b_to_g = $(OBJ_DIR)/sub_module_psi_B_TO_G.o
mod_psi_op = $(OBJ_DIR)/sub_module_psi_Op.o
mod_psi_op_mpi = $(OBJ_DIR)/sub_module_psi_Op_MPI.o
mod_psi_io = $(OBJ_DIR)/sub_module_psi_io.o
mod_psi_set_alloc = $(OBJ_DIR)/sub_module_psi_set_alloc.o
mod_type_ana_psi = $(OBJ_DIR)/sub_module_type_ana_psi.o
mod_set_pararph = $(OBJ_DIR)/sub_paraRPH.o
mod_fullanalysis = $(OBJ_DIR)/sub_analyse.o
mod_analysis = $(OBJ_DIR)/sub_module_analysis.o
mod_evr = $(OBJ_DIR)/EVR_Module.o
mod_ndgridfit = $(OBJ_DIR)/sub_main_nDfit.o
mod_hmax_mpi = $(OBJ_DIR)/sub_Hmax_MPI.o
mod_fullcontrol = $(OBJ_DIR)/sub_control.o
mod_arpack = $(OBJ_DIR)/sub_module_Arpack.o
mod_davidson = $(OBJ_DIR)/sub_module_Davidson.o
mod_davidson_mpi = $(OBJ_DIR)/sub_module_Davidson_MPI.o
mod_exactfact = $(OBJ_DIR)/sub_module_ExactFact.o
mod_filter = $(OBJ_DIR)/sub_module_Filter.o
mod_field = $(OBJ_DIR)/sub_module_field.o
mod_march = $(OBJ_DIR)/sub_module_propa_march.o
mod_march_mpi = $(OBJ_DIR)/sub_module_propa_march_MPI.o
mod_march_sg4 = $(OBJ_DIR)/sub_module_propa_march_SG4.o
mod_propa = $(OBJ_DIR)/sub_module_propagation.o
mod_propa_mpi = $(OBJ_DIR)/sub_module_propagation_MPI.o
mod_fullpropa = $(OBJ_DIR)/sub_propagation.o
#===============================================
$(OBJ_DIR)/sub_Auto_Basis.o : \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_mpi) \
          $(mod_system) \
          $(basismakegrid) \
          $(mod_constant) \
          $(mod_analysis) \
          $(mod_propa) \
          $(mod_psi) \
          $(mod_fullanalysis)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux) \
          $(mod_basis_btog_gtob_sgtype4_mpi)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4_MPI.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_basis_RCVec_SG4.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_param_RD.o : \
          $(mod_system) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_SymAbelian.o : \
          $(mod_system)
$(OBJ_DIR)/sub_SymAbelian_OF_Basis.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis)
$(OBJ_DIR)/sub_basis_El.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_module_BasisMakeGrid.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_ndindex) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_module_Basis_LTO_n.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_RotBasis.o : \
          $(mod_system) \
          $(mod_ndindex)
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
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_basis_btog_gtob_sgtype4)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB.o : \
          $(mod_system) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype2) \
          $(mod_param_sgtype2) \
          $(mod_basis_btog_gtob_mpi) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_MPI.o : \
          $(mod_system) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SGType2.o : \
          $(mod_system) \
          $(mod_basis_set_alloc) \
          $(mod_module_dind)
$(OBJ_DIR)/sub_module_basis_Grid_Param.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_basis_set_alloc.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(mod_rotbasis_param) \
          $(mod_basis_grid_param) \
          $(mod_symabelian) \
          $(mod_param_sgtype2) \
          $(mod_basis_l_to_n) \
          $(mod_param_rd) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_param_SGType2.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_quadra_DirProd.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_quadra_SincDVR.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_SparseBasis.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_module_dind)
$(OBJ_DIR)/sub_quadra_SparseBasis2n.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_Wigner.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_Ylm.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_box.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_dfst.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_fourier.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_herm.o : \
          $(mod_system) \
          $(mod_basis) \
          $(addnsvm_m) \
          $(basismakegrid) \
          $(mod_ndindex) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_quadra_inact.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_laguerre.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_legendre.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_read_data.o : \
          $(mod_system) \
          $(mod_basis) \
          $(mod_coord_keo) \
          $(mod_constant)
$(OBJ_DIR)/sub_CRP.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_psi) \
          $(mod_realwithunit) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(mod_mpi)
$(OBJ_DIR)/sub_MatOp.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_setop) \
          $(mod_psi) \
          $(mod_oppsi) \
          $(mod_ndindex) \
          $(mod_ana_psi) \
          $(mod_basis_btog_gtob_sgtype4)
$(OBJ_DIR)/sub_OpPsi.o : \
          $(mod_oppsi_sg4) \
          $(mod_oppsi_sg4_mpi) \
          $(mod_system) \
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
          $(mod_param_sgtype2)
$(OBJ_DIR)/sub_OpPsi_MPI.o : \
          $(mod_system) \
          $(mod_psi)
$(OBJ_DIR)/sub_OpPsi_SG4.o : \
          $(mod_system) \
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
$(OBJ_DIR)/sub_OpPsi_SG4_MPI.o : \
          $(mod_system) \
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
$(OBJ_DIR)/sub_lib_Op.o : \
          $(mod_system) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_setop)
$(OBJ_DIR)/sub_module_Op.o : \
          $(mod_setop) \
          $(mod_readop) \
          $(mod_matop) \
          $(mod_oppsi)
$(OBJ_DIR)/sub_module_OpGrid.o : \
          $(mod_system) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_ReadOp.o : \
          $(mod_system) \
          $(mod_opgrid) \
          $(mod_primop)
$(OBJ_DIR)/sub_module_SetOp.o : \
          $(mod_system) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_opgrid) \
          $(mod_readop) \
          $(mod_mpi) \
          $(mod_param_sgtype2) \
          $(mod_psi)
$(OBJ_DIR)/sub_main_Optimization.o : \
          $(mod_system) \
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
$(OBJ_DIR)/sub_module_BFGS.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis)
$(OBJ_DIR)/sub_module_Optimization.o : \
          $(mod_system) \
          $(mod_simulatedannealing) \
          $(mod_bfgs) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(basismakegrid) \
          $(mod_dnsvm) \
          $(mod_primop) \
          $(mod_op) \
          $(mod_auto_basis)
$(OBJ_DIR)/sub_module_SimulatedAnnealing.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis)
$(OBJ_DIR)/sub_Smolyak_DInd.o : \
          $(mod_system)
$(OBJ_DIR)/sub_Smolyak_RDP.o : \
          $(mod_system) \
          $(mod_smolyak_ba) \
          $(mod_smolyak_dind)
$(OBJ_DIR)/sub_Smolyak_ba.o : \
          $(mod_smolyak_dind) \
          $(mod_system)
$(OBJ_DIR)/sub_Smolyak_module.o : \
          $(mod_smolyak_dind) \
          $(mod_smolyak_rdp) \
          $(mod_smolyak_ba) \
          $(mod_system)
$(OBJ_DIR)/sub_Smolyak_test.o : \
          $(mod_system) \
          $(mod_smolyak_dind) \
          $(mod_smolyak_rdp) \
          $(mod_smolyak_ba) \
          $(mod_smolyak_test)
$(OBJ_DIR)/mod_psi.o : \
          $(mod_param_wp0) \
          $(mod_type_ana_psi) \
          $(mod_psi_set_alloc) \
          $(mod_ana_psi) \
          $(mod_psi_io) \
          $(mod_psi_b_to_g) \
          $(mod_psi_op) \
          $(mod_ana_psi_mpi) \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_wp0) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_WP0.o : \
          $(mod_system) \
          $(mod_psi_set_alloc) \
          $(mod_ana_psi) \
          $(mod_psi_io) \
          $(mod_psi_op) \
          $(mod_param_wp0)
$(OBJ_DIR)/sub_module_ana_psi.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_type_ana_psi) \
          $(mod_basis) \
          $(mod_psi_set_alloc) \
          $(mod_psi_io) \
          $(mod_psi_b_to_g) \
          $(mod_dnsvm) \
          $(mod_param_sgtype2) \
          $(mod_param_rd) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_ana_psi_MPI.o : \
          $(mod_system) \
          $(mod_psi_set_alloc) \
          $(mod_basis) \
          $(mod_param_sgtype2) \
          $(mod_ana_psi) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_param_WP0.o : \
          $(mod_system) \
          $(mod_coord_keo)
$(OBJ_DIR)/sub_module_psi_B_TO_G.o : \
          $(mod_basis) \
          $(mod_system) \
          $(mod_psi_set_alloc)
$(OBJ_DIR)/sub_module_psi_Op.o : \
          $(mod_basis) \
          $(mod_system) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_ana_psi)
$(OBJ_DIR)/sub_module_psi_Op_MPI.o : \
          $(mod_basis) \
          $(mod_system) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_param_sgtype2) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_psi_op)
$(OBJ_DIR)/sub_module_psi_io.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_psi_set_alloc)
$(OBJ_DIR)/sub_module_psi_set_alloc.o : \
          $(mod_system) \
          $(mod_basis) \
          $(mod_type_ana_psi) \
          $(mod_mpi_aux) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_type_ana_psi.o : \
          $(mod_system)
$(OBJ_DIR)/sub_Grid_SG4.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_op)
$(OBJ_DIR)/sub_diago_H.o : \
          $(mod_system) \
          $(mod_constant)
$(OBJ_DIR)/sub_ini_act_harm.o : \
          $(mod_system) \
          $(mod_op) \
          $(mod_primop)
$(OBJ_DIR)/sub_lib_act.o : \
          $(mod_system) \
          $(mod_primop) \
          $(mod_op) \
          $(mod_coord_keo) \
          $(mod_basis)
$(OBJ_DIR)/sub_paraRPH.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_primop) \
          $(mod_basis)
$(OBJ_DIR)/sub_NLO.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_analysis)
$(OBJ_DIR)/sub_VibRot.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_analysis)
$(OBJ_DIR)/sub_analyse.o : \
          $(mod_constant) \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_param_rd) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_intensity.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_analysis)
$(OBJ_DIR)/sub_module_analysis.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_crp) \
          $(mod_coord_keo)
$(OBJ_DIR)/ini_data.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_cap) \
          $(mod_basis) \
          $(mod_set_pararph) \
          $(mod_readop) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_propa) \
          $(mod_auto_basis) \
          $(mod_mpi_aux)
$(OBJ_DIR)/nb_harm.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis)
$(OBJ_DIR)/sub_namelist.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_primop) \
          $(mod_cap) \
          $(mod_hstep)
$(OBJ_DIR)/sub_HST_harm.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_primop) \
          $(mod_constant) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_ana_HS.o : \
          $(mod_system)
$(OBJ_DIR)/sub_changement_de_var.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_inactive_harmo.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_primop) \
          $(mod_basis)
$(OBJ_DIR)/EVR-T.o : \
          $(mod_system) \
          $(mod_ndgridfit)
$(OBJ_DIR)/EVR_Module.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_op) \
          $(mod_analysis)
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
$(OBJ_DIR)/cart.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_cart)
$(OBJ_DIR)/sub_main_nDfit.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop)
$(OBJ_DIR)/vib.o : \
          $(mod_system) \
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
$(OBJ_DIR)/sub_Hmax.o : \
          $(mod_system) \
          $(mod_op) \
          $(mod_psi) \
          $(mod_ana_psi_mpi) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_davidson) \
          $(mod_hmax_mpi) \
          $(mod_mpi_aux) \
          $(mod_march)
$(OBJ_DIR)/sub_Hmax_MPI.o : \
          $(mod_system) \
          $(mod_setop) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_TF_autocorr.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_propa)
$(OBJ_DIR)/sub_control.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_field) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_wp0) \
          $(mod_march)
$(OBJ_DIR)/sub_module_Arpack.o : \
          $(mod_constant) \
          $(mod_system) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_mpi_aux) \
          $(mod_wp0)
$(OBJ_DIR)/sub_module_Davidson.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(mod_system) \
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
$(OBJ_DIR)/sub_module_Davidson_MPI.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(mod_system) \
          $(mod_ana_psi_mpi) \
          $(mod_psi) \
          $(mod_psi_op_mpi) \
          $(mod_propa) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_ExactFact.o : \
          $(mod_system) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_field)
$(OBJ_DIR)/sub_module_Filter.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa)
$(OBJ_DIR)/sub_module_field.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_propa_march.o : \
          $(mod_system) \
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
$(OBJ_DIR)/sub_module_propa_march_MPI.o : \
          $(mod_system) \
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
$(OBJ_DIR)/sub_module_propa_march_SG4.o : \
          $(mod_system) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_op) \
          $(mod_oppsi_sg4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_propagation.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_field) \
          $(mod_op) \
          $(mod_exactfact) \
          $(mod_mpi_aux) \
          $(mod_realwithunit) \
          $(mod_type_ana_psi) \
          $(mod_coord_keo)
$(OBJ_DIR)/sub_module_propagation_MPI.o : \
          $(mod_propa) \
          $(mod_mpi_aux) \
          $(mod_system) \
          $(mod_op) \
          $(mod_psi_set_alloc) \
          $(mod_psi_op_mpi)
$(OBJ_DIR)/sub_propagation.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(mod_system) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_psi) \
          $(mod_field) \
          $(mod_march)
