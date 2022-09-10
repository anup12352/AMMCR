all:
	g++ main.cpp admixture_value.cpp admixture_value_p.cpp Aminus.cpp analytical_fitting.cpp Aplus.cpp beta.cpp betaminus.cpp betaplus.cpp conduction_dispersion.cpp components_BTE.cpp decompose_wave.cpp decompose_wave_p.cpp dedk.cpp df0dk.cpp df0dz.cpp DOS_value.cpp fitting_band.cpp f0.cpp find_cbm_vbm.cpp find_fermi.cpp FindMinInd.cpp J.cpp generate_required_data.cpp kminus.cpp kplus.cpp lambda_e_minus.cpp lambda_e_plus.cpp lambda_i_minus.cpp lambda_i_plus.cpp lambda_o_minus.cpp lambda_o_plus.cpp linspace.cpp make_band.cpp mu_elastic.cpp mu_overall.cpp mu_po.cpp mu_elasticH.cpp mu_overallH.cpp mu_poH.cpp n_DOS.cpp N_poph.cpp nu_alloy1.cpp nu_de.cpp nu_dis.cpp nu_ii.cpp nu_im.cpp nu_pe.cpp nu_npop_n.cpp nu_iv_n.cpp nu_ir.cpp nu_de_2D.cpp nu_npop_2D.cpp polyfit1.cpp polyval.cpp read_fermi.cpp read_input_file.cpp read_OUTCAR.cpp read_procar.cpp read_ispin.cpp rsquare.cpp copyright.cpp save_scattering_rate.cpp save_perturbation.cpp solve_g.cpp calculate_mobility.cpp save_results.cpp initialize_array.cpp conductivity_time.cpp mu_overall_time.cpp nu_pop.cpp nu_ii_p_funct.cpp nu_pop_p_funct.cpp nu_npop_p_funct.cpp nu_de_p_funct.cpp tight_binding_band_structure.cpp tight_binding_band_DOS.cpp polarizability.cpp nu_pop_2D.cpp nu_so_pop_2D.cpp  nu_rim_2D.cpp sigma_tensor.cpp sigma_freq.cpp save_scattering_rate_contribution.cpp nu_skew.cpp -o AMMCR

clean:
	rm -f *.o mobility.dat conductivity.dat themopower.dat scattering_rate.dat g.dat  conduction_band.dat valence_band.dat mobility_up_spin.dat conductivity_up_spin.dat thermopower_up_spin.dat scattering_rate_up_spin.dat g_up_spin.dat  conduction_band_up_spin.dat valence_band_up_spin.dat  mobility_down_spin.dat conductivity_down_spin.dat thermopower_down_spin.dat scattering_rate_down_spin.dat g_down_spin.dat  conduction_band_down_spin.dat valence_band_down_spin.dat EIGENVAL OUTCAR PROCAR DOSCAR *.dat *.o  -f *.txt *.o	

clean1: 
	rm OUTCAR EIGENVAL PROCAR DOSCAR
		
	
