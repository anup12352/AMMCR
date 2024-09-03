########################################################################
# The CFLAGS variable sets compile flags for gcc:
#  -g        compile with debug information
#  -Wall     give verbose compiler warnings
#  -O0       do not optimize generated code
#  -O3       
#  -std=c++14  use the C++14 
##########################################################################

# -MMD: Generates a dependency file alongside each object file.
# -MP: Adds phony targets for dependencies to avoid errors if files are deleted.
# DEPS = $(CPPS:.cpp=.d): This variable lists all the .d files corresponding to the .cpp files.
# -include $(DEPS): This directive tells make to include the dependency files. 
# The - ensures make ignores errors if the dependency files do not exist yet.

# %.o: %.cpp: A pattern rule that tells make how to build .o files from .cpp files.

# With this setup, make will regenerate the dependency files whenever the source files are compiled. 
# If a header file changes, the corresponding .d file will indicate which .cpp files need to be recompiled.

EXE = AMMCR

CPPS = main.cpp admixture_value.cpp admixture_value_p.cpp Aminus.cpp analytical_fitting.cpp Aplus.cpp beta.cpp betaminus.cpp \
	betaplus.cpp conduction_dispersion.cpp components_BTE.cpp decompose_wave.cpp decompose_wave_p.cpp dedk.cpp df0dk.cpp \
	df0dz.cpp DOS_value.cpp fitting_band.cpp f0.cpp find_cbm_vbm.cpp find_fermi.cpp FindMinInd.cpp J.cpp generate_required_data.cpp \
	kminus.cpp kplus.cpp lambda_e_minus.cpp lambda_e_plus.cpp lambda_i_minus.cpp lambda_i_plus.cpp lambda_o_minus.cpp lambda_o_plus.cpp \
	linspace.cpp make_band.cpp mu_elastic.cpp mu_overall.cpp mu_po.cpp mu_elasticH.cpp mu_overallH.cpp mu_poH.cpp n_DOS.cpp N_poph.cpp \
	nu_alloy1.cpp nu_de.cpp nu_dis.cpp nu_ii.cpp nu_im.cpp nu_pe.cpp nu_npop_n.cpp nu_iv_n.cpp nu_ir.cpp nu_de_2D.cpp nu_npop_2D.cpp polyfit1.cpp \
	polyval.cpp read_fermi.cpp read_input_file.cpp read_OUTCAR.cpp read_procar.cpp read_ispin.cpp rsquare.cpp copyright.cpp \
	save_scattering_rate.cpp save_perturbation.cpp solve_g.cpp calculate_mobility.cpp save_results.cpp initialize_array.cpp \
	conductivity_time.cpp mu_overall_time.cpp nu_pop.cpp nu_ii_p_funct.cpp nu_pop_p_funct.cpp nu_npop_p_funct.cpp nu_de_p_funct.cpp \
	tight_binding_band_structure.cpp tight_binding_band_DOS.cpp polarizability.cpp nu_pop_2D.cpp nu_so_pop_2D.cpp  nu_rim_2D.cpp \
	sigma_tensor.cpp sigma_freq.cpp save_scattering_rate_contribution.cpp nu_skew.cpp


OBJS =  $(CPPS:.cpp=.o)

DEPS = $(CPPS:.cpp=.d)

DebugFlag = -g
#DebugFlag = 

Compiler = g++ 

Flags = -std=c++14 -MMD -MP
#Flags =  -MMD -MP

#LIB = -I ../helper/eigen-3.4.0/eigen-3.4.0/

$(EXE): $(OBJS)
	$(Compiler) $(DebugFlag) $(Flags) -o $(EXE) $(OBJS) $(LIB)

-include $(DEPS)

%.o: %.cpp
	$(Compiler) $(DebugFlag) $(Flags) -c $< -o $@

clean_all_ubuntu:
	rm -f *.o *.d  AMMCR
	touch *

clean_all_window:
	del /Q *.o *.d AMMCR.exe
	@echo > NUL 

#	$(EXE)	
#clean:
#	$(EXE)	

#clean_time:		
#	touch *
	
#clean_obj:
#	rm -f *.o ../data_set/*.o ../model/*.o ../helper/*.o $(EXE)
