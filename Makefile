all: prepare_profiles_2D.cpp prepare_profiles.cpp prepare_profiles_no_corr.cpp libtc93.a
		g++ prepare_profiles.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles.out -Wall
		g++ prepare_profiles_no_corr.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles_no_corr.out -Wall
		g++ prepare_profiles_2D.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles_2D.out -Wall
		g++ prepare_profiles_yri.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles_yri.out -Wall
analysis:
		g++ analyse_model_1D_maxw_.cpp -c
		g++ analyse_model_2D_maxw_.cpp -c
		g++ analyse_model_1D_norm_.cpp -c   -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ analyse_model_2D_norm_.cpp -c   -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ analyse_model_1D_pach_.cpp -c
		g++ analyse_model_2D_pach_.cpp -c
		g++ analyse_model_1D_sum_maxw_.cpp -c
		g++ analyse_model_2D_sum_maxw_.cpp -c
		g++ analyse_model_1D_uni_.cpp -c
		g++ analyse_model_2D_uni_.cpp -c
		gfortran -c pikaia.f
		gfortran -c pikaia_1D.f
		gfortran analyse_model_1D_maxw_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_maxw.out
		gfortran analyse_model_2D_maxw_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_maxw.out
		gfortran analyse_model_1D_norm_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_norm.out  -lstdc++ -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas 
		gfortran analyse_model_2D_norm_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_norm.out  -lstdc++ -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		gfortran analyse_model_1D_pach_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_pach.out
		gfortran analyse_model_2D_pach_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_pach.out
		gfortran analyse_model_1D_sum_maxw_.o pikaia.o -lc++ -lstdc++ -o analyse_model_1D_sum_maxw.out
		gfortran analyse_model_2D_sum_maxw_.o pikaia.o -lc++ -lstdc++ -o analyse_model_2D_sum_maxw.out
		gfortran analyse_model_1D_uni_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_uni.out
		gfortran analyse_model_2D_uni_.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_uni.out
profiles: prepare_profiles.cpp
		g++ prepare_profiles.cpp -c
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles.o -o prepare_profiles_1D.out 
		g++ prepare_profiles_joint.cpp -c
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles_joint.o -o prepare_profiles_1D_joint.out 
		g++ prepare_profiles_2D.cpp -c 
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles_2D.o -o prepare_profiles_2D.out 
		g++ prepare_profiles_yri.cpp -c
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles_yri.o -o prepare_profiles_yri.out 
		g++ prepare_profiles_yri_joint.cpp -c
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles_yri_joint.o -o prepare_profiles_yri_joint.out 
		g++ prepare_profiles_2D_yri.cpp -c 
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles_2D_yri.o -o prepare_profiles_2D_yri.out 
		g++ prepare_profiles_sch.cpp -c
		g++ -L./develop/tc93thick/ -lTC93thick -L/opt/local/lib/ -lf95 -lm prepare_profiles_sch.o -o prepare_profiles_1D_sch.out 
		g++ prepare_profiles_joint_sch.cpp -c
		g++ -L./develop/tc93thick/ -lTC93thick -L/opt/local/lib/ -lf95 -lm prepare_profiles_joint_sch.o -o prepare_profiles_1D_joint_sch.out 
		g++ prepare_profiles_2D_sch.cpp -c 
		g++ -L./develop/tc93thick/ -lTC93thick -L/opt/local/lib/ -lf95 -lm prepare_profiles_2D_sch.o -o prepare_profiles_2D_sch.out 
		g++ prepare_profiles_yri_sch.cpp -c
		g++ -L./develop/tc93thick/ -lTC93thick -L/opt/local/lib/ -lf95 -lm prepare_profiles_yri_sch.o -o prepare_profiles_yri_sch.out 
		g++ prepare_profiles_yri_joint_sch.cpp -c
		g++ -L./develop/tc93thick/ -lTC93thick -L/opt/local/lib/ -lf95 -lm prepare_profiles_yri_joint_sch.o -o prepare_profiles_yri_joint_sch.out 
		g++ prepare_profiles_2D_yri_sch.cpp -c 
		g++ -L./develop/tc93thick/ -lTC93thick -L/opt/local/lib/ -lf95 -lm prepare_profiles_2D_yri_sch.o -o prepare_profiles_2D_yri_sch.out 


evidence: 
		g++ evidence.cpp analyse_model_1D_maxw.cpp fun_apriory_energy_maxw.cpp -o evidence_1D_maxw_energy.out
		g++ evidence.cpp analyse_model_1D_maxw.cpp fun_apriory_flat.cpp        -o evidence_1D_maxw_flat.out
		g++ evidence.cpp analyse_model_2D_maxw.cpp fun_apriory_energy_maxw.cpp -o evidence_2D_maxw_energy.out
		g++ evidence.cpp analyse_model_2D_maxw.cpp fun_apriory_flat.cpp        -o evidence_2D_maxw_flat.out
		g++ evidence.cpp analyse_model_1D_norm.cpp fun_apriory_energy_norm.cpp -o evidence_1D_norm_energy.out  -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ evidence.cpp analyse_model_1D_norm.cpp fun_apriory_flat.cpp        -o evidence_1D_norm_flat.out    -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ evidence.cpp analyse_model_2D_norm.cpp fun_apriory_energy_norm.cpp -o evidence_2D_norm_energy.out   -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ evidence.cpp analyse_model_2D_norm.cpp fun_apriory_flat.cpp        -o evidence_2D_norm_flat.out  -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ evidence.cpp analyse_model_1D_pach.cpp fun_apriory_energy_pach.cpp -o evidence_1D_pach_energy.out
		g++ evidence.cpp analyse_model_1D_pach.cpp fun_apriory_flat.cpp        -o evidence_1D_pach_flat.out
		g++ evidence.cpp analyse_model_2D_pach.cpp fun_apriory_energy_pach.cpp -o evidence_2D_pach_energy.out
		g++ evidence.cpp analyse_model_2D_pach.cpp fun_apriory_flat.cpp        -o evidence_2D_pach_flat.out
		g++ evidence.cpp analyse_model_1D_uni.cpp  fun_apriory_energy_uni.cpp  -o evidence_1D_uni_energy.out
		g++ evidence.cpp analyse_model_1D_uni.cpp  fun_apriory_flat.cpp        -o evidence_1D_uni_flat.out
		g++ evidence.cpp analyse_model_2D_uni.cpp  fun_apriory_energy_uni.cpp  -o evidence_2D_uni_energy.out
		g++ evidence.cpp analyse_model_2D_uni.cpp  fun_apriory_flat.cpp        -o evidence_2D_uni_flat.out
		g++ evidence_sum_maxw.cpp fun_apriory_flat_sum.cpp   analyse_model_1D_sum_maxw.cpp   -o evidence_1D_sum_maxw_flat.out
		g++ evidence_sum_maxw.cpp fun_apriory_energy_sum.cpp   analyse_model_1D_sum_maxw.cpp -o evidence_1D_sum_maxw_energy.out
		g++ evidence_sum_maxw.cpp fun_apriory_flat_sum.cpp analyse_model_2D_sum_maxw.cpp     -o evidence_2D_sum_maxw_flat.out
		g++ evidence_sum_maxw.cpp fun_apriory_energy_sum.cpp analyse_model_2D_sum_maxw.cpp -o evidence_2D_sum_maxw_energy.out
credential:
		g++ cred_intervals.cpp analyse_model_1D_maxw.cpp -o credential_1D_maxw.out 
		g++ cred_intervals.cpp analyse_model_2D_maxw.cpp -o credential_2D_maxw.out 
		g++ cred_intervals.cpp analyse_model_1D_norm.cpp -o credential_1D_norm.out   -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas  
		g++ cred_intervals.cpp analyse_model_2D_norm.cpp -o credential_2D_norm.out   -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ cred_intervals.cpp analyse_model_1D_pach.cpp -o credential_1D_pach.out 
		g++ cred_intervals.cpp analyse_model_2D_pach.cpp -o credential_2D_pach.out 
		g++ cred_intervals_sum_maxw.cpp analyse_model_1D_sum_maxw.cpp -o credential_1D_sum_maxw.out
		g++ cred_intervals_sum_maxw.cpp analyse_model_2D_sum_maxw.cpp -o credential_2D_sum_maxw.out









