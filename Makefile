all: prepare_profiles_2D.cpp prepare_profiles.cpp prepare_profiles_no_corr.cpp libtc93.a analyse_model.cpp
		g++ prepare_profiles.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles.out -Wall
		g++ prepare_profiles_no_corr.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles_no_corr.out -Wall
		g++ prepare_profiles_2D.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles_2D.out -Wall
		g++ analyse_model.cpp -c
		g++ analyse_model_pach.cpp -c
		gfortran -c pikaia.f
		gfortran analyse_model.o pikaia.o -lc++ -lstdc++ -o analyse_model.out 
		gfortran -c pikaia_1D.f
		gfortran analyse_model_pach.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_pach.out 
		rm profiles/profile*.dat
analysis:
		g++ analyse_model_1D_maxw.cpp -c
		g++ analyse_model_2D_maxw.cpp -c
		g++ analyse_model_1D_norm.cpp -c   -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ analyse_model_2D_norm.cpp -c   -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		g++ analyse_model_1D_pach.cpp -c
		g++ analyse_model_2D_pach.cpp -c
		g++ analyse_model_1D_sum_maxw.cpp -c
		g++ analyse_model_2D_sum_maxw.cpp -c
		g++ analyse_model_1D_uni.cpp -c
		g++ analyse_model_2D_uni.cpp -c
		gfortran -c pikaia.f
		gfortran -c pikaia_1D.f
		gfortran analyse_model_1D_maxw.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_maxw.out
		gfortran analyse_model_2D_maxw.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_maxw.out
		gfortran analyse_model_1D_norm.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_norm.out  -lstdc++ -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas 
		gfortran analyse_model_2D_norm.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_norm.out  -lstdc++ -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		gfortran analyse_model_1D_pach.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_pach.out
		gfortran analyse_model_2D_pach.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_pach.out
		gfortran analyse_model_1D_sum_maxw.o pikaia.o -lc++ -lstdc++ -o analyse_model_1D_sum_maxw.out
		gfortran analyse_model_2D_sum_maxw.o pikaia.o -lc++ -lstdc++ -o analyse_model_2D_sum_maxw.out
		gfortran analyse_model_1D_uni.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_1D_uni.out
		gfortran analyse_model_2D_uni.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_2D_uni.out
profiles: prepare_profiles.cpp
		g++ prepare_profiles.cpp -c
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles.o -o prepare_profiles_1D.out 
		g++ prepare_profiles_joint.cpp -c
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles_joint.o -o prepare_profiles_1D_joint.out 
		g++ prepare_profiles_2D.cpp -c 
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles_2D.o -o prepare_profiles_2D.out 
	
