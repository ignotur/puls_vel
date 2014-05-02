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
		g++ analyse_model_pach.cpp -c
		g++ analyse_model.cpp -c
		g++ analyse_model_igf.cpp -c -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas
		gfortran -c pikaia.f
		gfortran -c pikaia_1D.f
		gfortran analyse_model.o pikaia.o -lc++ -lstdc++ -o analyse_model.out 
		gfortran analyse_model_pach.o pikaia_1D.o -lc++ -lstdc++ -o analyse_model_pach.out 
		gfortran analyse_model_igf.o pikaia_1D.o -lc++ -lstdc++ -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -o analyse_model_igf.out
profiles: prepare_profiles.cpp
		g++ prepare_profiles.cpp -c
		g++ -L./develop/NE2001/ -lNE2001 -L/opt/local/lib/ -lf95 -lm prepare_profiles.o -o prepare_profiles.out 
		
