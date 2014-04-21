
all: prepare_profiles.cpp prepare_profiles_no_corr.cpp libtc93.a analyse_model.cpp
		g++ prepare_profiles.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles.out -Wall
		g++ prepare_profiles_no_corr.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles_no_corr.out -Wall
		g++ analyse_model.cpp -c
		gfortran -c pikaia.f
		gfortran analyse_model.o pikaia.o -lc++ -lstdc++ -o analyse_model.out 
		rm profiles/profile*.dat
analysis:
		g++ analyse_model.cpp -c
		gfortran -c pikaia.f
		gfortran analyse_model.o pikaia.o -lc++ -lstdc++ -o analyse_model.out 
		
