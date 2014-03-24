
all: prepare_profiles.cpp libtc93.a analyse_model.cpp
		g++ prepare_profiles.cpp -L. -ltc93 -L/opt/local/lib/ -lf95 -o prepare_profiles.out
		g++ analyse_model.cpp -o analyse_model.out
