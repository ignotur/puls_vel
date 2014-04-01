#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

double norm_distr ();

int main (int argv, char * argc[]) {

srand(time(0));

ofstream out ("rand_vel.txt");

double const rand_high_board = 2.14665e+9;
double const pi = 3.1415926;

double v, phi, theta, real_prob, v_z, v_x, dist, parral, parral_err;
double prob, which_err_prmotion, prmotion_err, dv;
double dv1, dv2, dv3;
	for (int i=0; i <100; i++)	{
			v = atof(argc[1]) * norm_distr();
			phi = rand() / rand_high_board;
			phi *= 2*pi;

		do {
			theta = rand() / rand_high_board;
			theta *=  pi;
			theta -= pi/2.;
			prob = rand() / rand_high_board;
			real_prob = cos(theta);
		} while (real_prob < prob);

		//	v_z = v*sin(theta);
		//	v_x = v*cos(phi)*cos(theta);

		do {
			dist = rand() / rand_high_board;
			dist *= 5;
			dist += 0.01;
			prob = rand() / rand_high_board;
			real_prob = 1./dist/100.;
		} while (real_prob < prob);

		parral = 1./dist;

		if (parral > 1.)
			parral_err = 0.15*norm_distr()*parral;
		else
			parral_err = norm_distr()*0.1;			
	
		which_err_prmotion = rand() / rand_high_board;
		which_err_prmotion*= 5;

		if (which_err_prmotion <= 4)	
			prmotion_err = 0.6 * norm_distr ();
		else
			prmotion_err = 2 + 0.5 * norm_distr ();

		dv = prmotion_err * 0.001 * dist / 206265.*9.4e8 + parral_err/(parral_err + 1./dist) * v + prmotion_err*0.001 /206265 * parral_err/(parral_err + 1./dist)*9.4e8;
//		dv1 = prmotion_err * 0.001 * dist / 206265.*9.4e8;
//		dv2 = parral_err/(parral_err + 1./dist) * v;
//		dv3 = prmotion_err*0.001 /206265 * parral_err/(parral_err + 1./dist)*9.4e8;

		v_z = (v+dv)*sin(theta);
	
		out<<abs(v)<<endl;

//	cout<<"v_z -- "<<v_z<<", dist -- "<<dist<<", parral -- "<<parral<<endl;
//	cout<<"parral_err -- "<<parral_err<<", prmotion_err -- "<<prmotion_err<<endl;
//	cout<<"v -- "<<v<<", dv -- "<<dv<<endl;
//	cout<<"dv1 -- "<<dv1<<", dv2 -- "<<dv2<<", dv3 -- "<<dv3<<endl;

	}


return 0;
}
