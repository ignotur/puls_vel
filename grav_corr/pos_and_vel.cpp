#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
//#include <gsl/gsl_linalg.h>
#include <cstring>
//#include <gsl/gsl_multifit.h>

using namespace std;

double const lcm     = 3.2407789e-22;           // 1 см  в кпк,   для внутренних переводов
double const lsec    = 3.1688955e-8;            // 1 сек в годах, для внутренних переводов
double const pi = 3.1415926;
double const rand_high_board = 2.14665e+9;


double dphi_dz  (double, double, double);
void   diff_equi(int, double *);
void Runge_Kutta (int, double, double *, void (*f)(int, double *));
double pos (double, double*, double);
double norm_distr();
double rho (double);

double move (double, double *, double *);
double phi (double, double, double);

int main (int argv, char * argc[]) {

double v[1000], z[1000], vUnc[1000], vel_corr[1000];
double x, x_, v_z, v_z_prev, v_z_copy;
double counter;

//ifstream in  ("v_b.txt");
ofstream out_pos ("pos.txt");
ofstream out_vel ("vel.txt");

double data [9][1700], pos[3], vel[3], t, proj[3];
double v_init;


//t = atof(argc[1]);

v_init = atof(argc[1]);

cout << "v_init -- "<<v_init << endl;

		pos[0] = 6.;
		pos[1] = 0;
		pos[2] = 0;
//		vel_generator(&vel[0]);
//		t = t_generator ();		
		vel[0] = 0;
		vel[1] = 225 *  1e5/lsec*lcm; 
		vel[2] = v_init *1e5/lsec*lcm;

	for (int i=1; i < 1000; i++)	{
		t = 1e5;

		data[0][i] = pos[0];
		data[1][i] = pos[1];
		data[2][i] = pos[2];
		data[3][i] = vel[0];
		data[4][i] = vel[1];
		data[5][i] = vel[2];
		move(t, &pos[0], &vel[0]);
		
		out_pos << 1e5*i << "\t" << pos[2] << endl; 
		out_vel << 1e5*i << "\t" << vel[2]*lsec/lcm/1e5 << endl;

	}



return 0;
}

double move (double t, double *pos, double *vel) {
double res[6];

res[0] = pos[0];
res[1] = pos[1];
res[2] = pos[2];
res[3] = vel[0];
res[4] = vel[1];
res[5] = vel[2];

Runge_Kutta (6, t, &res[0], &diff_equi);

pos[0]  = res[0]; 
pos[1]  = res[1];
pos[2]  = res[2];
vel[0]  = res[3];
vel[1]  = res[4];
vel[2]  = res[5];

return 0;
}


