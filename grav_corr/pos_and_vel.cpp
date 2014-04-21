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

double t_generator ();
void vel_generator (double * );
void pos_generator (double * );
double move (double, double *, double *);
double v_proj (double * , double * , double * );
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


t = atof(argc[1]);

v_init = atof(argc[2]);

cout << "Time is -- " << t <<", v_init -- "<<v_init << endl;

		pos[0] = 6.;
		pos[1] = 0;
		pos[2] = 0;
//		vel_generator(&vel[0]);
//		t = t_generator ();		
		vel[0] = 0;
		vel[1] = 225 *  1e5/lsec*lcm; 
		vel[2] = v_init *1e5/lsec*lcm;

	for (int i=1; i < 1000; i++)	{
//		pos_generator(&pos[0]);
//		pos[0] = 1. + i * 0.15;
///		pos[1] = 0;
//		pos[2] = 0;
//		vel_generator(&vel[0]);
//		t = t_generator ();		
//		vel[0] = 0;
//		vel[1] = 225 *  1e5/lsec*lcm; 
//		vel[2] = v_init *1e5/lsec*lcm;

		//t = 1e7;
		t = 1e5;

		data[0][i] = pos[0];
		data[1][i] = pos[1];
		data[2][i] = pos[2];
		data[3][i] = vel[0];
		data[4][i] = vel[1];
		data[5][i] = vel[2];
		move(t, &pos[0], &vel[0]);
//		v_proj(&vel[0], &pos[0], &proj[0]);
//		data[6][i] = proj[0];
//		data[7][i] = proj[1];
//		data[8][i] = proj[2];
	
//		cout <<"Original line -- "<< data[5][i]*lsec/lcm/1.e5<< "\t" << vel [2]*lsec/lcm/1.e5<<"\t" << pos[2] << "\t" << t <<endl;
//		cout <<"Second line   -- "<<pow(data[5][i], 2) << "\t" << pow(vel[2],2) <<"\t" << 2*phi (pos[0], pos[1], pos[2]) - 2*phi (pos[0], pos[1], 0) << endl;
//		cout <<"Third line    -- "<<  2*phi (pos[0], pos[1], pos[2]) - 2*phi (pos[0], pos[1], 0) + pow(vel[2], 2) - pow(data[5][i], 2) <<endl;
//		cout <<"Fourth line   -- "<< sqrt(6000*abs(pos[2])-2000*log(1+3.0*abs(pos[2]))) << "\t" << (vel[2] - data[5][i])*lsec/lcm/1e5 << endl;
		
		out_pos << 1e5*i << "\t" << pos[2] << endl; 
		out_vel << 1e5*i << "\t" << vel[2]*lsec/lcm/1e5 << endl;

		//out << data[0][i] << "\t" <<  abs((vel[2] - data[5][i]))*lsec/lcm/1e5 << endl;
	}



return 0;
}

double t_generator () {
double res;

		res = rand() / rand_high_board;
		res *=  5e6;
		res += 1e3;

return res;
}
	

void pos_generator (double * pos) {
double r, real_prob, prob, phi;

	do {
		r = rand() / rand_high_board;
		r *=  30;
		r -= 15;
		prob = rand() / rand_high_board;
		real_prob = rho(r);
	} while (real_prob < prob);

	phi = rand() / rand_high_board;
	phi *= 2*pi;
	
	pos[0] = r*sin(phi);
	pos[1] = r*cos(phi);
	pos[2] = 0;

}


void vel_generator (double * vel) {

double v, phi, theta, real_prob, prob;

	v = 450e5 * norm_distr() / lsec * lcm;
	
	phi = rand() / rand_high_board;
	phi *= 2*pi;

	do {
		theta = rand() / rand_high_board;
		theta *=  pi;
		theta -= pi/2.;
		prob = rand() / rand_high_board;
		real_prob = cos(theta);
	} while (real_prob < prob);

//	vel [0] = v * cos(phi)*cos(theta);
//	vel [1] = v * sin(phi)*cos(theta);
	vel [0] = 0;
	vel [1] = 0;
	vel [2] = v * sin(theta);
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

double v_proj (double * vel, double * pos, double * diff_proj) {
double scal_prod, length_PS;
double PS[3];
double length_before, length_after;

	PS[0] = pos[0] - 8.;
	PS[1] = pos[1];
	PS[2] = pos[2];

	length_PS = sqrt(PS[0]*PS[0] + PS[1]*PS[1] + PS[2]*PS[2]);

	length_before = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

	scal_prod = (vel[0]*PS[0] + vel[1]*PS[1] + vel[2]*PS[2]) / length_PS;

	vel[0] -= scal_prod * PS[0] / length_PS;
	vel[1] -= scal_prod * PS[1] / length_PS;
	vel[2] -= scal_prod * PS[2] / length_PS;

	diff_proj[0] = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
	diff_proj[1] = sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
	diff_proj[2] = vel[2];

//	cout<<length_before *lsec/lcm/1e5 <<"\t"<<diff_proj[0] *lsec/lcm/1e5 <<endl;

}



