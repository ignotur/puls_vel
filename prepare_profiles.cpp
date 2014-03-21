#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <cstdlib>

using namespace std;

double const pi = 3.1415926;

extern "C" {
void dmdsm_ (double *l, double *b, int *ndir, double *dmpsr, double *dist, char *limit, double *sm, double *smtau, double *smtheta);
}

void profile  (double * dist, double * prmot, double * res);
void pdf_dist (double *, double);
void pdf_prmot(double, double, double);
double delta_vl (double *, double);
double prob_vl  (...);
double f(double *, double, double, double);

int main (int argv, char * argc[]) 	{

ifstream in_dist  (argc[1]);
ifstream in_prmot (argc[2]);

ofstream out_prof;

double dist[10][1000], prmot[6][1000], profile[1000];
double entry_prmot[3], entry_dist[3]; 
double trash;
int n_dist, n_prmot;
n_dist  = 0;
n_prmot = 0;

string name_str;
stringstream name;

			

// Read all data (distance and proper motion)

	do {
		for (int i=0; i < 10; i++)
			in_dist >> dist[i][n_dist];
			//--------------------------------------------------------
			// first value distance (dist[1] and dist[2] = -1)
			// or parralax (dist[1] and dist[2] = err(parral), 
			// then four -1, then luminosity, then l and b in degrees)
			//--------------------------------------------------------
		n_dist++;
	} while (!in_dist.eof());
n_dist--;

	do {
		for (int i=0; i < 6; i++)
			in_prmot >> dist[i][n_prmot];
			//--------------------------------------------------------
			// first value mu_l, then err_mu_l, err_mu_l
			// then mu_b, err_mu_b, err_mu_b
			//--------------------------------------------------------
		n_prmot++;
	} while (!in_prmot.eof());
n_prmot--;

	if (n_dist != n_prmot)	{
		cout<<"List sizes are incompatible!"<<endl;
		return 0;
	}
	
// The main loop. Calculate profiles for pulsars one by one.

	for (int i=0; i < n_dist; i++)	{

		for (int j=0; j < 10; j++)		
			entry_dist [j] = dist [j][i];
		for (int j=0; j < 6; j++)	
			entry_prmot[j] = prmot[j][i];	

		profile(&entry_dist[0], &entry_prmot[0], &res[0]);	 // Call a function to compute profile	
	
		// Open file with appropriate name
		name << "profiles/profile_"<<i<<".dat";
		name_str = name.str();
		char *basic_name = new char [name_str.size()];
		memcpy(basic_name, name_str.c_str(), name_str.size());
		name.clear();
		name.str(string());	
		
		cout << basic_name <<endl;
		out_prof.open(basic_name);

		for (int j=0; j < 1000; j++)
			out_prof << j+1 << "\t" << res[j] <<endl;	// profile writing 

		out_prof.close();

	}


return 0;
}

//-----------------------------------------------------------
// This function computes probability that lat. proper motion
// is mu_l. The parameters of the Gaussian distribution are
// mu_c - center and mu_s - sigma
//----------------------------------------------------------- 

double pdf_prmot (double mu_l, double mu_c, double mu_s)	{
double res;

res = 1./(mu_s*sqrt(pi*2)) * exp (-pow(mu_l - mu_c, 2)/(2.*pow(mu_s, 2.)));

return res;
}

//-----------------------------------------------------------
// This function computes probability that distance is
// D. The parameters of the Gaussian distribution are
// dist[0] - center and dist[1] - sigma if given,
// otherwise use TC93 model of electron density
//----------------------------------------------------------- 

double prf_dist (double * dist, double D)	{
double res, parallax;
double DM1, DM2;
double sm, smtau, smtheta;
double dist1, dist2, z, Dcompar;
double l,b;
char limit;
int ndir;
ndir = -1;
z = 1.77;

	if (dist[1] != -1)  		{	// We base on the analytical form of the Gaussian
		parallax = 1./D;
		res = 1./(dist[1]*sqrt(pi*2)) * exp (-pow(parralax - dist[0], 2)/(2.*pow(dist[1], 2.))) 
	}
	else		{	// We find appropriate for D DM and then compare it with the actual DM 
		l = dist[8]/180.*pi;
		b = dist[9]/180.*pi;
		dist1 = dist[0];
		dist2 = D;
		Dcompar = z / sin(b);
		dmdsm_ (&l, &b, &ndir, &DM1, &dist1, &limit, &sm, &smtau, &smtheta);	
		if (D_compar > D)	{
			D = Dcompar;
			dmdsm_ (&l, &b, &ndir, &DM2, &dist2, &limit, &sm, &smtau, &smtheta);
			if (limit != ">")	{			
				cout<<"Something went wrong!!!"<<endl;
				cout<<"Exactly: D is "<< D <<", Dcompar is "<<Dcompar<<endl;
				cout<<"limit is "<<limit<<", and DM2 is "<<DM2<<endl;
				exit(2);
			}
		}
		else
			dmdsm_ (&l, &b, &ndir, &DM2, &dist2, &limit, &sm, &smtau, &smtheta);

		res = 1./(30.*sqrt(pi*2)) * exp (-pow(DM1 - DM2, 2)/(2.*pow(30., 2.)));
	
	}

return res;
}

//-------------------------------------------------
// This function calculates correction of velocity
// for differential rotation of the Galaxy
//-------------------------------------------------

double delta_vl (double * dist, double D)	{
double res;
double theta, R, R0, l;
R = 8.5;
R0= 8.5;
theta = 225;

l = dist[8]/180.*pi;

res = (theta/R) * (R0 * cos(l) - D) - theta*cos(l);

return res;
}

double prob_vl (double * entry_dist, double * entry_prmot, double vl)	{
double res;
double D, Dmin, Dmax;
double fl, fr, h;
double mu_c, mu_s;

mu_c = entry_prmot[0];
mu_s = entry_prmot[1];

Dmin = 2.5;

	do {

	D = Dmin;	
	fl = f(&entry_dist[0], mu_c - 3*mu_s, D, vl);
	fr = f(&entry_dist[0], mu_c - 3*mu_s, D+0.1, vl);
	h = 0.1;
	Dmin = D - h * fl / (fr-fl);

//	cout<<Dmin<<endl;

	} while (abs(f(&entry_dist[0], mu_c - 3*mu_s, Dmin, vl))> 0.0001);

	cout << "Dmin is "<<Dmin<<endl;

Dmax = 2.5;

	do {

	D = Dmax;	
	fl = f(&entry_dist[0], mu_c + 3*mu_s, D, vl);
	fr = f(&entry_dist[0], mu_c + 3*mu_s, D+0.1, vl);
	h = 0.1;
	Dmax = D - h * fl / (fr-fl);

//	cout<<Dmax<<endl;

	} while (abs(f(&entry_dist[0], mu_c + 3*mu_s, Dmax, vl))> 0.0001);

	cout << "Dmax is "<<Dmax<<endl;

double sum; 
int counter;

sum = 0.;

counter=0;

x_left  = Dmax;
x_right = Dmin;

eps=0.018;

h_init = (Dmin-Dmax)/20.;
x = Dmax; 
h = h_init;
	do {
		
		h = eps * h / abs(pdf_prmot(x, mu_c, mu_s) - pdf_prmot(x+h, mu_c, mu_s));
		if (h > 6.*mu_s/5.)
			h = 6*mu_s/5.;
		if (h + x > x_right)
			h = x_right - x;
		x_next = x + h;
		k1 = h * pdf_prmot(x, mu_c, mu_s);
		k2 = h * pdf_prmot(x+h/4., mu_c, mu_s);
		k3 = h * pdf_prmot(x+3.*h/4., mu_c, mu_s);
		k4 = h * pdf_prmot(x+h, mu_c, mu_s);
		sum += (k1 + 2*k2 + 2*k3 + k4) / 6.;
		counter++;
	
//		cout<<counter<<"\t"<<x<<"\t"<<x_next<<"\t"<<h<<"\t"<<sum<<endl;
		x = x_next;

	} while (x<x_right);


return res;
}

double f(double * dist, double mu, double D, double vl)	{
double res, b;

b = dist[9]/180.*pi;

res = (vl + delta_vl(dist, D))/(D*cos(b))*9.51e5/206265. - mu;

return res;
}
