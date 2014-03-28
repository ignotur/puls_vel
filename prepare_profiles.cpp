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
double pdf_dist (double *, double);
double pdf_prmot(double, double, double);
double delta_vl (double *, double);
double prob_vl (double *, double *, double);	
double prob_vl_special (double *, double *, double);
double f(double *, double, double, double);
double pdf_dist_fast (double *, double);
bool   compare       (double *, double *);
void   copy          (double *, double *);


int main (int argv, char * argc[]) 	{

ifstream in_dist  (argc[1]);
ifstream in_prmot (argc[2]);

ofstream out_prof;

double dist[10][1000], prmot[6][1000], res[1000];
double entry_prmot[6], entry_dist[10]; 
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
			in_prmot >> prmot[i][n_prmot];
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
	else 
		cout<<"Data have been read"<<endl;
	
// The main loop. Calculate profiles for pulsars one by one.

	for (int i=0; i < n_dist; i++)	{

		for (int j=0; j < 10; j++)		
			entry_dist [j] = dist [j][i];
		for (int j=0; j < 6; j++)	
			entry_prmot[j] = abs(prmot[j][i]);	
		
		cout<<"Working on profile -- "<<i<<endl;	

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

		for (int j=10; j < 1000; j++)
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

double pdf_dist (double * dist, double D)	{
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
		res = 1./(dist[1]*sqrt(pi*2)) * exp (-pow(parallax - dist[0], 2)/(2.*pow(dist[1], 2.))); 
	}
	else		{	// We find appropriate for D DM and then compare it with the actual DM 
		res = pdf_dist_fast (dist, D);	

	}

return res;
}


double pdf_dist_fast (double * dist, double D) {
double res, z;
double DM1, DM2, dist1, dist2;
double sm, smtau, smtheta;
char limit;
double static dist_in_use [10];
double static pdf_at_dist [300];
double intpart, fractpart;
double probl, probr;
double l,b, Dcompar;
int ndir;
ndir = -1;
z = 1.77;

	if (compare(dist, &dist_in_use[0]))	{
		fractpart = modf (D/0.05, &intpart);
	
//		cout<<intpart << "\t" <<pdf_at_dist[(int) intpart]<<endl;

		if (intpart >= 300)
			res = 0;
		else 				{
			probl = pdf_at_dist[(int) intpart];
			probr = pdf_at_dist[(int) intpart +1];
			res = (1-fractpart) * probl + fractpart * probr;
	//		res = probl;
		}
	}
	else	{
	cout<<"It is a different pulsar!"<<endl;	
		copy(&dist_in_use[0], dist);
		l = dist[8]/180.*pi;
		b = dist[9]/180.*pi;

		dist1 = dist[0];
		Dcompar = abs(z / sin(b));

		for (int i=0; i < 300; i++)	{
			dist2 = i*0.05;
			
			if (i==0)
				dist2 = 0.001;
				
			dmdsm_ (&l, &b, &ndir, &DM1, &dist1, &limit, &sm, &smtau, &smtheta);	
			dmdsm_ (&l, &b, &ndir, &DM2, &dist2, &limit, &sm, &smtau, &smtheta);
	
			res = 1./(15.*sqrt(pi*2)) * exp (-pow(DM1 - DM2, 2)/(2.*pow(15., 2.)));
		
			if (i*0.05*cos(b)>15.)
				res=0.;

			pdf_at_dist[i] = res;
		}

		dist2 = D;
		dmdsm_ (&l, &b, &ndir, &DM2, &dist2, &limit, &sm, &smtau, &smtheta);

		res = 1./(15.*sqrt(pi*2)) * exp (-pow(DM1 - DM2, 2)/(2.*pow(15., 2.)));
		
		if (D*cos(b)>15.)
			res=0.;

	}


return res;
}

bool   compare       (double * a, double * b)	{
bool flag;
flag = true;

	for (int i=0; i < 10; i++)	{	
//		cout<<a[i]<<"\t"<<b[i]<<endl;
		if (a[i] != b[i])
			flag=false;
	}

return flag;
}

void   copy          (double * a, double * b) 	{

cout<<"We start copying!"<<endl;

	for (int i=0; i < 10; i++)	{		
		a[i] = b[i];
		cout<<a[i]<<"\t"<<b[i]<<endl;
	}

}



//-------------------------------------------------
// This function calculates correction of velocity
// for differential rotation of the Galaxy
//-------------------------------------------------

double delta_vl (double * dist, double D)	{
double res;
double theta, R, R0, l, b;
R = 8.5;
R0= 8.5;
theta = 225;

l = dist[8]/180.*pi;
b = dist[9]/180.*pi;

res = (theta/R) * (R0 * cos(l) - D*cos(b)) - theta*cos(l);

//res = 0;

return res;
}

//-------------------------------------------------------------------------
// The main function of the program - calculates integral eq. (1)
// from the article by Brisken et al. (2003) (taking into account
// correction D*cos(b)) for particular velocity vl.
// It was tested for velocities more than 40 km/s
// The integrator analyze P(mu_l) profile of probability density function
// which is usually much sharper than P(D).
//--------------------------------------------------------------------------


double prob_vl (double * entry_dist, double * entry_prmot, double vl)	{
double res;
double mu_c, mu_s;
double x_left, x_right;
double h_init, h;
double sum, x, x_next, eps;
double k1, k2, k3, k4;
double h_D, h_D_next;
double Dmin, Dmax, D_next, Dinter;
double prob_l, prob_c, prob_r;
double D, fl, fr;
double D_prev, b;
double h_newton;
int emergence;

mu_c    = entry_prmot[0];
mu_s    = entry_prmot[1];

b       = entry_dist[9]/180.*pi; 

x_left  = mu_c - 3*mu_s;  // value of P(mu_l) is much larger 
x_right = mu_c + 3*mu_s;  // than 0 only from -3 to 3

eps = 1./(95.*5.+1.);     // this value is a result of modelling

h_init = 6*mu_s/20.;
x = mu_c + mu_s*3; 
h = h_init;

sum = 0.;


h=0.1;

Dmin = 2.5;
emergence=0;

// Here we search for Dmin (distance from which
// we are going to integrate)

	if (((int)vl)%100==0)
		cout<<"vl -- "<<vl<<endl;

do {
	D = Dmin;	
	fl = f(&entry_dist[0], mu_c - 3*mu_s, D, vl);
	fr = f(&entry_dist[0], mu_c - 3*mu_s, D+h, vl);
	if (D-h*fl/(fr-fl) <= 0. && emergence < 70)	{	
		Dmin/=2.;	
		h = 0.1*Dmin;				}
	else if (D-h*fl/(fr-fl) <= 0.)
		Dmin= (D + Dmin/2.)/2.;	
	else
		Dmin = D - h * fl / (fr-fl);

	emergence++;

//	cout<<Dmin<<endl;
	
	//	if (emergence==100)
	//		h /= 2.;	

		if (emergence>150)	{
			cout<<"Conditions for Dmin are not satisfied!"<<endl;
			exit(3);
		}

} while (abs(f(&entry_dist[0], mu_c - 3*mu_s, Dmin, vl))> 0.0001);

//cout<<Dmin<<endl;

Dmax = 2.5;
emergence=0;

// Here we search for Dmax (distance until that
// we are going to integrate)

h=0.1;

do {

	D = Dmax;	
	fl = f(&entry_dist[0], mu_c + 3*mu_s, D, vl);
	fr = f(&entry_dist[0], mu_c + 3*mu_s, D+h, vl);

	if (D-h*fl/(fr-fl) <= 0. && emergence < 70)	{
		Dmax/=2;
		h = 0.1 * Dmax;		}
	else
		Dmax = D - h * fl / (fr-fl);

//	cout<<Dmax<<endl;
	
	emergence++;

		if (emergence>150)	{
			cout<<"Conditions for Dmax are not satisfied!"<<endl;
			exit(4);
		}


} while (abs(f(&entry_dist[0], mu_c + 3*mu_s, Dmax, vl))> 0.0001);

h_newton = h;

//cout<<Dmax<<endl;

	D = Dmax;
	
	// So, starting integration process.
	// The step is adaptive and its value dependes on derivative
	// This integrator has two different steps. One - h is actual step
	// for sharpest function P(mu_l), which we do not need because
	// we integrate on distance, not proper motion.
	// So, based on h we calculate h_D


	// Let us check first should we integrate at all? It may happen 
	// that the PDF for distances is too small (say less than 1e-6 or 1e-7)

	Dinter = (Dmax + Dmin) / 2.;
	
	prob_l = pdf_dist(&entry_dist[0], Dmin);
	prob_c = pdf_dist(&entry_dist[0], Dinter);
	prob_r = pdf_dist(&entry_dist[0], Dmax);

//cout << prob_l<<"\t"<<prob_c<<"\t"<<prob_r<<endl; 

	if (prob_l < 1e-8 && prob_c < 1e-8 && prob_r < 1e-8)	{
	//	cout << prob_l<<"\t"<<prob_c<<"\t"<<prob_r<<endl; 
		return 0;
	}


	//cout<<Dmax<<"\t"<<Dmin<<endl;

	do {
		h = eps * h / abs(pdf_prmot(x, mu_c, mu_s) - pdf_prmot(x-h, mu_c, mu_s));
		if (h > 6.*mu_s/5.)
			h = 6*mu_s/5.;
		if (x-h < x_left)
			h = x - x_left;
		x_next = x - h;
		
		// Compute h_D by means on Newton algorithm
		
		D_prev = D;

		do {
			D_next = D;
			fl = f(&entry_dist[0], x_next, D_next, vl);
			fr = f(&entry_dist[0], x_next, D_next+h_newton, vl);
			D = D_next - h_newton * fl/(fr-fl);
		} while (abs(f(&entry_dist[0], x_next, D, vl))> 0.001);	
	
		h_D = D - D_prev;	
//cout<<D<<endl;

		k1 = h_D * pdf_prmot(x, mu_c, mu_s)         *  pdf_dist(&entry_dist[0], D);
		k2 = h_D * pdf_prmot(x+h/4., mu_c, mu_s)    *  pdf_dist(&entry_dist[0], D + h_D/4.);
		k3 = h_D * pdf_prmot(x+3.*h/4., mu_c, mu_s) *  pdf_dist(&entry_dist[0], D + 3.*h_D/4.);
		k4 = h_D * pdf_prmot(x+h, mu_c, mu_s)       *  pdf_dist(&entry_dist[0], D + h_D);
		sum += (k1 + 2*k2 + 2*k3 + k4) / 6.;
//		counter++;

		x = x_next;

	} while (x>x_left);

res = sum;
return res;
}

//----------------------------------------------------------------------
// This function computes residuals for proper motion
//----------------------------------------------------------------------

double f(double * dist, double mu, double D, double vl)	{
double res, b;

b = dist[9]/180.*pi;

res = (vl + delta_vl(dist, D))/(D)/9.51e5*206265. - mu;

return res;
}

//----------------------------------------------------------------------
// This function should call prob_vl multiple times and keep result
// in the massive res. As far as all calls are finished it should
// normalize result profile.
//----------------------------------------------------------------------


void profile(double * entry_dist, double * entry_prmot, double * res)	{
double sum;

sum = 0;


	if (entry_prmot[0] - 3.*entry_prmot[1] < 0.)		{
		cout << "We use standard (slow) scheme here."<<endl;
		for (int i=0; i < 1000; i++)				{	
			if (i==0 || res[i-1] > (sum / 1e6))	{
				res[i] = prob_vl_special (entry_dist, entry_prmot, i);
				sum += res[i];
			}
			else
				res[i] = 0;
		}
	}
	else							{
		cout << "We use fast scheme here."<<endl;
		for (int i=0; i < 25; i++)
			res[i] = 0.;
		for (int i=20; i < 1000; i++)				{
				res[i] = prob_vl (entry_dist, entry_prmot, i);
				sum += res[i];
		}
	}	

	
	// Here we normalise the profile

	for (int i=0; i < 1000; i++)
		if (sum != 0)
			res[i] /= sum;
}

//----------------------------------------------------------------------
// This is a spetial case for our integration procedure when we have 
// proper motion which is extremely small and less than its error.
// In this case we intend to integrate the profile as it is with
// some high-order numerical scheme.
//----------------------------------------------------------------------
double prob_vl_special (double * entry_dist, double * entry_prmot, double vl)	{
double sum;
double h, prob_c;
double b, mu_c, mu_s, lf, lc, lr, D;

sum = 0;
h   = 0.03;

	if (((int)vl)%100==0)
		cout<<"vl -- "<<vl<<endl;


b = entry_dist[9] * pi / 180.;
mu_c = entry_prmot[0];
mu_s = entry_prmot[1];


	for (int i=1; i < 500; i++)	{
//cout<<i<<endl;
		D = (double) i * h;
		lf = h * pdf_dist(entry_dist, D)         * pdf_prmot( (vl + delta_vl(entry_dist, D))        / (D )       * 206265/9.51e5, mu_c, mu_s ); 
		lc = h * pdf_dist(entry_dist, D + 0.5*h) * pdf_prmot( (vl + delta_vl(entry_dist, D + 0.5*h))/((D+0.5*h)) * 206265/9.51e5, mu_c, mu_s );
		lr = h * pdf_dist(entry_dist, D + 1.0*h) * pdf_prmot( (vl + delta_vl(entry_dist, D + 1.0*h))/((D+1.0*h)) * 206265/9.51e5, mu_c, mu_s );
		sum += (lf + 4 * lc + lr) / 6.;
	}

return sum;
}		
