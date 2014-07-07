#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <cstdlib>

using namespace std;

double const pi = 3.1415926;
double sigma_DM = 15.;


extern "C" {
void dmdsm_ (float *l, float *b, int *ndir, float *dmpsr, float *dist, char *limit, float *sm, float *smtau, float *smtheta, float *smiso);
}

void profile  (double * dist, double * prmot, double * res);
double pdf_dist (double *, double);
double pdf_prmot(double, double, double);
double delta_vl (double *, double);
double prob_vl_special (double *, double *, double);
double f(double *, double, double, double);
double pdf_dist_fast (double *, double);
bool   compare       (double *, double *);
void   copy          (double *, double *);


int main (int argv, char * argc[]) 	{

ifstream in_dist  (argc[1]);
ifstream in_prmot (argc[2]);

ofstream out_prof;

double dist[10][1000], prmot[6][1000], res[2000];
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
			entry_prmot[j] = prmot[j][i];	
		
		cout<<"Working on profile -- "<<i<<endl;	

		profile(&entry_dist[0], &entry_prmot[0], &res[0]);	 // Call a function to compute profile	
	
		// Open file with appropriate name
		name << "profiles/profile_"<< i + n_dist <<".dat";
		name_str = name.str();
		char *basic_name = new char [name_str.size()];
		memcpy(basic_name, name_str.c_str(), name_str.size());
		name.clear();
		name.str(string());	
		
		cout << basic_name <<endl;
		out_prof.open(basic_name);

		for (int j=10; j < 1500; j++)
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
double DM2, R, l, b;

res = 1;

	if 	(dist[1] != -1 && dist[2] != -1 && D > 1./dist[0])  
		res = 1./pow(D,2.) * exp(-0.5 * pow(dist[0] - 1/D, 2) / pow (dist[1], 2));
	else if (dist[1] != -1 && dist[2] != -1 && D < 1./dist[0])  
		res = 1./pow(D,2.) * exp(-0.5 * pow(dist[0] - 1/D, 2) / pow (dist[2], 2));
	else 
		res = pdf_dist_fast (dist, D);

//	if (dist[3] != -1 && dist[4] != -1 && dist[5] != -1 && dist[6] != -1) 
//		res *= 0.5 * (erf(dist[3] / sqrt(2) / dist[4]) - erf((dist[3] - D)/sqrt(2)/dist[4])) * 0.5 * (1. + erf((dist[5]-D)/sqrt(2)/dist[6])); 
//	else if (dist[3] != -1 && dist[4] != -1)
//		res *= 0.5 * (erf(dist[3] / sqrt(2) / dist[4]) - erf((dist[3] - D)/sqrt(2)/dist[4]));
//	else if (dist[5] != -1 && dist[6] != -1)
//		res *=  0.5 * (1. + erf((dist[5]-D)/sqrt(2)/dist[6])); 

	l = dist[8] * pi/180.;
	b = dist[9] * pi/180.;

//	R = sqrt(8.5*8.5 + pow(D*cos(b), 2) - 2*8.5*D*cos(b)*cos(l));

//	res *= pow(R/8.5, 1.9) * pow(D, 2) * exp(-abs(D*sin(b))/0.330 - 5*(R - 8.5)/8.5);

//	if (dist[7] != -1)
//		res *= 1./D * exp(-0.5*pow((log10(dist[7])+2*log10(D)+1.1)/0.9 ,2)); 



	if (isnan(D) || isinf(D))	{
		cout << "The function pdf_dist obtained D which is "<<D<<endl;
		exit(6);	
	}


return res;
}


double pdf_dist_fast (double * dist, double D) {
double res, z;
float DM1, DM2, dist1, dist2;
float sm, smtau, smtheta, smiso;
char limit;
double static dist_in_use [10];
double static pdf_at_dist [300];
double intpart, fractpart;
double probl, probr;
float l,b; 
double Dcompar;
int ndir;
ndir = -1;
z = 1.77;

	if (compare(dist, &dist_in_use[0]))	{
		fractpart = modf (D/0.05, &intpart);
	
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
		
	
	//	cout << i<< "\t" << DM1 << "\t" << DM2 << endl;

	
			if (i==0)
				dist2 = 0.001;
				
			dmdsm_ (&l, &b, &ndir, &DM1, &dist1, &limit, &sm, &smtau, &smtheta, &smiso);	
			dmdsm_ (&l, &b, &ndir, &DM2, &dist2, &limit, &sm, &smtau, &smtheta, &smiso);

			sigma_DM = 0.4 * DM1;
		
			res = 1./(sigma_DM*sqrt(pi*2)) * exp (-pow(DM1 - DM2, 2)/(2.*pow(sigma_DM, 2.)));
		
			if (i*0.05*cos(b)>15.)
				res=0.;

			pdf_at_dist[i] = res;
		}

		dist2 = D;
//cout<< "Here" << "\t" << dist2 <<endl;
//		dmdsm_ (&l, &b, &ndir, &DM2, &dist2, &limit, &sm, &smtau, &smtheta);
		dmdsm_ (&l, &b, &ndir, &DM2, &dist2, &limit, &sm, &smtau, &smtheta, &smiso);
//cout<< "Here" << endl;	

		sigma_DM = 0.4 * DM1;

		res = 1./(sigma_DM*sqrt(pi*2)) * exp (-pow(DM1 - DM2, 2)/(2.*pow(sigma_DM, 2.)));
		
		if (D*cos(b)>15.)
			res=0.;
	cout<<"Initialisation has finished"<<endl;
	}

//cout<<"Here"<<endl;

return res;
}



bool   compare       (double * a, double * b)	{
bool flag;
flag = true;

	for (int i=0; i < 10; i++)	{	
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

/*  Scalar product of two 3-vectors  (double precision) */
double product_ecliptic(double *va,double *vb)
{
  return va[0]*vb[0]+va[1]*vb[1]+va[2]*vb[2];
}

//-------------------------------------------------
// This function calculates correction of velocity
// for differential rotation of the Galaxy
//-------------------------------------------------

double delta_vl (double * dist, double D)	{
double res;
double theta, R, R0, l, b;
double sin_alpha, cos_alpha;
double apex_l, apex_b;
double cos_lambda, sin_lambda;
double cos_alpha_add;
double sign_pec, v0gal;
double vsun[4], b2000, l2000, vsl, vsb;
double p[3], q[3], r[3], vgal[3], galcart[3];

v0gal = 225;

vsun[0]=9.2; vsun[1]=10.5; vsun[2]=6.9;

apex_l = 58.87/180.*pi;
apex_b = 17.72/180.*pi;

R = 8.5;
R0= 8.5;
theta = 225; //225;

l = dist[8]/180.*pi;
b = dist[9]/180.*pi;

l2000 = l;
b2000 = b;

p[0] = -sin(l2000);
p[1] =  cos(l2000);
p[2] =  0.0;

q[0] = -sin(b2000)*cos(l2000);
q[1] = -sin(b2000)*sin(l2000);
q[2] =  cos(b2000);

r[0] =  cos(b2000)*cos(l2000);
r[1] =  cos(b2000)*sin(l2000);
r[2] =  sin(b2000);

/* Find the proper motions expected from galactic rotation using rotation */
/* curve model */
double vrgal, pmlrot, pmbrot, thetagal;

galcart[0] = cos(l2000)*cos(b2000);
galcart[1] = sin(l2000)*cos(b2000);
galcart[2] = sin(b2000);

for (int i=0;i < 3;i++)
    galcart[i]=galcart[i]*D; 

thetagal = atan2(galcart[1],R0-galcart[0]);
vrgal = v0gal; /* Use flat rotation curve */
vgal[0] = vrgal*sin(thetagal);
vgal[1] = vrgal*cos(thetagal)-v0gal;
vgal[2] = 0.0;

pmlrot=product_ecliptic(p,vgal);
pmbrot=product_ecliptic(q,vgal); 

vsl = -product_ecliptic(p,vsun);
vsb = -product_ecliptic(q,vsun); 

res = pmlrot+vsl;

return res;
}
//----------------------------------------------------------------------
// This function computes residuals for proper motion
//----------------------------------------------------------------------

double f(double * dist, double mu, double D, double vl)	{
double res, b;

b = dist[9]/180.*pi;

//res = abs(vl + delta_vl(dist, D))/(D)/9.51e5*206265. - abs(mu);

res = abs(mu + delta_vl(dist, D)/D/9.51e5*206265.) - vl/D/9.51e5*206265.;

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


		cout << "We use standard (slow) scheme here."<<endl;
		for (int i=0; i < 1500; i++)				{	
				res[i] = prob_vl_special (entry_dist, entry_prmot, i);
				sum += res[i];
		}
	
	// Here we normalise the profile

	for (int i=0; i < 1500; i++)
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
double sign_v;
double dang2vel;

dang2vel = 4.610573776;

sum = 0;
h   = 0.001;

	if (((int)vl)%100==0)
		cout<<"vl -- "<<vl<<endl;


b = entry_dist[9] * pi / 180.;
mu_c = entry_prmot[3];
mu_s = entry_prmot[4];

sign_v = mu_c / abs(mu_c);

vl = sign_v * vl;

	for (int i=1; i < 15000; i++)	{
		D = (double) i * h;
		lf = h * pdf_dist(entry_dist, D)         /D         * pdf_prmot( (vl /*+ delta_vl(entry_dist, D)*/)        / (D )      /dang2vel , mu_c, mu_s ); 
		lc = h * pdf_dist(entry_dist, D + 0.5*h) /(D+0.5*h) * pdf_prmot( (vl /*+ delta_vl(entry_dist, D + 0.5*h)*/)/((D+0.5*h))/dang2vel , mu_c, mu_s );
		lr = h * pdf_dist(entry_dist, D + 1.0*h) /(D+1.0*h) * pdf_prmot( (vl /*+ delta_vl(entry_dist, D + 1.0*h)*/)/((D+1.0*h))/dang2vel , mu_c, mu_s );
		sum += (lf + 4 * lc + lr) / 6.;
	}

return sum;
}		
