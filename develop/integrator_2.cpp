#include <iostream>
#include <cmath>

using namespace std;

double pdf_dist (double *, double);
double pdf_prmot (double, double, double);
double delta_vl (double *, double);
double f(double *, double, double, double);


double const pi = 3.1415926;

int main () {

double mu_c = abs(-2.92), mu_s = 5.6;
double entry_dist[10];
double x_left, x_right;
double h_init, h;
double sum, x, x_next, eps;
double k1, k2, k3, k4;
double h_D, h_D_next;
double Dmin, Dmax, D_next;
double D, fl, fr;
double vl, b;
double D_prev;
int counter;

//for (int i=20; i < 500; i++)	{

sum = 0.;

counter=0;

x_left  = mu_c - 3*mu_s;

	if (x_left < 0.)
		x_left = 0.;

x_right = mu_c + 3*mu_s;


//eps = 1./(i*5.+1.);
eps = 1./(95.*5.+1.);

//eps=0.018;

h_init = 6*mu_s/20.;
x = mu_c + mu_s*3; 
h = h_init;
D = Dmin;
vl = 40.;
// Let's find first Dmin, Dmax

entry_dist[0] = 3.72;
entry_dist[1] = 0.4;
entry_dist[8] = 234.5;
entry_dist[9] = 7.22;

b = entry_dist[9]/180.*pi;

Dmin = 2.5;

h=0.1;

int emergence=0;

	do {

	D = Dmin;
	cout<< entry_dist[0]<<endl;	
	fl = f(&entry_dist[0], mu_c - 3*mu_s, D, vl);
	fr = f(&entry_dist[0], mu_c - 3*mu_s, D+h, vl);

        cout<<entry_dist[8]<<"\t"<<entry_dist[9]<<endl;
	cout<<fl<<"\t"<<fr<<endl;
	cout<<vl<<"\t"<<D<<endl;
	cout<<"Here mu_c, mu_s -- "<<mu_c<<"\t"<<mu_s<<endl;
//	h = 0.1;
	if (D-h*fl/(fr-fl) <= 0.)	{
//		h*=2.;
		Dmin/=2;		}
	else
		Dmin = D - 0.1 * fl / (fr-fl);

	cout<< Dmin<<endl;

	emergence++;
	if (emergence>100)	{
		cout<<"Conditions for Dmin are not satisfied!"<<endl;
		exit(3);
	}

//	cout<<fl<<"\t"<<fr<<"\t"<< D-h*fl/(fr-fl) <<"\t" << D+h*fl/(fr-fl) <<"\t" <<Dmin<<endl;

	} while (abs(f(&entry_dist[0], mu_c - 3*mu_s, Dmin, vl))> 0.0001);

//	cout << "Dmin is "<<Dmin<<endl;

Dmax = 2.5;
emergence=0;

	do {

	D = Dmax;	
	fl = f(&entry_dist[0], mu_c + 3*mu_s, D, vl);
	fr = f(&entry_dist[0], mu_c + 3*mu_s, D+0.1, vl);

	if (D-h*fl/(fr-fl) <= 0.)	{
		Dmax/=2;		}
	else
		Dmax = D - 0.1 * fl / (fr-fl);

//	Dmax = D - 0.1 * fl / (fr-fl);
	if (emergence>100)	{
		cout<<"Conditions for Dmin are not satisfied!"<<endl;
		exit(4);
	}

//	cout<<Dmax<<endl;

	} while (abs(f(&entry_dist[0], mu_c + 3*mu_s, Dmax, vl))> 0.0001);

//	cout << "Dmax is "<<Dmax<<endl;

// Then we go from Dmin to Dmax

	D = Dmax;

	do {

		h = eps * h / abs(pdf_prmot(x, mu_c, mu_s) - pdf_prmot(x-h, mu_c, mu_s));
		if (h > 6.*mu_s/5.)
			h = 6*mu_s/5.;
		if (x-h < x_left)
			h = x - x_left;
		x_next = x - h;
		
		// Compute h_D
		//D = vl * 9.51e5 / (x_next * cos(b)) / 206265; 
//		cout<<"First estimate on D -- "<<D<<endl;

		D_prev = D;

		do {
		
		D_next = D;

		fl = f(&entry_dist[0], x_next, D_next, vl);
		fr = f(&entry_dist[0], x_next, D_next+0.1, vl);
		

		D = D_next - 0.1 * fl/(fr-fl);
//		cout<<"!"<<endl;
		} while (abs(f(&entry_dist[0], x_next, D, vl))> 0.001);	
	
		h_D = D - D_prev;	

		k1 = h_D * pdf_prmot(x, mu_c, mu_s)         *  pdf_dist(&entry_dist[0], D);
		k2 = h_D * pdf_prmot(x+h/4., mu_c, mu_s)    *  pdf_dist(&entry_dist[0], D + h_D/4.);
		k3 = h_D * pdf_prmot(x+3.*h/4., mu_c, mu_s) *  pdf_dist(&entry_dist[0], D + 3.*h_D/4.);
		k4 = h_D * pdf_prmot(x+h, mu_c, mu_s)       *  pdf_dist(&entry_dist[0], D + h_D);
		sum += (k1 + 2*k2 + 2*k3 + k4) / 6.;
		counter++;

//		cout<<"h_D is "<<h_D<<endl;
	
//		cout<<counter<<"\t"<<x<<"\t"<<x_next<<"\t"<<h<<"\t"<<sum<<endl;
		x = x_next;

//		cout<<x<<"\t"<<D<<"\t"<<D_next<<"\t"<<h_D<<endl;

	} while (x>x_left);
//if (counter<10)
//	cout<<"Too small number of steps!"<<endl;
//cout<<eps<<"\t"<<counter<<"\t"<<abs(sum - erf(3))<<endl;
//	cout<<i<<"\t"<<sum<<endl;

//cout<<erf(-3)<<"\t"<<erf(3)<<endl;
//}
return 0;
}

double pdf_prmot (double mu_l, double mu_c, double mu_s)	{
double res;

res = 1./(mu_s*sqrt(pi*2)) * exp (-pow(mu_l - mu_c, 2)/(2.*pow(mu_s, 2.)));

//res = exp(mu_l)/1000.;

return res;
}

double pdf_dist (double * dist, double D)	{
double parallax;
double res;

	if (dist[1] != -1)  		{	// We base on the analytical form of the Gaussian
		parallax = 1./D;
		res = 1./(dist[1]*sqrt(pi*2)) * exp (-pow(parallax - dist[0], 2)/(2.*pow(dist[1], 2.))); 
	}

return res;
}

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

double f(double * dist, double mu, double D, double vl)	{
double res, b;

b = dist[9]/180.*pi;

res = (vl + delta_vl(dist, D))/(D*cos(b))/9.51e5*206265. - mu;

return res;
}

