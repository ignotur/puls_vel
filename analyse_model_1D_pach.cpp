#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <sstream>

//extern "C" float dwod_ (int *, float*);

using namespace std;

double const pi = 3.1415926;


double model (double, double);
double integ (double *, double);

//---------------------------------------------------------------
// This is a function which models our distribution.
// As in the original work by Brisken we use sum of two Gaussian
// However, here we may use any model.
//---------------------------------------------------------------

double model (double vl, double v_)	{
double res;

//	res = ((1+(vl/v_))*log((1+pow(vl/v_, 2.))/pow(vl/v_, 2.))-1.)/(1+pow(vl/v_, 2))/pi;

	res = -4./pi/v_ *( log( (vl/v_) / sqrt(1+ vl*vl/v_/v_) ) + 0.5/(1.+ vl*vl/v_/v_) );


	if (isnan(res))	{
		cout << "vl -- "<< vl << ", v_ -- "<< v_ << endl;
		exit(0);
	}
	
//	res = abs(res);
	
	if (res < 0)	{
		cout << "Pachynski is negative" << endl;
		cout << "vl -- "<< vl << ", v_ -- "<< v_ << endl;
		cout << res << endl;
		cout << log ((1+pow(vl/v_, 2))/pow(vl/v_, 2))<<endl;
		exit(2);
	}

//	res = w / sqrt(2*pi*sigma_1*sigma_1) * exp (-pow(vl, 2.)/2./pow(sigma_1, 2.)) +
//	(1. - w) / sqrt(2*pi*sigma_2*sigma_2) * exp (-pow(vl, 2.)/2./pow(sigma_2, 2.));

return res;
}

//-----------------------------------------------------------------
// Here we find an integral eq. (2) using the simplest method of
// right rectangles on the grid 1 km/s.
//----------------------------------------------------------------- 

double integ (double * profile, double v_)	{
double sum;
sum = 0;

	for (int i=11; i < 1500; i++)	
		sum += profile[i] * model((double) i, v_);

return sum;
}


float dwod_ (int * n, float * x) 	{
float res;


ifstream in;

double static profile[400][2000];
bool static flag = false;
double entry_profile[2000], trash;
double sigma_1, sigma_2, w;
double v_;
long double L;

stringstream name;
string name_str;
int static i=0, counter;

	v_ = x[0] * 1000. + 10;
//	sigma_2 = x[1] * 1000.;
//	w = x[2];

//cout<< sigma_1<<"\t"<<sigma_2<<"\t"<<w<<endl;
cout <<v_ <<"\t"<<x[0]<< endl;

if (!flag)	{

	cout<<"We are reading files..."<<endl;

	do {
		name << "profiles/profile_"<<i<<".dat";
		name_str = name.str();
		char *basic_name = new char [name_str.size()];
		memcpy(basic_name, name_str.c_str(), name_str.size());
		name.clear();
		name.str(string());
		
		cout<<"Opening file "<<	basic_name << endl;

		in.open(basic_name);
		if (in.good())		{  // if a file exists, read it.
			counter = 11;
			do {
				in >> trash;
				in >> profile[i][counter];
				
				if (isnan(profile[i][counter]))	{
					cout << "The file -- "<<i<<", contains nan!"<<endl;
					exit(2);
				}

				if (isinf(profile[i][counter]))	{
					cout << "The file -- "<<i<<", contains inf!"<<endl;
					exit(3);
				}

				if ((profile[i][counter]) < 0)	{
					cout << "The file -- "<<i<<", contains pdf < 0!"<<endl;
					exit(5);
				}


				counter++;
			} while (!in.eof());
		}
		else {
			break;
		}
		in.close();
		i++;
	} while (1);
i--;

flag = true;

}
	// So, it is a right moment to start calculation of L from eq. (1)

L = 1.;
	for (int j=0; j < i; j++)	{
		for (int k = 11; k < 1500; k++)
			entry_profile[k] = profile[j][k];	
		L *= integ (&entry_profile[0], v_);
	
		if (isinf(L))	{
			cout<<"L is inf at this stage!"<<endl;
			cout<<"j = "<<j<<endl;
			exit(4);
		}
	
		if (isnan(L))	{
		cout << "L is nan" <<endl;
		exit(1);
		}
	
//		if (L == 0)	{
//			cout<<"L is zero!"<<endl;
//			cout<<"j = "<<j<<endl;
//			exit(6);
//		}

	}
//L *= 1e110;

cout<<"So, L is "<<log10(L)<<endl;

if (L!=0.)
res = log10(L);
else
res = -900;

return res;
}


