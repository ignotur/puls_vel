#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <sstream>

extern "C" float dwod_ (int *, float*);

using namespace std;

double const pi = 3.1415926;


double model (double, double, double, double);
double integ (double *, double, double, double);

//---------------------------------------------------------------
// This is a function which models our distribution.
// As in the original work by Brisken we use sum of two Gaussian
// However, here we may use any model.
//---------------------------------------------------------------

double model (double vl, double sigma_1, double sigma_2, double w)	{
double res;

	res = w / sqrt(2*pi*sigma_1*sigma_1) * exp (-pow(vl, 2.)/2./pow(sigma_1, 2.)) +
	(1. - w) / sqrt(2*pi*sigma_2*sigma_2) * exp (-pow(vl, 2.)/2./pow(sigma_2, 2.));

return res;
}

//-----------------------------------------------------------------
// Here we find an integral eq. (2) using the simplest method of
// right rectangles on the grid 1 km/s.
//----------------------------------------------------------------- 

double integ (double * profile, double sigma_1, double sigma_2, double w)	{
double sum;
sum = 0;

	for (int i=11; i < 1500; i++)	
		sum += profile[i] * model((double) i, sigma_1, sigma_2, w);

return sum;
}


float dwod_ (int * n, float * x) 	{
float res;


ifstream in;

double static profile[100][2000];
bool static flag = false;
double entry_profile[2000], trash;
double sigma_1, sigma_2, w;
double L;

stringstream name;
string name_str;
int static i=0, counter;

	sigma_1 = x[0] * 1000.;
	sigma_2 = x[1] * 1000.;
	w = x[2];

cout<< sigma_1<<"\t"<<sigma_2<<"\t"<<w<<endl;

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
		L *= integ (&entry_profile[0], sigma_1, sigma_2, w);

	}
//L *= 1e110;

cout<<"So, L is "<<log10(L)<<endl;

res = log10(L);

return res;
}


