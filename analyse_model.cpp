#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <sstream>

using namespace std;

double const pi = 3.1415926;

double model (double, double, double, double);
double integ (double *, double, double, double);

int main (int argv, char * argc[]) {

ifstream in;

double profile[100][1000], entry_profile[1000], trash;
double sigma_1, sigma_2, w;
double L;

stringstream name;
string name_str;
int i=0, counter;

	if (argv != 3)	{
		cout << "You have entered not enough parameters!"<<endl;
		cout << "We need sigma_1, sigma_2, w"<<endl;
		exit(1);
	}

	sigma_1 = atof(argc[1]);
	sigma_2 = atof(argc[2]);
	w = atof(argc[3]);

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
			counter = 40;
			do {
				in >> trash;
				in >> profile[i][counter];
				counter++;
			} while (!in.eof());
		}
		else {
			cout<<"We have read "<< i+1 << "files."<<endl;
			break;
		}
		i++;
	} while (1);

	// So, it is a right moment to start calculation of L from eq. (1)

L = 1.;

	for (int j=0; j < i; j++)	{
		for (int k = 40; k < 1000; k++)
			entry_profile[k] = profile[j][k];	
		L *= integ (&entry_profile[0], sigma_1, sigma_2, w);
	}

	cout << sigma_1 << "\t" << sigma_2 << "\t" << w << "\t" << L << endl;

return 0;
}

//---------------------------------------------------------------
// This is a function which models our distribution.
// As in the original work by Brisken we use sum of two Gaussian
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

	for (int i=40; i < 1000; i++)	
		sum += profile[i] * model((double) i, sigma_1, sigma_2, w);

return sum;
}

