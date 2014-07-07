#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

float dwod_ (int * n, float * x); 
double apriory (double, double, double);
	
int main () {
double h1, h2, h3;
long double sum;
int n;
float val[3];

ofstream out ("model.dat", ios::app);

h1 = 50.;
h2 = 50.;
h3 = 0.1;

sum = 0;

	for (int i=1; i < 20; i++)	
		for (int j=1; j < 20; j++)
			for (int k=1; k < 10; k++)	{
				val[0] =  (i)*h1/1000.;
				val[1] =  (j)*h2/1000.;
				val[2] =  (k-1)*h3; 
		n   = 3;
		sum += h1*h2*h3 * (long double)  apriory(i*h1, j*h2, k*h3) * pow ((long double)10., dwod_ (&n, &val[0]));
	}

	cout << "So, here is the result -- " << sum << endl;

out << sum << endl;


return 0;
}
