#include <iostream>
#include <cmath>

using namespace std;

float dwod_ (int * n, float * x); 
double apriory (double, double, double);
	
int main () {
double h1, h2, h3, sum;
int n;
float val[3];

h1 = 20.;
h2 = 20.;
h3 = 0.1;

sum = 0;

	for (int i=1; i < 50; i++)	
		for (int j=1; j < 50; j++)
			for (int k=1; k < 10; k++)	{
				val[0] =  (i)*h1/1000.;
				val[1] =  (j)*h2/1000.;
				val[2] =  (k-1)*h3; 
		n   = 3;
		sum += h1*h2*h3 * apriory(i*h1, j*h2, k*h3) * pow (10., dwod_ (&n, &val[0]));
	}

	cout << "So, here is the result -- " << sum << endl;


return 0;
}
