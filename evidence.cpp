#include <iostream>
#include <cmath>

using namespace std;

float dwod_ (int * n, float * x); 
double apriory (double v);
	
int main () {
double h, sum;
int n;
float val;

h = 5.;

sum = 0;

	for (int i=1; i < 200; i++)	{
		val =  (i-1)*h/1000.;
		n   = 1;
		sum += h * apriory(i*h) * pow (10., dwod_ (&n, &val));
	}

	cout << "So, here is the result -- " << sum << endl;


return 0;
}
