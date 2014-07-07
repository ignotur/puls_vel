#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

float dwod_ (int * n, float * x); 
double apriory (double v);
	
int main () {
double h;
long double sum;
int n;
float val;

ofstream out ("model.dat", ios::app);


h = 5.;

sum = 0;

	for (int i=1; i < 200; i++)	{
		val =  (i-1)*h/1000.;
		n   = 1;
		sum += h * (long double)  apriory(i*h) * pow ((long double)10., dwod_ (&n, &val));
	}

	cout << "So, here is the result -- " << sum << endl;

out << sum << endl;

return 0;
}
