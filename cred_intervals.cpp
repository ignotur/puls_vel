#include <iostream>
#include <cmath>

float dwod_ (int * n, float * x); 

using namespace std;

int main (int argc, char * argv[]) {

long double h, sum, sum_l, sum_r;
float val;
int n;
double v_max;

v_max = atof(argv[1]);

sum   = 0;
sum_l = 0;
sum_r = 0;
h=5;
n=1;

	for (int i=1; i < 200; i++)	{
		val = (i-1)* h/1000.;
		sum += h * pow((long double) 10., dwod_ (&n, &val));
	}

//	sum *= 0.3415;

	for (int i=0; i < 200; i++)	{
		val = (v_max + i * h)/1000.;
		sum_l += h * pow((long double) 10., dwod_(&n, &val));
		
		if (abs(sum_l/sum - 0.3415) < 0.03)	{
			cout << 1000*(val-v_max/1000.) << "\t" << sum <<"\t" << sum_l << "\t" << sum_l/sum <<endl;
			break;
		}
	}	

	for (int i=0; i < 200; i++)	{
		val = (v_max - i * h)/1000.;
		sum_r += h * pow((long double) 10., dwod_(&n, &val));
		
		if (val < 0)	{
			cout << "Alert! val < 0" << endl;
			exit(1);
		}

		if (abs(sum_r/sum - 0.3415) < 0.03)	{
			cout << 1000*(v_max/1000. - val) << "\t" <<sum <<"\t" << sum_r << "\t" << sum_r/sum <<endl;
			break;
		}
	}

return 0;
}
