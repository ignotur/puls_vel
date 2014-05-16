#include <iostream>
#include <cmath>

float dwod_ (int * n, float * x); 

using namespace std;

int main (int argc, char * argv[]) {

long double h, sum, sum_l, sum_r;
float val[3];
int n;
double v_max1, v_max2, v_max3;

v_max1 = atof(argv[1]);
v_max2 = atof(argv[2]);
v_max3 = atof(argv[3]);

double res_s1l,  res_s1r, res_s2l, res_s2r, res_s3l, res_s3r;

sum   = 0;
sum_l = 0;
sum_r = 0;
h=5;
n=3;

	//
	// The first sigma
	//
	val[0] = v_max1/1000.;
	val[1] = v_max2/1000.;
	val[2] = v_max3;


	for (int i=2; i < 200; i++)	{
		val[0] = (i-1)* h/1000.;
		sum += h * pow((long double) 10., dwod_ (&n, &val[0]));
	}

	cout << "Sum is " << sum << endl;

//	sum *= 0.3415;

	for (int i=0; i < 200; i++)	{
		val[0] = (v_max1 + i * h)/1000.;
		sum_l += h * pow((long double) 10., dwod_(&n, &val[0]));
		
		if (abs(sum_l/sum - 0.3415) < 0.03)	{
			cout << 1000*(val[0]-v_max1/1000.) << "\t" << sum <<"\t" << sum_l << "\t" << sum_l/sum <<endl;
			break;
		}
	}	

	res_s1l = 1000*(val[0]-v_max1/1000.); 

	for (int i=0; i < 200; i++)	{
		val[0] = (v_max1 - i * h)/1000.;
		sum_r += h * pow((long double) 10., dwod_(&n, &val[0]));
		
		if (val[0] < 0)	{
			cout << "Alert! val < 0" << endl;
			exit(1);
		}

		cout << sum_r/sum << endl;

		if (abs(sum_r/sum - 0.3415) < 0.03)	{
			cout << 1000*(v_max1/1000. - val[0]) << "\t" <<sum <<"\t" << sum_r << "\t" << sum_r/sum <<endl;
			break;
		}
	}

	res_s1r =  1000*(v_max1/1000. - val[0]);	

	//
	// The second sigma
	//
sum   = 0;
sum_l = 0;
sum_r = 0;


	val[0] = v_max1/1000.;
	val[1] = v_max2/1000.;
	val[2] = v_max3;


	for (int i=2; i < 200; i++)	{
		val[1] = (i-1)* h/1000.;
		sum += h * pow((long double) 10., dwod_ (&n, &val[0]));
	}


	for (int i=0; i < 200; i++)	{
		val[1] = (v_max2 + i * h)/1000.;
		sum_l += h * pow((long double) 10., dwod_(&n, &val[0]));
		
		if (abs(sum_l/sum - 0.3415) < 0.03)	{
			cout << 1000*(val[1]-v_max2/1000.) << "\t" << sum <<"\t" << sum_l << "\t" << sum_l/sum <<endl;
			break;
		}
	}	

	res_s2l = 1000*(val[1]-v_max2/1000.); 

	for (int i=0; i < 200; i++)	{
		val[1] = (v_max2 - i * h)/1000.;
		sum_r += h * pow((long double) 10., dwod_(&n, &val[0]));
		
		if (val[1] < 0)	{
			cout << "Alert! val < 0" << endl;
			exit(1);
		}

		if (abs(sum_r/sum - 0.3415) < 0.03)	{
			cout << 1000*(v_max2/1000. - val[1]) << "\t" <<sum <<"\t" << sum_r << "\t" << sum_r/sum <<endl;
			break;
		}
	}

	res_s2r = 1000*(v_max2/1000. - val[1]);

	//
	// w
	//
sum   = 0;
sum_l = 0;
sum_r = 0;


	val[0] = v_max1/1000.;
	val[1] = v_max2/1000.;
	val[2] = v_max3;

	for (int i=2; i < 200; i++)	{
		val[2] = i*h/1000.;
		sum += h * pow((long double) 10., dwod_ (&n, &val[0]));
	}


	for (int i=0; i < 200; i++)	{
		val[2] = v_max3 + i * h/1000.;
		sum_l += h * pow((long double) 10., dwod_(&n, &val[0]));
		
		if (abs(sum_l/sum - 0.3415) < 0.03)	{
			cout << 1000*(val[2]-v_max3) << "\t" << sum <<"\t" << sum_l << "\t" << sum_l/sum <<endl;
			break;
		}
	}	

	res_s3l = val[2]-v_max3;

	for (int i=0; i < 200; i++)	{
		val[2] = v_max3 - i * h/1000.;
		sum_r += h * pow((long double) 10., dwod_(&n, &val[0]));
		
		if (val[2] < 0)	{
			cout << "Alert! val < 0" << endl;
			exit(1);
		}

		if (abs(sum_r/sum - 0.3415) < 0.03)	{
			cout << 1000*(v_max3 - val[2]) << "\t" <<sum <<"\t" << sum_r << "\t" << sum_r/sum <<endl;
			break;
		}
	}

	res_s3r = v_max3 - val[2];

	cout << "----------- Summary ------------" << endl;
	cout << res_s1l  << "\t" << res_s1r << endl;
	cout << res_s2l  << "\t" << res_s2r << endl;
	cout << res_s3l  << "\t" << res_s3r << endl;
	cout << "----------- Summary ------------" << endl;

return 0;
}
