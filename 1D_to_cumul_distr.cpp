#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

int main (int argc, char * argv[]) {

ifstream in  (argv[1]);
ofstream out (argv[2]);

double data[2][10000];
int n=0;
double sum=0, total_sum=0;

	do {
		in >> data[0][n];
		in >> data[1][n];
		n++;
	} while (!in.eof());
n--;

	for (int i=0; i < n; i++)	
		total_sum += data[1][i];


	for (int i=0; i < n; i++)	{
		sum += data[1][i];
		out<<data[0][i]<<"\t"<<sum/total_sum<<endl;	}
	

return 0;
}
