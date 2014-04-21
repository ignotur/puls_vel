#include <iostream>
#include <fstream>

using namespace std;

int main ()	{

ifstream in_prmot ("prmot.txt");
ifstream in_dist  ("dist.txt");
ifstream in_ages  ("ages.txt");

ofstream out_prmot ("prmot_short.txt");
ofstream out_dist  ("dist_short.txt");

int n_dist=0, n_prmot=0, n_ages=0;
double prmot[6][1000], dist[10][1000], ages[1000];

	do {
		for (int i=0; i < 6; i++)
			in_prmot >> prmot[i][n_prmot];
		n_prmot++;
	} while (!in_prmot.eof());

n_prmot--;

	do {
		for (int i=0; i < 10; i++)
			in_dist >> dist[i][n_dist];
		n_dist++;
	} while (!in_dist.eof());
n_dist--;

	do {
		in_ages >> ages[n_ages];
		n_ages++;
	} while (!in_ages.eof());

n_ages--;
	if (n_ages != n_prmot || n_dist != n_prmot)	{
		cout << "The list sizes are incompatible!"<<endl;
		return 1;
	}
	
	for (int i=0; i < n_dist; i++)	
		if (ages[i] < 5e7)					{
			for (int j=0; j < 5; j++)
				out_prmot << prmot[j][i] << "\t";
			out_prmot << prmot[5][i] << endl;
			for (int j = 0; j < 9; j++)
				out_dist << dist[j][i] << "\t";
			out_dist << dist[9][i] << endl;
		}

return 0;
}
