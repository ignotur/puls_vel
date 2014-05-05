#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main (int argv, char * argc[]) {

ifstream in         (argc[1]);
ofstream out_smooth (argc[2]);
ofstream out_distr  (argc[3]);


double data[10000], dist[10000], integ[10000], sum;
int n = 0, n_of_flag, n_of_second_flag;
bool flag, second_flag;

	do {
		in >> dist[n];
		in >> data[n];
		n++;
	} while (!in.eof());
n--;

	for (int i=n-1; i >= 0; i--) 	{

		flag = false;
		second_flag = false;

		for (int j=i-1; j >=0; j--)	
			if (data[j] <= data[i])	{	
				flag = true;
				n_of_flag = j;	}


		if (flag)	
			cout<<i<<"\t"<<data[i]<<"\t"<<n_of_flag<<"\t"<<data[n_of_flag]<<endl;

/*		if (flag)							{
			for (int j=i-1; j >=0; j--)	
				data[j] += -(data[i+1] - data[i])/dist[i];
			data[i] = data[i+1] -(data[i+1] - data[i])/dist[i];					
		}
*/

		if (flag)							{
			sum = 0;
			for (int j=i; j >= n_of_flag; j--)	
				sum += data[j];
			
			for (int j=n_of_flag-1; j >=0; j--)
				if (data[j] <= sum/(dist[i] - dist[n_of_flag]+1))	{
					second_flag = true;
					n_of_second_flag = j-1;				}

			if (second_flag)	{
				sum = 0;
				for (int j=i; j >= n_of_second_flag; j--)	
					sum += data[j];
				n_of_flag = n_of_second_flag; 				}


			for (int j=i; j >= n_of_flag; j--) 				{	
				data[j] = sum/(dist[i] - dist[n_of_flag]+1);
				out_smooth << dist[j] << "\t" << data[j] << endl;	}
		
		i = n_of_flag;
		}
		else
			out_smooth << dist[i] << "\t" << data[i] << endl;

	}

	for (int i=0; i < n; i++) 	{
		integ[i] = (data[i] - data[i+1]) * dist[i];
		out_distr << dist[i] << "\t" << integ[i] <<endl;
	}
		


return 0;
}
