#include <cmath>

using namespace std;

double const pi = 3.1415926;

double apriory (double v1, double v2, double w)	{
double res, v_1, v_2;

v_1 = sqrt(8./pi) * v1*w ;
v_2 = sqrt(8./pi) * v2 * (1. - w);

res = 33./20.*33./20. / pow(15. + v_1, 2.) / pow(15. + v_2, 2.);

return res;
}
