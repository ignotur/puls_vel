#include <cmath>

using namespace std;

double const pi = 3.1415926;

double apriory (double v)	{
double res, v_;

v_ = v * sqrt(8./pi);

res = 33./20. / pow(15. + v_, 2.);

return res;
}
