#include <cmath>

using namespace std;

double const pi = 3.1415926;

double apriory (double v1, double v2, double w)	{
double res, v_;

v_ = sqrt(8./pi) * (v1*w + v2*(1.-w));

res = 33./20. / pow(15. + v_, 2.);

return res;
}
