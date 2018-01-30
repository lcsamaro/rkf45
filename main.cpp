// Runge-Kutta-Fehlberg method (RKF45)

#include <cmath>
#include <functional>
#include <iostream>

typedef double Real;
typedef int Index;

/*
Butcher Tableau
time | c
-----+-----
     | 5 order
     | 4 order
*/
const Real tableau[8][7]{
  {       0.0,           0.0,            0.0,            0.0,             0.0,        0.0,      0.0 },
  {   1.0/4.0,       1.0/4.0,            0.0,            0.0,             0.0,        0.0,      0.0 },
  {   3.0/8.0,      3.0/32.0,       9.0/32.0,            0.0,             0.0,        0.0,      0.0 },
  { 12.0/13.0, 1932.0/2197.0, -7200.0/2197.0,  7296.0/2197.0,             0.0,        0.0,      0.0 },
  {       1.0,   439.0/216.0,           -8.0,   3680.0/513.0,   -845.0/4104.0,        0.0,      0.0 },
  {   1.0/2.0,     -8.0/27.0,            2.0, -3544.0/2565.0,   1859.0/4104.0, -11.0/40.0,      0.0 },
  {       0.0,    16.0/135.0,            0.0, 6656.0/12825.0, 28561.0/56430.0,  -9.0/50.0, 2.0/55.0 },
  {       0.0,    25.0/216.0,            0.0,  1408.0/2565.0,   2197.0/4104.0,   -1.0/5.0,      0.0 }
};

#define B(row,col) (tableau[row-1][col-1])
#define FOREACH for (Index j = 0; j < N; j++)
#define TK (x[0])
#define XK (x[j+1])
#define YK (y[j])
#define K1 (k[j][0])
#define K2 (k[j][1])
#define K3 (k[j][2])
#define K4 (k[j][3])
#define K5 (k[j][4])
#define K6 (k[j][5])

template <Index N>
void srkf45(std::function<Real(Real*)> f[N],
	Real a, Real b, Real h, Real y[N], Real tol,
	std::function<void(Real*, Real*, Real)> callback) {

	Real k[N][6], x[N + 1], dx[N], r[N + 1];
	
	r[0] = a;
	for (Index i = 0; i < N; i++) r[i + 1] = y[i];
	callback(r, r, 0.0);

	while (a <= b) {
		a += h;
		r[0] = a;

		TK = a;
		FOREACH XK = YK;
		FOREACH K1 = f[j](x);

		TK = a + B(2,1)*h;
		FOREACH XK = YK + h*(B(2,2)*K1);
		FOREACH K2 = f[j](x);

		TK = a + B(3,1)*h;
		FOREACH XK = y[j] + h*(B(3,2)*K1 + B(3,3)*K2);
		FOREACH K3 = f[j](x);

		TK = a + B(4,1)*h;
		FOREACH XK = y[j] + h*(B(4,2)*K1 + B(4,3)*K2 + B(4,4)*K3);
		FOREACH K4 = f[j](x);

		TK = a + B(5,1)*h;
		FOREACH XK = y[j] + h*(B(5,2)*K1 + B(5,3)*K2 + B(5,4)*K3 + B(5,5)*K4);
		FOREACH K5 = f[j](x);

		TK = a + B(6,1)*h;
		FOREACH XK = y[j] + h*(B(6,2)*K1 + B(6,3)*K2 + B(6,4)*K3 + B(6,5)*K4 + B(6,6)*K5);
		FOREACH K6 = f[j](x);

		Real mrel = 0.0, s = 0.0;
		for (Index j = 0; j < N; j++) {
			auto rk4 = YK + h*(B(8,2)*K1 + B(8,4)*K3 + B(8,5)*K4 + B(8,6)*K5);

			dx[j] = (B(7, 2)*K1 + B(7, 4)*K3 + B(7, 5)*K4 + B(7, 6)*K5 + B(7, 7)*K6);
			auto rk5 = YK + h*dx[j];

			auto rel = abs(rk5 - rk4);
			if (rel >= mrel) {
				s = 0.84 * pow(tol*h / rel, 0.25);
				mrel = rel;
			}

			r[j + 1] = rk5;
		}

		//std::cout << "rel: " << mrel << "\ts: " << s << '\n';

		if (mrel > tol * 10.0) {
			std::cout << " ** rel > tol*10, rel " << mrel << "  tol " << tol;
			a -= h;

			std::cout << "  h: " << h;

			h *= s;
			//h *= 0.1;

			std::cout << "  adjusted h: " << h << '\n';

			continue;
		}

		for (Index j = 0; j < N; j++) y[j] = r[j + 1];

		h *= s;

		callback(r, dx, mrel);
	}
}


template <Index N>
Real InitialStepSize(std::function<Real(Real*)> f[N], Real a, Real b, Real y[N], Real tol) {
	Real x[N + 1], k[N][6];

	auto rel = [&](Real h) -> Real {
		TK = a;
		FOREACH XK = YK;
		FOREACH K1 = h*f[j](x);

		TK = a + B(2, 1)*h;
		FOREACH XK = YK + B(2, 2)*K1;
		FOREACH K2 = h*f[j](x);

		TK = a + B(3, 1)*h;
		FOREACH XK = y[j] + B(3, 2)*K1 + B(3, 3)*K2;
		FOREACH K3 = h*f[j](x);

		TK = a + B(4, 1)*h;
		FOREACH XK = y[j] + B(4, 2)*K1 + B(4, 3)*K2 + B(4, 4)*K3;
		FOREACH K4 = h*f[j](x);

		TK = a + B(5, 1)*h;
		FOREACH XK = y[j] + B(5, 2)*K1 + B(5, 3)*K2 + B(5, 4)*K3 + B(5, 5)*K4;
		FOREACH K5 = h*f[j](x);

		TK = a + B(6, 1)*h;
		FOREACH XK = y[j] + B(6, 2)*K1 + B(6, 3)*K2 + B(6, 4)*K3 + B(6, 5)*K4 + B(6, 6)*K5;
		FOREACH K6 = h*f[j](x);

		Real rel = 0.0;
		for (int j = 0; j < N; j++) {
			auto rk4 = YK + B(8, 2)*K1 + B(8, 4)*K3 + B(8, 5)*K4 + B(8, 6)*K5;
			auto rk5 = YK + B(7, 2)*K1 + B(7, 4)*K3 + B(7, 5)*K4 + B(7, 6)*K5 + B(7, 7)*K6;
			rel = fmax(rel, abs(rk5 - rk4));
		}
		return rel;
	};

	Real l = tol, c, r = (b - a);

	for (int i = 0; i < 32; i++) {
		c = (r + l)*0.5;
		if (rel(c) < tol) l = c;
		else r = c;
	}

	return c;
}

