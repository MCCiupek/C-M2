#include "stdafx.h"
#include "Option.h"

#include <iostream>

using namespace std;

Option::Option()
{
	cout << "Je suis le constructeur d'Option par défaut..." << endl;
}

Option::Option(double St, double K, double r, double sigma, double T, double t)
{
	cout << "Je suis le constructeur d'Option avec args..." << endl;
	this->St = St;
	this->K = K;
	this->r = r;
	this->sigma = sigma;
	this->T = T;
	this->t = t;

	this->isCall = true;
	this->mult = 1;

	this->areDsCalculated = false;
}

Option::~Option()
{
	cout << "Je suis le destructeur d'Option" << endl;
}

double Option::N(double x)
{
	double a1 =  0.319381530;
	double a2 = -0.356563782;
	double a3 =  1.781477937;
	double a4 = -1.821255978;
	double a5 =  1.330274429;
	double k;

	if (x >= 0.0)
	{
		k = 1 / (1 + 0.2316419*x);
		return (1 - NP(x) * ((a1*k) + (a2*k*k) + (a3*k*k*k) + (a4*k*k*k*k) + (a5*k*k*k*k*k)));
	}
	else
	{
		return (1 - N(-x));
	}
}

double Option::NP(double x)
{
	return (1.0 / sqrt(2.0 * 3.1415) * exp(-x*x*0.5));
}

double Option::calcul()
{
	// call = ( 1)*St*N( 1*d1) - ( 1)*K*exp(-r*(T-t))*N( 1*d2);
	// put  = (-1)*St*N(-1*d1) - (-1)*K*exp(-r*(T-t))*N(-1*d2);

	// call = ( 1)*(St*N( 1*d1) - K*exp(-r*(T-t))*N( 1*d2));
	// put  = (-1)*(St*N(-1*d1) - K*exp(-r*(T-t))*N(-1*d2));

	calcul_d1_d2();
	return (mult*(St*N(mult*d1) - K*exp(-r*(T - t))*N(mult*d2)));
}

void Option::calcul_d1_d2()
{
	if (this->areDsCalculated == false) {
		d1 = (log(St / K) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
		d2 = d1 - sigma*sqrt(T - t);
		this->areDsCalculated = true;
	}
}

void Option::display_values() {
	cout << "St   : " << St    << endl;
	cout << "K    : " << K     << endl;
	cout << "r    : " << r     << endl;
	cout << "sigma: " << sigma << endl;
	cout << "T    : " << T     << endl;
	cout << "t    : " << t     << endl;

	cout << "N(0) = " << N(0) << endl;
	cout << "N(1.96) = " << N(1.96) << endl;
	cout << "N(100000000) = " << N(100000000) << endl;

	calcul_d1_d2();
	cout << "d1 = " << d1 << endl;
	cout << "d2 = " << d2 << endl;
}

bool Option::isParityOk(bool display)
{
	setCall();
	double call = calcul();
	setPut();
	double put = calcul();

	double epsilon = 0.00001;

	// C - P = St - K*exp(-r(T-t))
	// C - P - St + K*exp(-r(T-t)) = 0

	if (display)
	{
		cout << "C-P     = " << call - put << endl;
		cout << "St-K... = " << St - K*exp(-r*(T - t)) << endl;
		cout << "Diff    = " << (call - put - St + K*exp(-r*(T - t))) << endl;
	}

	return ((call - put - St +K*exp(-r*(T-t))) < epsilon);
}

void Option::calcul_mult() {
	mult = -1;
	if (isCall) {
		mult = 1;
	}
}

vector<double> Option::greeks() {
	vector<double> v(5);
	// Il faut recalculer d1 et d2...
	calcul_d1_d2();

	// delta_c = N(d1)
	// delta_p = N(d1) - 1
	// delta   = N(d1) + (mult>0?0:-1)
	double delta = N(d1) + (mult > 0 ? 0 : -1);

	// rho_c   = K*(T-t)*exp(-r(T-t))*N(d2)
	// rho_p   = -K*(T - t)*exp(-r(T - t))*N(-d2)
	// rho     = mult*K*(T - t)*exp(-r(T - t))*N(mult*d2)
	double rho = mult*K*(T - t)*exp(-r*(T - t))*N(mult*d2);

	// gamma   = NP(d1)/(St*sigma*sqrt(T-t))
	double gamma = NP(d1) / (St*sigma*sqrt(T - t));

	// vega    = St * sqrt(T-t) * NP(d1)
	double vega = St * sqrt(T - t) * NP(d1);

	// theta_c = -St*NP(d1)*sigma/(2*sqrt(T-t) - r*K*exp(-r*(T-t))*N(d2)
	// theta_p = -St*NP(d1)*sigma/(2*sqrt(T-t) + r*K*exp(-r*(T-t))*N(-d2)
	// theta   = -St*NP(d1)*sigma/(2*sqrt(T-t) - mult*r*K*exp(-r(T-t))*N(mult*d2)
	double theta = ((-1*St*NP(d1)*sigma)/(2*sqrt(T-t))) - (mult*r*K*exp(-r*(T - t))*N(mult*d2));

	cout << "delta = " << delta << endl;
	cout << "gamma = " << gamma << endl;
	cout << "theta = " << theta << endl;
	cout << "vega  = " << vega << endl;
	cout << "rho   = " << rho << endl;

	v.push_back(delta);
	v.push_back(gamma);
	v.push_back(theta);
	v.push_back(vega);
	v.push_back(rho);

	// Return ???
	return v;
}