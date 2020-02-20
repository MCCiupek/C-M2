#ifndef OPTION_272
#define OPTION_272

#include <vector>

using namespace std;

class Option {
	private:
		double d1;
		double d2;

		double St;
		double K;
		double r;
		double sigma;
		double T;
		double t;

		bool isCall;
		double mult;

		bool areDsCalculated;

	public:
		// Constructeurs
		Option();
		Option(double St, double K, double r, double sigma,
			   double T, double t);
		// Destructeur
		~Option();

		// ~~~ Setters ~~~
		void setCall() { isCall = true; calcul_mult(); }
		void setPut() { isCall = false; calcul_mult(); }

		void setSt(double inSt) { this->St = inSt; this->areDsCalculated = false; }
		void setK(double inK) { this->K = inK; this->areDsCalculated = false; }
		void setr(double inr) { this->r = inr; this->areDsCalculated = false; }
		void setT(double inT) { this->T = inT; this->areDsCalculated = false; }
		void sett(double inT) { this->t = inT; this->areDsCalculated = false; }
		void setSigma(double inSigma) { this->sigma = inSigma; this->areDsCalculated = false; }
		// ~~~ Setters ~~~

		double N (double x);
		double NP(double x);

		double calcul();
		void calcul_d1_d2();

		void display_values();

		bool isParityOk(bool display = false);

		vector<double> greeks();

	private:
		void calcul_mult();
};

#endif // !OPTION_272

