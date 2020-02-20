// Princing.cpp : définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include <iostream>

#include "Option.h"

using namespace std;

int main()
{
	cout << "Hello World" << endl;

	Option * myOption = new Option(100,110,0.05,0.15,1,0);

	myOption->display_values();

	myOption->setCall();
	cout << "Call : " << myOption->calcul() << endl;
	myOption->greeks();
	myOption->setPut();
	cout << "Put  : " << myOption->calcul() << endl;
	myOption->greeks();

	// cout << myOption.isParityOk(true) << endl;

	delete myOption;

	system("PAUSE");

	return 0;
}
