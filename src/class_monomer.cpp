/*
 * monomer.cpp
 *
 *  Created on: Dec 21, 2020
 *      Author: iharden
 */


#include "class_monomer.hpp"
#include "functions.hpp"

using namespace std;

Monomer::Monomer(string argv) {
	ifstream monomer(argv);
	if(!monomer) {
		cerr << "Could not find file " << argv << endl;
	    throw runtime_error("File not found");
	}
	name=argv;

	vector<string> file{};
	string s{};
	vector<string> res{};
	size_t found{};

	while(getline(monomer, s))
		file.push_back(s);

	// NUMBER OF ELECTRONS
	s = "Number of Electrons";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line);
			nel=stoi(res[5]);
			break;
		}
	}

	// NUMBER OF BASIS FUNCTIONS
	s = "Number of basis functions";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line);
			nbasis=stoi(res[5]);
			break;
		}
	}

	// EHF
	s = "E(0)";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line);
			ehf=stod(res[2]);
			break;
		}
	}

	s="DLPNO BASED TRIPLES CORRECTION";
	bool triples=false;
	for(string line:file) {
		found=line.find(s);
		if(found!=string::npos) {
			triples=true;
			break;
		}
	}

	if(triples==true) {
	// E(CCSD)
		s = "E(CCSD)";
		for(string line:file) {
			found = line.find(s);
			if(found!=string::npos) {
				split(res, line);
				eccsd=stod(res[2]);
				break;
			}
		}

		//ECORR
		ecorr=eccsd-ehf;

		// E(CCSD(T))
		s = "E(CCSD(T))";
		for(string line:file) {
			found = line.find(s);
			if(found!=string::npos) {
				split(res, line);
				eccsdt=stod(res[2]);
				break;
			}
		}

		//ETRIPLES
		et = eccsdt-eccsd;

		ecorrt = ecorr + et;
	}

	else {
		// GET E(CCSD) IF NO TRIPLES ARE CALCULATED
		s="E(TOT)";
		for(string line:file) {
			found=line.find(s);
			if(found!=string::npos) {
				split(res, line);
				eccsd=stod(res[2]);
				break;
			}
		}
		ecorr=eccsd-ehf;
		eccsdt=eccsd;
		et=eccsdt-eccsd;
		ecorrt=ecorr+et;
	}
}
