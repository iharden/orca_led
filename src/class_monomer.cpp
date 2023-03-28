/*
 * monomer.cpp
 *
 *  Created on: Dec 21, 2020
 *      Author: iharden
 */


#include "class_monomer.hpp"
#include "functions.hpp"
#include "fmt/format.h"
#include "fmt/core.h"
#include "fmt/format-inl.h"
#include "fmt/printf.h"
#include "fmt/ostream.h"

using namespace std;

Monomer::Monomer(string_view argv) {
	ifstream monomer(string{argv});
	if(!monomer) {
		cerr << "Could not find file " << argv << endl;
	    throw runtime_error("File not found");
	}
	name=argv;

	vector<string> file{};
	string s{};
	vector<string_view> res{};
	string_view sv{};

	while(getline(monomer, s))
		file.push_back(s);

	// NUMBER OF ELECTRONS
	sv = "Number of Electrons";
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			split(res, line);
			nel=stoi(string{res[5]});
			break;
		}
	}

	// NUMBER OF BASIS FUNCTIONS
	sv = "Number of basis functions";
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			split(res, line);
			nbasis=stoi(string{res[5]});
			break;
		}
	}

	// EHF
	sv = "E(0)";
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			split(res, line);
			ehf=stod(string{res[2]});
			break;
		}
	}

	sv="DLPNO BASED TRIPLES CORRECTION";
	bool triples=false;
	for(string line:file) {
		if(line.find(sv)!=string::npos) {
			triples=true;
			break;
		}
	}

	if(triples==true) {
	// E(CCSD)
		sv = "E(CCSD)";
		for(const string& line:file) {
			if(line.find(sv)!=string::npos) {
				split(res, line);
				eccsd=stod(string{res[2]});
				break;
			}
		}

		//ECORR
		ecorr=eccsd-ehf;

		// E(CCSD(T))
		sv = "E(CCSD(T))";
		for(const string& line:file) {
			if(line.find(sv)!=string::npos) {
				split(res, line);
				eccsdt=stod(string{res[2]});
				break;
			}
		}

		//ETRIPLES
		et = eccsdt-eccsd;

		ecorrt = ecorr + et;
	}

	else {
		// GET E(CCSD) IF NO TRIPLES ARE CALCULATED
		sv="E(TOT)";
		for(const string& line:file) {
			if(line.find(sv)!=string::npos) {
				split(res, line);
				eccsd=stod(string{res[2]});
				break;
			}
		}
		ecorr=eccsd-ehf;
		eccsdt=eccsd;
		et=eccsdt-eccsd;
		ecorrt=ecorr+et;
	}
}
