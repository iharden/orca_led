/*
 * monomer.cpp
 *
 *  Created on: Dec 21, 2020
 *      Author: iharden
 */


#include "class_monomer.hpp"

using namespace std;
using namespace boost;

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
			split(res, line, is_any_of(" "), token_compress_on);
			nel=stoi(res[6]);
			break;
		}
	}

	// NUMBER OF BASIS FUNCTIONS
	s = "Number of basis functions";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line, is_any_of(" "), token_compress_on);
			nbasis=stoi(res[5]);
			break;
		}
	}

	// EHF
	s = "E(0)";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line, is_any_of(" "), token_compress_on);
			ehf=stod(res[2]);
			break;
		}
	}

	// E(CCSD)
	s = "E(CCSD)";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line, is_any_of(" "), token_compress_on);
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
			split(res, line, is_any_of(" "), token_compress_on);
			eccsdt=stod(res[2]);
			break;
		}
	}

	//ETRIPLES
	et = eccsdt-eccsd;
}
