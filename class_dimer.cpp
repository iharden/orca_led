/*
 * class_dimer.cpp
 *
 *  Created on: Dec 22, 2020
 *      Author: iharden
 */

#include "class_dimer.hpp"

using namespace std;
using namespace boost;

Dimer::Dimer(string argv) {
	ifstream dimer(argv);
	if(!dimer) {
		cerr << "Could not find file " << argv << endl;
		throw runtime_error("File not found");
	}

	// GET FILE NAME
	name=argv;

	vector<string> file{};
	string s{}, line{};
	vector<string> res{};
	size_t found{};
	int ifrag{};

	while(getline(dimer, s))
		file.push_back(s);

	// IS LED PRESENT?
	s = "LOCAL ENERGY DECOMPOSITION FOR DLPNO-CC METHODS";
	bool led=false;
	for(string line:file) {
		found=line.find(s);
		if(found!=string::npos) {
			led=true;
			break;
		}
	}
	if(led==false) {
		cout << "File corrupted! No LED section in Dimer file found. Please try again \n";
		throw runtime_error("File corrupted");
	}

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
	s = "of contracted basis functions";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line, is_any_of(" "), token_compress_on);
			nbasis=stoi(res[7]);
			break;
		}
	}

	// HFTYPE
	s="Hartree-Fock type";
	for(string line:file) {
		found=line.find(s);
		if(found!=string::npos) {
			split(res, line, is_any_of(" "), token_compress_on);
			hftype=res[5];
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

	//EHF USED FOR LED ANALYSIS
	s="REFERENCE ENERGY E(0) DECOMPOSITION (Eh)";
	for(size_t i=0;i<file.size();++i) {
		line=file[i];
		found = line.find(s);
		if(found!=string::npos) {
			line=file[i+7];
			split(res, line, is_any_of(" "), token_compress_on);
			eref=stod(res[3]);
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

	// NUMBER OF FRAGMENTS
	s = "Number of fragments";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line, is_any_of(" "), token_compress_on);
			nfrag=stoi(res[4]);
		}
	}

	// FRAGEHF
	for(int i=0;i<nfrag;++i) {
		format frmt("INTRA-FRAGMENT REF. ENERGY FOR FRAGMENT   %1%");
		ifrag=i+1;
		frmt % ifrag;
		s = frmt.str();
		for(size_t j=0;j<file.size();++j) {
			line=file[j];
			found=line.find(s);
			if(found!=string::npos) {
				line = file[j+7];
				split(res, line, is_any_of(" "), token_compress_on);
				fragehf.push_back(stod(res[3]));
				break;
			}
		}
	}

	//FRAGHFJ, FRAGHFK, FRAGCCDISP
	for(int i=2;i<=nfrag;++i) {
		for(int j=1;j<i;++j) {
			s = str(format("Interaction of fragments  %1% and  %2%") % i % j );
			for(size_t k=0;k<file.size();++k) {
				line=file[k];
				found=line.find(s);
				if(found!=string::npos) {
					line=file[k+1];
					split(res, line, is_any_of(" "), token_compress_on);
					fraghfJ.push_back(stod(res[2]));

					line=file[k+2];
					split(res, line, is_any_of(" "), token_compress_on);
					fraghfK.push_back(stod(res[2]));

					line=file[k+3];
					split(res, line, is_any_of(" "), token_compress_on);
					fragccDispstrong.push_back(stod(res[3]));

					line=file[k+4];
					split(res, line, is_any_of(" "), token_compress_on);
					fragccDispweak.push_back(stod(res[3]));
					break;
				}
			}
		}
	}

	// SUM OF NONDISPERSION (STRONG PAIRS AND WEAK PAIRS)
	s = "Sum of non dispersive correlation terms:";
	for(size_t i=0;i<file.size();++i) {
		line=file[i];
		found=line.find(s);
		if(found!=string::npos) {
			line=file[i+1];
			split(res, line, is_any_of(" "), token_compress_on);
			nondispStrong = stod(res[4]);

			line=file[i+2];
			split(res, line, is_any_of(" "), token_compress_on);
			nondispWeak = stod(res[4]);
			break;
		}
	}

	// INTER-STRONG PAIRS WEAK PAIRS AND INTER TRIPLES
	for(int i=2;i<=nfrag;++i) {
		for(int j=1;j<i;++j) {
			if(hftype=="RHF")
				s = str(format("Interaction correlation for Fragments   %1% and   %2%:") % i % j );
			else if(hftype=="UHF")
				s = str(format("Interaction correlation for Fragments   %1% and  %2%:") % i % j );S
			else
				throw runtime_error("Unknown HFType");

			for(size_t k=0;k<file.size();++k) {
				line=file[k];
				found=line.find(s);
				if(found!=string::npos) {
					line=file[k+2];
					split(res, line, is_any_of(" "), token_compress_on);
					fragInterstrong.push_back(stod(res[3]));

					line=file[k+3];
					split(res, line, is_any_of(" "), token_compress_on);
					fragIntertriples.push_back(stod(res[2]));

					line=file[k+4];
					split(res, line, is_any_of(" "), token_compress_on);
					fragInterweak.push_back(stod(res[3]));
					break;
				}
			}
		}
	}

	// CALCULATE GAMMA
	if(fragccDispstrong.size()!=fragInterstrong.size()) {
		cerr << "Number of Dispersion Strong pairs does not match number of Inter Strong pairs \n";
		throw logic_error("Vectors do not have the same size");
	}

	for(size_t i=0;i<fragccDispstrong.size();++i) {
		if(fragInterstrong[i] != 0)
			gamma.push_back(fragccDispstrong[i]/fragInterstrong[i]);
		else
			gamma.push_back(0.0);
	}

	// CALCULATE DISPERSION TRIPLES
	for(size_t i=0;i<fragIntertriples.size();++i) {
		fragccDisptriples.push_back(gamma[i]*fragIntertriples[i]);
	}

	// CALCULATE TOTAL DISPERSION CONTRIBUTION
	for(size_t i=0;i<fragccDisptriples.size();++i) {
		fragccDispges.push_back(fragccDispstrong[i]+fragccDispweak[i]+fragccDisptriples[i]);
	}

	// GET ALL THE INTRA CONTRIBUTIONS (STRONG PAIRS, WEAK PAIRS, TRIPLES AND SINGLES)
	s="INTER- vs INTRA-FRAGMENT CORRELATION ENERGIES (Eh)";
	for(size_t j=0;j<file.size();++j) {
		line=file[j];
		found=line.find(s);
		if(found!=string::npos) {
			line=file[j+5];
			split(res, line, is_any_of(" "), token_compress_on);
			for(int i=0;i<nfrag;++i) {
				fragIntrastrong.push_back(stod(res[3+i]));
			}

			line=file[j+6];
			split(res, line, is_any_of(" "), token_compress_on);
			for(int i=0;i<nfrag;++i) {
				fragIntratriples.push_back(stod(res[2+i]));
			}

			line=file[j+7];
			split(res, line, is_any_of(" "), token_compress_on);
			for(int i=0;i<nfrag;++i) {
				fragIntraweak.push_back(stod(res[3+i]));
			}

			line=file[j+8];
			split(res, line, is_any_of(" "), token_compress_on);
			for(int i=0;i<nfrag;++i) {
				fragIntrasingles.push_back(stod(res[2+i]));
			}
			break;
		}
	}
}
