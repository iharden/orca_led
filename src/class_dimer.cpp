/*
 * class_dimer.cpp
 *
 *  Created on: Dec 22, 2020
 *      Author: iharden
 */

#include "class_dimer.hpp"
#include "functions.hpp"
//#include <format>
using namespace std;

Dimer::Dimer(string_view argv) {
	ifstream dimer(string{argv});
	if(!dimer) {
		cerr << "Could not find file " << argv << endl;
		throw runtime_error("File not found");
	}

	// GET FILE NAME
	name=argv;

	vector<string> file{};
	string line{}, s{};
	vector<string_view> res{};

	while(getline(dimer, line))
		file.push_back(line);

	string_view sv;
	// IS LED PRESENT?
	sv = "LOCAL ENERGY DECOMPOSITION FOR DLPNO-CC METHODS";
	bool led=false;
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			led=true;
			break;
		}
	}
	if(led==false) {
		cerr << "File corrupted! No LED section in Dimer file found. Please try again \n";
		throw runtime_error("File corrupted");
	}

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

	// HFTYPE
	sv="Hartree-Fock type";
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			split(res, line);
			hftype=string{res[4]};
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

	//EHF USED FOR LED ANALYSIS
	sv="REFERENCE ENERGY E(0) DECOMPOSITION (Eh)";
	for(size_t i=0;i<file.size();++i) {
		if(file[i].find(sv)!=string::npos) {
			line=file[i+7];
			split(res, line);
			eref=stod(string{res[3]});
			break;
		}
	}

	sv="DLPNO BASED TRIPLES CORRECTION";
	bool triples=false;
	for(const string& line:file) {
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

	// NUMBER OF FRAGMENTS
	sv = "Number of fragments";
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			split(res, line);
			nfrag=stoi(string{res[4]});
		}
	}

	// FRAGEHF
	for(int i=0;i<nfrag;++i) {
		s=fmt::format("INTRA-FRAGMENT REF. ENERGY FOR FRAGMENT   {}", i+1);
		for(size_t j=0;j<file.size();++j) {
			if(file[j].find(s)!=string::npos) {
				line = file[j+7];
				split(res, line);
				fragehf.push_back(stod(string{res[3]}));
				break;
			}
		}
	}

	//FRAGHFJ, FRAGHFK, FRAGCCDISP
	for(int i=2;i<=nfrag;++i) {
		for(int j=1;j<i;++j) {
			s = fmt::format("Interaction of fragments  {} and  {}", i, j);
			for(size_t k=0;k<file.size();++k) {
				if(file[k].find(s)!=string::npos) {
					line=file[k+1];
					split(res, line);
					fraghfJ.push_back(stod(string{res[2]}));

					line=file[k+2];
					split(res, line);
					fraghfK.push_back(stod(string{res[2]}));

					line=file[k+3];
					split(res, line);
					fragccDispstrong.push_back(stod(string{res[3]}));

					line=file[k+4];
					split(res, line);
					fragccDispweak.push_back(stod(string{res[3]}));
					break;
				}
			}
		}
	}

	// Calculate HF interaction
	for(size_t i=0;i<fraghfK.size();++i) {
		hfint.push_back(fraghfJ[i]+fraghfK[i]);
	}

	// SUM OF NONDISPERSION (STRONG PAIRS AND WEAK PAIRS)
	sv = "Sum of non dispersive correlation terms:";
	for(size_t i=0;i<file.size();++i) {
		if(file[i].find(sv)!=string::npos) {
			line=file[i+1];
			split(res, line);
			nondispStrong = stod(string{res[4]});

			line=file[i+2];
			split(res, line);
			nondispWeak = stod(string{res[4]});
			break;
		}
	}

	// INTER-STRONG PAIRS WEAK PAIRS AND INTER TRIPLES
	for(int i=2;i<=nfrag;++i) {
		for(int j=1;j<i;++j) {
			if(hftype=="RHF")
				s = fmt::format("Interaction correlation for Fragments   {} and   {}:",i, j);
			else if(hftype=="UHF")
				s = fmt::format("Interaction correlation for Fragments   {} and  {}:",i, j);
			else
				throw runtime_error("Unknown HFType");

			for(size_t k=0;k<file.size();++k) {
				if(file[k].find(s)!=string::npos) {
					line=file[k+2];
					split(res, line);
					fragInterstrong.push_back(stod(string{res[3]}));

					int test=-1;
					if(triples) {
						test=0;
						line=file[k+3+test];
						split(res, line);
						fragIntertriples.push_back(stod(string{res[2]}));
					}
					else {
						fragIntertriples.push_back(0.0);
					}


					line=file[k+4+test];
					split(res, line);
					fragInterweak.push_back(stod(string{res[3]}));
					break;
				}
			}
		}
	}

	// CALCULATE TOTAL CC-INTERACTION
	for(size_t i=0;i<fragInterstrong.size();++i) {
		ccint.push_back(fragInterstrong[i]+fragInterweak[i]+fragIntertriples[i]);
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
	sv="INTER- vs INTRA-FRAGMENT CORRELATION ENERGIES (Eh)";
	for(size_t j=0;j<file.size();++j) {
		if(file[j].find(sv)!=string::npos) {
			line=file[j+5];
			split(res, line);
			for(int i=0;i<nfrag;++i) {
				fragIntrastrong.push_back(stod(string{res[3+i]}));
			}

			// ONLY PRESENT IF TRIPLES ARE AROUND!!!
			int test=-1;
			if(triples) {
				test=0;
				line=file[j+6];
				split(res, line);
				for(int i=0;i<nfrag;++i) {
					fragIntratriples.push_back(stod(string{res[2+i]}));
				}
			}
			else {
				for(int i=0;i<nfrag;++i) {
					fragIntratriples.push_back(0.0);
				}
			}

			line=file[j+7+test];
			split(res, line);
			for(int i=0;i<nfrag;++i) {
				fragIntraweak.push_back(stod(string{res[3+i]}));
			}

			line=file[j+8+test];
			split(res, line);
			for(int i=0;i<nfrag;++i) {
				fragIntrasingles.push_back(stod(string{res[2+i]}));
			}
			break;
		}
	}

	for(int i=0;i<nfrag;++i) {
		fragIntrages.push_back(fragIntrastrong[i]+fragIntratriples[i]+fragIntraweak[i]+fragIntrasingles[i]);
	}

	// CALCULATE NON-DISPERSIVE TRIPLES
	for(int i=0;i<nfrag;++i) {
		if(eccsdt !=0)
			fragccNonDisptriples.push_back((1-gamma[i])*fragIntertriples[i]);
		else
			fragccNonDisptriples.push_back(0.0);
	}

	// CALCULATE PAIRWISE NONDISP CONTRIBUTIONS
	for(size_t i=0;i<fragccDispges.size();++i) {
		fragccNonDispges.push_back(ccint[i]-fragccDispges[i]);
	}

	// DELOCALIZED TRIPLES
	sv="Delocalized correction (triples)";
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			split(res, line);
			delocalizedtriples=stod(string{res[3]});
			break;
		}
	}

	// DELOCALIZED Strong-pairs
	sv="Delocalized correlation";
	for(const string& line:file) {
		if(line.find(sv)!=string::npos) {
			split(res, line);
			delocalizedstrongpairs=stod(string{res[2]});
			break;
		}
	}
}
