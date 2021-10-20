/*
 * class_dimer.cpp
 *
 *  Created on: Dec 22, 2020
 *      Author: iharden
 */

#include "class_dimer.hpp"
#include "functions.hpp"

using namespace std;

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

	// HFTYPE
	s="Hartree-Fock type";
	for(string line:file) {
		found=line.find(s);
		if(found!=string::npos) {
			split(res, line);
			hftype=res[4];
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

	//EHF USED FOR LED ANALYSIS
	s="REFERENCE ENERGY E(0) DECOMPOSITION (Eh)";
	for(size_t i=0;i<file.size();++i) {
		line=file[i];
		found = line.find(s);
		if(found!=string::npos) {
			line=file[i+7];
			split(res, line);
			eref=stod(res[3]);
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

	// NUMBER OF FRAGMENTS
	s = "Number of fragments";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line);
			nfrag=stoi(res[4]);
		}
	}

	// FRAGEHF
	for(int i=0;i<nfrag;++i) {
		s=fmt::format("INTRA-FRAGMENT REF. ENERGY FOR FRAGMENT   {}", i+1);
		for(size_t j=0;j<file.size();++j) {
			line=file[j];
			found=line.find(s);
			if(found!=string::npos) {
				line = file[j+7];
				split(res, line);
				fragehf.push_back(stod(res[3]));
				break;
			}
		}
	}

	//FRAGHFJ, FRAGHFK, FRAGCCDISP
	for(int i=2;i<=nfrag;++i) {
		for(int j=1;j<i;++j) {
			s = fmt::format("Interaction of fragments  {} and  {}", i, j);
			for(size_t k=0;k<file.size();++k) {
				line=file[k];
				found=line.find(s);
				if(found!=string::npos) {
					line=file[k+1];
					split(res, line);
					fraghfJ.push_back(stod(res[2]));

					line=file[k+2];
					split(res, line);
					fraghfK.push_back(stod(res[2]));

					line=file[k+3];
					split(res, line);
					fragccDispstrong.push_back(stod(res[3]));

					line=file[k+4];
					split(res, line);
					fragccDispweak.push_back(stod(res[3]));
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
	s = "Sum of non dispersive correlation terms:";
	for(size_t i=0;i<file.size();++i) {
		line=file[i];
		found=line.find(s);
		if(found!=string::npos) {
			line=file[i+1];
			split(res, line);
			nondispStrong = stod(res[4]);

			line=file[i+2];
			split(res, line);
			nondispWeak = stod(res[4]);
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
				line=file[k];
				found=line.find(s);
				if(found!=string::npos) {
					line=file[k+2];
					split(res, line);
					fragInterstrong.push_back(stod(res[3]));

					int test=-1;
					if(triples) {
						test=0;
						line=file[k+3+test];
						split(res, line);
						fragIntertriples.push_back(stod(res[2]));
					}
					else {
						fragIntertriples.push_back(0.0);
					}


					line=file[k+4+test];
					split(res, line);
					fragInterweak.push_back(stod(res[3]));
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
	s="INTER- vs INTRA-FRAGMENT CORRELATION ENERGIES (Eh)";
	for(size_t j=0;j<file.size();++j) {
		line=file[j];
		found=line.find(s);
		if(found!=string::npos) {
			line=file[j+5];
			split(res, line);
			for(int i=0;i<nfrag;++i) {
				fragIntrastrong.push_back(stod(res[3+i]));
			}

			// ONLY PRESENT IF TRIPLES ARE AROUND!!!
			int test=-1;
			if(triples) {
				test=0;
				line=file[j+6];
				split(res, line);
				for(int i=0;i<nfrag;++i) {
					fragIntratriples.push_back(stod(res[2+i]));
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
				fragIntraweak.push_back(stod(res[3+i]));
			}

			line=file[j+8+test];
			split(res, line);
			for(int i=0;i<nfrag;++i) {
				fragIntrasingles.push_back(stod(res[2+i]));
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
	s="Delocalized correction (triples)";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line);
			delocalizedtriples=stod(res[3]);
			break;
		}
	}

	// DELOCALIZED Strong-pairs
	s="Delocalized correlation";
	for(string line:file) {
		found = line.find(s);
		if(found!=string::npos) {
			split(res, line);
			delocalizedstrongpairs=stod(res[2]);
			break;
		}
	}
}
