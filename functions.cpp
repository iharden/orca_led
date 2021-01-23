/*
 * functions.cpp
 *
 *  Created on: Dec 24, 2020
 *      Author: iharden
 */

#include "functions.hpp"

using namespace std;
using namespace boost;

namespace {
	constexpr double CF=627.509608; //Eh to kcal/mol
	string fs1="%-15s %5d \n";
	string fs3="%-15s %5s \n";
	string fs2="%s %s \n";

	string fs4="%-15s %5.5f \n";  // general infos
	string fs5="%-70s %5.5f \n";  // mostly used for OS

	string fs6="%-20s %5.5f \n";  // print EInt for CSV

	string fs7="%-50s %12.5f %12.5f %12.5f \n"; // compare OS
	string fs8="%-50s, %12.5f, %12.5f, %12.5f \n"; // compare CSV
}

double get_time(chrono::time_point<chrono::high_resolution_clock>  start, chrono::time_point<chrono::high_resolution_clock> finish) {
	chrono::duration<double, milli> elapsed = finish - start;
	return elapsed.count();
}

void do_led(const Dimer& dim, ostream& os, ostream& csv) {

	print_header(os);
	do_generalinfo(dim, os, csv);
	do_hartreefock(dim, os, csv);
	do_intra(dim, os, csv);
	do_ccsddisp(dim, os, csv);
	do_triplesdisp(dim, os, csv);
	do_ccsdnondisp(dim, os, csv);
}

void do_led(const Dimer& dim, vector<Monomer>& mons, ostream& os, ostream& csv) {

	map<string, double> summary{};
	vector<string> insertOrder{};

	do_led(dim, os, csv);
	do_generalinfo(dim, mons, os);
	do_geoprep(dim, mons, os, summary, insertOrder, csv);
	do_hfint(dim, mons, os, summary, insertOrder, csv);
	do_ccsdint(dim, mons, os, summary, insertOrder, csv);
	do_triplesint(dim, mons, os, summary, insertOrder, csv);
	do_consistency(dim, mons, os, summary, insertOrder, csv);
	do_summary(os, summary, insertOrder, csv);
}

void do_compare(const vector<string> comps, ostream& os, ostream& csv) {
	if(comps.size()!=2) {
		cout << "You have to specify exactly two .led files \n";
		throw runtime_error("Number of arguments is wrong");
	}

	ifstream f1{comps[0]};
	ifstream f2{comps[1]};

	if(!f1 || !f2)
		throw runtime_error("do_compare():File not found");

	print_header(os);

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "Name of files \n");
	for(size_t i=0;i<comps.size();++i) {
		fmt::fprintf(os, "   %-s \n", comps[i]);
	}
	fmt::fprintf(os, "\n");


	vector<string> v1{}, v2{}, res{};
	vector<string> insertOrder{};
	string s{}, line{};

	while(getline(f1, s))
		v1.push_back(s);

	while(getline(f2, s))
		v2.push_back(s);


	size_t found{};

	s="E(Geo-prep)";
	vector<bool> geoprep{false,false};
	for(string line:v1) {
		found=line.find(s);
		if(found!=string::npos) {
			geoprep[0]=true;
			break;
		}
	}

	for(string line:v2) {
		found=line.find(s);
		if(found!=string::npos) {
			geoprep[1]=true;
			break;
		}
	}

	if(geoprep[0] != geoprep[1]) {
		cerr << "*** \n";
		cerr << "WARNING: E(GEO-PREP) is present in one file but not in the other \n";
		cerr << "Please repeat the LED analysis and be consistent this time \n";
		cerr << "I will abbort \n";
		cerr << "*** \n";
		cerr << " \n";
		throw runtime_error("do_compare():Inconsistency found");
	}


	int nvalue{};
	if(geoprep[0]==false)
		nvalue=5;
	else
		nvalue=6;

	map<string, double> sum1{}, sum2{};
	s="Summary";
	for(size_t i=0;i<v1.size();++i) {
		line=v1[i];
		found=line.find(s);
		if(found!=string::npos) {
			for(int j=0;j<nvalue;++j) {
				line=v1[i+3+j];
				split(res, line, is_any_of(" "), token_compress_on);
				sum1[res[0]]=stod(res[1]);
				insertOrder.push_back(res[0]);

			}
			break;
		}
	}

	for(size_t i=0;i<v2.size();++i) {
		line=v2[i];
		found=line.find(s);
		if(found!=string::npos) {
			for(int j=0;j<nvalue;++j) {
				line=v2[i+3+j];
				split(res, line, is_any_of(" "), token_compress_on);
				sum2[res[0]]=stod(res[1]);
			}
			break;
		}
	}

	fmt::fprintf(os, "%-50s %12s %12s %12s \n", " ", "File 1", "File 2", "Difference");
	fmt::fprintf(csv, "%-50s %12s %12s %12s \n", " ,", "File 1,", "File 2,", "Difference");
	fmt::fprintf(os, "\n");
	for(size_t i=0;i<insertOrder.size();++i) {
		line=insertOrder[i];
		fmt::fprintf(os, fs7, line, sum1[line], sum2[line], sum2[line]-sum1[line]);
		fmt::fprintf(csv, fs8, line, sum1[line], sum2[line], sum2[line]-sum1[line]);
	}

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
}

void print_header(ostream& os) {
	time_t t = time(0);
	tm* now=localtime(&t);

	fmt::fprintf(os, "                  ************************************************************** \n");
	fmt::fprintf(os, "                  ***                                                        *** \n");
	fmt::fprintf(os, "                  ***                        ORCA_LED                        *** \n");
	fmt::fprintf(os, "                  ***                      iharden 01/21                     *** \n");
	fmt::fprintf(os, "                  ***                                                        *** \n");
	fmt::fprintf(os, "                  ***                                                        *** \n");
	fmt::fprintf(os, "                  ************************************************************** \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "*** \n");
	fmt::fprintf(os, "The purpose of this little program is to facilitate the LED analysis as implemented in ORCA \n");
	fmt::fprintf(os, "ORCA_LED does not help you with the calculations. It works only with ORCA output files \n");
	fmt::fprintf(os, "Its main purpose is to print the essential results of LED + it performs the decomposition of the interaction energy \n");
	fmt::fprintf(os, "However, it might be useful when you want to compare the outcomes of many LED analyzes \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "*** \n");
	fmt::fprintf(os, "Some general remarks:");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "Without using the compare-mode (--compare) I will only calculate the interaction energies \n");
	fmt::fprintf(os, "To include the electronic preparation I need DLPNO-CCSD(T) calculations of the fragments in the dimer geometry \n");
	fmt::fprintf(os, "To include the geometric preparation I need DLPNO-CCSD(T) calculations of the fragments in their equilibrium geometries \n");
	fmt::fprintf(os, "Type 'orca_led --help' to get a list of options and some command line examples \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");


	os << "Program was started at " << put_time(now, "%c %Z") << "\n";
	fmt::fprintf(os, "\n");
}

void do_generalinfo(const Dimer& dim, ostream& os, ostream& csv) {

	fmt::fprintf(os, "Name of Dimer-file \n");
	fmt::fprintf(os, "   %-s \n", dim.name);
	fmt::fprintf(os, "\n");

	fmt::fprintf(os, "******* GENERAL INFORMATION ABOUT THE DIMER SYSTEM ******* \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, fs1, "nel:", dim.nel);
	fmt::fprintf(os, fs1, "nbasis:", dim.nbasis);
	fmt::fprintf(os, fs1, "nfragment:", dim.nfrag);
	fmt::fprintf(os, fs3, "HFType:", dim.hftype);
	fmt::fprintf(os, "\n");

	fmt::fprintf(csv, fs2, "name:,", dim.name);
	fmt::fprintf(csv, fs1, "nel:,", dim.nel);
	fmt::fprintf(csv, fs1, "nbasis:,", dim.nbasis);
	fmt::fprintf(csv, fs1, "nfragment:,", dim.nfrag);
	fmt::fprintf(csv, fs2, "HFType:,", dim.hftype);
	fmt::fprintf(csv, "\n");
}

void do_hartreefock(const Dimer& dim, ostream& os, ostream& csv) {

	string tmp{};
	int N=0;

	//CHECK IF E(HF) AND E(REF) ARE EQUAL
	if(abs(dim.ehf-dim.eref) > 0.0001) {
		fmt::fprintf(os, "WARNING: The Hartree-Fock energy (E(HF)) differs from the HF-Reference energy E(Ref) used in LED! \n");
		fmt::fprintf(os, "Typically, this happens because LED makes use of RI-JK while the initial HF calculation might not \n");
		fmt::fprintf(os, "I will print both energies but for the consistency check I am using E(Ref) \n");
		fmt::fprintf(os, "\n");
	}

	fmt::fprintf(os, "******* OVERALL ENERGIES ******* \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, fs4, "E(HF):", dim.ehf);
	fmt::fprintf(os, fs4, "E(Ref):", dim.eref);
	fmt::fprintf(os, fs4, "E(CCSD):", dim.eccsd);
	fmt::fprintf(os, fs4, "E(CCSD(T)):", dim.eccsdt);
	fmt::fprintf(os, fs4, "E(CORR):", dim.ecorr);
	fmt::fprintf(os, fs4, "E(Triples):", dim.et);
	fmt::fprintf(os, "\n");

	fmt::fprintf(csv, fs4, "E(HF):,", dim.ehf);
	fmt::fprintf(csv, fs4, "E(Ref):,", dim.eref);
	fmt::fprintf(csv, fs4, "E(CCSD):,", dim.eccsd);
	fmt::fprintf(csv, fs4, "E(CCSD(T)):,", dim.eccsdt);
	fmt::fprintf(csv, fs4, "E(CORR):,", dim.ecorr);
	fmt::fprintf(csv, fs4, "E(Triples):,", dim.et);
	fmt::fprintf(csv, "\n");

	fmt::fprintf(os, "******* DECOMPOSITION OF HARTREE-FOCK REFERENCE ENERGY ******* \n");
	fmt::fprintf(os, "\n");
	for(int i=0;i<dim.nfrag;++i) {
		tmp=fmt::format("E(HF) for Fragment {}:", i+1);
		fmt::fprintf(os, fs5, tmp, dim.fragehf[i]);

		tmp=fmt::format("E(HF)_{}:,", i+1);
		fmt::fprintf(csv, fs4, tmp, dim.fragehf[i]);
	}
	fmt::fprintf(os, "\n");
	fmt::fprintf(csv, "\n");

	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("Electrostatic Interaction between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.fraghfJ[N]);

			tmp=fmt::format("ElStat_{}_{}:,", i, j);
			fmt::fprintf(csv, fs4, tmp, dim.fraghfJ[N]);
			++N;
		}
	}
	fmt::fprintf(os, "\n");

	N=0;
	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("Exchange Interaction between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.fraghfK[N]);

			tmp=fmt::format("Exchange_{}_{}:,", i, j);
			fmt::fprintf(csv, fs4, tmp, dim.fraghfK[N]);
			++N;
		}
	}
	fmt::fprintf(csv, "\n");

	double sum=0.0;
	for(size_t i=0;i<dim.fragehf.size();++i)
		sum += dim.fragehf[i];
	for(size_t i=0;i<dim.fraghfJ.size();++i) {
		sum += dim.fraghfJ[i];
		sum += dim.fraghfK[i];
	}
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "Consistency check:\n");
	fmt::fprintf(os, fs5, "Total Reference energy:", dim.eref);
	fmt::fprintf(os, fs5, "Sum of the fragments HF reference and interaction energies:", sum);
	fmt::fprintf(os, "\n");
}

void do_intra(const Dimer& dim, std::ostream& os, ostream& csv) {

	string tmp{};

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "******* DECOMPOSITION OF INTRA-EXCITATIONS ******* \n");
	fmt::fprintf(os, "\n");

	os << fixed;
	os << left << setw(40) << " ";
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(9) << "Fragment" << setw(5) << i+1;
	}
	os << "\n";
	os << "\n";

	os << left << setw(40) << "Intra strong pairs";
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << dim.fragIntrastrong[i];
	}
	os << "\n";

	os << left << setw(40) << "Intra triples";
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << dim.fragIntratriples[i];
	}
	os << "\n";

	os << left << setw(40) << "Intra weak pairs";
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << dim.fragIntraweak[i];
	}
	os << "\n";

	os << left << setw(40) << "Intra singles";
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << dim.fragIntrasingles[i];
	}
	os << "\n";
	os << "\n";

	for(int i=0;i<dim.nfrag;++i) {
		tmp=fmt::format("IntraSP_{}:,", i+1);
		fmt::fprintf(csv, fs4, tmp, dim.fragIntrastrong[i]);
		tmp=fmt::format("IntraTrip_{}:,", i+1);
		fmt::fprintf(csv, fs4, tmp, dim.fragIntratriples[i]);
		tmp=fmt::format("IntraWP_{}:,", i+1);
		fmt::fprintf(csv, fs4, tmp, dim.fragIntraweak[i]);
		tmp=fmt::format("IntraSing_{}:,", i+1);
		fmt::fprintf(csv, fs4, tmp, dim.fragIntrasingles[i]);
	}
	csv << "\n";
}

void do_ccsddisp(const Dimer& dim, std::ostream& os, ostream& csv) {
	string tmp{};
	int N=0;

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "******* DECOMPOSITION OF CCSD-Dispersion ******* \n");
	fmt::fprintf(os, "\n");

	N=0;
	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("Dispersion (strong pairs) between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.fragccDispstrong[N]);

			tmp=fmt::format("DispSP_{}_{}:,", i, j);
			fmt::fprintf(csv, fs4, tmp, dim.fragccDispstrong[N]);
			++N;
		}
	}
	fmt::fprintf(os, "\n");

	N=0;
	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("Dispersion (weak pairs) between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.fragccDispweak[N]);

			tmp=fmt::format("DispWP_{}_{}:,", i, j);
			fmt::fprintf(csv, fs4, tmp, dim.fragccDispweak[N]);
			++N;
		}
	}
	fmt::fprintf(os, "\n");

	double sum=0.0;
	for(size_t i=0;i<dim.fragccDispstrong.size();++i) {
		sum += dim.fragccDispstrong[i];
		sum += dim.fragccDispweak[i];
	}
	fmt::fprintf(os, fs5, "Sum of CCSD dispersion terms:", sum);
	fmt::fprintf(os, "\n");

	fmt::fprintf(csv, fs4, "SumDisp:,", sum);
	fmt::fprintf(csv, "\n");
}

void do_triplesdisp(const Dimer& dim, std::ostream& os, ostream& csv) {
	string tmp{};
	int N=0;

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "******* Triples correction to Dispersion ******* \n");
	fmt::fprintf(os, "\n");

	N=0;
	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("gamma (Disp(SP)/Inter(SP)) between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.gamma[N]);
			++N;
		}
	}
	fmt::fprintf(os, "\n");

	N=0;
	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("Triples correction to Dispersion between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.fragccDisptriples[N]);

			tmp=fmt::format("DispTrip_{}_{}:,", i, j);
			fmt::fprintf(csv, fs4, tmp, dim.fragccDisptriples[N]);
			++N;
		}
	}
	fmt::fprintf(os, "\n");
	fmt::fprintf(csv, "\n");
}

void do_ccsdnondisp(const Dimer& dim, std::ostream& os, ostream& csv) {

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "******* DECOMPOSITION OF CCSD-Nondispersion ******* \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "Warning: NonDisp (Strong Pairs) includes Singles excitations! \n");
	fmt::fprintf(os, fs5, "Sum of Nondispersion(strong pairs):", dim.nondispStrong);
	fmt::fprintf(os, fs5, "Sum of Nondispersion(weak pairs):", dim.nondispWeak);
	fmt::fprintf(os, fs5, "Sum of Nondispersion terms:", dim.nondispStrong+dim.nondispWeak);

	fmt::fprintf(csv, fs4, "NonDispSP:,", dim.nondispStrong);
	fmt::fprintf(csv, fs4, "NonDispWP:,", dim.nondispWeak);
	fmt::fprintf(csv, fs4, "NonDispSum:,", dim.nondispStrong+dim.nondispWeak);
	fmt::fprintf(csv, "\n");
}

void do_generalinfo(const Dimer& dim, vector<Monomer>& mons, ostream& os) {

	string tmp{}, tmpp{};
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "************************************************************** \n");
	fmt::fprintf(os, "******* I will start now with the Interaction Energies ******* \n");
	fmt::fprintf(os, "************************************************************** \n");
	fmt::fprintf(os, "\n");

	if(dim.nfrag ==2)
		tmp = "Dimer";
	else if(dim.nfrag==3)
		tmp = "Trimer";
	else if (dim.nfrag==4)
		tmp = "Quatromer";
	else
		tmp = "Oligomer";

	tmpp=fmt::format("You gave me {} monomer files and {} fragments are present in the {}-file! \n", mons.size(), dim.nfrag, tmp);
	fmt::fprintf(os, tmpp);

	tmpp=fmt::format("Name of {}-file \n", tmp);
	fmt::fprintf(os, tmpp);
	fmt::fprintf(os, "   %-s \n", dim.name);

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "Name of Monomer-files \n");
	for(size_t i=0;i<mons.size();++i) {
		fmt::fprintf(os, "   %-s \n", mons[i].name);
	}

	os << fixed;
	os << left << setw(30) << " ";
	os << left << setw(14) << tmp;
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(9) << "Monomer" << setw(5) << i+1;
	}
	os << "\n";
	os << "\n";

	os << left << setw(30) << "E(HF):";
	os << left << setw(14) << setprecision(5) << dim.ehf;
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << mons[i].ehf;
	}
	os << "\n";

	os << left << setw(30) << "E(CCSD):";
	os << left << setw(14) << setprecision(5) << dim.eccsd;
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << mons[i].eccsd;
	}
	os << "\n";

	os << left << setw(30) << "E(CCSD(T)):";
	os << left << setw(14) << setprecision(5) << dim.eccsdt;
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << mons[i].eccsdt;
	}
	os << "\n";

	os << left << setw(30) << "E(CORR):";
	os << left << setw(14) << setprecision(5) << dim.ecorr;
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << mons[i].ecorr;
	}
	os << "\n";

	os << left << setw(30) << "E(Triples):";
	os << left << setw(14) << setprecision(5) << dim.et;
	for(int i=0;i<dim.nfrag;++i) {
		os << left << setw(14) <<  setprecision(5) << mons[i].et;
	}
	os << "\n";
	os << "\n";

	if(mons.size() == dim.nfrag) {
		fmt::fprintf(os, "I will calculate the electronic preparation \n");
		fmt::fprintf(os, "\n");
	}

	if(mons.size() == 2*dim.nfrag) {
		fmt::fprintf(os, "I will calculate the geometric and the electronic preparation \n");
		fmt::fprintf(os, "\n");
		// CHECK IF YOU MESSED UP ORDER
		for(int i=0;i<dim.nfrag;++i) {
			if(mons[i].eccsdt < mons[dim.nfrag+i].eccsdt) {
				tmp=fmt::format("Fragment {} in dimer has lower energy than in equilibrium geometry. I will therefore swap them \n", i+1);
				fmt::fprintf(os, tmp);
				swap(mons[i],mons[i+dim.nfrag]);
			}
		}
	}

	fmt::fprintf(os, "If not stated otherwise, energies from now on are in kcal/mol !!! \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
}


void do_geoprep(const Dimer& dim, vector<Monomer>& mons, ostream& os, map<string, double>& summary, vector<string>& insertOrder, ostream& csv) {

	double sum=0.0;
	string tmp{};

	if(mons.size() == 2*dim.nfrag) {
		fmt::fprintf(os, "******* Geometric preparation energy ******* \n");
		fmt::fprintf(os, "\n");
		sum=0.0;
		for(size_t i=0;i<dim.nfrag;++i) {
			sum+=(mons[i].eccsdt-mons[dim.nfrag+i].eccsdt)*CF;
			tmp=fmt::format("Geometric preparation for Fragment {}:", i+1);
			fmt::fprintf(os, fs5, tmp, (mons[i].eccsdt-mons[dim.nfrag+i].eccsdt)*CF);

			tmp=fmt::format("GeoPrep_{}:,", i+1);
			fmt::fprintf(csv, fs6, tmp, (mons[i].eccsdt-mons[dim.nfrag+i].eccsdt)*CF);
		}
		summary["E(Geo-prep)"]=sum;
		insertOrder.push_back("E(Geo-prep)");

		fmt::fprintf(os, "\n");
		fmt::fprintf(os, "\n");
		fmt::fprintf(os, "\n");
		fmt::fprintf(csv, "\n");
	}
}

void do_hfint(const Dimer& dim, vector<Monomer>& mons, ostream& os, map<string, double>& summary, vector<string>& insertOrder, ostream& csv) {
	string tmp{};

	int N=0;
	double sum=0.0;

	fmt::fprintf(os, "******* HF-Interaction energy ******* \n");
	fmt::fprintf(os, "\n");
	for(int i=0;i<dim.nfrag;++i) {
		tmp=fmt::format("Electronic preparation for Fragment {}:", i+1);
		fmt::fprintf(os, fs5, tmp, (dim.fragehf[i]-mons[i].ehf)*CF);

		tmp=fmt::format("ElPrep_{}:,", i+1);
		fmt::fprintf(csv, fs6, tmp, (dim.fragehf[i]-mons[i].ehf)*CF);
	}
	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("Electrostatic Interaction between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.fraghfJ[N]*CF);

			tmp=fmt::format("ElStat_{}_{}:,", i, j);
			fmt::fprintf(csv, fs6, tmp, dim.fraghfJ[N]*CF);
			++N;
		}
	}
	N=0;
	for(int i=2;i<=dim.nfrag;++i) {
		for(int j=1;j<i;++j) {
			tmp=fmt::format("Exchange Interaction between fragments {} and {}:", i, j);
			fmt::fprintf(os, fs5, tmp, dim.fraghfK[N]*CF);

			tmp=fmt::format("Exchange_{}_{}:,", i, j);
			fmt::fprintf(csv, fs6, tmp, dim.fraghfK[N]*CF);
			++N;
		}
	}

	sum=0.0;
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "Consistency check:\n");
	for(int i=0;i<dim.nfrag;++i)
		sum += mons[i].ehf;

	fmt::fprintf(os, fs5, "Total HF interaction energy:", (dim.eref-sum)*CF);
	fmt::fprintf(csv, fs6, "HFint:;", (dim.eref-sum)*CF);

	sum=0.0;
	for(int i=0;i<dim.nfrag;++i)
		sum += dim.fragehf[i] - mons[i].ehf;
	for(size_t i=0;i<dim.fraghfJ.size();++i) {
		sum+= dim.fraghfJ[i];
		sum+= dim.fraghfK[i];
	}

	summary["E(HF)_int"]=sum*CF;
	insertOrder.push_back("E(HF)_int");
	fmt::fprintf(os, fs5, "Sum of elec.preparation/electrostatic- and exchange Interaction:", sum*CF);
	fmt::fprintf(os, "\n");
	fmt::fprintf(csv, "\n");
}
void do_ccsdint(const Dimer& dim, vector<Monomer>& mons, ostream& os, map<string, double>& summary, vector<string>& insertOrder, ostream& csv) {

	double sum=0.0;

	fmt::fprintf(os, "******* CCSD Interaction energy ******* \n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "Warning: NonDisp (Strong Pairs) includes Singles excitations! \n");

	sum=0.0;
	for(int i=0;i<dim.nfrag;++i)
		sum += mons[i].ecorr;

	summary["E(NonDisp)_int"]=((dim.nondispStrong+dim.nondispWeak)-sum)*CF;
	insertOrder.push_back("E(NonDisp)_int");
	fmt::fprintf(os, fs5, "Total Nondispersion interaction energy:", ((dim.nondispStrong+dim.nondispWeak)-sum)*CF);
	fmt::fprintf(csv, fs6, "NonDisp_int:,", ((dim.nondispStrong+dim.nondispWeak)-sum)*CF);

	sum=0.0;
	for(size_t i=0;i<dim.fragccDispstrong.size();++i) {
		sum += dim.fragccDispstrong[i];
		sum += dim.fragccDispweak[i];
	}

	summary["E(Disp)"]=sum*CF;
	insertOrder.push_back("E(Disp)");
	fmt::fprintf(os, fs5, "Total Dispersion interaction energy:", sum*CF);
	fmt::fprintf(csv, fs6, "Disp_int:,", sum*CF);
	fmt::fprintf(os, "\n");
	fmt::fprintf(csv, "\n");
}

void do_triplesint(const Dimer& dim, vector<Monomer>& mons, ostream& os, map<string, double>& summary, vector<string>& insertOrder, ostream& csv) {

	fmt::fprintf(os, "******* Triples Interaction energy ******* \n");
	fmt::fprintf(os, "\n");

	double sum=0.0;
	for(int i=0;i<dim.nfrag;++i)
		sum+=mons[i].et;

	summary["E(Triples)_int"]=(dim.et-sum)*CF;
	insertOrder.push_back("E(Triples)_int");
	fmt::fprintf(os, fs5, "Total Triples interaction energy:", (dim.et-sum)*CF);
	fmt::fprintf(csv, fs6, "Triples_int:,", (dim.et-sum)*CF);
	fmt::fprintf(os, "\n");
	fmt::fprintf(csv, "\n");
}

void do_consistency(const Dimer& dim, vector<Monomer>& mons, ostream& os, map<string, double>& summary, vector<string>& insertOrder, ostream& csv) {

	fmt::fprintf(os, "Consistency check:\n");
	double sum=0.0;
	if(mons.size() == dim.nfrag) {
		for(int i=0;i<dim.nfrag;++i)
				sum+=mons[i].eccsdt;
		fmt::fprintf(os, fs5, "Total interaction energy (without geometric preparation):", (dim.eccsdt-sum)*CF);
		fmt::fprintf(csv, fs6, "Total_int:,", (dim.eccsdt-sum)*CF);
	}

	if(mons.size() == 2*dim.nfrag) {
		for(int i=dim.nfrag;i<mons.size();++i)
			sum+=mons[i].eccsdt;
		fmt::fprintf(os, fs5, "Total interaction energy (with geometric preparation):", (dim.eccsdt-sum)*CF);
		fmt::fprintf(csv, fs6, "Total_int:,", (dim.eccsdt-sum)*CF);
	}

	sum=0.0;
	for(auto it=summary.begin();it!=summary.end();++it)
		sum+=it->second;

	summary["E(Ges)_int"]=sum;
	insertOrder.push_back("E(Ges)_int");
	fmt::fprintf(os, fs5, "Sum of all mentioned interaction terms:", sum);
}

void do_summary(ostream& os, map<string, double>& summary, vector<string>& insertOrder, ostream& csv) {

	string tmp{};

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "************************************************************** \n");
	fmt::fprintf(os, "*******                    Summary                     ******* \n");
	fmt::fprintf(os, "************************************************************** \n");
	fmt::fprintf(os, "\n");

	for(size_t i=0;i<insertOrder.size();++i) {
		tmp=insertOrder[i];
		fmt::fprintf(os, fs5, tmp, summary[tmp]);
		fmt::fprintf(csv, fs6, tmp+",", summary[tmp]);
	}

	fmt::fprintf(os, "\n");
	fmt::fprintf(os, "\n");
}
