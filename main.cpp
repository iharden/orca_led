/*
 * main.cpp
 *
 *  Created on: Dec 21, 2020
 *      Author: iharden
 */

#include <boost/program_options.hpp>

#include "class_monomer.hpp"
#include "class_dimer.hpp"
#include "functions.hpp"

using namespace std;
using namespace boost::program_options;

int main(int argc, char* argv[]) {
	chrono::time_point<chrono::high_resolution_clock> tstart, tend, ttotalstart, ttotalend;
	ttotalstart = chrono::high_resolution_clock::now();
	string dimername{};
	string tmp{};
	vector<string> monomers{};
	vector<string> comps{};

	// BOOST PROGRAM_OPTIONS DESCRIPTIONS
	options_description desc("Options");
	desc.add_options()
		("help,h", "produce help message")
		("dimer,d", value<string>(&dimername), "set Dimer-file")
		("monomers,m", value<vector<string>>(&monomers)->multitoken(), "set Monomer-files")
		("compare,c",value<vector<string>>(&comps)->multitoken(), "set files for comparison");

	variables_map vm;
	store(parse_command_line(argc, argv, desc), vm);
	notify(vm);

	if(vm.count("help") || argc == 1) {
		cout << desc << "\n";
		cout << "How you should use this program: \n";
		cout << "Start with the dimer file. Only one dimer file can be specified \n";
		cout << "Continue with the monomer files from fragment 1 to fragment N in Dimer geometry (for electronic preparation) \n";
		cout << "Continue with the monomer files from fragment 1 to fragment N in equilibrium geometry (for geometric preparation) \n";
		cout << "The program will produce a DIMERNAME.led file which contains the results \n";
		cout << "It will also write a DIMERNAME.csv file which can be opened by Excel or LibreOffice \n";
		cout << "If you mess something up with the ordering, the program might help you \n";
		cout << "\n";
		cout << "Command line example: \n";
		cout << "orca_led --dimer dimer.out --monomers frag1.out frag2.out frag1_opt.out frag2_opt.out \n";
		cout << "\n";
		cout << "For the compare mode you need to specify exactly two .led files \n";
		cout << "Command line example: \n";
		cout << "orca_led --compare file1.led file2.led \n";
		cout << "\n";
		return 0;
	}

	if(vm.count("compare")) {
		ofstream oss{"comparison.led"};
		ofstream csv{"comparison.csv"};
		do_compare(comps, oss, csv);
		ttotalend = chrono::high_resolution_clock::now();
		oss << "Total computation time: " << get_time(ttotalstart, ttotalend) << " ms\n";
		oss.close();
		csv.close();
		return 0;
	}

	// CALL CONSTRUCTORS FOR INPUT FILES
	tstart = chrono::high_resolution_clock::now();
	Dimer dimer(dimername);

	vector<Monomer> mons{};
	if(monomers.size()!=0) {
		for(string s : monomers) {
			mons.push_back(Monomer(s));
		}
	}

	tend = chrono::high_resolution_clock::now();
	cout << "Time for building data structures: " << get_time(tstart, tend) << "ms \n";

	cout << "Your Dimer has " << dimer.nfrag << " fragments \n";
	cout << "You gave me " << mons.size() << " monomer files \n";

	// NOW CALL THE CORRECT do_led() FUNCTION AND WRITE OUTPUT TO FILE
	tmp=dimername+".led";
	ofstream os(tmp);
	string tmpp=dimername+".csv";
	ofstream csv(tmpp);

	if(mons.size()<dimer.nfrag) {
		cout << "More fragments then monomer files present! I don't know what to do and therefore I will just print the LED of the dimer \n";
		do_led(dimer, os, csv);
	}
	else if(mons.size()==dimer.nfrag) {
		cout << "Perfect! You gave me as many monomer files as fragments are present. Therefore I will calculate the electronic preparation \n";
		do_led(dimer, mons, os, csv);
	}
	else if(mons.size()==2*dimer.nfrag) {
		cout << "Perfect! You gave me twice as many monomer files as fragments are present. Therefore I will calculate the electronic and the geometric preparation \n";
		do_led(dimer, mons, os, csv);
	}
	else {
		cout << "Wrong number of monomer files \n";
		cout << "Number of monomer files should equal nfrag or should be 2*nfrag \n";
		throw runtime_error("Wrong number of monomer files");
	}

	cout << "I will write output to " << tmp << "\n";

	ttotalend = chrono::high_resolution_clock::now();
	os << "Total computation time: " << get_time(ttotalstart, ttotalend) << " ms\n";
	os.close();
	csv.close();
}

