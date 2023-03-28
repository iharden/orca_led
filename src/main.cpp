/*
 * main.cpp
 *
 *  Created on: Dec 21, 2020
 *      Author: iharden
 */

#include <future>
#include <algorithm>
#include <string_view>
#include "class_monomer.hpp"
#include "class_dimer.hpp"
#include "functions.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	chrono::time_point<chrono::high_resolution_clock> tstart, tend, ttotalstart, ttotalend;
	ttotalstart = chrono::high_resolution_clock::now();

	string dimername{};
	string tmp{};
	vector<string_view> monomers{};
	vector<string_view> comps{};

	vector<string_view> cmd_arg{};
	for(int i=1;i<argc;++i)
		cmd_arg.emplace_back(argv[i]);

	if((find(cmd_arg.begin(), cmd_arg.end(), "-h") != cmd_arg.end()) || (find(cmd_arg.begin(), cmd_arg.end(), "--help") != cmd_arg.end())){
		cout << "______________________________________________________________________________________________________________________________ \n";
		cout << "\n";
		cout << "  -h, --help:       Print this message \n";
		cout << "  -d, --dimer:      Specify ORCA output file that contains LED results \n";
		cout << "  -m, --monomers:   Specify list of ORCA output files of monomers for electronic/geometric preparation \n";
		cout << "  -c, --compare:    Specify two .led files obtained from previous orca_led runs to compare two binding energies. \n"
				"                    Works only if the decomposition of binding energies was carried out! \n \n";
		cout << "______________________________________________________________________________________________________________________________ \n \n";

		cout << "How you should use this program: \n";
		cout << "Start with the dimer file. Only one dimer file can be specified \n";
		cout << "Continue with the monomer files from fragment 1 to fragment N in dimer geometry (for electronic preparation) \n";
		cout << "Continue with the monomer files from fragment 1 to fragment N in equilibrium geometry (for geometric preparation) \n";
		cout << "The program will produce a DIMERNAME.led file which contains the results \n";
		cout << "It will also write a DIMERNAME.csv file which can be opened by Excel or LibreOffice \n";
		cout << "If you mess something up with the ordering, the program might help you \n";
		cout << "\n";
		cout << "Command line example: \n";
		cout << "orca_led --dimer dimer.out --monomers frag1.out frag2.out frag1_opt.out frag2_opt.out \n";
		cout << "\n";
		cout << "For the compare mode you have to specify exactly two .led files \n";
		cout << "Command line example: \n";
		cout << "orca_led --compare file1.led file2.led \n";
		cout << "\n";
		return 0;
	}

	if((find(cmd_arg.begin(), cmd_arg.end(), "-c") != cmd_arg.end()) || (find(cmd_arg.begin(), cmd_arg.end(), "--compare") != cmd_arg.end())) {

		auto found=find(cmd_arg.begin(), cmd_arg.end(), "-c");
		if(found==cmd_arg.end())
			found=find(cmd_arg.begin(), cmd_arg.end(), "--compare");

		int idx=distance(cmd_arg.begin(),found);
		comps.push_back(cmd_arg[idx+1]);
		comps.push_back(cmd_arg[idx+2]);

		ofstream oss{"comparison.led"};
		ofstream csv{"comparison.csv"};
		//do_compare(comps, oss, csv);
		ttotalend = chrono::high_resolution_clock::now();
		oss << "Total computation time: " << get_time(ttotalstart, ttotalend) << " ms\n";
		oss.close();
		csv.close();
		return 0;
	}

	if((find(cmd_arg.begin(), cmd_arg.end(), "-d") != cmd_arg.end()) || (find(cmd_arg.begin(), cmd_arg.end(), "--dimer") != cmd_arg.end())) {

		auto found=find(cmd_arg.begin(), cmd_arg.end(), "-d");
		if(found==cmd_arg.end())
			found=find(cmd_arg.begin(), cmd_arg.end(), "--dimer");

		int idx=distance(cmd_arg.begin(),found);
		dimername=cmd_arg[idx+1];
	}

	if((find(cmd_arg.begin(), cmd_arg.end(), "-m") != cmd_arg.end()) || (find(cmd_arg.begin(), cmd_arg.end(), "--monomers") != cmd_arg.end())) {

		auto found=find(cmd_arg.begin(), cmd_arg.end(), "-m");
		if(found==cmd_arg.end())
			found=find(cmd_arg.begin(), cmd_arg.end(), "--monomers");

		int idx=distance(cmd_arg.begin(),found);

		while(cmd_arg[idx] != "-d" && cmd_arg[idx] != "--dimer" && idx < cmd_arg.size()-1) {
			monomers.push_back(cmd_arg[idx+1]);
			++idx;
		}
	}

	// CALL CONSTRUCTORS FOR INPUT FILES
	tstart = chrono::high_resolution_clock::now();

	vector<Monomer> mons{};
	/*future<void> fut = async([&mons, &monomers] () {
		for(string_view s:monomers)
			mons.emplace_back(s);
	});*/

	vector<future<Monomer>> futs;
	for(const auto& monomer : monomers) {
		futs.push_back(async([](const auto& mon) {return Monomer(mon);},monomer));
	}

	Dimer dimer(dimername);
	for(auto& fut:futs)
		mons.push_back(fut.get());
	//fut.get();

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

