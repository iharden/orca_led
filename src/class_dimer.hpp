/*
 * class_dimer.hpp
 *
 *  Created on: Dec 22, 2020
 *      Author: iharden
 */

#ifndef CLASS_DIMER_HPP_
#define CLASS_DIMER_HPP_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <exception>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

class Dimer {
public:
	int nel = 0;
	int nbasis = 0;
	int nfrag = 0;
	std::string hftype{};
	std::string name{};

	double ehf = 0;
	double eref = 0;
	double eccsd = 0;
	double eccsdt = 0;
	double ecorr = 0;
	double et = 0;
	double ecorrt = 0;

	double nondispStrong = 0;
	double nondispWeak = 0;

	std::vector<double> fragehf{};

	std::vector<double> fraghfJ{};
	std::vector<double> fraghfK{};
	std::vector<double> hfint{};

	std::vector<double> fragccDispstrong{};
	std::vector<double> fragccDispweak{};
	std::vector<double> fragccDisptriples{};
	std::vector<double> fragccDispges{};

	std::vector<double> fragInterstrong{};
	std::vector<double> fragInterweak{};
	std::vector<double> fragIntertriples{};
	std::vector<double> ccint{};
	std::vector<double> gamma{};

	std::vector<double> fragIntrastrong{};
	std::vector<double> fragIntraweak{};
	std::vector<double> fragIntratriples{};
	std::vector<double> fragIntrasingles{};
	std::vector<double> fragIntrages{};

	std::vector<double> fragccNonDisptriples{};
	std::vector<double> fragccNonDispges{};

	double delocalizedstrongpairs{};
	double delocalizedtriples{};

	Dimer();
	Dimer(std::string argv);
};


#endif /* CLASS_DIMER_HPP_ */
