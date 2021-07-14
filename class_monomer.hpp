/*
 * monomer.hpp
 *
 *  Created on: Dec 21, 2020
 *      Author: iharden
 */

#ifndef CLASS_MONOMER_HPP_
#define CLASS_MONOMER_HPP_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <exception>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

class Monomer {
public:
	int nel = 0;
	int nbasis = 0;
	std::string name{};

	double ehf = 0;
	double eccsd = 0;
	double eccsdt = 0;
	double ecorr = 0;
	double et = 0;
	double ecorrt = 0;

	Monomer() = delete;
	Monomer(std::string argv);
};



#endif /* CLASS_MONOMER_HPP_ */
