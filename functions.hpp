/*
 * functions.hpp
 *
 *  Created on: Dec 24, 2020
 *      Author: iharden
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

#include <iostream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <map>
#include <chrono>

#include <boost/algorithm/string.hpp>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/ostream.h>

#include "class_monomer.hpp"
#include "class_dimer.hpp"

// FUNCTIONS CALLED IN MAIN()
double get_time(std::chrono::time_point<std::chrono::high_resolution_clock>  start, std::chrono::time_point<std::chrono::high_resolution_clock> finish);
void do_led(const Dimer& dim, std::ostream& os, std::ostream& csv);
void do_led(const Dimer& dim, std::vector<Monomer>& mons, std::ostream& os, std::ostream& csv);
void do_compare(const std::vector<std::string> comps, std::ostream& os, std::ostream& csv);

// FUNCTIONS CALLED IN DO_LED(DIM)
void print_header(std::ostream& os);
void do_generalinfo(const Dimer& dim, std::ostream& os, std::ostream& csv);
void do_hartreefock(const Dimer& dim, std::ostream& os, std::ostream& csv);
void do_intra(const Dimer& dim, std::ostream& os, std::ostream& csv);
void do_ccsddisp(const Dimer& dim, std::ostream& os, std::ostream& csv);
void do_triplesdisp(const Dimer& dim, std::ostream& os, std::ostream& csv);
void do_ccsdnondisp(const Dimer& dim, std::ostream& os, std::ostream& csv);

// FUNCTIONS CALLED IN DO_LED(DIM,MONS)
void do_generalinfo(const Dimer& dim, std::vector<Monomer>& mon, std::ostream& os);
void do_geoprep(const Dimer& dim, std::vector<Monomer>& mon, std::ostream& os, std::map<std::string,double>& summary, std::vector<std::string>& insertOrder, std::ostream& csv);
void do_hfint(const Dimer& dim, std::vector<Monomer>& mon, std::ostream& os, std::map<std::string,double>& summary, std::vector<std::string>& insertOrder, std::ostream& csv);
void do_ccsdint(const Dimer& dim, std::vector<Monomer>& mon, std::ostream& os, std::map<std::string,double>& summary, std::vector<std::string>& insertOrder, std::ostream& csv);
void do_triplesint(const Dimer& dim, std::vector<Monomer>& mon, std::ostream& os, std::map<std::string,double>& summary, std::vector<std::string>& insertOrder, std::ostream& csv);
void do_consistency(const Dimer& dim, std::vector<Monomer>& mon, std::ostream& os, std::map<std::string,double>& summary, std::vector<std::string>& insertOrder, std::ostream& csv);
void do_summary(std::ostream& os, std::map<std::string,double>& summary, std::vector<std::string>& insertOrder, std::ostream& csv);
#endif /* FUNCTIONS_HPP_ */
