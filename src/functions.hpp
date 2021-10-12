/*
 * coordinate_transformations.hpp
 *
 *  Created on: Oct 5, 2021
 *      Author: iharden
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_


#include "class_atom.hpp"

#include<cmath>
#include<vector>
#include<tuple>
#include<array>
#include<exception>
#include<unordered_map>
#include<algorithm>
#include<chrono>
#include<omp.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

double get_time(std::chrono::time_point<std::chrono::high_resolution_clock>  start, std::chrono::time_point<std::chrono::high_resolution_clock> finish);
std::vector<double> get_masses(const std::vector<Atom>& atoms);
std::tuple<double, double, double> get_centerofmass(const std::vector<Atom>& atoms, const std::vector<double>& mass);

std::vector<Atom> quaternion(const std::vector<Atom>& sys_one, const std::vector<Atom>& sys_two);
std::vector<Atom> Kabsch(const std::vector<Atom>& sys_one, const std::vector<Atom>& sys_two);
double calc_rmsd(const std::vector<Atom>& sys_one_aligned, const std::vector<Atom>& sys_two);


#endif /* FUNCTIONS_HPP_ */
