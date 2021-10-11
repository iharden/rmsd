
#include "class_atom.hpp"
#include "functions.hpp"

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<tuple>
#include<omp.h>
#include<exception>
#include<unordered_map>
#include<algorithm>

using namespace std;

int main(int argc, char *argv[]) {

	double wertx  = 0.0, werty = 0.0, wertz = 0.0;
	string sym, dummy;
	int natoms = 0;
	bool Quaternion=false;

	// ********************************************************************
	// CHECKING THE INPUT FILES AND READING THEM IN                 *******
	// ********************************************************************
	for(int i=1;i<argc;++i) {
		string s=argv[i];
		if(s.find("-q") != string::npos || s.find("--quarternion") != string::npos)
			Quaternion=true;
	}
	if(argc > 4)
	  throw runtime_error("To many arguments given!");

	ifstream input(argv[1]);
	if(!input)
	  throw runtime_error(string ("File could not be found: ") + argv[1]);

	ifstream inputzwei(argv[2]);
	if(!inputzwei)
	  throw runtime_error(string ("File could not be found: ") + argv[2]);

	getline(input,dummy);
	natoms = stoi(dummy);
	getline(input,dummy);

	vector<Atom> sys_one, sys_two;

	while(input >> sym >> wertx >> werty >> wertz) {
		sys_one.push_back(Atom(sym, wertx, werty, wertz));
	}

	getline(inputzwei,dummy);
	getline(inputzwei,dummy);

	while(inputzwei >> sym >> wertx >> werty >> wertz) {
		sys_two.push_back(Atom(sym, wertx, werty, wertz));
	}

	input.close();
	inputzwei.close();

	// ********************************************************************
	// CHECK ATOM ORDER				                              *******
	// ********************************************************************

	if(sys_one.size() != sys_two.size() || sys_one.size() != natoms || sys_two.size() != natoms)
		throw logic_error("Number of atoms does not match or wrong number of atoms was specified in one of the xyz file!");

	for(int i=0;i<sys_one.size();++i) {
		if(sys_one[i].symbx != sys_two[i].symbx) {
		  string a,b;
		  a=sys_one[i].symbx;
		  b=sys_two[i].symbx;
		  transform(a.begin(),a.end(), a.begin(), ::toupper);
		  transform(b.begin(),b.end(), b.begin(), ::toupper);
		  if(a != b)
			  throw logic_error("Element order does not match!");
		}
	}

	// ********************************************************************
	// CALCULATING THE CENTER OF MASS                               *******
	// ********************************************************************

	vector<double> mass = get_masses(sys_one);

	tuple<double, double, double> com_one = get_centerofmass(sys_one, mass);
	tuple<double, double, double> com_two = get_centerofmass(sys_two, mass);


	// ********************************************************************
	// SHIFTING CENTER OF MASS TO THE ORIGIN                        *******
	// ********************************************************************

	#pragma omp parallel for
	for(int i=0;i<natoms;++i) {
		sys_one[i].pos[0] -= get<0>(com_one);
		sys_one[i].pos[1] -= get<1>(com_one);
		sys_one[i].pos[2] -= get<2>(com_one);

		sys_two[i].pos[0] -= get<0>(com_two);
		sys_two[i].pos[1] -= get<1>(com_two);
		sys_two[i].pos[2] -= get<2>(com_two);
	}

	// ********************************************************************
	// GET ALIGNED SYSTEM FROM COORDINATE TRANSFORMATION            *******
	// ********************************************************************

	//auto start = chrono::high_resolution_clock::now();
	vector<Atom> sys_one_aligned;
	if(Quaternion) {
		sys_one_aligned = quaternion(sys_one, sys_two);
	}
	else {
		sys_one_aligned = Kabsch(sys_one, sys_two);
	}
	//auto end = chrono::high_resolution_clock::now();

	// ********************************************************************
	// ACTUAL RMSD-CALCULATION                                      *******
	// ********************************************************************

	double rmsd = calc_rmsd(sys_one_aligned, sys_two);
	cout << argv[1] << "  " << argv[2] << "   " << rmsd << "   " << endl;

	ofstream output("system_one_aligned.xyz");
	output << natoms << "\n \n";
	for(int i=0;i<sys_one_aligned.size();++i) {
		output << sys_one_aligned[i].symbx << " " << sys_one_aligned[i].pos[0] << " " << sys_one_aligned[i].pos[1] << " " << sys_one_aligned[i].pos[2] << " \n";
	}
	output.close();
}





