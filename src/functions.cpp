/*
 * coordinate_transformations.cpp
 *
 *  Created on: Oct 5, 2021
 *      Author: iharden
 */


/* Here the actual coordinate transformations take place. Calculating the RMSD once the molecules are aligend is trivial */

#include "class_atom.hpp"
#include "functions.hpp"

using namespace std;

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double get_time(chrono::time_point<chrono::high_resolution_clock>  start, chrono::time_point<chrono::high_resolution_clock> finish) {

	chrono::duration<double> elapsed = finish - start;
	return elapsed.count();
}

vector<double> get_masses(const vector<Atom>& atoms) {

	unordered_map<string,double> masses{
		{"H",1.008}, {"HE",4.0026}, {"LI",6.941}, {"BE",9.0122},
		{"B",10.81}, {"C",12.011}, {"N",14.007}, {"O",15.999},
		{"F",18.998}, {"NE",20.180}, {"NA",22.9897}, {"MG",24.305},
		{"AL",26.9815},  {"SI",28.0855}, {"P",30.9738},  {"S",32.065},
		{"CL",35.453},  {"K",39.0983}, {"AR",39.948},  {"CA",40.078},
		{"SC",44.9559},  {"TI",47.867}, {"V",50.9415},  {"CR",51.9961},
		{"MN",54.938},  {"FE",55.845}, {"NI",58.6934},  {"CO",58.9332},
		{"CU",63.546}, {"ZN",65.39}, {"GA",69.723},  {"GE",72.64},
		{"AS",74.9216},  {"SE",78.96},  {"BR",79.904},  {"KR",83.8},
		{"RB",85.4678},  {"SR",87.62},  {"Y",88.9059}, {"ZR",91.224},
		{"NB",92.9064},  {"MO",95.94},  {"TC",98},  {"RU",101.07},
		{"RH",102.9055},  {"PD",106.42},  {"AG",107.8682}, {"CD",112.411},
		{"IN",114.818},  {"SN",118.71},  {"SB",121.76},  {"I",126.9045},
		{"TE",127.6},  {"XE",131.293},  {"CS",132.9055},  {"BA",137.327},
		{"LA",138.9055},  {"CE",140.116},  {"PR",140.9077},  {"ND",144.24},
		{"PM",145},  {"SM",150.36}, {"EU",151.964},  {"GD",157.25},
		{"TB",158.9253},  {"DY",162.5}, {"HO",164.9303},  {"ER",167.259},
		{"TM",168.9342},  {"YB",173.04}, {"LU",174.967},  {"HF",178.49},
		{"TA",180.9479},  {"W",183.84}, {"RE",186.207},  {"OS",190.23},
		{"IR",192.217},  {"PT",195.078}, {"AU",196.9665},  {"HG",200.59},
		{"TL",204.3833},  {"PB",207.2}, {"BI",208.9804},  {"PO",209},
		{"AT",210},  {"RN",222}, {"FR",223},  {"RA",226},
		{"AC",227},  {"PA",231.0359}, {"TH",232.0381},  {"NP",237},
		{"U",238.0289},  {"AM",243}, {"PU",244},  {"CM",247},
		{"BK",247},  {"CF",251},  {"ES",252},  {"FM",257},
		{"MD",258},  {"NO",259}, {"RF",261},  {"LR",262},
		{"DB",262},  {"BH",264}, {"SG",266},  {"MT",268},
		{"RG",272},  {"HS",277},
		{"h",1.008}, {"He",4.0026}, {"Li",6.941}, {"Be",9.0122},
		{"b",10.81}, {"c",12.011}, {"n",14.007}, {"o",15.999},
		{"f",18.998}, {"Ne",20.180}, {"Na",22.9897}, {"Mg",24.305},
		{"Al",26.9815},  {"Si",28.0855}, {"p",30.9738},  {"s",32.065},
		{"Cl",35.453},  {"k",39.0983}, {"Ar",39.948},  {"Ca",40.078},
		{"Sc",44.9559},  {"Ti",47.867}, {"v",50.9415},  {"Cr",51.9961},
		{"Mn",54.938},  {"Fe",55.845}, {"Ni",58.6934},  {"Co",58.9332},
		{"Cu",63.546}, {"Zn",65.39}, {"Ga",69.723},  {"Ge",72.64},
		{"As",74.9216},  {"Se",78.96},  {"Br",79.904},  {"Kr",83.8},
		{"Rb",85.4678},  {"Sr",87.62},  {"y",88.9059}, {"Zr",91.224},
		{"Nb",92.9064},  {"Mo",95.94},  {"Tc",98},  {"Ru",101.07},
		{"Rh",102.9055},  {"Pd",106.42},  {"Ag",107.8682}, {"Cd",112.411},
		{"In",114.818},  {"Sn",118.71},  {"Sb",121.76},  {"i",126.9045},
		{"Te",127.6},  {"Xe",131.293},  {"Cs",132.9055},  {"Ba",137.327},
		{"La",138.9055},  {"Ce",140.116},  {"Pr",140.9077},  {"Nd",144.24},
		{"Pm",145},  {"Sm",150.36}, {"Eu",151.964},  {"Gd",157.25},
		{"Tb",158.9253},  {"Dy",162.5}, {"Ho",164.9303},  {"Er",167.259},
		{"Tm",168.9342},  {"Yb",173.04}, {"Lu",174.967},  {"Hf",178.49},
		{"Ta",180.9479},  {"w",183.84}, {"Re",186.207},  {"Os",190.23},
		{"Ir",192.217},  {"Pt",195.078}, {"Au",196.9665},  {"Hg",200.59},
		{"Tl",204.3833},  {"Pb",207.2}, {"Bi",208.9804},  {"Po",209},
		{"At",210},  {"Rn",222}, {"Fr",223},  {"Ra",226},
		{"Ac",227},  {"Pa",231.0359}, {"Th",232.0381},  {"Np",237},
		{"u",238.0289},  {"Am",243}, {"Pu",244},  {"Cm",247},
		{"Bk",247},  {"Cf",251},  {"Es",252},  {"Fm",257},
		{"Md",258},  {"No",259}, {"Rf",261},  {"Lr",262},
		{"Db",262},  {"Bh",264}, {"Sg",266},  {"Mt",268},
		{"Rg",272},  {"Hs",277}
	};

	vector<double> mass;
	for(int i=0;i<atoms.size();++i) {
	  mass.push_back(masses[atoms[i].symbx]);
	}
	return mass;
}

tuple<double, double, double> get_centerofmass(const vector<Atom>& atoms, const vector<double>& mass) {

	double massges=0;
	#pragma omp parallel for reduction(+:massges)
	for(int i=0; i<atoms.size();++i) {
		massges += mass[i];
	}

	double rx = 0.0;
	double ry = 0.0;
	double rz = 0.0;

	#pragma omp parallel for reduction(+:rx,ry,rz)
	for(int i=0; i<atoms.size(); ++i) {
		rx += mass[i] * atoms[i].pos[0];
		ry += mass[i] * atoms[i].pos[1];
		rz += mass[i] * atoms[i].pos[2];
	}

	rx /= massges;
	ry /= massges;
	rz /= massges;

	return make_tuple(rx, ry, rz);
}

std::vector<Atom> quaternion(const std::vector<Atom>& sys_one, const std::vector<Atom>& sys_two) {

	// ********************************************************************
	// DEFINING AND WRITING CORRELATION MATRIX                      *******
	// ********************************************************************

	int natoms=sys_one.size();
	array<array<double,3>,3> R;

	for(int i=0;i<3;++i) {
	  for(int j=0;j<3;++j) {
		  R[i][j] = 0.0;
		  for(int k=0;k<natoms;++k) {
			  R[i][j] += sys_one[k].pos[i] * sys_two[k].pos[j];
		  }
	  }
	}

	// ********************************************************************
	// DEFINING AND WRITING F-MATRIX                                *******
	// ********************************************************************


	Eigen::Matrix4d F;
	F(0,0) = R[0][0] + R[1][1] + R[2][2];
	F(0,1) = R[1][2] - R[2][1];
	F(0,2) = R[2][0] - R[0][2];
	F(0,3) = R[0][1] - R[1][0];

	F(1,0) = F(0,1);
	F(1,1) = R[0][0] - R[1][1] - R[2][2];
	F(1,2) = R[0][1] + R[1][0];
	F(1,3) = R[0][2] + R[2][0];

	F(2,0) = F(0,2);
	F(2,1) = F(1,2);
	F(2,2) = -R[0][0] + R[1][1] - R[2][2];
	F(2,3) = R[1][2] + R[2][1];

	F(3,0) = F(0,3);
	F(3,1) = F(1,3);
	F(3,2) = F(2,3);
	F(3,3) = -R[0][0] - R[1][1] + R[2][2];

	// ********************************************************************
	// DIAGONALIZING THE F-MATRIX                                   *******
	// ********************************************************************

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(4);
	solver.compute(F);
	auto evecs = solver.eigenvectors();

	array<double,4> q;
	for(int i=0;i<4;++i) {
	  q[i] = evecs(i,3);
	}

	// ********************************************************************
	// DEFINING AND WRITING THE TRANSFORMATION MATRIX               *******
	// ROTATION MATRIX TAKEN FROM                                   *******
	// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
	// ********************************************************************

	Eigen::Matrix3d U;
	//U(0,0) = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
	U(0,0) = 1 - 2 * (q[2] * q[2] + q[3] * q[3]);
	U(0,1) = 2 * (q[1] * q[2] - q[0] * q[3]);
	U(0,2) = 2 * (q[1] * q[3] + q[0] * q[2]);

	U(1,0) = 2 * (q[1] * q[2] + q[0] * q[3]);
	//U(1,1) = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	U(1,1) = 1 - 2 * (q[1] * q[1] + q[3] * q[3]);
	U(1,2) = 2 * (q[2] * q[3] - q[0] * q[1]);

	U(2,0) = 2 * (q[1] * q[3] - q[0] * q[2]);
	U(2,1) = 2 * (q[2] * q[3] + q[0] * q[1]);
	//U(2,2) = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
	U(2,2) = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);

	//cout << "Transformation Matrix U" << endl;
	//cout << U << endl;
	//cout << "The determinant of U is: " << U.determinant() << endl;

	// ********************************************************************
	// TRANSFORMING COORDINATES OF x-SET                            *******
	// ********************************************************************

	vector<Atom> sys_one_aligned(natoms);
	#pragma omp parallel for
	for(int i=0;i<natoms;++i) {
	  sys_one_aligned[i].symbx = sys_one[i].symbx;
	  for(int k=0;k<3;++k) {
		  for(int l=0;l<3;++l) {
			  sys_one_aligned[i].pos[k] += U(k,l) * sys_one[i].pos[l];
		  }
	  }
	}
	return sys_one_aligned;
}

std::vector<Atom> Kabsch(const std::vector<Atom>& sys_one, const std::vector<Atom>& sys_two) {

	// ********************************************************************
	// SET UP H-MATRIX												*******
	// ********************************************************************
	int natoms = sys_one.size();
	Eigen::Matrix3d H;

	for(int i=0;i<3;++i) {
		for(int j=0;j<3;++j) {
			H(i,j) = 0.0;
			for(int k=0;k<natoms;++k) {
				H(i,j) += sys_one[k].pos[i] * sys_two[k].pos[j];
			}
		}
	}

	// ********************************************************************
	// Build Rotation Matrix									    *******
	// ********************************************************************

	Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
	auto d=(svd.matrixV()*(svd.matrixU().transpose())).determinant();
	d=sgn(d);

	Eigen::Matrix3d tmp;
	tmp << 1,0,0,0,1,0,0,0,d;

	Eigen::Matrix3d R = svd.matrixV()*tmp*(svd.matrixU().transpose());

	// ********************************************************************
	// TRANSFORMING COORDINATES OF x-SET                            *******
	// ********************************************************************

	vector<Atom> sys_one_aligned(natoms);
	#pragma omp parallel for
	for(int i=0;i<natoms;++i) {
	  sys_one_aligned[i].symbx = sys_one[i].symbx;
	  for(int k=0;k<3;++k) {
		  for(int l=0;l<3;++l) {
			  sys_one_aligned[i].pos[k] += R(k,l) * sys_one[i].pos[l];
		  }
	  }
	}
	return sys_one_aligned;
}

double calc_rmsd(const std::vector<Atom>& sys_one_aligned, const std::vector<Atom>& sys_two) {

	int natoms=sys_one_aligned.size();
	vector<double> delta(natoms,0.0);
	double sumdelta = 0.0;
	double rmsd = 0.0;

	#pragma omp parallel for
	for(int i=0;i<natoms;++i) {
	  for(int k=0;k<3;++k)
		  delta[i] += (sys_one_aligned[i].pos[k] - sys_two[i].pos[k]) * (sys_one_aligned[i].pos[k] - sys_two[i].pos[k]); // distance^2
	}

	#pragma omp parallel for reduction(+:sumdelta)
	for(int i=0; i<natoms; i++) {
	  sumdelta +=delta[i];
	}

	sumdelta /= natoms;
	rmsd = sqrt(sumdelta);
	return rmsd;
}


















