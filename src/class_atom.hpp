#ifndef CLASS_ATOM_HPP_
#define CLASS_ATOM_HPP_

#include <iostream>
#include <array>
#include <string>

class Atom
{
public:
  std::string symbx;
  std::array<double,3> pos{0.0,0.0,0.0};
  Atom();
  Atom(std::string &, double &, double &, double &);
};

#endif /* CLASS_ATOM_HPP_ */
