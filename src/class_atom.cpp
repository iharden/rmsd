#include "class_atom.hpp"

Atom::Atom() {

}

Atom::Atom(std::string &a, double &b, double &c, double &d)
{
  symbx = a;
  pos[0] = b;
  pos[1] = c;
  pos[2] = d;
}

