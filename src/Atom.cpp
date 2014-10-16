#include "Atom.hpp"


Atom::Atom()
{
  type="A";
  property=0;
}

Atom::Atom(std::string new_type, double prop_val)
{
  type=new_type;
  property = prop_val;
}

void Atom::setProperty(double val)
{
  property=val;
}

void Atom::setType(std::string new_type)
{
  type = new_type;
}

std::string Atom::getType()
{
  return type;
}

double Atom::getProperty()
{
  return property;
}
