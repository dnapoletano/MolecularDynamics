#include "Atom.hpp"

Atom::Atom()
: Position{0.0,0.0,0.0}, Velocity{0.0,0.0,0.0}
{}

Atom::Atom(const Vec3d& pos, const Vec3d& vel)
: Position{pos}, Velocity{vel}
{}

Atom::~Atom()
{
}