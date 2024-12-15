#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <vector>
#include <fstream>
#include "World.h"
#include "Species.h"
#include "funkc.h"

namespace Output {

	enum modes{none, all, screen, fields, particles, diagnostics, convergence};

	void fields(World& world, std::vector<Species>& species, std::string name1 = "");
	void screenOutput(World& world, std::vector<Species>& species);
	void diagOutput(World& world, std::vector<Species>& species);
	void particles(World& world, std::vector<Species>& species, int num_parts_to_output_base, std::string name1 = "");

	void convergence(type_calc L2, int it, int ts); //TODO disable?
	void convergence(int NR_it, int ts);
    //namespace
}
std::istream& operator>>(std::istream& in, Output::modes& type);
std::ostream& operator<<(std::ostream& out, Output::modes& type);

#endif