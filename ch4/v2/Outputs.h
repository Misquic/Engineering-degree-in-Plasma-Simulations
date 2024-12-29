#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <vector>
#include <fstream>
#include "World.h"
#include "Species.h"
#include "funkc.h"

namespace Output {

	enum modes{none, all, screen, fields, particles, diagnostics, convergence};

	void fieldsOutput(World& world, std::vector<Species>& species, std::string name1 = "");
	void screenOutput(World& world, std::vector<Species>& species);
	void diagOutput(World& world, std::vector<Species>& species);
	void particlesOutput(World& world, std::vector<Species>& species, int num_parts_to_output_base, std::string name1 = "");

	void convergenceOutput(type_calc L2, int it, int ts); //TODO disable?
	void convergenceOutput(int NR_it, int ts);
    //namespace
std::ostream& operator<<(std::ostream& out, Output::modes& type);
std::istream& operator>>(std::istream& in, Output::modes& type);
}

#endif