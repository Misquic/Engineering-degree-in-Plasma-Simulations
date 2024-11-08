#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <vector>
#include <fstream>
#include "World.h"
#include "Species.h"

namespace Output {
	void fields(World &world, std::vector<Species> &species, std::string name1 = "");
	void screenOutput(World &world, std::vector<Species> &species);
	void diagOutput(World &world, std::vector<Species> &species);

	void convergence(type_calc L2, int it, int ts);
	void convergence(int NR_it, int ts);
    //namespace
}

#endif