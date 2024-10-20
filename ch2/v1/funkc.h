#ifndef FUNKC_H
#define FUNKC_H

#include <vector>
#include <string>
#include <algorithm>

int calc_potential();
int calc_electric_field();
int calc_charge_density();

int update_particles_velocity();
int update_particles_positions();

int load_particles();

int perform_collisions();
int inject_particles();
int output_diagnostics();

int output_results();



#endif
