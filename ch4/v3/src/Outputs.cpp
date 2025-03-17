#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Outputs.h"


void Output::fieldsOutput(World& world, std::vector<Species>& species, std::string name1){

	#ifndef DEBUG
	for(Species& sp: species){
		sp.computeGasProperties();
	}

    std::stringstream name;
	name<<"results/" << name1 << "fields_"<<std::setfill('0')<<std::setw(5)<<world.getTs()<<".vti";
    //name << "results/fields.vti";

    /*open output file*/
    std::ofstream out(name.str());
   	if(!out.is_open()){
        std::cerr<<"Could not open "<<name.str()<<std::endl;
        return;
    }

	/*ImageData is vtk format for structured Cartesian meshes*/
	out<<"<VTKFile type=\"ImageData\">\n";
	type_calc3 x0 = world.getX0();
	type_calc3 dx = world.getDx();
	out<<"<ImageData Origin=\""<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<"\" ";
	out<<"Spacing=\""<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<"\" ";
	out<<"WholeExtent=\"0 "<<world.ni-1<<" 0 "<<world.nj-1<<" 0 "<<world.nk-1<<"\">\n";
	
	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";

	/*node volumes, scalar*/
	out<<"<DataArray Name=\"NodeVol\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.node_vol;
	out<<"</DataArray>\n";

	/*objects ids, scalar*/
	out<<"<DataArray Name=\"ObjectID\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.object_id;
	out<<"</DataArray>\n";

	/*node_type, scalar*/
	out<<"<DataArray Name=\"NodeType\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.node_type;
	out<<"</DataArray>\n";

	/*potential, scalar*/
	out<<"<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.phi;
	out<<"</DataArray>\n";

	/*charge density, scalar*/
	out<<"<DataArray Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.rho;
	out<<"</DataArray>\n";

	/*species number densities*/
	for (Species &sp : species)
	{
		out<<"<DataArray Name=\"nd."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den;
		out<<"</DataArray>\n";
	}

	/*species average number densities*/
	for (Species &sp : species)
	{
		out<<"<DataArray Name=\"avg_nd."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den_avg;
		out<<"</DataArray>\n";
	}

	/*species velocity stream*/
	for (Species &sp : species)
	{
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.vel;
		out<<"</DataArray>\n";
	}

	/*species temperature*/
	for (Species &sp : species)
	{
		out<<"<DataArray Name=\"T."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.T;
		out<<"</DataArray>\n";
	}
	
	/*electric field, 3 component vector*/
	out<<"<DataArray Name=\"ef\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.ef;
	out<<"</DataArray>\n";
	
	/*close out tags*/
	out<<"</PointData>\n";

	/*cell data*/
	out<<"<CellData>\n";

	/*species macroparticle count per cell*/
	for (Species &sp : species)
	{
		out<<"<DataArray Name=\"mpc."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.macro_part_count;
		out<<"</DataArray>\n";
	}

	out<<"</CellData>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();

	/*clear samples if not steady state*/
 	if (!world.steadyState())
 		for (Species &sp:species) sp.clearSamples();
	#endif
};


//writes information to the screen
void Output::screenOutput(World &world, std::vector<Species> &species)
{
	std::cout<<"\r                                                                  \rts: "<<world.getTs();
	for (Species &sp:species)
		std::cout<<std::setprecision(3)<<"\t "<<sp.name<<":"<<sp.getNumParticles()<< " ";
	//std::cout<<std::endl;
	std::cout << std::flush;

}

//file stream handle
namespace Output {
std::ofstream f_diag;
}

/*save runtime diagnostics to a file*/
void Output::diagOutput(World &world, std::vector<Species> &species)
{
	using namespace Output;	//to get access to f_diag

	//is the file open?
	if (!f_diag.is_open())
	{
		f_diag.open("results/runtime_diags.csv");
		f_diag<<"ts,time,wall_time";
		for (Species &sp:species)
			f_diag<<",mp_count."<<sp.name<<",real_count."<<sp.name
				  <<",px."<<sp.name<<",py."<<sp.name<<",pz."<<sp.name
			      <<",KE."<<sp.name;
		f_diag<<",PE,E_total"<<std::endl;
	}

	f_diag<<world.getTs()<<","<<world.getTime();
	f_diag<<","<<world.getWallTime();

	double tot_KE = 0;
	for (Species &sp:species)
	{
		double KE = sp.getKE();	//species kinetic energy
		tot_KE += KE;		//increment total energy
		double3 mom = sp.getMomentum();

		f_diag<<","<<sp.getNumParticles()<<","<<sp.getMicroCount()
			  <<","<<mom[0]<<","<<mom[1]<<","<<mom[2]<<","<<KE;
	}

	//write out system potential and total energy
	double PE = world.getPE();
	f_diag<<","<<PE<<","<<(tot_KE+PE);

	f_diag<<"\n";	//use \n to avoid flush to disc
	if (world.getTs()%25==0) f_diag.flush();
}

namespace Output{
	std::ofstream f_convergence;
}

void Output::convergenceOutput(type_calc L2, int it, int ts){
	using namespace Output;	//to get access to f_diag
	//is the file open?
	if (!f_convergence.is_open())
	{
		f_convergence.open("convergence.csv");
		f_convergence<<"L2,PCG_it,ts,NR_it\n";
	}
	f_convergence<<L2<<","<<it<<","<<ts<<",\n";
};
void Output::convergenceOutput(int NR_it, int ts){
	using namespace Output;	//to get access to f_diag

	//is the file open?
	if (!f_convergence.is_open())
	{
		f_convergence.open("convergence.csv");
		f_convergence<<"L2,PCG_it,ts,NR_it\n";
	}
	f_convergence<<","<<","<< ts << "," << NR_it << "\n";
};

void Output::particlesOutput(World& world, std::vector<Species>& species, int num_parts_to_output_base, std::string name1){
		/*loop over all species*/
	for (Species &sp:species) {
		int num_parts_to_output = num_parts_to_output_base;
		//open a phase_sp_it.vtp
		std::stringstream name;
		name<<"results/" << name1 <<"parts_"<<sp.name<<"_"<<std::setfill('0')<<std::setw(5)<<world.getTs()<<".vtp";

		/*open output file*/
		std::ofstream out(name.str());
		if (!out.is_open()){
			std::cerr<<"Could not open "<<name.str()<<std::endl;
			return;
		}

		/*build a list of particles to output*/
		int num_of_parts = sp.getNumParticles();
		std::vector<std::unique_ptr<Particle>> to_output;

		if(num_of_parts != 0){
			int step = num_of_parts/num_parts_to_output; //always rounding down, so num_of_outputted_parts will achive num_parts_to_output, before num_of_outputted_parts*step will achive sp.particles.size()-1
			
			if(num_parts_to_output > num_of_parts){
				num_parts_to_output = num_of_parts;
				step = 1;
			}
			/*list of particles to output*/ //doing my way, because in book is weird
			
			to_output.reserve(num_parts_to_output);
			int num_of_outputted_parts = 0;

			while(num_of_outputted_parts < num_parts_to_output){
				if(num_of_outputted_parts*step<num_of_parts){ //sp.getNumParticles returns sp.particles.size(), additional protection
					to_output.emplace_back(std::make_unique<Particle>(sp.getConstPartRef(num_of_outputted_parts*step)));
					//out << num_of_outputted_parts*step << "\n";
				}
				num_of_outputted_parts++;
			}	
		}
		// double dp = num_parts_to_output/(double)sp.getNumParticles();
		// std::cout << sp.name << " " << dp << " " << sp.getNumParticles() << "\n";
		// double counter = 0;
		// std::vector<const Particle*> to_output;
		// for (const Particle &part : sp.getConstPartRef())	{
		// 	counter+=dp;
		// 	if (counter>1){ //save particle
		// 		to_output.emplace_back(&part);
		// 		counter=-1;}
		// }

		/*header*/
		out<<"<?xml version=\"1.0\"?>\n";
		out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		out<<"<PolyData>\n";
		out<<"<Piece NumberOfPoints=\""<<to_output.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
		out<<"NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

		/*points*/
		out<<"<Points>\n";
		out<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		// for (const Particle *part: to_output)
		for (const std::unique_ptr<Particle>& part: to_output)
			out<<part->pos<<"\n";
		out<<"</DataArray>\n";
		out<<"</Points>\n";

		/*velocities*/
		out<<"<PointData>\n";
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		// for (const Particle *part: to_output)
		for (const std::unique_ptr<Particle>& part: to_output)
			out<<part->vel<<"\n";
		out<<"</DataArray>\n";
		out<<"</PointData>\n";

		out<<"</Piece>\n";
		out<<"</PolyData>\n";
		out<<"</VTKFile>\n";

		out.close();
	}

	// for(Species& sp: species){
	// 	std::stringstream name;
	// 	name<<"results/" << name1 << "particles_" << sp.name << std::setfill('0')<<std::setw(5)<<world.getTs()<<".vti";

	// 	std::ofstream out(name.str());
	// 	if(!out.is_open()){
	// 		std::cerr << "Could not open " << name.str() << std::endl;
	// 		return;
	// 	}
	// 	int num_of_parts = sp.getNumParticles();
	// 	if(num_of_parts != 0 && num_parts_to_output <= num_of_parts){
	// 		/*list of particles to output*/ //doing my way, because in book is weird
	// 		std::vector<std::unique_ptr<Particle>> to_output;
	// 		to_output.reserve(num_parts_to_output);
	// 		int step = num_of_parts/num_parts_to_output; //always rounding down, so num_of_outputted_parts will achive num_parts_to_output, before num_of_outputted_parts*step will achive sp.particles.size()-1
	// 		int num_of_outputted_parts = 0;
	// 		while(num_of_outputted_parts < num_parts_to_output){
	// 			if(num_of_outputted_parts*step<num_of_parts){ //sp.getNumParticles returns sp.particles.size(), additional protection
	// 				to_output.emplace_back(std::make_unique<Particle>(sp.getConstPartRef(num_of_outputted_parts*step)));
	// 				out << num_of_outputted_parts*step << "\n";
	// 			}
	// 			num_of_outputted_parts++;
	// 		}		
			

	// 		out << "step: " << step << " num_of_outputted_parts: " << num_of_outputted_parts;
	// 	}
	// 	out << " particles.size(): " << sp.getNumParticles();

	// }


};

std::istream& Output::operator>>(std::istream& in, Output::modes& type){
    std::string input;
    in >> input;
	lower(input);

    if (input == "none" || input == "0") {
        type = Output::modes::none;
    } else if (input == "all" || input == "1") {
        type = Output::modes::all;
    } else if (input == "screen" || input == "2") {
        type = Output::modes::screen;
    } else if (input == "fields" || input == "3") {
        type = Output::modes::fields;
    } else if (input == "particles" || input == "4") {
        type = Output::modes::particles;
    } else if (input == "diagnostics" || input == "5") {
        type = Output::modes::diagnostics;
	} else if (input == "convergence" || input == "6") {
        type = Output::modes::convergence;
    } else {
		try{
        	throw(std::invalid_argument("Wrong Output::modes value"));
		}
        catch(std::invalid_argument& e){
			std::cerr << e.what() << ", setting Output::modes::fields \n";
		}
		type = Output::modes::fields;  // DEFAULT
    }

    return in;

};

std::ostream& Output::operator<<(std::ostream& out, Output::modes& type){
    switch(type){
    case Output::modes::none:
        out << "none";
        break;
    case Output::modes::all:
        out << "all";
        break;
    case Output::modes::screen:
        out << "screen";
        break;
    case Output::modes::fields:
        out << "fields";
        break;
    case Output::modes::particles:
        out << "particles";
        break;
    case Output::modes::diagnostics:
        out << "diagnostics";
        break;
    case Output::modes::convergence:
        out << "convergence";
        break;
    default:
        break;
    }
    return out;
};
