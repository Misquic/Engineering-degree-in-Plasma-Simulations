#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Outputs.h"


void Output::fields(World& world, std::vector<Species>& species, std::string name1){
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
	
	/*electric field, 3 component vector*/
	out<<"<DataArray Name=\"ef\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.ef;
	out<<"</DataArray>\n";
	
	/*close out tags*/
	out<<"</PointData>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();
};


//writes information to the screen
void Output::screenOutput(World &world, std::vector<Species> &species)
{
	std::cout<<"\r                                                                                \rts: "<<world.getTs();
	for (Species &sp:species)
		std::cout<<std::setprecision(3)<<"\t "<<sp.name<<":"<<sp.getNumParticles()<< " ";
	//std::cout<<std::endl;
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
		f_diag.open("runtime_diags.csv");
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

void Output::convergence(type_calc L2, int it, int ts){
	using namespace Output;	//to get access to f_diag

	//is the file open?
	if (!f_convergence.is_open())
	{
		f_convergence.open("convergence.csv");
		f_convergence<<"L2,PCG_it,ts,NR_it\n";
	}
	f_convergence<<L2<<","<<it<<","<<ts<<",\n";
};
void Output::convergence(int NR_it, int ts){
	using namespace Output;	//to get access to f_diag

	//is the file open?
	if (!f_convergence.is_open())
	{
		f_convergence.open("convergence.csv");
		f_convergence<<"L2,PCG_it,ts,NR_it\n";
	}
	f_convergence<<","<<","<< ts << "," << NR_it << "\n";
};