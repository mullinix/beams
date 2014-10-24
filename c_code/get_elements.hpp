#ifndef GET_ELEMENTS_HPP
#define GET_ELEMENTS_HPP

// basic io
#include <iostream>
// enhanced file io
#include <fstream>
// c++ strings
#include <string>
// c strings
#include <cstring>
// enhanced string manipulation
#include <sstream>
// dynamic arrays via vectors
#include <vector>
// enhanced vector/array algorithms
#include <algorithm>
// for mainpulation of stdio
#include <iomanip>

/*includes and macros for later*/
//#include "Eigen/Dense"
//#include <cmath>
//#define PI 3.141592653589793238462643

/*Space of variable and function names in which we shall operate*/
using namespace std;

//global structures for manipulation
struct input { //struct for the inputs file data
	int composite_number;
	vector<string> material_files;
	string loads_file;
	string nodes_file;
	string bcs_file;
	vector <int> displacement_dofs;
	double total_length;
	double a;
	double Omega;
} inputs;

struct node { //struct for the nodes with index and location
	int index;
	vector <double> location;
};

struct beam { //struct for the beams with index, width, height, node connections
	int index;
	double width;
	double height;
	vector<int> node_connections;
};

struct material { //struct for the material
	double E;
	double nu;
	double rho;
	double stress_yield;
	double stress_ultimate;
	string geometry_file;
	vector<beam> beams;
};

struct load { //struct for the applied force with node index, force, direction
	int node_index;
	double force;
	vector<double> direction;
};

struct bc {
	int node_id;
	vector <bool> clamp_dir;
};

//declare helper functions
vector<string> parse_csv(string input_txt);
int read_inputs(string inputs_filename);
int read_material(string material_filename, int matl_idx);
int read_geometry(string geometry_filename, int matl_idx);
int read_nodes(string nodes_filename);
int read_loads(string loads_filename);
int read_bcs(string bcs_filename);
void print_inputs(input inputs);
void print_materials(vector<material> materials);
void print_loads(vector<load> loads);
void print_nodes(vector<node> nodes);
void print_bcs(vector<bc> bcs);
void print_structs(void);
double get_Machine_Epsilon(void);

//construct global variables
vector<material> materials;
vector<load> loads;
vector<node> nodes;
vector<bc> bcs;

//helper functions
vector<string> parse_csv(string input_txt) {
	string tmp_txt;
	vector<string> parsed_txt;
	char tmp_char;
	for(int idx=0; idx<input_txt.length(); idx++) {
		tmp_char = input_txt[idx];
		if (tmp_char==',') {
			parsed_txt.push_back(tmp_txt);
			tmp_txt.clear();
		}
		else {
			tmp_txt+=tmp_char;
		}
	}
	if(tmp_txt.compare("")!=0) {
		parsed_txt.push_back(tmp_txt);
	}
	return parsed_txt;
}

int read_inputs(string inputs_filename) {
	ifstream inputs_fptr;
	inputs_fptr.open(inputs_filename.c_str());
	if (!inputs_fptr.good()){
		cout << "Error: input file " << inputs_filename << " could not be opened!" << endl; 
		return 100;
	}
	string varname,varvalue,dlm;
	vector <string> inputs_array;
	inputs_array.push_back("composite_number");
	inputs_array.push_back("material_files");
	inputs_array.push_back("loads_file");
	inputs_array.push_back("nodes_file");
	inputs_array.push_back("bcs_file");
	inputs_array.push_back("displacement_dofs");
	inputs_array.push_back("total_length");
	inputs_array.push_back("a");
	inputs_array.push_back("Omega");
	while(inputs_fptr >> varname >> dlm >> varvalue) {
		vector<string>::const_iterator array_location = find(inputs_array.begin(), inputs_array.end(), varname);
		if(array_location == inputs_array.end()) {
			cout << "Error: invalid entry '" << varname; 
			cout << "' in inputs file '" << inputs_filename;
			cout << "'!" << endl;
			return 101;
		}
		if(varname=="composite_number") {
			istringstream (varvalue) >> inputs.composite_number;
		}
		else if(varname == "loads_file") {
			inputs.loads_file = varvalue;
		}
		else if(varname == "nodes_file") {
			inputs.nodes_file = varvalue;
		}
		else if(varname=="bcs_file") {
			inputs.bcs_file = varvalue;
		}
		else if(varname=="displacement_dofs") {
			vector<string> csv = parse_csv(varvalue);
			for(int i=0;i<csv.size(); i++) {
				int tmp;
				istringstream (csv[i]) >> tmp;
				inputs.displacement_dofs.push_back(tmp);
			}
		}
		else if(varname == "material_files") {
			vector<string> csv = parse_csv(varvalue);
			if(csv.size() != inputs.composite_number) {
				cout << "Error: object " << varname; 
				cout << " has a different size than composite_number!" << endl;
				return 102;
			}
			else {
				for(int idx=0;idx<csv.size();idx++) {
					inputs.material_files.push_back(csv[idx]);
				}
			}
		}
		else if(varname== "total_length") {
			istringstream (varvalue) >> inputs.total_length;
		}
		else if(varname== "a") {
			istringstream (varvalue) >> inputs.a;
		}
		else if(varname== "Omega") {
			istringstream (varvalue) >> inputs.Omega;
		}
		else {
			return 103;
		}
	}
	return 0;
}

int read_material(string material_filename, int matl_idx) {
	ifstream material_fptr;
	material_fptr.open(material_filename.c_str());
	if (!material_fptr.good()){
		cout << "Error: material file " << material_filename << " could not be opened!" << endl; 
		return 200;
	}
	vector <string> material_array;
	material_array.push_back("E");
	material_array.push_back("nu");
	material_array.push_back("rho");
	material_array.push_back("stress_yield");
	material_array.push_back("stress_ultimate");
	material_array.push_back("geometry_file");
	string varname,varvalue,dlm;
	while(material_fptr >> varname >> dlm >> varvalue) {
		// error checking
		vector<string>::const_iterator array_location = find(material_array.begin(), material_array.end(), varname);
		if(array_location == material_array.end()) {
			cout << "Error: invalid entry '" << varname; 
			cout << "' in material file '" << material_filename;
			cout << "'!" << endl;
			return 201;
		}
		else {
			if(varname=="E"){
				istringstream (varvalue) >> materials[matl_idx].E;
			}
			else if(varname=="nu"){
				istringstream (varvalue) >> materials[matl_idx].nu;
			}
			else if(varname=="rho"){
				istringstream (varvalue) >> materials[matl_idx].rho;
			}
			else if(varname=="stress_yield"){
				istringstream (varvalue) >> materials[matl_idx].stress_yield;
			}
			else if(varname=="stress_ultimate"){
				istringstream (varvalue) >> materials[matl_idx].stress_ultimate;
			}
			else if(varname=="geometry_file"){
				materials[matl_idx].geometry_file = varvalue;
			}
			else {
				return 202;
			}
		}
	}
	return 0;
}

int read_geometry(string geometry_filename, int matl_idx){
	ifstream geometry_fptr;
	geometry_fptr.open(geometry_filename.c_str());
	if (!geometry_fptr.good()){
		cout << "Error: geometry file " << geometry_filename << " could not be opened!" << endl; 
		return 300;
	}
	double index,width,height;
	string node_connections;
	int beam_idx = 0;
	while(geometry_fptr >> index >> width >> height >> node_connections) {
		materials[matl_idx].beams.push_back(beam());
		vector<string> csv = parse_csv(node_connections);
		for(int idx=0;idx<csv.size();idx++) {
			int tmp;
			istringstream(csv[idx]) >> tmp;
			materials[matl_idx].beams[beam_idx].node_connections.push_back(tmp);
		}
		materials[matl_idx].beams[beam_idx].index = index;
		materials[matl_idx].beams[beam_idx].width = width;
		materials[matl_idx].beams[beam_idx].height = height;
		beam_idx++;
	}
	return 0;
}

int read_nodes(string nodes_filename){
	ifstream nodes_fptr;
	nodes_fptr.open(nodes_filename.c_str());
	if (!nodes_fptr.good()){
		cout << "Error: nodes file " << nodes_filename << " could not be opened!" << endl; 
		return 400;
	}
	double index;
	string location;
	int node_idx = 0;
	while(nodes_fptr >> index >> location) {
		nodes.push_back(node());
		nodes[node_idx].index = index;
		vector<string> csv = parse_csv(location);
		if( csv.size() > 3) {
			cout << "Error: Parsed more than 3 spatial locations from ";
			cout << "nodes file '" << nodes_filename;
			cout << "' at line: " << node_idx+1 << endl;
			return 401;
		}
		for(int idx=0; idx<csv.size(); idx++) {
			double tmp;
			istringstream (csv[idx]) >> tmp;
			nodes[node_idx].location.push_back(tmp);
		}
		node_idx++;
	}
	return 0;
}

int read_loads(string loads_filename) {
	double machEps = get_Machine_Epsilon();
	ifstream loads_fptr;
	loads_fptr.open(loads_filename.c_str());
	if (!loads_fptr.good()){
		cout << "Error: loads file " << loads_filename << " could not be opened!" << endl; 
		return 500;
	}
	int load_node;
	double force;
	string direction;
	int load_idx=0;
	while(loads_fptr >> load_node >> force >> direction) {
		loads.push_back(load());
		loads[load_idx].node_index = load_node;
		loads[load_idx].force = force;
		vector<string> csv = parse_csv(direction);
		double load_norm=0;//this is the "norm^2", but we don't want interference
		for(int idx=0;idx<csv.size();idx++) {
			double tmp;
			istringstream(csv[idx]) >> tmp;
			loads[load_idx].direction.push_back(tmp);
			load_norm += tmp*tmp;
		}
		if(load_norm > 1.0 + 20*machEps || load_norm < 1.0 - 20*machEps) {
			cout << "Error: The load vector is not normalized!" << endl;
			cout << "File: " << loads_filename << endl;
			cout << "\tLine (Load)#: " << load_idx+1 << endl;
			cout << setprecision(18) << "\tNorm^2: " << load_norm << endl; 
			return 501;
		}
		load_idx++;
	}
	return 0;
}

int read_bcs(string bcs_filename) {
	ifstream bcs_fptr;
	bcs_fptr.open(bcs_filename.c_str());
	if (!bcs_fptr.good()){
		cout << "Error: bc file " << bcs_filename << " could not be opened!" << endl; 
		return 600;
	}
	int bc_node;
	string clamp_str;
	int bc_idx=0;
	while(bcs_fptr >> bc_node >> clamp_str) {
		bcs.push_back(bc());
		bcs[bc_idx].node_id = bc_node;
		vector<string> csv = parse_csv(clamp_str);
		for(int idx=0;idx<csv.size();idx++) {
			double tmp;
			istringstream(csv[idx]) >> tmp;
			bcs[bc_idx].clamp_dir.push_back(tmp);
		}
		bc_idx++;
	}
	return 0;
}

void print_structs(void) {
	print_inputs(inputs);
	print_materials(materials);
	print_loads(loads);
	print_nodes(nodes);
	print_bcs(bcs);
}

void print_inputs(input inputs) {
	cout << "inputs: " <<endl;
	cout << "\tcomposite number: " << inputs.composite_number << endl;
	cout << "\tmaterials files: " << endl;
	for(int i=0; i<inputs.material_files.size(); i++) {
		cout << "\t\t" << inputs.material_files[i] << endl;
	}
	cout << "\tloads file: " << inputs.loads_file << endl;
	cout << "\tnodes file: " << inputs.nodes_file << endl;
	cout << "\tbcs file: " << inputs.bcs_file << endl;
	cout << "\tdisplacement_dofs: ";
	int last_elt = inputs.displacement_dofs.size()-1;
	for(int i=0; i<last_elt; i++) {
		cout << inputs.displacement_dofs[i] << ", ";
	}
	cout << inputs.displacement_dofs[last_elt] << endl;
}

void print_materials(vector<material> materials) {
	for(int i=0; i<materials.size(); i++) {
		cout << "material[" << i << "]:" << endl;
		cout << "\tE: " << materials[i].E*1e-9 << "GPa"<< endl;
		cout << "\tnu: " << materials[i].nu << endl;
		cout << "\tstress_yield: " << materials[i].stress_yield << endl;
		cout << "\tstress_ultimate: " << materials[i].stress_ultimate << endl;
		cout << "\tgeometry_file: " << materials[i].geometry_file << endl;
		cout << "\tbeams: " << materials[i].beams.size() << endl;
	}
}

void print_loads(vector<load> loads) {
	for(int i=0; i<loads.size(); i++) {
		cout << "load[" << i << "]:" << endl;
		cout << "\tnode: " << loads[i].node_index << endl;
		cout << "\tforce: " << loads[i].force << "N"<< endl;
		cout << "\tdirection: ";
		cout << "< " << loads[i].direction[0] << ", ";
		cout << loads[i].direction[1] << ", ";
		cout << loads[i].direction[2] << " >" << endl;
	}
}

void print_bcs(vector<bc> bcs) {
	for(int i=0; i<bcs.size(); i++) {
		cout << "bc[" << i << "]:" << endl;
		cout << "\tnode to clamp: " << bcs[i].node_id << endl;
		cout << "\tclamp dof?: < ";
		int last_elt = bcs[i].clamp_dir.size()-1;
		for (int j=0; j<last_elt; j++) {
			cout << bcs[i].clamp_dir[j] << ", ";
		}
		cout << bcs[i].clamp_dir[last_elt] << " >"<< endl;
	}
}

void print_nodes(vector<node> nodes) {
	cout << "total nodes: " << nodes.size() << endl;
	
}

double get_Machine_Epsilon(void) {
        double machEps = 1.0;
        do
           machEps /= 2.0;
        while ((double) (1.0 + (machEps / 2.0)) != 1.0);
        return machEps;
    }

#endif
