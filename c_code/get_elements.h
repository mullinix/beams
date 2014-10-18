#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>

/*includes and macros for later*/
//#include "Eigen/Dense"
//#include <cmath>
//#define PI 3.141592653589793238462643

/*Space of variable and function names in which we shall operate*/
using namespace std;
//using namespace Eigen;

//global structures for manipulation
struct input { //struct for the inputs file data
	int composite_number;
	vector<string> material_files;
	string loads_file;
	string nodes_file;
	string bcs_file;
	vector <int> displacement_dofs;
} inputs;

struct node { //struct for the nodes with index and location
	int index;
	double x;
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
int read_nodes(string nodes_filename, int matl_idx);
int read_loads(string loads_filename);
void print_inputs(input inputs);
void print_materials(vector<material> materials);
void print_loads(vector<load> loads);
void print_structs(input inputs, vector<material> materials, vector<load> loads);

//construct global variables
vector<material> materials;
vector<load> loads;
vector<node> nodes;
vector<bc> bcs;
