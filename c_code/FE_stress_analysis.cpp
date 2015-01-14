/* 
 * The following library holds all of the logic and tech required to 
 * read in the beam for analysis
*/
#include "get_elements.hpp"
#include <Eigen/Dense>

#include <cstdio>
#include <cstdlib>
#include <sys/time.h>
#include <iostream>

// specify namespace usage
using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* warn compiler about user fns  */

// shape functions
void build_shapes(double xe, double xe1);
void poly_diff(double poly[], int n, double poly_return[]);
void poly_val(double poly[], int n, double x, double poly_return[]);
void poly_shift(double poly[], int n);
double poly_int(double poly[], int n, double a, double b);
void poly_mult(double poly1[],double poly2[], int m, int n, double poly_return[]);
void mult_vecs(int i, int j,vector<VectorXd> input_mx,double arry_out[]);
void build_p(double xe, double xe1, double rho, double A);
void build_m(double xe, double xe1, double rho, double A);
void build_g(double xe, double xe1, double rho, double A);
void build_k(double xe, double xe1, double E, double A, double Iyy, double Izz);
void build_sigma(double xe, double xe1, double len, double rho, double A, double a, double Omega);
void build_global_matrices(void);
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

// rotation parameters
// TODO: need Omega function, here Omega = 0;			
double Omega;
double a;

// builders
void build_stiffness_matrix(void);


/*  end warn compiler */

//declare global variables
int elt_dofs=0;
int nodal_dofs = 0;
int global_dofs = 0;
int u_elt_dofs=0;
int v_elt_dofs=0;
int w_elt_dofs=0;
vector <int> global_delete_index;

// shape vectors
MatrixXd N_u;
MatrixXd N_v;
MatrixXd N_w;
vector<VectorXd> N_uvw;
vector<VectorXd> N_uvw_x;
vector<VectorXd> N_uvw_xx;
// element matrices
MatrixXd p;
MatrixXd m;
MatrixXd g;
MatrixXd k;
MatrixXd sigma;
//MatrixXd f;
// global matrices
MatrixXd P;
MatrixXd M;
MatrixXd G;
MatrixXd K;
MatrixXd Sigma;
// file output strs
const char* P_file = "P.dat";
const char* M_file = "M.dat";
const char* G_file = "G.dat";
const char* K_file = "K.dat";
const char* Sigma_file = "Sigma.dat";
FILE * p_fp;
FILE * m_fp;
FILE * g_fp;
FILE * k_fp;
FILE * sigma_fp;

//params


int main(int argc, char** argv) {
	//check input; if no file specified use a default
	string inputs_file;
	if(argc > 1) {
		inputs_file += argv[1];
	}
	else {
		inputs_file += "inputs.txt";
	}
	double thyme;
	long usecs;
	struct timeval tstart,tend;
	gettimeofday(&tstart,NULL);
	// try to open the file. terminate upon failure.
	int return_error;
	return_error = read_inputs(inputs_file);
	if(return_error != 0)
		return return_error;
	// once the file is open, try to open the data descriptors for each 
	// material (allows for a composite structure).
	for(int i=0; i<inputs.composite_number; i++) {
		// add a material
		materials.push_back(material());
		// read in the material properties
		return_error = read_material(inputs.material_files[i],i);
		if(return_error != 0) 
			return return_error;
		// read in the geometric properties
		return_error = read_geometry(materials[i].geometry_file, i);
		if(return_error != 0) 
			return return_error;	
	}
	// read in the applied loads.
	return_error = read_loads(inputs.loads_file);
	if(return_error != 0)
		return return_error;
	//read in the node locations
	return_error = read_nodes(inputs.nodes_file);
	if(return_error != 0) 
		return return_error;
	return_error = read_bcs(inputs.bcs_file);
	if(return_error != 0) 
		return return_error;
	gettimeofday(&tend,NULL);
	// print what we have read in.
	//print_structs();
	usecs = (tend.tv_sec - tstart.tv_sec) *1000000 + (tend.tv_usec - tstart.tv_usec);
	thyme = usecs*1e-6;
	printf("Time to read inputs: %.5e\n",thyme);
	/*
	 * Here we will begin doing work.
	//test element matrices
	//build_shapes(1.0);
	//build_p(1.0,1.0,420.0);
	//build_m(1.0,1.0,420.0);
	//build_g(1.0,1.0,60.0);
	//build_k(1.0,1.0,1.0,1.0,1.0);
	//build_sigma(1.0,1.0,420.0,1.0,1.0,1.0);
	*/
	a = inputs.a;
	Omega = inputs.Omega;
	
	for(int i=0;i<inputs.displacement_dofs.size();i++) {
		nodal_dofs += inputs.displacement_dofs[i];
	}
	
	global_dofs=nodal_dofs*nodes.size();
	//cout << nodal_dofs << " " <<  nodes.size() << " " << global_dofs << endl;
	
	P = MatrixXd::Zero(global_dofs,global_dofs);
	M = MatrixXd::Zero(global_dofs,global_dofs);
	G = MatrixXd::Zero(global_dofs,global_dofs);
	K = MatrixXd::Zero(global_dofs,global_dofs);
	Sigma = MatrixXd::Zero(global_dofs,global_dofs);

	gettimeofday(&tstart,NULL);
	build_global_matrices();
	gettimeofday(&tend,NULL);

	usecs = (tend.tv_sec - tstart.tv_sec) *1000000 + (tend.tv_usec - tstart.tv_usec);
	thyme = usecs*1e-6;
	printf("Time to build and assemble: %.5e\n",thyme);
	
	gettimeofday(&tstart,NULL);
	
	p_fp = fopen(P_file,"w");
	m_fp = fopen(M_file,"w");
	g_fp = fopen(G_file,"w");
	k_fp = fopen(K_file,"w");
	sigma_fp = fopen(Sigma_file,"w");
	
	if(p_fp==NULL){
		printf("Error: Couldn't open %s for printing!\n",P_file);
		return 1;
	}
	if(m_fp==NULL){
		printf("Error: Couldn't open %s for printing!\n",M_file);
		return 1;
	}
	if(g_fp==NULL){
		printf("Error: Couldn't open %s for printing!\n",G_file);
		return 1;
	}
	if(k_fp==NULL){
		printf("Error: Couldn't open %s for printing!\n",K_file);
		return 1;
	}
	if(sigma_fp==NULL){
		printf("Error: Couldn't open %s for printing!\n",Sigma_file);
		return 1;
	}
	
	for(int i=0; i<P.rows(); i++) {
		for(int j=0; j<P.cols(); j++) {
			fprintf(p_fp,"%.15e\t",P(i,j));
			fprintf(m_fp,"%.15e\t",M(i,j));
			fprintf(g_fp,"%.15e\t",G(i,j));
			fprintf(k_fp,"%.15e\t",K(i,j));
			fprintf(sigma_fp,"%.15e\t",Sigma(i,j));
		}
		fprintf(p_fp,"\n");
		fprintf(m_fp,"\n");
		fprintf(g_fp,"\n");
		fprintf(k_fp,"\n");
		fprintf(sigma_fp,"\n");
	}
	
	fclose(p_fp);
	fclose(m_fp);
	fclose(g_fp);
	fclose(k_fp);
	fclose(sigma_fp);

	
	/*ofstream data_stream_write(P_file);
	data_stream_write.precision(16);
	data_stream_write.setf(ios::scientific);
	data_stream_write << P << endl;
	data_stream_write.close();
	
	data_stream_write.open(M_file);
	data_stream_write.precision(16);
	data_stream_write.setf(ios::scientific);
	data_stream_write << M << endl;
	data_stream_write.close();
	
	data_stream_write.open(G_file);
	data_stream_write.precision(16);
	data_stream_write.setf(ios::scientific);
	data_stream_write << G << endl;
	data_stream_write.close();
	
	data_stream_write.open(K_file);
	data_stream_write.precision(16);
	data_stream_write.setf(ios::scientific);
	data_stream_write << K << endl;
	data_stream_write.close();

	data_stream_write.open(Sigma_file);
	data_stream_write.precision(16);
	data_stream_write.setf(ios::scientific);
	data_stream_write << Sigma << endl;
	data_stream_write.close();*/
	
	gettimeofday(&tend,NULL);

	usecs = (tend.tv_sec - tstart.tv_sec) *1000000 + (tend.tv_usec - tstart.tv_usec);
	thyme = usecs*1e-6;
	printf("Time to write files: %.5e\n",thyme);
	
	return 0;
}

void build_shapes(double xe, double xe1) {
	//build u shape vectors
	double a=0;//xe
	double b=xe1-xe;//xe1;//
	int u_dofs = inputs.displacement_dofs[0];
	int dofs = u_dofs*2;
	u_elt_dofs=dofs;
	double poly[dofs],poly_der[dofs];
	double poly0[dofs],poly_der0[dofs];
	double polyL[dofs],poly_derL[dofs];
	double Array_MX[dofs][dofs];
	MatrixXd U(dofs,dofs);
	MatrixXd Ux(dofs,dofs);
	MatrixXd Uxx(dofs,dofs);
	for(int i=0; i<dofs; i++) {
		poly[i]=1.0;
	}
	poly_val(poly,dofs,a, poly0);//0.0,poly0);
	poly_val(poly,dofs,b,polyL);//l,polyL);
	if(dofs==2) {
		for(int i=0; i<dofs; i++) {
			U(0,i) = poly0[i];
			U(1,i) = polyL[i];
		}
	}
	else {
		poly_diff(poly,dofs,poly_der);
		poly_val(poly_der,dofs,a,poly_der0);//0.0,poly_der0);
		poly_shift(poly_der0,dofs);
		poly_val(poly_der,dofs,b,poly_derL);//l,poly_derL);
		poly_shift(poly_derL,dofs);
		for(int i=0; i<dofs; i++) {
			U(0,i) = poly0[i];
			U(1,i) = poly_der0[i];
			U(2,i) = polyL[i];
			U(3,i) = poly_derL[i];		
		}
	}
	MatrixXd Uinv = U.inverse();

	// build v shape vectors
	int v_dofs = inputs.displacement_dofs[1];
	dofs = v_dofs*2;
	v_elt_dofs=dofs;
	MatrixXd V(dofs,dofs);
	for(int i=0; i<dofs; i++) {
		poly[i]=1.0;
	}
	poly_val(poly,dofs,a,poly0);
	poly_val(poly,dofs,b,polyL);
	if(dofs==2) {
		for(int i=0; i<dofs; i++) {
			V(0,i) = poly0[i];
			V(1,i) = polyL[i];
		}
	}
	else {
		poly_diff(poly,dofs,poly_der);
		poly_val(poly_der,dofs,a,poly_der0);
		poly_shift(poly_der0,dofs);
		poly_val(poly_der,dofs,b,poly_derL);
		poly_shift(poly_derL,dofs);
		for(int i=0; i<dofs; i++) {
			V(0,i) = poly0[i];
			V(1,i) = poly_der0[i];
			V(2,i) = polyL[i];
			V(3,i) = poly_derL[i];		
		}
	}
	MatrixXd Vinv = V.inverse();
	// build w shape vectors
	int w_dofs = inputs.displacement_dofs[2];
	dofs = w_dofs*2;
	MatrixXd W(dofs,dofs);
	for(int i=0; i<dofs; i++) {
		poly[i]=1.0;
	}
	poly_val(poly,dofs,a,poly0);
	poly_val(poly,dofs,b,polyL);
	if(dofs==2) {
		for(int i=0; i<dofs; i++) {
			W(0,i) = poly0[i];
			W(1,i) = polyL[i];
		}
	}
	else {
		poly_diff(poly,dofs,poly_der);
		poly_val(poly_der,dofs,a,poly_der0);
		poly_shift(poly_der0,dofs);
		poly_val(poly_der,dofs,b,poly_derL);
		poly_shift(poly_derL,dofs);
		for(int i=0; i<dofs; i++) {
			W(0,i) = poly0[i];
			W(1,i) = poly_der0[i];
			W(2,i) = polyL[i];
			W(3,i) = poly_derL[i];		
		}
	}
	MatrixXd Winv = W.inverse();
	dofs = u_dofs+v_dofs+w_dofs;
	elt_dofs = 2*dofs;
	w_elt_dofs=dofs;
	int tmp,j;
	for( int i=0; i<2; i++ ) {
		tmp = u_dofs;
		j=0;
		while(tmp>0) {
			MatrixXd val_vec = Uinv.col(i+j);
			MatrixXd diff_vec = MatrixXd::Zero(val_vec.rows(),val_vec.cols());
			MatrixXd diff2_vec = diff_vec;
			int num_vals = val_vec.rows();
			double poly[num_vals], poly_x[num_vals], poly_xx[num_vals];
			for (int i_idx=0;i_idx<num_vals;i_idx++) {
				poly[i_idx] = val_vec(i_idx);
				poly_x[i_idx] = 0;
				poly_xx[i_idx] = 0;
			}
			poly_diff(poly,num_vals,poly_x);
			poly_diff(poly_x,num_vals,poly_xx);
			for (int i_idx=0;i_idx<num_vals;i_idx++) {
				diff_vec(i_idx) = poly_x[i_idx];
				diff2_vec(i_idx) = poly_xx[i_idx];
			}
			N_uvw.push_back(val_vec);
			N_uvw_x.push_back(diff_vec);
			N_uvw_xx.push_back(diff2_vec);
			tmp--;
			j++;
		}
		tmp = v_dofs;
		j=i;
		while(tmp>0) {
			MatrixXd val_vec = Vinv.col(i+j);
			MatrixXd diff_vec = MatrixXd::Zero(val_vec.rows(),val_vec.cols());
			MatrixXd diff2_vec = diff_vec;
			int num_vals = val_vec.rows();
			double poly[num_vals], poly_x[num_vals], poly_xx[num_vals];
			for (int i_idx=0;i_idx<num_vals;i_idx++) {
				poly[i_idx] = val_vec(i_idx);
				poly_x[i_idx] = 0;
				poly_xx[i_idx] = 0;
			}
			poly_diff(poly,num_vals,poly_x);
			poly_diff(poly_x,num_vals,poly_xx);
			for (int i_idx=0;i_idx<num_vals;i_idx++) {
				diff_vec(i_idx) = poly_x[i_idx];
				diff2_vec(i_idx) = poly_xx[i_idx];
			}
			N_uvw.push_back(val_vec);
			N_uvw_x.push_back(diff_vec);
			N_uvw_xx.push_back(diff2_vec);
			tmp--;
			j++;
		}
		tmp = w_dofs;
		j=i;
		while(tmp>0) {
			MatrixXd val_vec = Winv.col(i+j);
			MatrixXd diff_vec = MatrixXd::Zero(val_vec.rows(),val_vec.cols());
			MatrixXd diff2_vec = diff_vec;
			int num_vals = val_vec.rows();
			double poly[num_vals], poly_x[num_vals], poly_xx[num_vals];
			for (int i_idx=0;i_idx<num_vals;i_idx++) {
				poly[i_idx] = val_vec(i_idx);
				poly_x[i_idx] = 0;
				poly_xx[i_idx] = 0;
			}
			poly_diff(poly,num_vals,poly_x);
			poly_diff(poly_x,num_vals,poly_xx);
			for (int i_idx=0;i_idx<num_vals;i_idx++) {
				diff_vec(i_idx) = poly_x[i_idx];
				diff2_vec(i_idx) = poly_xx[i_idx];
			}
			N_uvw.push_back(val_vec);
			N_uvw_x.push_back(diff_vec);
			N_uvw_xx.push_back(diff2_vec);
			tmp--;
			j++;
		}
	}
	N_u = MatrixXd::Zero(elt_dofs,1);
	N_v = MatrixXd::Zero(elt_dofs,1);
	N_w = MatrixXd::Zero(elt_dofs,1);
	
	for (int i=0; i<u_dofs; i++) {
		N_u(i) = 1.0;
		N_u(dofs+i) = 1.0;
	}
	for (int i=0; i<v_dofs; i++) {
		N_v(i+u_dofs) = 1.0;
		N_v(dofs+i+u_dofs) = 1.0;
	}
	for (int i=0; i<w_dofs; i++) {
		N_w(i+u_dofs+v_dofs) = 1.0;
		N_w(dofs+i+(u_dofs+v_dofs)) = 1.0;
	}
	
	/* Print results. Uncomment the following for testing. */
	//TEST PRINTING - SHAPES
	//for (int i=0; i<2*dofs; i++) {
		//cout << "N_uvw[" << i << "]: " << N_uvw[i].transpose() << endl;
	//}
	//for (int i=0; i<2*dofs; i++) {
		//cout << "N_uvw_x[" << i << "]: " << N_uvw_x[i].transpose() << endl;
	//}
	//for (int i=0; i<2*dofs; i++) {
		//cout << "N_uvw_xx[" << i << "]: " << N_uvw_xx[i].transpose() << endl;
	//}
	//cout << "N_u: " << N_u.transpose() << endl;
	//cout << "N_v: " << N_v.transpose() << endl;
	//cout << "N_w: " << N_w.transpose() << endl;
	
	//cout << "Products test [p]: " << endl;
	//cout << N_u*N_u.transpose() + N_v*N_v.transpose() << endl;
	
	return;
}

void build_p(double xe, double xe1, double rho, double A) {
	p = MatrixXd::Zero(elt_dofs,elt_dofs);
	/*
	MatrixXd p_mx = N_u*N_u.transpose() + N_v*N_v.transpose();
	double val;
	for (int i=0; i<elt_dofs; i++ ){
		for (int j=0; j<elt_dofs; j++ ){
			if(p_mx(i,j)) {
				int m = N_uvw[i].rows();
				int n = N_uvw[j].rows();
				double vec[m+n-1];
				for (int idx=0;idx<m+n-1;idx++) {
					vec[idx]=0;
				}
				mult_vecs(i,j,N_uvw,vec);
				val = poly_int(vec, m+n-1, xe, xe1);
				p(i,j) = rho*A*val;
			}
		}
	}
	*/
	
	double l = xe1-xe;
	p << 140,     0,      0, 0, 0,  70,     0,      0, 0, 0,
	       0,   156,   22*l, 0, 0,   0,    54,  -13*l, 0, 0,
	       0,  22*l,  4*l*l, 0, 0,   0,  13*l, -3*l*l, 0, 0,
	       0,     0,      0, 0, 0,   0,     0,      0, 0, 0,
	       0,     0,      0, 0, 0,   0,     0,      0, 0, 0,
	      70,     0,      0, 0, 0, 140,     0,      0, 0, 0,
	       0,    54,   13*l, 0, 0,   0,   156,  -22*l, 0, 0,
	       0, -13*l, -3*l*l, 0, 0,   0, -22*l,  4*l*l, 0, 0,
	       0,     0,      0, 0, 0,   0,     0,      0, 0, 0,
	       0,     0,      0, 0, 0,   0,     0,      0, 0, 0;
	p *= l*rho*A/420;
	
	
}

void build_m(double xe, double xe1, double rho, double A) {
	m = MatrixXd::Zero(elt_dofs,elt_dofs);
	
	/*
	MatrixXd m_mx = N_w*N_w.transpose();
	double val;
	for (int i=0; i<elt_dofs; i++ ){
		for (int j=0; j<elt_dofs; j++ ){
			if(m_mx(i,j)) {
				int m_sz = N_uvw[i].rows();
				int n_sz = N_uvw[j].rows();
				double vec[m_sz+n_sz-1];
				for (int idx=0;idx<m_sz+n_sz-1;idx++) {
					vec[idx]=0;
				}
				mult_vecs(i,j,N_uvw,vec);
				val = poly_int(vec, m_sz+n_sz-1, xe, xe1);//0.0, l);
				m(i,j) = rho*A*val;
			}
		}
	}
	m += p;
	*/
	
	double l = xe1-xe;
	m << 140,     0,      0,     0,      0,  70,     0,      0,     0,      0,
	       0,   156,   22*l,     0,      0,   0,    54,  -13*l,     0,      0,
	       0,  22*l,  4*l*l,     0,      0,   0,  13*l, -3*l*l,     0,      0,
	       0,     0,      0,   156,   22*l,   0,     0,      0,    54,  -13*l,
	       0,     0,      0,  22*l,  4*l*l,   0,     0,      0,  13*l, -3*l*l,
	      70,     0,      0,     0,      0, 140,     0,      0,     0,      0,
	       0,    54,   13*l,     0,      0,   0,   156,  -22*l,     0,      0,
	       0, -13*l, -3*l*l,     0,      0,   0, -22*l,  4*l*l,     0,      0,
	       0,     0,      0,    54,   13*l,   0,     0,      0,   156,  -22*l,
	       0,     0,      0, -13*l, -3*l*l,   0,     0,      0, -22*l,  4*l*l;
	       
	m *= l*rho*A/420;
	
	
}

void build_g(double xe, double xe1, double rho, double A) {
	g = MatrixXd::Zero(elt_dofs,elt_dofs);
	/*
	MatrixXd g_mx = N_v*N_u.transpose()-N_u*N_v.transpose();
	double val;
	for (int i=0; i<elt_dofs; i++ ){
		for (int j=0; j<elt_dofs; j++ ){
			// instead of logic, we need the (-) values...
			int m_sz = N_uvw[i].rows();
			int n_sz = N_uvw[j].rows();
			double vec[m_sz+n_sz-1];
			for (int idx=0;idx<m_sz+n_sz-1;idx++) {
				vec[idx]=0;
			}
			mult_vecs(i,j,N_uvw,vec);
			val = poly_int(vec, m_sz+n_sz-1, xe, xe1);
			g(i,j) = g_mx(i,j)*rho*A*val;
		}
	}
	*/
		
	double l = xe1-xe;
	double aa = 21;
	double bb = 3*l;
	double cc = 2*l;
	double dd = 9;
	g <<    0, -aa, -bb,   0,   0,   0, -dd,  cc,   0,   0,
		   aa,   0,   0,   0,   0,  dd,   0,   0,   0,   0,
		   bb,   0,   0,   0,   0,  cc,   0,   0,   0,   0,
		    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
		    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
		    0, -dd, -cc,   0,   0,   0, -aa,  bb,   0,   0,
		   dd,   0,   0,   0,   0,  aa,   0,   0,   0,   0,
		  -cc,   0,   0,   0,   0, -bb,   0,   0,   0,   0,
		    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
		    0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
		    
	g*=(l*rho*A/60.0);
	
}

void build_k(double xe, double xe1, double E, double A, double Iyy, double Izz){
	k = MatrixXd::Zero(elt_dofs,elt_dofs);
	/*
	MatrixXd u_mx = N_u*N_u.transpose();
	MatrixXd v_mx = N_v*N_v.transpose();
	MatrixXd w_mx = N_w*N_w.transpose();
	double val;
	for (int i=0; i<elt_dofs; i++ ){
		for (int j=0; j<elt_dofs; j++ ){
			// instead of logic, we need the (-) values...
			int m_sz = N_uvw_x[i].rows();
			int n_sz = N_uvw_x[j].rows();
			double vec[m_sz+n_sz-1];
			for (int idx=0;idx<m_sz+n_sz-1;idx++) {
				vec[idx]=0;
			}
			mult_vecs(i,j,N_uvw_x,vec);
			val = poly_int(vec, m_sz+n_sz-1, xe,xe1);//0.0, l);
			k(i,j) += u_mx(i,j)*E*A*val;
			for (int idx=0;idx<m_sz+n_sz-1;idx++) {
				vec[idx]=0;
			}
			mult_vecs(i,j,N_uvw_xx,vec);
			val = poly_int(vec, m_sz+n_sz-1, xe,xe1);//0.0, l);
			k(i,j) += v_mx(i,j)*E*Izz*val;
			k(i,j) += w_mx(i,j)*E*Iyy*val;
		}
	}
	*/
	
	double l = xe1-xe;
	cout << "xe1: " << xe1 << " xe: " << xe << " l: " << l << endl;
	k <<  A*l*l,       0,         0,       0,         0,  -A*l*l,        0,         0,        0,         0,
	          0,  12*Izz,   6*Izz*l,       0,         0,       0,  -12*Izz,   6*Izz*l,        0,         0,
	          0, 6*Izz*l, 4*Izz*l*l,       0,         0,       0, -6*Izz*l, 2*Izz*l*l,        0,         0,
	          0,       0,         0,  12*Iyy,   6*Iyy*l,       0,        0,         0,  -12*Iyy,   6*Iyy*l,
	          0,       0,         0, 6*Iyy*l, 4*Iyy*l*l,       0,        0,         0, -6*Iyy*l, 2*Iyy*l*l,
	     -A*l*l,       0,         0,       0,         0,   A*l*l,        0,         0,        0,         0, 
	          0, -12*Izz,  -6*Izz*l,       0,         0,       0,   12*Izz,  -6*Izz*l,        0,         0,
	          0, 6*Izz*l, 2*Izz*l*l,       0,         0,       0, -6*Izz*l, 4*Izz*l*l,        0,         0,
	          0,       0,         0, -12*Iyy,  -6*Iyy*l,       0,        0,         0,   12*Iyy,  -6*Iyy*l,
	          0,       0,         0, 6*Iyy*l, 2*Iyy*l*l,       0,        0,         0, -6*Iyy*l, 4*Iyy*l*l;
	     
	k *= E/(l*l*l);
	
	         
}

void build_sigma(double xe, double xe1, double len, double rho, double A, double a, double Omega){
	sigma = MatrixXd::Zero(elt_dofs,elt_dofs);
	
	/*
	MatrixXd v_mx = N_v*N_v.transpose();
	MatrixXd w_mx = N_w*N_w.transpose();
	double L = len;
	double l = xe1-xe;
	double val;
	double mult_vec[3] = { a*(L-xe)+0.5*(L*L-xe*xe), -a-xe, -0.5};
	for (int i=0; i<elt_dofs; i++ ){
		for (int j=0; j<elt_dofs; j++ ){
			// instead of logic, we need the (-) values...
			int m_sz = N_uvw_x[i].rows();
			int n_sz = N_uvw_x[j].rows();
			double vec[m_sz+n_sz-1],vec2[m_sz+n_sz-1];
			for (int idx=0;idx<m_sz+n_sz-1;idx++) {
				vec[idx]=0;
				vec2[idx]=0;
			}
			mult_vecs(i,j,N_uvw_x,vec);
			poly_mult(vec,mult_vec,m_sz+n_sz-1,3,vec2);
			val = poly_int(vec2, m_sz+n_sz-1, 0.0, l);
			sigma(i,j) += v_mx(i,j)*val*A*rho;
			sigma(i,j) += w_mx(i,j)*val*A*rho;
		}
	}
	*/
	
	/*
	double L = len;
	double l = xe1-xe;
	double a1 = 420*(3*L*L/5+6*a*L/5)/l-72*l-252*a;
	double a2 = 21*L*L+42*a*L-15*l*l-42*a*l;
	double a3 = 2*l*(14*L*L+28*a*L-2*l*l-7*a*l);
	double b1 = 21*L*L+42*a*L+6*l*l;
	double b2 = -l*(7*L*L+14*a*L-3*l*l-7*a*l);
	double b3 = 2*l*(14*L*L+28*a*L-9*l*l-21*a*l);
	*/
	double L = len;
	double l = xe1-xe;
	double aa = L*L-xe*xe;
	double bb = 2*a*(L-xe);
	double cc = l*(a+xe);
	// (2,2) = (4,4) = (7,7) = (9,9) = -(2,7) = -(7,2) = -(4,9) = -(9,4)
	double a22 = -72*pow(l,2)+252*(aa+bb-cc);
	// (2,3) = (3,2) = (4,5) = (5,4) = -(3,7) = -(7,3) = -(5,9) = -(9,5)
	double b23 = -15*pow(l,3)+21*l*(aa+bb-2*cc);
	// (2,8) = (8,2) = (4,10) = (10,4) = -(7,8) = -(8,7) = -(9,10) = -(10,9)
	double c28 = 6*pow(l,3)+21*l*(aa+bb);
	// (3,3) = (5,5)
	double d33 = -4*pow(l,4)+14*l*l*(2*aa+2*bb-cc);
	// (3,8) = (8,3) = (5,10) = (10,5)
	double e38 = 3*pow(l,4)+7*l*l*(-aa-bb+cc);
	// (8,8) = (10,10)             
	double f88 = -18*pow(l,4)+14*l*l*(2*aa+2*bb-3*cc);
	
	
	sigma << 0,   0,    0,    0,    0,   0,    0,    0,    0,    0,
			 0,  a22,  b23,   0,    0,   0, -a22,  c28,    0,    0,
			 0,  b23,  d33,   0,    0,   0, -b23,  e38,    0,    0,
			 0,   0,    0,  a22,  b23,   0,    0,    0, -a22,  c28,
			 0,   0,    0,  b23,  d33,   0,    0,    0, -b23,  e38,
			 0,   0,    0,    0,    0,   0,    0,    0,    0,    0,
			 0, -a22, -b23,   0,    0,   0,  a22, -c28,    0,    0,
			 0,  c28,  e38,   0,    0,   0, -c28,  f88,    0,    0,
			 0,   0,    0, -a22, -b23,   0,    0,    0,  a22, -c28,
			 0,   0,    0,  c28,  e38,   0,    0,    0, -c28,  f88;
	
			 
	sigma *= (rho*A)/(420*l);//Omega*Omega*
	
}

void poly_diff(double poly[], int n, double poly_return[]) {
	for(int i=0; i<n-1; i++) {
		poly_return[i] = (i+1.0)*poly[i+1];
	}
	poly_return[n-1]=0.0;
	return;
}

double poly_int(double poly[], int n, double a, double b) {
	double I=0.0;
	for(int i=0; i<n; i++) {
		I += (poly[i]/(i+1.0))*(pow(b,i+1.0)-pow(a,i+1.0));
	}
	return I;
}

void poly_val(double poly[], int n, double x, double poly_return[]) {
	for(int i=1; i<n; i++) {
		poly_return[i] = pow(x,i)*poly[i];
	}
	poly_return[0] = poly[0];
	return;
}

void poly_shift(double poly[], int n) {
		/*  We now wish to shift values in a polynomial.
		 *  Why? So that the values match with the corresponding
		 *  coefficients for a linear system.
		 *  This is typically necessary for derivatives.
		 *                                                      */
	for(int i=n-1; i>0; i--) {
		poly[i] = poly[i-1];
	}
	poly[0]=0.0;
	return;
}

void poly_mult(double poly1[],double poly2[], int m, int n, double poly_return[]) {
	//My implementation
	for(int i=0; i<m; i++) {
		for(int j=0;j<n; j++) {
			poly_return[i+j] += poly1[i]*poly2[j];
		}
	} 
	// Matlab's implementation
	/* for(int i=0; i<m+n-1; i++) {
		int minj = (0>i-n+1)?0:i-n+1;
		int maxj = (i<m)?i:m;
		double tmp_sum=0;
		for(int j=minj; j<=maxj; j++) {
			tmp_sum += poly1[j]*poly2[i-j+1];
		}
		poly_return[i] = tmp_sum;
	}
	*/
	return;
}

void mult_vecs(int i, int j,vector<VectorXd> input_mx,double arry_out[]) {
	int m = input_mx[i].rows();
	int n = input_mx[j].rows();
	double arry1[m];
	double arry2[n];
	MatrixXd vec1 = input_mx[i];
	MatrixXd vec2 = input_mx[j];
	for( int idx=0; idx<m; idx++) {
		arry1[idx] = vec1(idx);
	}
	for( int idx=0; idx<n; idx++) {
		arry2[idx] = vec2(idx);
	}
	poly_mult(arry1,arry2,m,n,arry_out);
	return;
}

void build_global_matrices(void){
	int global_index[2*nodal_dofs];
	for (int clear_idx=0; clear_idx<2*nodal_dofs; clear_idx++) {
		global_index[clear_idx]=0;
	}
	int dims = nodes[0].location.size();
	double pos_array[dims][2];
	double diff_array[dims];
	double cos_array[dims];
	for(int mat_idx=0; mat_idx<materials.size(); mat_idx++) {
		for(int bar_idx=0; bar_idx<materials[mat_idx].beams.size(); bar_idx++) {
			for(int node_idx=0; node_idx<2; node_idx++) {
				int nid =  materials[mat_idx].beams[bar_idx].node_connections[node_idx];
				for (int n_idx = 0; n_idx<nodal_dofs; n_idx++) {
					global_index[nodal_dofs*node_idx+n_idx] = nodal_dofs*(nid-1)+n_idx;
				}
				for(int pos_idx=0; pos_idx<dims; pos_idx++) {
					pos_array[pos_idx][node_idx] = nodes[nid-1].location[pos_idx];
				}
			}

			//calculate element length
			double l=0;
			for(int diff_idx=0; diff_idx<dims; diff_idx++) {
				diff_array[diff_idx] = pos_array[diff_idx][1]-pos_array[diff_idx][0];
				l += diff_array[diff_idx]*diff_array[diff_idx];
			}
			l = sqrt(l);
			for(int diff_idx=0; diff_idx<dims; diff_idx++) {
				cos_array[diff_idx] = diff_array[diff_idx]/l;
			}
			Eigen::MatrixXd T = Eigen::MatrixXd::Zero(2,6);
			T << cos_array[0], cos_array[1], cos_array[2],            0,            0,            0,
					        0,            0,            0, cos_array[0], cos_array[1], cos_array[2];
			Eigen::MatrixXd u_global = Eigen::MatrixXd::Zero(6,1);
			u_global << pos_array[0][0], pos_array[1][0], pos_array[2][0],
			            pos_array[0][1], pos_array[1][1], pos_array[2][1];
			Eigen::MatrixXd u_local = T*u_global;
			double xe = u_local(0);//pos_array[0][0];//
			double xe1 = u_local(1);//pos_array[0][1];//
			//build element matrices
			build_shapes(xe,xe1);
			double breadth = materials[mat_idx].beams[bar_idx].width;
			double depth = materials[mat_idx].beams[bar_idx].height;
			double A = breadth*depth;
			double Izz = breadth*pow(depth,3.0)/12.0;
			double Iyy = depth*pow(breadth,3.0)/12.0;
			build_p(xe,xe1,materials[mat_idx].rho,A);
			build_m(xe,xe1,materials[mat_idx].rho,A);
			build_g(xe,xe1,materials[mat_idx].rho,A);			
			build_k(xe,xe1,materials[mat_idx].E,A,Iyy,Izz);	
			// a and Omega are currently globals	
			build_sigma(xe,xe1,inputs.total_length,materials[mat_idx].rho,A,a,Omega);
			// Set global values
			for(int i=0;i<2*nodal_dofs;i++) {
				for(int j=0; j<2*nodal_dofs; j++) {
					//cout << "global_index(" << i << "," << j << "): (" << global_index[i] << "," << global_index[j] << ")" << endl;
					P(global_index[i],global_index[j]) += p(i,j);
					M(global_index[i],global_index[j]) += m(i,j);
					G(global_index[i],global_index[j]) += g(i,j);
					K(global_index[i],global_index[j]) += k(i,j);
					Sigma(global_index[i],global_index[j]) += sigma(i,j);
				}
			}
		}
	}
	// use input BCs to decide which rows/cols to delete
	global_delete_index.clear();
	for(int i=0; i<bcs.size(); i++) {
		for(int j=0; j<bcs[i].clamp_dir.size(); j++) {
			if(bcs[i].clamp_dir[j]) {
				global_delete_index.push_back(nodal_dofs*(bcs[i].node_id-1)+j);
			}
		}
	}
	// sort BCs so we don't have to reindex
	sort(global_delete_index.rbegin(), global_delete_index.rend());
	// apply BCs
	for (int i=0; i<global_delete_index.size(); i++) {
		removeColumn(P,global_delete_index[i]);
		removeRow(P,global_delete_index[i]);
		removeColumn(M,global_delete_index[i]);
		removeRow(M,global_delete_index[i]);
		removeColumn(G,global_delete_index[i]);
		removeRow(G,global_delete_index[i]);
		removeColumn(K,global_delete_index[i]);
		removeRow(K,global_delete_index[i]);
		removeColumn(Sigma,global_delete_index[i]);
		removeRow(Sigma,global_delete_index[i]);
	}
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

//eof
