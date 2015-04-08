#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <time.h>
#include <math.h>
#include <mex.h>
#define intPts 10
#define PI 3.141592653589793238462643
#define EPS 1e-14

/*Space of variable and function names in which we shall operate*/
using namespace std;
using namespace Eigen;

/* These are the functions which perform the gaussian quadrature using 
 * legendre polynomial functions
 * */
void lege_coef();
double lege_eval(int n, double x);
double lege_diff(int n, double x);
void lege_roots();
double lege_inte(double (*f)(double,Vector4d), double a, double b, Vector4d params);

/*Eigen helper functions*/
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
int getParams();

/* System parameters */
double L = 0.1524; // length of beam [m] (0.45)
clock_t init, dt, tfinal;//ttemp,
int N = getParams(); // number of elements [1]

//~ N = getParams();
const double g = 9.81; // accel due to gravity ()
const double rho = 800;//7850;// // density of steel [kg/m^3] (steel: 7850, Alum: 2700)
// m = L*w*h*rho; // m: total beam mass [kg]
const double E = 2.07e11; // E: Young's modulus of steel [N/m^2] (steel: 21e10, alum: 7.2e10)
const double G = E/2.6;//19.3e9; // Shear modulus of Pascals (steel: 19.3e9)
// nu = 3/10; // nu: Poisson's ratio (steel: 0.3, alum: 0.35)
const double mu = 0.833;// 10*(1+nu)/(12+11*nu); // "shear coefficient" for rectangular cross section (approximate) http://en.wikipedia.org/wiki/Timoshenko_beam_theory#Shear_coefficient

// // rectangular beam
const double w = 0.0254;//0.0254;// // w: beam width (lateral) [m] (0.02)
const double h = 0.0046;//w*.25;// // h: beam height (vertical) [m] (0.003)

// R&G element parameters
// allow for twist
// nodal angular displacement wrt origin
VectorXd theta = VectorXd::LinSpaced(N,0.0,45*PI/180.0);//linspace(0,45*pi/180,N); 
VectorXd z = VectorXd::LinSpaced(N,0.0,L);//linspace(0,L,N); // physical locations of nodes from origin
// allow for varying cross section
const double breadth_taper = 2.56;//1;//
const double depth_taper = 2.29;//1;//
VectorXd breadth = VectorXd::LinSpaced(N,w,w/breadth_taper);// distribution - use for "b" in RG elemnts
VectorXd depth = VectorXd::LinSpaced(N,h,h/depth_taper);// distribution - use for "h" in RG elements
const double Omega = 100;//6000*60;// // rotation speed, rps - rpm*(secs/min)
const double offset = 0.0; // distance from center of rotation

// params for calculating cross sections and moments of inertia
MatrixXd c = MatrixXd::Zero(N,3);
MatrixXd a = MatrixXd::Zero(N,5);
MatrixXd d = MatrixXd::Zero(N,5);
// lengths at each node
MatrixXd l = MatrixXd::Zero(N,1);
// coupled moments of inertia
MatrixXd Ixx = MatrixXd::Zero(N,1); // full actualized moment
MatrixXd Ixxp = MatrixXd::Zero(N,1); // moment before twist 
MatrixXd Iyy = MatrixXd::Zero(N,1);
MatrixXd Iyyp = MatrixXd::Zero(N,1); 
MatrixXd Ixy = MatrixXd::Zero(N,1);

// tabulated parameters, tables B1 & B2
Matrix4d H;
MatrixXd R(12,4);
MatrixXd Q(20,4);
MatrixXd P(28,4);
/* end system parameters */


/* These are the supporting functions for the sub-matrix assembly */
double Ui(int i,double l);
double Li(int i,double l);
double ViFun(double x,Vector4d params);
double Vi(int i,double l,double theta1,double theta2);
double SiFun(double x,Vector4d params);
double Si(int i,double l,double theta1,double theta2);

/* These are the generating functions for the 4x4 matrices in the 
 * element matrix
 * */
Matrix4d AK(int idx);
Matrix4d BK(int idx);
Matrix4d CK(int idx);	
Matrix4d DK(int idx);
Matrix4d AM(int idx);
Matrix4d BM(int idx);
Matrix4d CM(int idx);
Matrix4d DM(int idx);
Matrix4d EK(int idx);
Matrix4d FK(int idx);

/* Declaring global variables for the gaussian quadrature so that 
 * we don't have to recompute the weights and positions for every element
 * */
double lroots[intPts];
double weight[intPts];
double lcoef[intPts + 1][intPts + 1] = {{0}};

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	//~ double multiplier;      /* input scalar */
	//~ double *inMatrix;       /* 1xN input matrix */
	//~ mwSize ncols;           /* size of matrix */
	//~ 
	//~ double *outMatrix;      /* output matrix */
/* variable declarations here */
	/* initialize integration parameters */	
	lege_coef();
	lege_roots();	
	//timing params
	dt=0;
	init = clock();
	// Ixyp - pretwisted moment is uncoupled, so this is zero.
	/* Moments of Inertia and Area
	* calculation of areas and inertias along nonuniform twisted beam
	* Rao & Gupta p.107
	* */
	double bb,dd;
	for(int idx = 1; idx < N; idx++ ) {//matlab index = 2:N, c++ = 0:N-1
		bb = breadth(idx)-breadth(idx-1);
		dd = depth(idx)-depth(idx-1);
		// c 
		c(idx,0) = bb*dd;
		c(idx,1) = breadth(idx-1)*dd+depth(idx-1)*bb;
		c(idx,2) = breadth(idx-1)*depth(idx-1);
		// a
		a(idx,0) = bb*pow(dd,3);
		a(idx,1) = breadth(idx-1)*pow(dd,3)+3*bb*pow(dd,2)*depth(idx-1);
		a(idx,2) = 3*(breadth(idx-1)*depth(idx-1)*pow(dd,3)+bb*dd*pow(depth(idx-1),2));
		a(idx,3) = 3*breadth(idx-1)*pow(depth(idx-1),2)*dd+bb*pow(depth(idx-1),3);
		a(idx,4) = breadth(idx-1)*pow(depth(idx-1),3);
		// d
		d(idx,0) = dd*pow(bb,3);
		d(idx,1) = depth(idx-1)*pow(bb,3)+3*dd*pow(bb,2)*breadth(idx-1);
		d(idx,2) = 3*(breadth(idx-1)*depth(idx-1)*pow(bb,3)+dd*bb*pow(breadth(idx-1),2));
		d(idx,3) = 3*depth(idx-1)*pow(breadth(idx-1),2)*bb+dd*pow(breadth(idx-1),3);
		d(idx,4) = depth(idx-1)*pow(breadth(idx-1),3);
		// element node legths
		l(idx) = z(idx)-z(idx-1);
		// inertial terms (commented should work for non-rectangular)
		Ixxp(idx) = (breadth(idx)*pow(depth(idx),3.0))/12.0; //1/(12*l(idx)^4)*(a(idx,1)*z(idx)^4+...
	//                 a(idx,2)*l(idx)*z(idx)^3+a(idx,3)*l(idx)^2*z(idx)^2+...
	//                 a(idx,4)*l(idx)^3*z(idx)+a(idx,5)*l(idx));//(breadth(idx)*depth(idx)^3)/12; //
		Iyyp(idx) = (depth(idx)*pow(breadth(idx),3.0))/12.0; //1/(12*l(idx)^4)*(d(idx,1)*z(idx)^4+...
	//                 d(idx,2)*l(idx)*z(idx)^3+d(idx,3)*l(idx)^2*z(idx)^2+...
	//                 d(idx,4)*l(idx)^3*z(idx)+d(idx,5)*l(idx));//(depth(idx)*breadth(idx)^3)/12; //
		Ixx(idx) = Ixxp(idx)*(pow(cos(theta(idx)),2.0))+Iyyp(idx)*(pow(sin(theta(idx)),2.0));
		Iyy(idx) = Iyyp(idx)*(pow(cos(theta(idx)),2.0))+Ixxp(idx)*(pow(sin(theta(idx)),2.0));
		Ixy(idx) = (Ixxp(idx)-Iyyp(idx))*(sin(2.0*theta(idx))/2.0);
	}
	//~ cout << "a:" << endl; cout << a << endl;
	//~ cout << "d:" << endl; cout << d << endl;
	//~ cout << "l:" << endl; cout << l << endl;
	//~ cout << "theta:" << endl; cout << theta << endl;

	//// Rao & Gupta parameters
	// element parameters

	// Table B1
	H <<           0.0, 0.0, 1.0, 1.0,
				   0.0, 0.0, 1.0, 1.0,
				   1.0, 1.0, 2.0, 2.0,
				   1.0, 1.0, 2.0, 2.0;
	  
	R <<         144.0,-144.0,-72.0,-72.0,
				-144.0, 144.0, 72.0, 72.0,
				 -72.0,  72.0, 36.0, 36.0,
				 -72.0,  72.0, 36.0, 36.0,//j=1
				-144.0, 144.0, 84.0, 60.0,
				 144.0,-144.0,-84.0,-60.0,
				  84.0, -84.0,-48.0,-36.0,
				  60.0, -60.0,-36.0,-24.0,//j=2
				  36.0, -36.0,-24.0,-12.0,
				 -36.0,  36.0, 24.0, 12.0,
				 -24.0,  24.0, 16.0,  8.0,
				 -12.0,  12.0,  8.0,  4.0;//j=3
			   
	Q <<          36.0, -36.0,-18.0,-18.0,
				 -36.0,  36.0, 18.0, 18.0,
				 -18.0,  18.0,  9.0,  9.0,
				 -18.0,  18.0,  9.0,  9.0,//j=1
				 -72.0,  72.0, 42.0, 30.0,
				  72.0, -72.0,-42.0,-30.0,
				  42.0, -42.0,-24.0,-18.0,
				  30.0, -30.0,-18.0,-12.0,//j=2
				  36.0, -36.0,-30.0,-12.0,
				 -36.0,  36.0, 30.0, 12.0,
				 -30.0,  30.0, 22.0, 11.0,
				 -12.0,  12.0, 11.0,  4.0,//j=3
				   0.0,   0.0,  6.0,  0.0,
				   0.0,   0.0, -6.0,  0.0,
				   6.0,  -6.0, -8.0, -2.0,
				   0.0,   0.0, -2.0,  0.0,//j=4
				   0.0,   0.0,  0.0,  0.0,
				   0.0,   0.0,  0.0,  0.0,
				   0.0,   0.0,  1.0,  0.0,
				   0.0,   0.0,  0.0,  0.0;//j=5
			   
	// Table B2
	P <<           4.0,  -4.0, -2.0, -2.0,
				  -4.0,   4.0,  2.0,  2.0,
				  -2.0,   2.0,  1.0,  1.0,
				  -2.0,   2.0,  1.0,  1.0,//j=1
				 -12.0,  12.0,  7.0,  5.0,
				  12.0, -12.0, -7.0, -5.0,
				   7.0,  -7.0, -4.0, -3.0,
				   5.0,  -5.0, -3.0, -2.0,//j=2
				   9.0,  -9.0, -8.0, -3.0,
				  -9.0,   9.0,  8.0,  3.0,
				  -8.0,   8.0,  6.0,  3.0,
				  -3.0,   3.0,  3.0,  1.0,//j=3
				   4.0,  -2.0,  2.0, -1.0,
				  -2.0,   0.0, -3.0,  0.0,
				   2.0,  -3.0, -4.0, -1.0,
				  -1.0,   0.0, -1.0,  0.0,//j=4
				  -6.0,   3.0,  2.0,  1.0,
				   3.0,   0.0,  0.0,  0.0,
				   2.0,   0.0,  1.0,  0.0,
				   1.0,   0.0,  0.0,  0.0,//j=5
				   0.0,   0.0, -1.0,  0.0,
				   0.0,   0.0,  0.0,  0.0,
				  -1.0,   0.0,  0.0,  0.0,
				   0.0,   0.0,  0.0,  0.0,//j=6
				   1.0,   0.0,  0.0,  0.0,
				   0.0,   0.0,  0.0,  0.0,
				   0.0,   0.0,  0.0,  0.0,
				   0.0,   0.0,  0.0,  0.0;//j=7
	//// Assemble matrices
	MatrixXd Vk,Wk,Rot,Coup,Shear,Vm,Wm,Am,Dm;
	Vk = MatrixXd::Zero(2*(N+1),2*(N+1));
	Wk = MatrixXd::Zero(2*(N+1),2*(N+1));
	Rot = MatrixXd::Zero(2*(N+1),2*(N+1));
	Coup = MatrixXd::Zero(2*(N+1),2*(N+1));
	Shear = MatrixXd::Zero(2*(N+1),2*(N+1));
	Vm = MatrixXd::Zero(2*(N+1),2*(N+1));
	Wm = MatrixXd::Zero(2*(N+1),2*(N+1));
	Am = MatrixXd::Zero(2*(N+1),2*(N+1));
	Dm = MatrixXd::Zero(2*(N+1),2*(N+1));
	Matrix4d R1,R2,R3,R123,R321;
	R1 <<     1,0,0,0,
			  0,0,1,0,
			  0,1,0,0,
			  0,0,0,1;
	R2 <<     0,1,0,0,
			  1,0,0,0,
			  0,0,1,0,
			  0,0,0,1;
	R3 <<     1,0,0,0,
			  0,1,0,0,
			  0,0,0,1,
			  0,0,1,0;
	R123 = R1*R2*R3;
	R321 = R3*R2*R1;
	int j_idx,k_idx;
	Matrix4d ak,bk,ck,dk,am,bm,cm,dm,ek,fk;
	for(int elt=1;elt<N;elt++) {
		j_idx = 2*(elt)-1;//Starting location of the node
		k_idx = 4;//Number to traverse in the matrix
	//     j_idx = 8*(elt)-7;
	// 	k_idx = j_idx+15;
		//~ otra_j = 4*elt-7;
		//~ otra_k = otra_j+3;
		ak = R123*AK(elt)*R321;
		bk = R123*BK(elt)*R321;
		ck = R123*CK(elt)*R321;
		dk = R123*DK(elt)*R321;
		am = R123*AM(elt)*R321;
		bm = R123*BM(elt)*R321;
		cm = R123*CM(elt)*R321;
		dm = R123*DM(elt)*R321;
		ek = R123*EK(elt)*R321;
		fk = R123*FK(elt)*R321;
	/*
	//     Kel = [ak+ek-fk,      ek-fk,       dk,  zeros(4);
	//               ek-fk,   ck+ek-fk, zeros(4),  zeros(4);
	//                  dk,   zeros(4), bk+ek-fk,     ek-fk;
	//            zeros(4),   zeros(4),    ek-fk,  ck+ek-fk];
	//        
	//     Mel = [   am+bm,         am,       dm,  zeros(4);
	//                  am,         am,       am,  zeros(4);
	//                  dm,         am,    am+cm,        am;
	//            zeros(4),   zeros(4),       am,        am];
	*/       
		if(elt==1){
			  Vk.block(0,0,4,4) = ak+ek-fk;
			  Wk.block(0,0,4,4) = bk+ek-fk;
			  Rot.block(0,0,4,4) = ek-fk;
			  Shear.block(0,0,4,4) = ck+ek-fk;
			  Coup.block(0,0,4,4) = dk;
			  Vm.block(0,0,4,4) = am+bm;
			  Wm.block(0,0,4,4) = am+cm;
			  Am.block(0,0,4,4) = am;
			  Dm.block(0,0,4,4) = dm;
		}
		  /* uncomment for testing: */
		  /*
	//     K_factor = E*Ixx(elt)/(l(elt)^2);// should be l^3 in denominator!
	//     JJ = mu*G*breadth(elt)*depth(elt)*l(elt)^2/(30*E*Ixx(elt));
	//     fprintf(1,'6l: //.5e\n',6*l(elt));
	//     fprintf(1,'2l^2: //.5e\n',2*l(elt)^2);
	//     fprintf(1,'36J: //.5e\n',36*JJ);
	//     M_factor = (rho*breadth(elt)*depth(elt))/(420*g)*l(elt)^2;// should be times Ixx, not l^2!
	//     PP = (14*Ixx(10;//elt))/(breadth(elt)*depth(elt)*l(elt)^2);
		  */
		 // assembly
		Vk.block(j_idx,j_idx,k_idx,k_idx) += ak+ek-fk;
		Wk.block(j_idx,j_idx,k_idx,k_idx) = bk+ek-fk;
		Rot.block(j_idx,j_idx,k_idx,k_idx) = ek-fk;
		Shear.block(j_idx,j_idx,k_idx,k_idx) = ck+ek-fk;
		Coup.block(j_idx,j_idx,k_idx,k_idx) = dk;
		Vm.block(j_idx,j_idx,k_idx,k_idx) = am+bm;
		Wm.block(j_idx,j_idx,k_idx,k_idx) = am+cm;
		Am.block(j_idx,j_idx,k_idx,k_idx) = am;
		Dm.block(j_idx,j_idx,k_idx,k_idx) = dm;
	}
	removeRow(Vk,1);removeRow(Vk,2);
	removeColumn(Vk,1);removeColumn(Vk,2);

	removeRow(Wk,1);removeRow(Wk,2);
	removeColumn(Wk,1);removeColumn(Wk,2);

	removeRow(Rot,1);removeRow(Rot,2);
	removeColumn(Rot,1);removeColumn(Rot,2);

	removeRow(Shear,1);removeRow(Shear,2);
	removeColumn(Shear,1);removeColumn(Shear,2);

	removeRow(Coup,1);removeRow(Coup,2);
	removeColumn(Coup,1);removeColumn(Coup,2);

	removeRow(Vm,1);removeRow(Vm,2);
	removeColumn(Vm,1);removeColumn(Vm,2);

	removeRow(Wm,1);removeRow(Wm,2);
	removeColumn(Wm,1);removeColumn(Wm,2);

	removeRow(Am,1);removeRow(Am,2);
	removeColumn(Am,1);removeColumn(Am,2);

	removeRow(Dm,1);removeRow(Dm,2);
	removeColumn(Dm,1);removeColumn(Dm,2);

	MatrixXd mt;
	mt = MatrixXd::Zero(Vk.rows(),Vk.cols());
	MatrixXd K(8*N,8*N);
	K <<     Vk,   Rot,  Coup,    mt,
			Rot, Shear,    mt,    mt,
		   Coup,    mt,    Wk,   Rot,
			 mt,    mt,   Rot, Shear;
		
	MatrixXd M(8*N,8*N);
	M <<  Vm, Am, Dm, mt,
		  Am, Am, Am, mt,
		  Dm, Am, Wm, Am,
		  mt, mt, Am, Am;

	tfinal=clock()-init;
//~ 
	cout << "Assembly time:" << endl;
	cout << (double)tfinal / ((double)CLOCKS_PER_SEC) << endl;
	init = clock();
	ofstream myfile;
	myfile.open ("k.dat");
	myfile << K;
	myfile.close();
	myfile.open ("m.dat");
	myfile << M;
	myfile.close();
	tfinal=clock()-init;
	cout << "File write time:" << endl;
	cout << (double)tfinal / ((double)CLOCKS_PER_SEC) << endl;
}///////// END MAIN //////////////////////

/* gauss legendre integration, rosettacode.org */
/* http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature */
void lege_coef()
{
	int n, i;
	lcoef[0][0] = lcoef[1][1] = 1;
	for (n = 2; n <= intPts; n++) {
		lcoef[n][0] = -(n - 1) * lcoef[n - 2][0] / n;
		for (i = 1; i <= n; i++)
			lcoef[n][i] = ((2 * n - 1) * lcoef[n - 1][i - 1]
					 - (n - 1) * lcoef[n - 2][i] ) / n;
	}
}
double lege_eval(int n, double x)
{
	int i;
	double s = lcoef[n][n];
	for (i = n; i; i--)
		s = s * x + lcoef[n][i - 1];
	return s;
}
double lege_diff(int n, double x)
{
	return n * (x * lege_eval(n, x) - lege_eval(n - 1, x)) / (x * x - 1);
}
void lege_roots()
{
	int i;
	double x, x1;
	for (i = 1; i <= intPts; i++) {
		x = cos(PI * (i - .25) / (intPts + .5));
		do {
			x1 = x;
			x -= lege_eval(intPts, x) / lege_diff(intPts, x);
		} while (fabs(x-x1)>EPS);
		/*  x != x1 is normally a no-no, but this task happens to be 
		 *  well behaved. */
		lroots[i - 1] = x;
 
		x1 = lege_diff(intPts, x);
		weight[i - 1] = 2 / ((1 - x * x) * x1 * x1);
	}
}
double lege_inte(double (*f)(double, Vector4d), double a, double b, Vector4d params)
{
	double c1 = (b - a) / 2, c2 = (b + a) / 2, sum = 0;
	int i;
	for (i = 0; i < intPts; i++)
		sum += weight[i] * f(c1 * lroots[i] + c2,params);
	return c1 * sum;
}
/* end gauss legendre integration, rosettacode.org */

/* Rao & Gupta Subfunctions */
double Ui(int i,double l) {
	return (1.0/(i))*pow(l,i);
}
double Li(int i,double l) {
	return pow(l,(i-1));
}
double ViFun(double x,Vector4d params){
	int i=params(0);
	double l=params(1);
	double theta1=params(2);
	double theta2=params(3);
	return pow(x,(i-1))*pow((cos((theta2-theta1)*(x/l)+theta1)),2);
}
double Vi(int i,double l,double theta1,double theta2){
	Vector4d params;
	params << i,l,theta1,theta2;
	return lege_inte(ViFun,0,l,params);
}
double SiFun(double x,Vector4d params){
	int i=params(0);
	double l=params(1);
	double theta1=params(2);
	double theta2=params(3);
	return pow(x,(i-1))*pow((sin((theta2-theta1)*(x/l)+theta1)),2);
}
double Si(int i,double l,double theta1,double theta2){
	Vector4d params;
	params << i,l,theta1,theta2;
	return lege_inte(SiFun,0,l,params);
}

/* Rao & Gupta 4x4 sub-matrices */
// stiffness
Matrix4d AK(int idx) {
	Matrix4d ak;
	ak = MatrixXd::Zero(4,4);
	double ak_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ){
			ak_sum=0;
			for(int i=1; i<=5; i++ ) {
				for(int j=1; j<=3; j++ ) {
					ak_sum += d(idx,i-1)*
							  Li(i+j+H(I,J),l(idx))*Ui(9-i-j,l(idx))*
							  R(I+4*(j-1),J)+(a(idx,i-1)-d(idx,i-1))*
							  (Li(i+j+H(I,J),l(idx))*Vi(9-i-j,l(idx),\
							  theta(idx-1),theta(idx))*R(I+4*(j-1),J));
				}
			}
			ak(I,J) = E/(12*pow(l(idx),10))*ak_sum;
		}
	}
	return ak;
}
Matrix4d BK(int idx) {
	Matrix4d bk;
	bk =  MatrixXd::Zero(4,4);
	double bk_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			bk_sum=0;
			for(int i=1; i<=5; i++ ) {
				for(int j=1; j<=3; j++ ) {
					bk_sum += a(idx,i-1)*
							  Li(i+j+H(I,J),l(idx))*Ui(9-i-j,l(idx))*
							  R(I+4*(j-1),J)+(d(idx,i-1)-a(idx,i-1))*
							  (Li(i+j+H(I,J),l(idx))*Vi(9-i-j,l(idx),\
							  theta(idx-1),theta(idx))*R(I+4*(j-1),J));
				}
			}
			bk(I,J) = E/(12*pow(l(idx),10))*bk_sum;
		}
	}
	return bk;
}
Matrix4d CK(int idx) {
	Matrix4d ck;
	ck =  MatrixXd::Zero(4,4);
	double ck_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			ck_sum=0;
			for(int i=1; i<=3; i++ ) {
				for(int j=1; j<=5; j++ ) {
					ck_sum += c(idx,i-1)*
							  Li(i+j+H(I,J),l(idx))*Ui(9-i-j,l(idx))*
							  Q(I+4*(j-1),J);
				}
			}
			ck(I,J) = mu*G/(pow(l(idx),8))*ck_sum;
		}
	}
	return ck;
}
Matrix4d DK(int idx) {
	Matrix4d dk;
	dk =  MatrixXd::Zero(4,4);
	double dk_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			dk_sum=0;
			for(int i=1; i<=5; i++ ) {
				for(int j=1; j<=3; j++ ) {
					dk_sum += (a(idx,i-1)-d(idx,i-1))*
							  (Li(i+j+H(I,J),l(idx))*Si(9-i-j,l(idx),\
							  theta(idx-1),theta(idx))*R(I+4*(j-1),J));
				}
			}
			dk(I,J) = E/(12*pow(l(idx),10))*dk_sum;
		}
	}
	return dk;
}
// mass
Matrix4d AM(int idx) {
	Matrix4d am;
	am =  MatrixXd::Zero(4,4);
	double am_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			am_sum=0;
			for(int i=1; i<=3; i++ ) {
				for(int j=1; j<=7; j++ ) {
					am_sum += c(idx,i-1)*
							  Li(i+j+H(I,J),l(idx))*
							  Ui(11-i-j,l(idx))*P(I+4*(j-1),J);
				}
			}
			am(I,J) = rho/(g*pow(l(idx),8))*am_sum;
		}
	}
	return am;
}
Matrix4d BM(int idx){
	Matrix4d bm;
	bm =  MatrixXd::Zero(4,4);
	double bm_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			bm_sum=0;
			for(int i=1; i<=5; i++ ) {
				for(int j=1; j<=5; j++ ) {
					bm_sum += d(idx,i-1)*Li(i+j+H(I,J),l(idx))*
							  Ui(11-i-j,l(idx))*Q(I+4*(j-1),J)+
							  (a(idx,i-1)-d(idx,i-1))*(Li(i+j+H(I,J),l(idx))*
							  Vi(11-i-j,l(idx),theta(idx-1),\
							  theta(idx))*Q(I+4*(j-1),J));
				}
			}
			bm(I,J) = rho/(12*g*pow(l(idx),10))*bm_sum;
		}
	}
	return bm;
}
Matrix4d CM(int idx){
	Matrix4d cm;
	cm =  MatrixXd::Zero(4,4);
	double cm_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			cm_sum=0;
			for(int i=1; i<=5; i++ ) {
				for(int j=1; j<=5; j++ ) {
					cm_sum += a(idx,i-1)*Li(i+j+H(I,J),l(idx))*
							  Ui(11-i-j,l(idx))*Q(I+4*(j-1),J)+
							  (d(idx,i-1)-a(idx,i-1))*(Li(i+j+H(I,J),l(idx))*
							  Vi(11-i-j,l(idx),theta(idx-1),\
							  theta(idx))*Q(I+4*(j-1),J));
				}
			}
			cm(I,J) = rho/(12*g*pow(l(idx),10))*cm_sum;
		}
	}
	return cm;
}
Matrix4d DM(int idx){
	Matrix4d dm;
	dm =  MatrixXd::Zero(4,4);
	double dm_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			dm_sum=0;
			for(int i=1; i<=5; i++ ) {
				for(int j=1; j<=5; j++ ) {
					dm_sum += (a(idx,i-1)-d(idx,i-1))*
							  (Li(i+j+H(I,J),l(idx))*Si(11-i-j,\
							  l(idx),theta(idx-1),theta(idx))*
							  Q(I+4*(j-1),J));
				}
			}
			dm(I,J) = rho/(12*g*pow(l(idx),10))*dm_sum;
		}
	}
	return dm;
}
// rotation induced stiffness
Matrix4d EK(int idx){
	Matrix4d ek;
	ek =  MatrixXd::Zero(4,4);
	double ek_sum1, ek_sum2, ek_sum3;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			ek_sum1=0;
			ek_sum2=0;
			ek_sum3=0;
			for(int i=1; i<=3; i++ ) {
				for(int j=1; j<=5; j++ ) {
					ek_sum1 += c(idx,i-1)*Li(i+j+H(I,J),l(idx))*
							   Ui(9-i-j,l(idx))*Q(I+4*(j-1),J);
					ek_sum2 += c(idx,i-1)*Li(i+j+H(I,J),l(idx))*
							   Ui(10-i-j,l(idx))*Q(I+4*(j-1),J);
					ek_sum3 += c(idx,i-1)*Li(i+j+H(I,J),l(idx))*
							   Ui(11-i-j,l(idx))*Q(I+4*(j-1),J);
				}
			}
			ek(I,J) = (rho*pow(Omega,2))/(g*pow(l(idx),8))*(
				(offset*L+0.5*L*L-offset*z(idx)-pow(z(idx),2))*ek_sum1+
				-(offset+z(idx))*ek_sum2-0.5*ek_sum3);
		}
	}
	return ek;
}
Matrix4d FK(int idx){
	Matrix4d fk;
	fk =  MatrixXd::Zero(4,4);
	double fk_sum;
	for(int  I=0; I<4; I++ ) {
		for(int J=0; J<4; J++ ) {
			fk_sum=0;
			for(int i=1; i<=3; i++ ) {
				for(int j=1; j<=7; j++ ) {
					fk_sum += c(idx,i-1)*Li(i+j+H(I,J),l(idx))*
							  Ui(11-i-j,l(idx))*P(I+4*(j-1),J);
				}
			}
			fk(I,J) = (2*rho*Omega*Omega)/(g*pow(l(idx),8))*fk_sum;
		}
	}
	return fk;
}

/* Eigen helper functions*/

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

int getParams() {
	int i;
	ifstream infile("inputs.txt");
	if( infile ) {
		while (infile >> i);
	}
	infile.close();
	return i;
}
