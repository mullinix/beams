#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include <time.h>
#include <math.h>
#define intPts 10 //number of legendre polynomial points
#define PI 3.141592653589793238462643
#define EPS 1e-14 //integration tolerance for gauss-legendre coeff calculation

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

/* Read parameters text file */
int getParams();

/* timing for efficiency */
clock_t init, dt, tfinal;

/* System parameters */
// number of elements
int N = 100;
// length of beam [m] (0.45)
double L = 0.1524; 
// accel due to gravity ()
double g = 9.81; 
// density of steel [kg/m^3] (steel: 7850, Alum: 2700)
double rho = 800;
// E: Young's modulus of steel [N/m^2] (steel: 21e10, alum: 7.2e10)
double E = 2.07e11; 
// Shear modulus of Pascals (steel: 19.3e9)
double G = E/2.6;
// "shear coefficient" for rectangular cross section (approximate) http://en.wikipedia.org/wiki/Timoshenko_beam_theory#Shear_coefficient
double mu = 0.833;// 10*(1+nu)/(12+11*nu); 
// w: beam width [m] (0.0254) 
double w = 0.0254;
// h: beam height [m] (0.0046)
double h = 0.0046;
// total twist angle from base to tip (45.0)
double twist_angle = 45.0;
// taper of width (w) (2.56)
double breadth_taper = 2.56;
// taper of height (h) (2.29)
double depth_taper = 2.29;
// rotation speed, rps - rpm*(secs/min) (100)
double Omega = 100; 
// distance from center of rotation to base of beam
double offset = 0.0; 

int gotParams = getParams();
/* other, unused parameters */
/* 
// m: total beam mass [kg]
// m = L*w*h*rho; 
// nu: Poisson's ratio (steel: 0.3, alum: 0.35)
// nu = 3/10; 
*/

// nodal angular displacement wrt origin
VectorXd theta = VectorXd::LinSpaced(N,0.0,twist_angle*PI/180.0); //linear twist
VectorXd z = VectorXd::LinSpaced(N,0.0,L); // physical locations of nodes from origin
// allow for varying cross section
VectorXd breadth = VectorXd::LinSpaced(N,w,w/breadth_taper);// distribution - use for "b" in RG elemnts
VectorXd depth = VectorXd::LinSpaced(N,h,h/depth_taper);// distribution - use for "h" in RG elements

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

int main() {
	FILE *k_dat;
	FILE *m_dat;
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
	Matrix4d ak,bk,ck,dk,am,bm,cm,dm,ek,fk,zz;
	zz = MatrixXd::Zero(4,4);
	MatrixXd Kel,Mel,Ktemp,Mtemp;
	Kel = MatrixXd::Zero(16,16);
	Mel = MatrixXd::Zero(16,16);
	Ktemp = MatrixXd::Zero(16,16);
	Mtemp = MatrixXd::Zero(16,16);
	//~ cout << Kel << endl;
	MatrixXd M,K;
	M = MatrixXd::Zero(8*(N+1),8*(N+1));
	K = MatrixXd::Zero(8*(N+1),8*(N+1));
	for(int elt=1;elt<N;elt++) {
		ak = AK(elt);
		bk = BK(elt);
		ck = CK(elt);
		dk = DK(elt);
		am = AM(elt);
		bm = BM(elt);
		cm = CM(elt);
		dm = DM(elt);
		ek = EK(elt);
		fk = FK(elt);
		
	    Kel << ak+ek-fk,      ek-fk,       dk,  	  zz,
	              ek-fk,   ck+ek-fk, 	   zz,  	  zz,
	                 dk,   		 zz, bk+ek-fk,     ek-fk,
					 zz,   		 zz,    ek-fk,  ck+ek-fk;

	    Mel <<    am+bm,         am,       dm,  	  zz,
	                 am,         am,       am,        zz,
	                 dm,         am,    am+cm,        am,
	                 zz,         zz,       am,        am;
		ofstream el_mx;
		el_mx.open("el_mx.dat");
		el_mx << Mel << endl;
		el_mx.close();
		for(int el_idx=0;el_idx<8;el_idx++) {
			Ktemp.row(el_idx) = Kel.row(2*el_idx);
			Ktemp.row(8+el_idx) = Kel.row(2*el_idx+1);
			Mtemp.row(el_idx) = Mel.row(2*el_idx);
			Mtemp.row(8+el_idx) = Mel.row(2*el_idx+1);
		}
		for(int el_idx=0;el_idx<8;el_idx++) {
			Kel.col(el_idx) = Ktemp.col(2*el_idx);
			Kel.col(8+el_idx) = Ktemp.col(2*el_idx+1);
			Mel.col(el_idx) = Mtemp.col(2*el_idx);
			Mel.col(8+el_idx) = Mtemp.col(2*el_idx+1);
		}
		if(elt==1){
			K.block(0,0,16,16) = Kel;
			M.block(0,0,16,16) = Mel;
		}
		K.block(8*elt,8*elt,16,16) += Kel;
		M.block(8*elt,8*elt,16,16) += Mel;
		  /* uncomment for testing: */
//	NOTE: Adjust AK, BK prefactor to divide by l^11, not l^10!
//  NOTE: Adjust CK prefactor to divide by l^9, not l^8!
	     double K_factor = E*Ixx(elt)/pow(l(elt),3);// should be l^3 in denominator!
	     double JJ = mu*G*breadth(elt)*depth(elt)*pow(l(elt),2)/(30*E*Ixx(elt));
		 double M_factor = (rho*breadth(elt)*depth(elt))/(420*g)*pow(l(elt),2);// should be times Ixx, not l^2!
	     double PP = (14*Ixx(elt))/(breadth(elt)*depth(elt)*pow(l(elt),2));
	     cout << "6l: " << 6*l(elt) << endl;
	     cout << "2l^2: " << 2*pow(l(elt),2) << endl;
	     cout << "4l^2: " << 4*pow(l(elt),2) << endl;
	     cout << "3*l*J: " << 3*l(elt)*JJ << endl;
	     cout << "36*J: " << 36*JJ << endl;
	     cout << "l^2*J: " << pow(l(elt),2)*JJ << endl;
	     cout << "4l^2*J: " << 4*pow(l(elt),2)*JJ << endl;
	     cout << "P: " << PP << endl;
	     cout << "K factor: " << K_factor << endl;
	     cout << "M Factor: " << M_factor << endl;
	     
	}

	tfinal=clock()-init;
	cout << "Assembly time: " << (double)tfinal / ((double)CLOCKS_PER_SEC) << endl;
	init = clock();
	k_dat = fopen ("k.dat", "w");
	m_dat = fopen ("m.dat", "w");
	N = K.rows();
	for( int i=0; i<N; i++ ) {
		for( int j=0; j<N; j++ ) {
			fprintf(k_dat, "%.15e ", K(i,j));
			fprintf(m_dat, "%.15e ",  M(i,j));
		}
		fprintf(k_dat,"\n");
		fprintf(m_dat,"\n");
	}
	fclose(k_dat);
	fclose(m_dat);
	tfinal=clock()-init;
	cout << "File write time: " << (double)tfinal / ((double)CLOCKS_PER_SEC) << endl;
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
			ak(I,J) = E/(12*pow(l(idx),11))*ak_sum;
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
			bk(I,J) = E/(12*pow(l(idx),11))*bk_sum;
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
			ck(I,J) = ((mu*G)/pow(l(idx),9))*ck_sum;
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
/* Read custom parameters from file */
int getParams(){	
	string line;
	ifstream infile;
	infile.open("inputs.txt");
	if( infile.is_open() ) {
		while ( getline(infile,line) ){
			stringstream stream(line); //convert line in to stream
			/* conversion temporary vars */
			stringstream name_stream;
			stringstream val_stream;
			/* temp holding vars */
			string name_read;
			string val_read;
			/* actual vals we are interested in */
			double value;
			string name;
			/* get the name and value from the text file */
			getline(stream,name_read,'=');
			getline(stream,val_read,'=');
			/* this converts the string value into a double */
			val_stream << val_read;
			val_stream >> value;
			/* this parses whitespace from the string */
			name_stream << name_read;
			name_stream >> name;
			/* determine which variable to change, change it */
			if(name=="N") N = (int)value;
			else if(name=="L") L = value;
			else if(name=="g") g = value;
			else if(name=="rho") rho = value;
			else if(name=="E") E = value;
			else if(name=="G") G = value;
			else if(name=="mu") mu = value;
			else if(name=="w") w = value;
			else if(name=="h") h = value;
			else if(name=="twist_angle") twist_angle = value;
			else if(name=="breadth_taper") breadth_taper = value;
			else if(name=="depth_taper") depth_taper = value;
			else if(name=="Omega") Omega = value;
			else if(name=="offset") offset = value;
			/* print the name and value, for troubleshooting */
		    //cout << "Name: " << name << ", Value: " << value << endl;
		}
	}
	else {
		cout << "Not reading inputs.txt!" << endl; 
	}
	infile.close();
	return 1;
}
