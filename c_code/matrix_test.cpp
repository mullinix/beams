#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#define n 256

using namespace std;
using namespace Eigen;

//~ typedef Eigen::Matrix<double, n, n> Matrixnd;
//~ typedef Eigen::Matrix<double, n, 1> Vectornd;
int main()
{
  Eigen::setNbThreads(4);
  clock_t init, ttemp, dt, tfinal;
  dt=0;
  MatrixXd A(n,n), ApAT(n,n);
  VectorXd b(n), x(n);
  cout << "Solve Ax=b, size="<< n <<  ":" << endl;
  A = MatrixXd::Random(n,n);
  b = VectorXd::Random(n);
  //~ ApAT.noalias()= A*A.transpose();
  init=clock();
  for(int i=0;i<1000;i++){
	  x = A.partialPivLu().solve(b);
	  ttemp=clock();
	  if(i%100==99){
		cout << "iteration: " << i+1 << endl;
		double error = (A*x - b).norm()/b.norm(); // norm() is L2 norm
		cout << "The relative error is:\n" << error << endl;
	  }
      dt += clock()-ttemp;
  }
  tfinal=clock()-dt-init;
  cout << "Computation time:" << endl;
  cout << (double)tfinal / ((double)CLOCKS_PER_SEC) << endl;
}
