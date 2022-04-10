#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

// double ***TENSOR(int nrl, int nrh, int ncl, int nch, int ndl, int ndh){
//   double ***m;
//   //*** ALLOCATE ROW POINTER ***//
//   m = (double ***) malloc ((nrh-nrl+1) * sizeof(double));
//   for (int i = nrl, i <= nrh, i++){
//     m[i] = (double **) malloc ((nch-ncl+1) * sizeof(double));
//     for (int j = ncl, j <= nch, j++){
//       m[i][j] = (double *) malloc ((ndh-ndl+1) * sizeof(double));
//     }
//   }
//   return m;
// }

//void FREE_RTENSOR(int ***m){
//  delete m;
//}
void startuvwpphi(double ***U, double ***V, double ***W, double ***F, double ***G, double ***H, double ***P, double ***RHS, double ***PHI, double ***NEW_PHI, double nx, double ny, double nz){
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	U[i][j][k] = 0;
	V[i][j][k] = 0;
	W[i][j][k] = 0;
	F[i][j][k] = 0;
	G[i][j][k] = 0;
	H[i][j][k] = 0;
	P[i][j][k] = 0;
	RHS[i][j][k] = 0;
	PHI[i][j][k] = 0;
	NEW_PHI[i][j][k] = 0;
      }
    }
  }
  cout << "UVP is start to run now\n";
}
void initialize(string filename, double ***U, double ***V, double ***W, double ***P, double ***PHI, int nx, int ny, int nz, int resolution){
  //OPEN THE DAT FILE
  //=============================================//
  // ===== Set the initial pressure ============ //
  //=============================================//
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	P[i][j][k] = 1.0;
      }
    }
  }
  cout << "We can initialize the initial pressure now" << endl;
  //=============================================//
  //======= Set the crossflow velocity ==========//
  //=============================================//
  for(int j = 1; j <= ny; j++){
    for(int k = 1; k <= nz; k++){
      U[0][j][k] = 0.5;
    }
  }
  cout << "Jet can be run now " << endl;
  //=============================================//
  // ======= Set the jet velocity ============== //
  //=============================================//
  for(int i = 4*resolution+1; i <= 5*resolution; i++){
    for(int j = 4*resolution+1; j <= 5*resolution; j++){
      W[i][j][0] = 1.0;
      PHI[i][j][0] = 1.0;
    }
  }
  cout << "Initialize is okay now\n";
}
void continuee(string filename, double ***U, double ***V, double ***W, double ***P, double ***PHI, int nx, int ny, int nz){
    //OPEN THE DAT FILE
  string word1, word2, word3, word4;
  ifstream myfileI;
  double t;
  double resolution;
  double Re;
  double dt;
  double ReadU, ReadV, ReadW, ReadP, ReadPHI;
  double period;
  myfileI.open(filename);
  myfileI >> word1 >> t >> word2 >> resolution >> word3 >> Re >> word4 >> dt >> period;
  for(int i = 0 ; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; j <= nz+1; k++){
	myfileI >> ReadU >> ReadV >> ReadW >> ReadP >> ReadPHI;
	U[i][j][k] = ReadU;
	V[i][j][k] = ReadV;
	W[i][j][k] = ReadW;
	P[i][j][k] = ReadP;
	PHI[i][j][k] = ReadPHI;
      }
    }
  }
  myfileI.close();  
  cout << "Continue is also okay!\n";
}
void comp_delta(double dt, int nx, int ny, int nz, double dx, double dy, double dz, double ***U, double ***V, double ***W, double Re, double tau){
  double umax, vmax, wmax;
  umax = U[0][0][0];
  vmax = V[0][0][0];
  wmax = W[0][0][0];
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	if (U[i][j][k] > umax){U[i][j][k] = umax;}
	if (V[i][j][k] > vmax){V[i][j][k] = vmax;}
	if (W[i][j][k] > wmax){W[i][j][k] = wmax;}
      }
    }
  }
  double DT[4];
  DT[0] = Re/(2.0*(1/(dx*dx)+1/(dy*dy)+1/(dz*dz)));
  DT[1] = dx/abs(umax);
  DT[2] = dy/abs(vmax);
  DT[3] = dx/abs(wmax);
  dt = DT[0];
  for(int i = 0; i < 4 ; i++){
    if( dt < DT[i]){
      dt = DT[i];
    }
  }
  if (dt > 0.1){dt = 0.1;}
  cout << "We can compile the value of delta" ;
}
void setbcond(double ***U, double ***V, double ***W, double  ***P, double ***Phi, int nx, int ny, int nz, int resolution){
  //==========================================================================================//
  //------------------NEUMANN BOUNDARY CONDITION----------------------------------------------//
  //==========================================================================================//
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      //=============AT THE TOP BOUNDARY CONDITION==============//
      U[i][j][nz+1]   = U[i][j][nz];
      V[i][j][nz+1]   = V[i][j][nz];
      W[i][j][nz+1]   = W[i][j][nz];
      P[i][j][nz+1]   = P[i][j][nz];
      Phi[i][j][nz+1] = Phi[i][j][nz];
    }
  }
  cout << "Neumann is possible to use" << endl;
  //=============AT x = 0 and x = nx+1 BOUNDARY CONDITION============//
  for(int j = 1; j <= ny; j++){
    for(int k = 1; k <= nz; k++){
      //=============AT x = 0 BOUNDARY CONDITION==============//
      U[0][j][k]   = U[1][j][k];
      V[0][j][k]   = V[1][j][k];
      W[0][j][k]   = W[1][j][k];
      P[0][j][k]   = P[1][j][k];
      Phi[0][j][k] = Phi[1][j][k];
      
      //=============AT x = n+1 BOUNDARY CONDITION==============//
      U[nx+1][j][k]   = U[nx][j][k];
      V[nx+1][j][k]   = V[nx][j][k];
      W[nx+1][j][k]   = W[nx][j][k];
      P[nx+1][j][k]   = P[nx][j][k];
      Phi[nx+1][j][k] = Phi[nx][j][k];
    }
  }
  cout << "x=0 and x= nx+1 is fine\n" ;
  //=============AT y = 0 and y = ny+1 BOUNDARY CONDITION============//
  for(int i = 1; i <= nx; i++){
    for(int k = 1; k <= nz; k++){
      //=============AT x = 0 BOUNDARY CONDITION==============//
      U[i][0][k]   = U[i][1][k];
      V[i][0][k]   = V[i][1][k];
      W[i][0][k]   = W[i][1][k];
      P[i][0][k]   = P[i][1][k];
      Phi[i][0][k] = Phi[i][1][k];
      
      //=============AT x = n+1 BOUNDARY CONDITION==============//
      U[i][ny+1][k]   = U[i][ny][k];
      V[i][ny+1][k]   = V[i][ny][k];
      W[i][ny+1][k]   = W[i][ny][k];
      P[i][ny+1][k]   = P[i][ny][k];
      Phi[i][ny+1][k] = Phi[i][ny][k];
    }
  }
  cout << "y = 0 abd y = ny+1 is also good\n";
  //==================================================================//
  //==DIRICHLET INFLOW BOUNDARY CONDITION AND NO SLIP-WALL CONDITION==//
  //==================================================================//
  //====NOSLIP====//
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      U[i][j][0] = 0;
      V[i][j][0] = 0;
      W[i][j][0] = 0;
      P[i][j][0] = 0;
      Phi[i][j][0] = 0; 
    }
  }
  cout << "Inflow is also good right now\n";
  //=======================================//
  //==DIRICHLET INFLOW BOUNDARY CONDITION==//
  //=======================================//
  int block = resolution;
  for(int i = 4*block+1; i <= 5*block; i++){
    for(int j = 4*block+1; j <= 5*block; j++){
      Phi[i][j][0] = 1.;
      W[i][j][0] = 1.;
    }
  }
  cout << "another inflow is also good\n";
  //================================================//
  //======== CROSSFLOW VELOCITY ====================//
  //================================================//
  for(int j = 1; j <= ny; j++){
    for(int k = 1; k <= nz; k++){
      U[0][j][k] = 0.5;
    }
  }
  cout << "BCON function is useable" << endl;
}
void comp_FGH(double ***U, double ***V, double ***W, double ***F, double ***G, double ***H, double dt, double dx, double dy, double dz, int nx, int ny, int nz, double Gx, double Gy, double Gz,
	      double gamma, double Re){
  double d2udx2, d2udy2, d2udz2, d2vdx2, d2vdy2, d2vdz2, d2wdx2, d2wdy2, d2wdz2;
  double du2dx, duvdx, duwdx, dvudy, dv2dy, dvwdy, dwudz, dwvdz, dw2dz;
  double du2dxp1, du2dxp2, duvdxp1, duvdxp2, duwdxp1, duwdxp2;
  double dvudyp1, dvudyp2, dv2dyp1, dv2dyp2, dvwdyp1, dvwdyp2;
  double dwudzp1, dwudzp2, dwvdzp1, dwvdzp2, dw2dzp1, dw2dzp2;
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      for(int k = 1; k <= nz; k++){
	//============================================//
	//--------------VISCOUS TERM------------------//
	//============================================//
	d2udx2 = (U[i+1][j][k]-2*U[i][j][k]+U[i-1][j][k])/(dx*dx);
	d2udy2 = (U[i][j+1][k]-2*U[i][j][k]+U[i][j-1][k])/(dy*dy);
	d2udz2 = (U[i][j][k+1]-2*U[i][j][k]+U[i][j][k-1])/(dz*dz);
	d2vdx2 = (V[i+1][j][k]-2*V[i][j][k]+V[i-1][j][k])/(dx*dx);
	d2vdy2 = (V[i][j+1][k]-2*V[i][j][k]+V[i][j-1][j])/(dy*dy);
	d2vdz2 = (V[i][j][k+1]-2*V[i][j][k]+V[i][j][k-1])/(dz*dz);
	d2wdx2 = (W[i+1][j][k]-2*W[i][j][k]+W[i-1][j][k])/(dx*dx);
	d2wdy2 = (W[i][j+1][k]-2*W[i][j][k]+W[i][j+1][k])/(dy*dy);
	d2wdz2 = (W[i][j][k+1]-2*W[i][j][k]+W[i][j][k-1])/(dz*dz);
	//========================================================//
	//----------------CONVECTION TERM-------------------------//
	//========================================================//
	//***                 PART1                            ***//
	du2dxp1 = 1/(4*dx)*((U[i+1][j][k]+U[i][j][k])*(U[i+1][j][k]+U[i][j][k])-(U[i-1][j][k]+U[i][j][k])*(U[i-1][j][k]+U[i][j][k]));
	duvdxp1 = 1/(4*dx)*((U[i][j+1][k]+U[i][j][k])*(V[i+1][j][k]+V[i][j][k])-(U[i-1][j][k]+U[i-1][j+1][k])*(V[i-1][j][k]+V[i][j][k]));
	duwdxp1 = 1/(4*dx)*((U[i][j][k+1]+U[i][j][k])*(W[i+1][j][k]+W[i][j][k])-(U[i-1][j][k]+U[i-1][j][k+1])*(W[i-1][j][k]+W[i][j][k]));
        dvudyp1 = 1/(4*dy)*((V[i+1][j][k]+V[i][j][k])*(U[i][j+1][k]+U[i][j][k])-(V[i][j-1][k]+V[i+1][j-1][k])*(U[i][j-1][k]+U[i][j][k]));
	dv2dyp1 = 1/(4*dy)*((V[i][j+1][k]+V[i][j][k])*(V[i][j+1][k]+V[i][j][k])-(V[i][j-1][k]+V[i][j][k])*(V[i][j-1][k]+V[i][j][k]));
	dvwdyp1 = 1/(4*dy)*((V[i][j][k+1]+V[i][j][k])*(W[i][j+1][k]+W[i][j][k])-(V[i][j-1][k]+V[i][j-1][k+1])*(W[i][j-1][k]+W[i][j][k]));
	dwudzp1 = 1/(4*dz)*((W[i+1][j][k]+W[i][j][k])*(U[i][j][k+1]+U[i][j][k])-(W[i][j][k-1]+W[i+1][j][k-1])*(U[i][j][k-1]+U[i][j][k]));
	dwvdzp1 = 1/(4*dz)*((W[i][j+1][k]+W[i][j][k])*(V[i][j][k+1]+V[i][j][k])-(W[i][j][k-1]+W[i][j+1][k-1])*(V[i][j][k-1]+V[i][j][k]));
	dw2dzp1 = 1/(4*dz)*((W[i][j][k+1]+W[i][j][k])*(W[i][j][k+1]+W[i][j][k])-(W[i][j][k-1]+W[i][j][k])*(W[i][j][k-1]+W[i][j][k]));
	//*** PART2 ****/
	du2dxp2 = gamma/(4*dx)*(abs(U[i+1][j][k]+U[i][j][k])*(U[i+1][j][k]-U[i][j][k])-abs(U[i-1][j][k]+U[i][j][k])*(U[i-1][j][k]-U[i][j][k]));
	duvdxp2 = gamma/(4*dx)*(abs(U[i][j+1][k]+U[i][j][k])*(V[i+1][j][k]-V[i][j][k])-abs(U[i-1][j][k]+U[i-1][j+1][k])*(V[i-1][j][k]-V[i][j][k]));
	duwdxp2 = gamma/(4*dx)*(abs(U[i][j][k+1]+U[i][j][k])*(W[i+1][j][k]-W[i][j][k])-abs(U[i-1][j][k]+U[i-1][j][k+1])*(W[i-1][j][k]-W[i][j][k]));
        dvudyp2 = gamma/(4*dy)*(abs(V[i+1][j][k]+V[i][j][k])*(U[i][j+1][k]-U[i][j][k])-abs(V[i][j-1][k]+V[i+1][j-1][k])*(U[i][j-1][k]-U[i][j][k]));
	dv2dyp2 = gamma/(4*dy)*(abs(V[i][j+1][k]+V[i][j][k])*(V[i][j+1][k]-V[i][j][k])-abs(V[i][j-1][k]+V[i][j][k])*(V[i][j-1][k]-V[i][j][k]));
	dvwdyp2 = gamma/(4*dy)*(abs(V[i][j][k+1]+V[i][j][k])*(W[i][j+1][k]-W[i][j][k])-abs(V[i][j-1][k]+V[i][j-1][k+1])*(W[i][j-1][k]-W[i][j][k]));
	dwudzp2 = gamma/(4*dz)*(abs(W[i+1][j][k]+W[i][j][k])*(U[i][j][k+1]-U[i][j][k])-abs(W[i][j][k-1]+W[i+1][j][k-1])*(U[i][j][k-1]-U[i][j][k]));
	dwvdzp2 = gamma/(4*dz)*(abs(W[i][j+1][k]+W[i][j][k])*(V[i][j][k+1]-V[i][j][k])-abs(W[i][j][k-1]+W[i][j+1][k-1])*(V[i][j][k-1]-V[i][j][k]));
	dw2dzp2 = gamma/(4*dz)*(abs(W[i][j][k+1]+W[i][j][k])*(W[i][j][k+1]-W[i][j][k])-abs(W[i][j][k-1]+W[i][j][k])*(W[i][j][k-1]-W[i][j][k]));
        //============SUMMATION OF CONVECTION TERM=============//
	du2dx = du2dxp1+du2dxp2;
	duvdx = duvdxp1+duvdxp2;
	duwdx = duwdxp1+duwdxp2;
	dvudy = dvudyp1+dvudyp2;
	dv2dy = dv2dyp1+dv2dyp2;
	dvwdy = dvwdyp1+dvwdyp2;
	dwudz = dwudzp1+dwudzp2;
	dwvdz = dwvdzp1+dwvdzp2;
	dw2dz = dw2dzp1+dw2dzp2;
	//=============================================//
	//=============COMPUTE F G H-------------------//
	//=============================================//
	F[i][j][k] = U[i][j][k]+dt*((1/Re)*(d2udx2+d2udy2+d2udz2)-(du2dx+dvudy+dwudz)+Gx);
	G[i][j][k] = V[i][j][k]+dt*((1/Re)*(d2vdx2+d2vdy2+d2vdz2)-(duvdx+dv2dy+dwvdz)+Gy);
	H[i][j][k] = W[i][j][k]+dt*((1/Re)*(d2wdx2+d2wdy2+d2wdz2)-(duwdx+dvwdy+dw2dz)+Gz);
      }
    }
  }
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      F[i][j][0] = F[i][j][1];
      G[i][j][0] = G[i][j][1];
      H[i][j][0] = H[i][j][1];
      F[i][j][nz+1] = F[i][j][nz];
      G[i][j][nz+1] = G[i][j][nz];
      H[i][j][nz+1] = H[i][j][nz];
    }
  }
  for(int i = 1; i <= nx; i++){
    for(int k = 1; k <= nz; k++){
      F[i][0][k] = F[i][1][k];
      G[i][0][k] = G[i][1][k];
      H[i][0][k] = H[i][1][k];
      F[i][ny+1][k] = F[i][ny][k];
      G[i][ny+1][k] = G[i][ny][k];
      H[i][ny+1][k] = H[i][ny][k];
    }
  }
  for(int j = 1; j <= ny; j++){
    for(int k = 1; k <= nz; k++){
      F[0][j][k] = F[1][j][k];
      G[0][j][k] = G[1][j][k];
      H[0][j][k] = H[1][j][k];
      F[nx+1][j][k] = F[nx][j][k];
      G[nx+1][j][k] = G[nx][j][k];
      H[nx+1][j][k] = H[nx][j][k];
    }
  }
  cout << "Comp FGH is usable" << endl;
}
void comp_RHS(double ***F, double ***G, double ***H, double ***RHS, int nx, int ny, int nz, double dx, double dy, double dz, double dt){
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      for(int k = 1; k <= nz; k++){
	RHS[i][j][k] = (1/dt)*((F[i+1][j][k]-F[i][j][k])/(dx)+(G[i][j+1][k]-G[i][j][k])/(dy)+(H[i][j][k+1]-H[i][j][k])/(dz));
      }
    }
  }
  cout << "CompRHS is able to use" << endl; 

}
int poisson(double ***P, double ***RHS, int nx, int ny, int nz,double dx, double dy, double dz, int itmax, double omg, double res){
    for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      P[i][j][0] = P[i][j][1];
      P[i][j][nz+1] = P[i][j][nz];
    }
  }
  for(int i = 1; i <= nx; i++){
    for(int k = 1; k <= nz; k++){
      P[i][0][k] = P[i][1][k];
      P[i][ny+1][k] = P[i][ny][k];
    }
  }
  for(int j = 1; j <= ny; j++){
    for(int k = 1; k <= nz; k++){
      P[0][j][k] = P[1][j][k];
      P[nx+1][j][k] = P[nx][j][k];
    }
  }
  int it = 0;
  while (true) {
    int N = 0;
    double Pstar;
    for(int i = 1; i <= nx; i++){
      for(int j = 1; j <= ny; j++){
	for(int k = 1; k <= nz; k++){
	  double Pprevious = P[i][j][k];
	  Pstar = (1/6.0)*(P[i+1][j][k]+P[i-1][j][k]+P[i][j+1][j]+P[i][j-1][k]+P[i][j][k+1]+P[i][j][k])+(dx*dx/6.0)*RHS[i][j][k];
	  P[i][j][k] = (1-omg)*P[i][j][k]+omg*Pstar;
	  if(abs(Pprevious-Pstar) <= pow(10, -8)){
	    N+=1;
	    cout << N <<" "<< it<< endl;
	  }
	}
      }
    }
    it+=1;
    if (it == 10){break;}
    if (N == nx*ny*nz){break;}
  }
  cout << "The POISSION is useable" << endl;
  return 0;
}
void adap_UVW(double ***U, double ***V, double ***W, double ***F, double ***G, double ***H, double ***P, int nx, int ny, int nz, double dt, double dx, double dy, double dz){
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      for(int k = 1; k <= nz; k++){
	U[i][j][k] = F[i][j][k] - (dt/dx)*(P[i+1][j][k]-P[i][j][k]);
	V[i][j][k] = G[i][j][k] - (dt/dy)*(P[i][j+1][k]-P[i][j][k]);
	W[i][j][k] = H[i][j][k] - (dt/dz)*(P[i][j][k+1]-P[i][j][k]);
      }
    }
  }
  cout << "Updated UVW!!!!!!!" << endl;
}
void simulation_phi(double ***PHI, double ***NEW_PHI, double ***U, double ***V, double ***W, double dt, double dx, double dy, double dz, int nx, int ny, int nz, double Re){
  double RHS1, RHS2, RHS;
  //simulation for phi
  for(int i = 1; i <= nx; i++){
    for(int j = 1; j <= ny; j++){
      for(int k = 1; k <= nz; k++){
	RHS1 = (1/(Re*dx*dx))*(PHI[i+1][j][k]+PHI[i-1][j][k]+PHI[i][j+1][k]+PHI[i][j-1][k]+PHI[i][j][k+1]+PHI[i][j][k-1]-6*PHI[i][j][k]);
	RHS2 = (1/dx)*(U[i][j][k]*(PHI[i+1][j][k]-PHI[i][j][k])+V[i][j][k]*(PHI[i][j+1][k]-PHI[i][j][k])+W[i][j][k]*(PHI[i][j][k+1]-PHI[i][j][k]));
	RHS = RHS1-RHS2;
	NEW_PHI[i][j][k] = PHI[i][j][k]+dt*RHS;
      }
    }
  }
  cout << "simulation of phi is successed !!" << endl;
}
void update(double ***PHI, double ***NEW_PHI, int nx, int ny, int nz){
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	PHI[i][j][k] = NEW_PHI[i][j][k]; 
      }
    }
  }
  cout << "Updated PHI !! " << endl;
}
void paraview(string fileName, double ***U, double ***V, double ***W, double ***P, double ***PHI, double dx, double dy, double dz, double nx, double ny, double nz){
  ofstream myfile;
  myfile.open(fileName);
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";
  //Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS" <<  " " << nx+2 << " " << nz+2 << " " << ny+2 << endl;
  myfile << "POINTS " << int((nx+2)*(ny+2)*(nz+2)) << " float\n";
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	myfile << dx*i << " " << dy*j << " " << dz*k << endl;
	cout << dx*i << " " << dy*j << " " << dz*k << endl;
      }
    }
  }
  myfile << "\n";
  myfile << "POINT_DATA" <<  " " << int((nx+2)*(ny+2)*(nz+2)) << endl;
  myfile << "\n";
  myfile << "SCALARS U float 1" << endl;
  myfile << "LOOKUP_TABLE default\n";
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	myfile << U[i][j][k] <<  endl;
      }
    }
  }
  myfile << "\n";
  myfile << "SCALARS V float 1" << endl;
  myfile << "LOOKUP_TABLE default\n";
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	myfile << V[i][j][k] <<  endl; 
      }
    }
  }
  myfile << "\n";
  myfile << "SCALARS W float 1" << endl;
  myfile << "LOOKUP_TABLE default\n";
  
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	myfile << W[i][j][k] <<  endl; 
      }
    }
  }
  myfile << "\n";
  myfile << "SCALARS P float 1" << endl;
  myfile << "LOOKUP_TABLE default\n";
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	myfile << P[i][j][k] <<  endl; 
      }
    }
  }
  myfile << "\n";
  myfile << "SCALARS PHI float 1" << endl;
  myfile << "LOOKUP_TABLE default\n";
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	myfile << PHI[i][j][k] <<  endl; 
      }
    }
  }
  myfile.close();
  cout << "Paraview file is written " <<  endl;
}
void save_restartfile(double ***U, double ***V, double ***W, double ***P, double ***PHI, double t, double Re, double dt, int Resolution, int nx, int ny, int nz, int period){
  ofstream MyfileO;
  MyfileO.open("backup.dat");
  MyfileO << "t=" << " " << t << " " << "Resolution=" << " " << Resolution << " " << "Re=" << " " << Re << " " << "dt=" << " " << dt << " " << period <<endl;
  for(int i = 0; i <= nx+1; i++){
    for(int j = 0; j <= ny+1; j++){
      for(int k = 0; k <= nz+1; k++){
	MyfileO << U[i][j][k] << " " << V[i][j][k] << " " << W[i][j][k] << " " << P[i][j][k] << " " << PHI[i][j][k] << endl;
      }
    }
  }
  MyfileO.close();
  cout << "Backup file is written " << endl;
}
int main(){
  int period = 0;
  int NUMBER;
  cout << "# Would you like to continue the simulation ? "     << endl;
  cout << "Press 1 then enter if you want to start the new simulation"  << endl;
  cout << "Press 2 then enter if you want to continue the simulation" << endl;
  cout << "Press the number: ";
  cin >> NUMBER ;
  if(NUMBER != 1 && NUMBER != 2){cerr << "Illegal value of the number!\n";exit(1);}
  int resolution;
  double t;
  double Re;
  double dt;
  if(NUMBER == 1){
    t = 0;
    cout << "# Enter the resolution R of the grid as the integer between 1-20 \n";
    cout << "# The domain of this simulation is x:y:z = 20:10:15 respectively \n";
    cout << "# The Number of grid in x = 20*R, y = 10*R and z = 15*R. The Total number of grid is 3000*N*N*N grids in the simulation \n"; 
    cout << "# For example, if you enter 20 in the input the total number of grids is 3000*8000 = 24000000 grids respectively \n";
    cout << "Enter the resolution R (R is an integer between 1-20";
    cin >> resolution;
    if(resolution <= 0 && resolution >= 21){cerr << "Rewrite the new resolution\n"; exit(1);}
    cout << "\n";
    cout << "The resoulution is " << resolution << endl;
    cout << "The total number of grid is = " << 3000*resolution*resolution*resolution << endl;
    cout << "#Enter the Reynold number R (1000 - 10000)";
    cin >> Re;
    cout << "\n";
    if(Re <= 1000 && Re >= 10000){cerr << "Inappropiate value of Reynold number Re\n"; exit(1);}
    cout << "\n";
    cout << "#Enter the Reynold number dt (0.00001 - 0.2)";
    cin>> dt;
    cout << "\n";
    if(dt <= 0.00001 && Re >= 0.2){cerr << "Inappropiate value of Reynold number Re\n"; exit(1);}
  }
  if(NUMBER == 2){
    ifstream readfile1;
    string word1, word2, word3, word4;
    int period;
    readfile1.open("backup.dat");
    readfile1 >> word1 >> t >> word2 >> resolution >> word3 >> Re >> word4 >> dt >> period;
  }
  int nx = 20*resolution;
  int ny = 10*resolution;
  int nz = 15*resolution;
  //======================STARTING TO DECLARE POINTER===============//
  double ***U;
  U = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    U[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      U[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***V;
  V = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    V[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      V[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***W;
  W = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    W[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      W[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***F;
  F = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    F[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      F[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***G;
  G = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    G[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      G[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***H;
  H = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    H[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      H[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***P;
  P = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    P[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      P[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***PHI;
  PHI = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    PHI[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      PHI[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***NEW_PHI;
  NEW_PHI = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    NEW_PHI[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      NEW_PHI[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  double ***RHS;
  RHS = (double ***) malloc ((nx+2) * sizeof(double));
  for(int i = 0; i < nx+2; i++){
    RHS[i] = (double **) malloc((ny+2) * sizeof(double));
    for(int j = 0; j < ny+2; j++){
      RHS[i][j] = (double *) malloc ((nz+2) * sizeof(double));
    }
  }
  startuvwpphi(U, V, W, F, G, H, P, RHS, PHI, NEW_PHI, nx, ny, nz);
  //=========================================================//
  //================END OF DECLARING POINTER=================//
  //=========================================================//
  //===============PROCEED THE SIMULATION====================//
  //=========================================================//
  double Gx = 0;
  double Gy = 0;
  double Gz = 0;
  double gamma = 1.0;
  double itmax = 1000;
  double omg = 1.7;
  double res = 1;
  double dx, dy, dz;
  dx = 1.0 ; dy = 1.0 ; dz =1.0;
  //==========SELECT THE FILE FOR SIMULATION=================//
  if(NUMBER == 1){
    int period = 0;
    initialize("DATA_0.vtk", U, V, W, P, PHI, nx, ny, nz, resolution);
    
    paraview("DATA_0.vtk", U, V,W, P, PHI, dx,dy, dz,nx,ny,nz);
  }
  if(NUMBER == 2){
    continuee("backup.dat", U, V, W, P, PHI, nx, ny, nz);
  }
  int PERIOD = 0;
  while (t < 50){
    setbcond(U, V, W, P, PHI,nx,ny, nz, resolution);
    comp_FGH(U,V,W,F,G,H,dt,dx,dy,dz,nx,ny,nz,Gx,Gy,Gz,gamma,Re);
    comp_RHS(F, G, H, RHS, nx, ny, nz, dx, dy, dz, dt);
    poisson(P, RHS, nx,ny,nz,dx, dy, dz, itmax, omg, res);
    adap_UVW(U, V,W, F, G, H, P, nx, ny, nz,dt, dx, dy, dz);
    simulation_phi(PHI,NEW_PHI,U,V,W, dt, dx, dy, dz, nx, ny, nz, Re);
    update(PHI,NEW_PHI,nx,ny,nz);
    //Write new paraview file
    string fileName;
    period +=1;
    fileName = "DATA_"+to_string(period)+".vtk";
    paraview(fileName, U, V,W, P, PHI, dx,dy, dz,nx,ny,nz);
    PERIOD += 1;
    cout << "The data number " << period << " is saved";
    if( PERIOD == 20 ){
      save_restartfile(U, V, W,  P, PHI, t, Re, dt, resolution, nx, ny, nz, period);
      cout << "You have already save the restartfile" << endl;
      PERIOD = 0;
    }
    t += dt;
    if( t == 50){cout << "Congratulation!!! You finished you final project absolutely!!" << endl; break;}
  }
  return 0;  
}
