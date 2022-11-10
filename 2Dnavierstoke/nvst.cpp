#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <iomanip>
using namespace std;
double **MATRIX(int x,int y){
  double **m;
  m = (double **) malloc ((x)* sizeof(double));
  for(int i = 0; i < y; i++){
    m[i] = (double *) malloc ((y) * sizeof(double));
  }
  return m;
}
void initialize(double **U, double **V, double **F, double **G, double **P, double **RHS, double **PHI, double **NEW_PHI, int imax, int jmax){
}
void setbcon(double **U, double **V, double **F, double **G, double **P, double **RHS, double **PHI, double **PHI_NEW, int imax, int jmax){
//set u_velocity
  // Inflow Condition
  for(int j = 0; j <= jmax-1; j++){
    U[0][j] = 1.0;
  }
  // Outflow Condition
  for(int j = 0; j <= jmax-1; j++){
    U[imax-2][j] = U[imax-3][j];
    U[imax-1][j] = U[imax-1][j];
  }
  // No-Slip Condition
  for(int i = 0; i <= imax-1; i++){
    U[i][jmax-1] = -1.0*U[i][jmax-2]; //Top
    U[i][0] = -1.0*U[i][1]; // Bottom
  }
  //set v_velocity
  // inlet condition
  for(int j = (jmax-1)/2.0; j <= jmax-1; j++){
    V[0][j] = -1.0*V[1][j];
  }
  // outlet
  for(int j = 0; j <= jmax-1; j++){
    V[imax-1][j] = V[imax-2][j];
  }
  // No-Slip Condition
  for(int i = 1; i <= imax-1; i++){
    V[i][jmax-2] = 0;
    V[i][0] = 0.0;
  }
  //set phi condition
  for(int j = 0; j <= jmax-1; j++){
    PHI[imax-1][j] = PHI[imax-2][j];
  }
  for(int i = 0; i <= imax-1; i++){
    PHI[i][0] = 0.0;
    PHI[i][jmax-1] = 0.0;
  }
  //Set Phi Condition;
  for(int j = 0; j <= (jmax-1)/2; j++){
    PHI[0][j] = 1.0;
  }
}
void comp_FG(double **U, double **V, double **F, double **G, int imax, int jmax,  double Re, double dx, double dy, double dt){
  double d2udx2, d2udy2, d2vdx2, d2vdy2;
  double du2dx, duvdx, dvudy, dv2dy;
  double du2dxp1, duvdxp1, dvudyp1, dv2dyp1;
  double du2dxp2, duvdxp2, dvudyp2, dv2dyp2;
  //Evaluate F
  double gx = 0; double gy = 0;
  for(int i = 1; i <= imax-3; i++){
    for(int j = 1; j <= jmax-2; j++){
      d2udx2 = (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
      d2udy2 = (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
      du2dxp1 = 1/(4*dx)*((U[i+1][j]+U[i][j])*(U[i+1][j]+U[i][j])   -   (U[i-1][j]+U[i][j])*(U[i-1][j]+U[i][j]));
      dvudyp1 = 1/(4*dy)*((V[i+1][i]+V[i][j])*(U[i][j+1]+U[i][j])   -   (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]));
      du2dxp2 = 0;//1/(4*dx)*(abs(U[i+1][j]+U[i][j])*(U[i][j]-U[i+1][j])-abs(U[i-1][j]+U[i][j])*(-U[i][j]+U[i-1][j]));
      dvudyp2 = 0;//1/(4*dy)*(abs(V[i+1][i]+V[i][j])*(U[i][j]-U[i][j+1])-abs(V[i][j-1]+V[i+1][j-1])*(-U[i][j]+U[i][j-1]));
      dvudy = dvudyp1;
      du2dx = du2dxp1;
      F[i][j] = U[i][j]+dt*((1/Re)*(d2udx2+d2udy2)-(du2dx+dvudy)+gx);
    }
  }
  //Evaluate G
  for(int i = 1; i <= imax-2; i++){
    for(int j = 1; j <= jmax-3; j++){
      d2vdx2 = (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
      //cout << "part 1 " << n << endl;
      d2vdy2 = (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
      //cout << "part 2 " << n << endl;
      duvdxp1 = ((U[i][j+1]+U[i][j])*(V[i+1][j]+V[i][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))/(4*dx);
      //cout << "part 3 " << n << endl;
      dv2dyp1 = ((V[i][j+1]+V[i][j])*(V[i][j+1]+V[i][j])-(V[i][j-1]+V[i][j])*(V[i][j-1]+V[i][j]))/(4*dy);
      //cout << "part 4 " << n << endl;
      duvdxp2 = 0;//1/(4*dx)*(abs(U[i][j+1]+U[i][j])*(V[i][j]-V[i+1][j])-abs(U[i-1][j]+U[i-1][j+1])*(-V[i][j]+V[i-1][j]));
      //cout << "part 5 " << n << endl;
      dv2dyp2 = 0;//1/(4*dy)*(abs(V[i][j+1]+V[i][j])*(V[i][j]-V[i][j+1])-abs(V[i][j-1]+V[i][j])*(-V[i][j]+V[i][j-1]));
      //cout << "part 6 " << n << endl;
      duvdx = duvdxp1;
      dv2dy = dv2dyp1;
      G[i][j] = V[i][j]+dt*((1/Re)*(d2vdx2+d2vdy2)-(duvdx+dv2dy)+gy);
      //cout << n << endl;
    }
  }
  //NEWMANN BOUNDARY CONDITION//
  for(int j = 1; j <= jmax-2; j++){
    F[0][j] = U[0][j];
    F[imax-2][j] = U[imax-2][j];
  }
  //NEWMANN BOUUNDARY CONDITION//
  for(int i = 1; i <= imax-2; i++){
    G[i][0] = V[i][0];
    G[i][jmax-2] = V[i][jmax-2];
  }
  double Fmax = F[0][0];
  double Gmax = G[0][0];
  double Umax = U[0][0];
  double Vmax = V[0][0];
  for(int i = 0; i <= imax-2; i++){
    for(int j = 0; j <= jmax-2; j++){
      if(Fmax < F[i][j]){
	Fmax = F[i][j];
      }
      if(Gmax < G[i][j]){
	Gmax = G[i][j];
      }
      if(Vmax < V[i][j]){
	Vmax = V[i][j];
      }
      if(Umax < U[i][j]){
	Umax = U[i][j];
      }
    }
  }
  cout << "Fmax " << "\t"<< Fmax << "\t" << "Gmax " << "\t" << Gmax << "\t" << "Umax " << "\t" << Umax << "\t" << "Vmax " << "\t" << Vmax <<  endl;
}
void comp_RHS(double **F, double **G, double **RHS, int imax, int jmax, double dx, double dy, double dt){
  for(int i = 1; i <= imax-1; i++){
    for(int j = 1; j <= jmax-1; j++){
      RHS[i][j] = (1/dt)*((F[i][j]-F[i-1][j])/(dx)+(G[i][j]-G[i][j-1])/(dy));
    }
  }
  double RHSmax = RHS[0][0];
  for(int i = 0; i <= imax-2; i++){
    for(int j = 0; j <= jmax-2; j++){
      if(RHSmax <= RHS[i][j]){
	RHSmax = RHS[i][j];
      }
    }
  }
  cout << "RHSmax " << "\t" << RHSmax << "\t";
}
void poisson(double **P, double **RHS, int imax, int jmax, double dx, double dy, double dt, int max_iter){
  double Pmax = P[0][0];
  for(int i =0 ; i <= imax-2; i++){
    for(int j = 0; j <= jmax-2; j++){
      if(Pmax <= P[i][j]){
	Pmax = P[i][j];
      }
    }
  }
  cout << "Pmax " << "\t" << Pmax << endl;
  // int stepsize = 0;
  // while(true){
  //   //Still not complete yet lol
  //   int N = 0;
  //   double Pstar;
  //   int it = 0;
  //   for(int i = 1; i <= imax-2; i++){
  //     for(int j = 1; j <= jmax-2; j++){
  //       double p = P[i][j] ;
  // 	P[i][j] = ((dx*dx*dy*dy)/(2*(dx*dx+dy*dy)))*((1/(dx*dx))*(P[i+1][j]+P[i-1][j])+(1/(dy*dy))*(P[i][j+1]+P[i][j-1]))-RHS[i][j];
  //     }
  //   }
  //   stepsize += 1;
  //   cout << stepsize << endl;
  //   if(stepsize == 30){break;}
  // }
  // for(int j = 1; j <= jmax-1; j++){
  //   P[0][j]=  P[1][j];// P[imax-1][j] = P[imax-2][j];
  // }
  // for(int i = 1; i <= imax-1; i++){
  //   P[i][0] = P[i][1]; P[i][jmax-1] = P[i][jmax-2];
  // }
  double SOR = 1.7;
  for(int iter = 1; iter <= max_iter; iter++){
    //Setting up my pressure
    for(int j = 0; j <= jmax-1; j++){
      P[0][j] = P[1][j];
      P[imax-1][j] = P[imax-2][j];
    }
    for(int i = 0; i <= imax-1; i++){
      P[i][0] = P[i][1];
      P[i][jmax-1] = P[i][jmax-2];
    }
    //Iteration zone narok
    for(int i = 1; i<= imax-2; i++){
      for(int j = 1; j <= jmax-2; j++){
	double eiw, eie, ejs, ejn;
	if(i == 1){eiw = 1;}else{eiw = 1;}
	if(i == imax-2){eie = 1;}else{eie = 1;}
	if(j == 1){ejs = 1;}else{ejs = 1;}
	if(j == jmax-2){ejn = 1;}else{ejn = 1;}
	P[i][j] = (1-SOR)*P[i][j]+(SOR/((eie+eiw)/(dx*dx)+(ejs+ejn)/(dy*dy)))*((eie*P[i+1][j]+eiw*P[i-1][j])/(dx*dx)+(ejn*P[i][j+1]+ejs*P[i][j-1])/(dy*dy)-RHS[i][j]);
      }
    }
    if(iter % 5 == 0){
      double summ = 0.0;
      double max_norm = 0.0;
      double res_ij = 0.0;
      for(int i = 1; i<= imax-2; i++){
	for(int j = 1; j <= jmax-2; j++){
	  double eiw, eie, ejs, ejn;
	  if(i == 1){eiw = 1;}else{eiw = 1;}
	  if(i == imax-2){eie = 1;}else{eie = 1;}
	  if(j == 1){ejs = 1;}else{ejs = 1;}
	  if(j == jmax-2){ejn = 1;}else{ejn = 1;}
	  res_ij = (eie*(P[i+1][j]-P[i][j])-eiw*(P[i][j]-P[i-1][j]))/(dx*dx)+(ejn*(P[i][j+1]-P[i][j])-ejs*(P[i][j]-P[i][j-1]))/(dy*dy)-RHS[i][j];
	  summ += res_ij*res_ij;
	  if(res_ij > max_norm){max_norm = res_ij;}
	}
      }
      double norm  = sqrt(summ/((imax-2)*(jmax-2)));
      cout << iter << "\t" << norm << "\t" << max_norm << "\n";
      if(iter  >= 10 && norm <= 0.1) break;    
    }
  }
  cout << "=================================" << endl;
}
void simulation_uv(double **F, double **G, double **P, double **U, double **V, double Re, int imax, int jmax, double dx, double dy, double dt){
  for(int i = 1; i <= imax-3; i++){
    for(int j = 1; j <= jmax-2; j++){
      U[i][j] = F[i][j]-(dt*(P[i+1][j]-P[i][j]))/dx;
    }
  }
  for(int i = 1; i <= imax-2; i++){
    for(int j = 1; j <= jmax-3; j++){
      V[i][j] = G[i][j]-(dt*(P[i][j+1]-P[i][j]))/dy;
    }
  }
  //set u_velocity
  // Inflow Condition
  for(int j = 0; j <= jmax-1; j++){
    U[0][j] = 1.0;
  }
  // Outflow Condition
  for(int j = 0; j <= jmax-1; j++){
    U[imax-2][j] = U[imax-3][j];
    U[imax-1][j] = U[imax-1][j];
  }
  // No-Slip Condition
  for(int i = 0; i <= imax-1; i++){
    U[i][jmax-1] = -1.0*U[i][jmax-2]; //Top
    U[i][0] = -1.0*U[i][1]; // Bottom
  }
  //set v_velocity
  // inlet condition
  for(int j = 0; j <= jmax-1; j++){
    V[0][j] = -1.0*V[1][j];
  }
  // outlet
  for(int j = 0; j <= jmax-1; j++){
    V[imax-1][j] = V[imax-2][j];
  }
  // No-Slip Condition
  for(int i = 1; i <= imax-1; i++){
    V[i][jmax-2] = 0.0;
    V[i][0] = 0.0;
  }
}
void simulation_phi(double **U, double **V, double **P, double **PHI, double **NEW_PHI, double Re, int imax, int jmax, double dx, double dy, double dt){
  double LHS;
  double RHS;
  for(int i = 1; i<= imax-2; i++){
    for(int j =1; j <= jmax-2; j++){
      LHS = (U[i][j]/(2*dx))*(PHI[i+1][j]-PHI[i-1][j])+(V[i][j]/(2*dy))*(PHI[i][j+1]-PHI[i][j-1]);
      RHS = (1/Re)*(PHI[i+1][j]-2*PHI[i][j]+P[i-1][j])*(1/(dx*dx))+(1/Re)*(PHI[i][j+1]-2*PHI[i][j]+PHI[i][j-1])*(1/(dy*dy));
      NEW_PHI[i][j] = (RHS-LHS)*dt+PHI[i][j];
    }
  }
  
  //set phi condition
  for(int j = 0; j <= jmax-1; j++){
    PHI[imax-1][j] = PHI[imax-2][j];
    PHI[0][j] = 1.0;
  }
  for(int i = 0; i <= imax-1; i++){
    PHI[i][0] = 0.0;
    PHI[i][jmax-1] = 0.0;
  }
}
void update_phi(double **PHI, double **NEW_PHI, int imax, int jmax){
  for(int i =0 ; i <= imax-1; i++){
    for(int j = 0; j <= jmax-1; j++){
      PHI[i][j] = NEW_PHI[i][j];
    }
  }
}
void paraview(string filename, double **U, double **V, double **P, double **PHI, int imax, int jmax, double dx, double dy){
  ofstream myfile;
  myfile.open(filename);
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << imax-2 << " " <<  1 << " " << jmax-2 << "\n";
  myfile << "POINTS " << (imax-2)*1*(jmax-2) << " float\n";
  for(int j = 0; j <= jmax-3; j++){
    for(int i = 0; i <= imax-3; i++){
      myfile << dx*i << " " << dy*j << " 0\n"; 
    }
  }
  //DATA
  myfile << "\n";
  myfile << "POINT_DATA" << " " << (imax-2)*1*(jmax-2) << "\n";
  myfile << "\n";
  myfile << "SCALARS U float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int j = 1; j <= jmax-2; j++){
    for(int i = 1; i <= imax-2; i++){
      myfile << U[i][j] << endl;
    }
  }
  myfile << "\n";
  myfile << "SCALARS V float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int j = 1; j <= jmax-2; j++){
    for(int i = 1; i <= imax-2; i++){
      myfile << V[i][j] << endl;
    }
  }
  myfile << "\n";
  myfile << "SCALARS P float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int j = 1; j <= jmax-2; j++){
    for(int i = 1; i <= imax-2; i++){
      myfile << P[i][j] << endl;
    }
  }
  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int j = 1; j <= jmax-2; j++){
    for(int i = 1; i <= imax-2; i++){
      myfile << PHI[i][j] << endl;
    }
  }
  myfile.close();
}
int main(){
  int nx = 400;
  int ny = 200;
  int imax = nx+2;
  int jmax = ny+2;
  //Unknown bruh bruh bruh
  double length = 8.0;
  double width = 1.0;
  double dx = length/imax;
  double dy = width/jmax;
  double t = 0;
  double dt = 0.00001;
  double Re = 300;
  //MATRIX
  double **U;U = (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax-1; i++){U[i] = (double *) malloc ((jmax)*sizeof(double));}
  double **V;V = (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax; i++){V[i] = (double *) malloc ((jmax)*sizeof(double));}  
  double **F;F = (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax-1; i++){F[i] = (double *) malloc ((jmax)*sizeof(double));}
  double **G;G = (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax; i++){G[i] = (double *) malloc ((jmax)*sizeof(double));}
  double **P;P =  (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax; i++){P[i] = (double *) malloc ((jmax)*sizeof(double));}
  double **RHS; RHS =  (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax; i++){RHS[i] = (double *) malloc ((jmax)*sizeof(double));}
  double **PHI; PHI = (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax-1; i++){PHI[i] = (double *) malloc ((jmax)*sizeof(double));}
  double **NEW_PHI; NEW_PHI = (double **) malloc ((imax)*sizeof(double));
  for(int i = 0; i <= imax-1; i++){NEW_PHI[i] = (double *) malloc ((jmax)*sizeof(double));}
  initialize(U, V, F, G, P, RHS, PHI, NEW_PHI, imax, jmax);
  int period = 0;
  int pa = 0;
  int paraviewperiod = 100;
  while(true){
    setbcon(U,V,F,G,P,RHS,PHI,NEW_PHI,imax,jmax);
    comp_FG(U, V, F, G, imax, jmax, Re, dx, dy, dt);
    comp_RHS(F,G,RHS,imax,jmax, dx, dy, dt);
    if(t == 0){poisson(P, RHS, imax, jmax, dx, dy, dt, 100);}
    poisson(P, RHS, imax, jmax, dx, dy, dt, 30);
    simulation_uv(F, G, P, U, V, Re, imax, jmax, dx, dy, dt);
    simulation_phi(U, V, P, PHI, NEW_PHI, Re, imax, jmax, dx, dy, dt);
    string filename;
    period += 1;
    filename = "DATA_"+to_string(period)+".vtk";
    pa += 1;
    if(pa == 100){paraview(filename,U,V,P,PHI,imax,jmax,dx,dy);pa = 0;}
    t += dt;
    if(t == 50){break;}
  }
  return 0;
}
