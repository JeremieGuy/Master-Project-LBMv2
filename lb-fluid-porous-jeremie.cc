// complie with g++ -O3 lf-fluid.cc
// use C++17 or more   g++ -O3 -std=c++17 lf-fluid.cc

// this code implements a D2Q9 LB fluid model. The lattice is defined
// in the class D2Q9Fluid as a 2D vector of object of the class
// D2Q9Cell (vectors of f_in and f_out, and a local collision method,
// either BGK or Bounce-Back).  Physical parameters (viscosity, dx,
// dt, F), and lattice parameters (c_i, w_i, cs2) are defined in
// stuct.
// The methods of D2Q9Fluid allows us to setup a flow: here a
// Poiseueil flow in a canal with a body force and a porous region.

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <random>

using namespace std;

mt19937 generator;
uniform_real_distribution<double> rng;

struct D2Q9Param{
  vector<double> w;
  vector<vector<int>> c;
  double cs2;
};


struct PhysicsParam{
  double dt;
  double dx;
  double tau;
  double rho_0;
  double canalHeight;
  double canalLength;
  vector<double> F;
};

class D2Q9Cell
{
public:
  vector<double> fin, fout; 
  double Kx;   // resistance to flow; Kx=0 for a non-porous cell
  
  void (*myCollision)(D2Q9Cell&, const D2Q9Param& latticeParam, const PhysicsParam& physicsParam);

  D2Q9Cell(){
    fin=vector<double>(9);
    fout=vector<double>(9);
    Kx=0.; // default value
  }

  void setCollision(void (*collision)(D2Q9Cell&, const D2Q9Param& latticeParam, const PhysicsParam& physicsParam)){
    myCollision=collision;
  }
  
};  // end class D2Q9Cell


// This function implements the BGK collision on a D2Q9 lattice
void bgk(D2Q9Cell& cell, const D2Q9Param& latticeParam, const PhysicsParam& physicsParam)
{
  double rho;
  auto u=array{0., 0.};
  auto feq=0.;
  double Q;
  double FF;
  
  // compute density and velocity
  rho=0.;
  for(int k=0;k<9;k++)rho+=cell.fin[k];

  u[0]=0.;   u[1]=0.;
  
  for(int k=1;k<9;k++){
    u[0]+=latticeParam.c[k][0]*cell.fin[k];
    u[1]+=latticeParam.c[k][1]*cell.fin[k];
  }
  u[0]=u[0]/rho;      u[1]=u[1]/rho; 

  for(int k=0;k<9;k++){

    Q=0.;
    for(int a=0; a<2; a++)
      for(int b=0; b<2; b++)Q+=latticeParam.c[k][a]*latticeParam.c[k][b]*u[a]*u[b];

    Q =  Q/(2*latticeParam.cs2*latticeParam.cs2);
    Q -= (u[0]*u[0] + u[1]*u[1])/(2*latticeParam.cs2);

    feq=rho*latticeParam.w[k]*(1 
			       + (latticeParam.c[k][0]*u[0] + latticeParam.c[k][1]*u[1])/latticeParam.cs2
                               + Q
			       );

    // body force term: here we make it prop to rho, but this is a choice
    FF  = rho*(latticeParam.c[k][0]*(physicsParam.F[0] - cell.Kx*u[0])
	       + latticeParam.c[k][1]*physicsParam.F[1] );
    FF  = FF*(latticeParam.w[k]*physicsParam.dt*physicsParam.dt)/(physicsParam.dx*latticeParam.cs2);
    cell.fout[k] = cell.fin[k] + (1./physicsParam.tau)*(feq - cell.fin[k]);
    cell.fout[k] += FF;

  }  // end loop over k


} // end bgk()

void bounceBack(D2Q9Cell& cell, const D2Q9Param& latticeParam, const PhysicsParam& physicsParam)
{
  cell.fout[0]=cell.fin[0];
  for(int k=1;k<9;k++)cell.fout[k]=cell.fin[(k+3)%8 + 1];  // this is the opposit of k
}  // end bounceBack()


class D2Q9Fluid
{
public:
  int nx, ny;
  int ii, jj;
  PhysicsParam physicsParam;
  vector<vector<D2Q9Cell>> lattice;
  D2Q9Param latticeParam;
  string text;
  
  D2Q9Fluid(int Nx, int Ny, double dt_){
    nx=Nx;    ny=Ny;

    // define already known physical parameters
    physicsParam.dt=dt_;  

    // define lattice parameters
    latticeParam.w=vector<double>{4./9., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36.};
    latticeParam.c=vector<vector<int>>{{0,0}, {1,0}, {1,1}, {0,1}, {-1, 1}, {-1,0},
				  {-1,-1}, {0,-1}, {1,-1}};
    latticeParam.cs2=1./3.;

    // create a 2D lattice of cells
    lattice=vector<vector<D2Q9Cell>>( nx, vector<D2Q9Cell>(ny) );

  }  //---- end constructor


  // init porous region of resistence Kx_ for all j and between i1 and i2
  void setPorousRegion(double Kx_, int i1, int i2){
    for(int i=i1;i<=i2;i++)
      for(int j=1; j<ny-1; j++)
	lattice[i][j].Kx=Kx_;
  }

  void setGridSpacing(double dx_){
    physicsParam.dx=dx_;
  }
  
  void init_point(int i, int j, int k){
    lattice[i][j].fin[k]=1;
  }

  // init the lattice with fin[i]=w[i]*rho_0
  void init(double rho_){
    physicsParam.rho_0=rho_;
    for(int i=0;i<nx;i++)
      for(int j=0; j<ny; j++)
	for(int k=0;k<9;k++)lattice[i][j].fin[k]=latticeParam.w[k]*physicsParam.rho_0;
  }


  // define the cells of the lattice as a canal with the upper and lower walls with bounce back rule
  void setBounceBackCanal(double canal_height, double visco, double Fx, double Fy){

    for(int i=0;i<nx;i++)lattice[i][0].setCollision(&bounceBack);
    //    for(int i=0;i<nx;i++)lattice[i][0].setCollision(&lowerWall); // this does not work
    for(int i=0;i<nx;i++)lattice[i][ny-1].setCollision(&bounceBack);


    for(int i=0;i<nx;i++)
         for(int j=1;j<ny-1;j++)lattice[i][j].setCollision(&bgk);

    physicsParam.dx=canal_height/(double)(ny-1);  // ny-1 because of bounce back in the y-axis
    physicsParam.canalHeight=canal_height;
    physicsParam.canalLength=nx*physicsParam.dx; // we have periodic boundary conditions along x 
    physicsParam.tau=0.5 + 3*visco*physicsParam.dt/(physicsParam.dx*physicsParam.dx);
    physicsParam.F=vector<double>{Fx, Fy};

    text="(Velocity profile, LB Poiseuielle flow, 2 bounceback walls)";

  }


    // define the cells of a periodic lattice with BGK collision
  void setPeriodicSystem(double dx, double visco){
    physicsParam.dx=dx;
    physicsParam.tau=0.5 + 3*visco*physicsParam.dt/(physicsParam.dx*physicsParam.dx);
    for(int i=0;i<nx;i++)
         for(int j=0;j<ny;j++)lattice[i][j].setCollision(&bgk);
  }

  // this implements a propagation with periodic boundaries
  void propagation(){
    for(int i=0;i<nx;i++)
      for(int j=0;j<ny;j++)
	for(int k=0;k<9;k++){
	  ii=(i+nx-latticeParam.c[k][0])%nx;
	  jj=(j+ny-latticeParam.c[k][1])%ny;
	  lattice[i][j].fin[k]=lattice[ii][jj].fout[k];
      }// end  loops
  } // end propagation
  
  // exectute the method collision in each lattce cell
  void collision(){
    for(int j=0;j<ny;j++)
      for(int i=0;i<nx;i++)lattice[i][j].myCollision(lattice[i][j], latticeParam, physicsParam);
  }

  void printTotalMassAndMomentum(){
    double rhoTot=0.;
    double Jx=0.;
    double Jy=0.;
    for(int i=0;i<nx;i++)
      for(int j=0;j<ny;j++)
	for(int k=0;k<9;k++){
	  rhoTot += lattice[i][j].fin[k];
	  Jx += latticeParam.c[k][0]*lattice[i][j].fin[k];
	  Jy += latticeParam.c[k][1]*lattice[i][j].fin[k];
	}
    cout<<"Total mass="<<rhoTot<<"  total Jx="<<Jx<<"  total Jy="<<Jy<<endl;
  }

  void printVelocities(int i, int j){

    auto u=vector{0.,0.};
    auto rho=0.;
    for(int k=0; k<9;k++)rho+=lattice[i][j].fin[k];
    
    for(int k=1; k<9;k++){
      u[0]+=latticeParam.c[k][0]*lattice[i][j].fin[k];
      u[1]+=latticeParam.c[k][1]*lattice[i][j].fin[k];
    }
    u[0]=u[0]/rho;        u[1]=u[1]/rho;
    
    cout<<rho<<"  "<< u[0]<<"  "<< u[1]<<endl;
  }  // end printVelocities


  // return lattice[i][j].fin[k]
  double getFin(int i, int j, int k){
    return lattice[i][j].fin[k];
  }
  
  // return the density rho at lattice point (i,j)
  double getDensity(int i, int j){
      auto rho=0.;
      for(int k=0; k<9;k++)rho+=lattice[i][j].fin[k];
      return rho;
  }
 
  // return the velocity in physical units at lattice point (i,j)
  vector<double> getVelocity(int i, int j){
      auto rho=0.;
      auto u=vector{0., 0.};

      for(int k=0; k<9;k++)rho+=lattice[i][j].fin[k];
    
      for(int k=1; k<9;k++){
	u[0]+=latticeParam.c[k][0]*lattice[i][j].fin[k];
	u[1]+=latticeParam.c[k][1]*lattice[i][j].fin[k];
      }
      u[0]=(physicsParam.dx/physicsParam.dt)*u[0]/rho;
      u[1]=(physicsParam.dx/physicsParam.dt)*u[1]/rho;
      return u;
  }

  
  void PSprintVelocityProfile(int i, int max_iter, ofstream& fout){

    auto u=vector{0.,0.};
    double rho;
    double dx=physicsParam.dx;
    double dt=physicsParam.dt;
    double v=dx/dt;

    fout<<"/nx "<<nx<<" def"<<endl;
    fout<<"/ny "<<ny<<" def"<<endl;

    fout<<"/dx "<<physicsParam.dx<<" def"<<endl;
    fout<<"/dt "<<physicsParam.dt<<" def"<<endl;
    fout<<"/Fx "<<physicsParam.F[0]<<" def"<<endl;
    fout<<"/Fy "<<physicsParam.F[1]<<" def"<<endl;
    fout<<"/tau "<<physicsParam.tau<<" def"<<endl;
    fout<<"/rho_0 "<<physicsParam.rho_0<<" def"<<endl;
    fout<<"/t_max "<<max_iter*dt<<" def"<<endl;
    fout<<"/H "<< physicsParam.canalHeight<<" def"<<endl;
    fout<<"/L "<< physicsParam.canalLength<<" def"<<endl;
    fout<<"/caption-text "<< text<<" def"<<endl;

    fout<<endl;
    
    fout<<"startReadData"<<endl;
    for(int j=0;j<ny;j++){
      rho=0.; u[0]=0.;  u[1]=0.;
      for(int k=0; k<9;k++)rho+=lattice[i][j].fin[k];
    
      for(int k=1; k<9;k++){
	u[0]+=latticeParam.c[k][0]*lattice[i][j].fin[k];
	u[1]+=latticeParam.c[k][1]*lattice[i][j].fin[k];
      }
      u[0]=u[0]/rho;        u[1]=u[1]/rho;

      //      cout<<j*dx<<"  "<<rho<<"  "<<v*u[0]<<"   "<<v*u[1]<<endl;
      fout<<j*dx<<"  "<<rho<<"  "<<v*u[0]<<"   "<<v*u[1]<<"  addData"<<endl;
    }
    fout<<"endReadData\n"<<endl;
  } // end PSprintVelocityProfile()

  void PSprintVelocityPressureAlongCanal(int max_iter, ofstream& fout){

    auto u=vector{0.,0.};
    double rho;
    double  rho_i,   ux_i,  uy_i;
    double dx=physicsParam.dx;
    double dt=physicsParam.dt;
    double v=dx/dt;

    fout<<"/nx "<<nx<<" def"<<endl;
    fout<<"/ny "<<ny<<" def"<<endl;

    fout<<"/dx "<<physicsParam.dx<<" def"<<endl;
    fout<<"/dt "<<physicsParam.dt<<" def"<<endl;
    fout<<"/Fx "<<physicsParam.F[0]<<" def"<<endl;
    fout<<"/Fy "<<physicsParam.F[1]<<" def"<<endl;
    fout<<"/tau "<<physicsParam.tau<<" def"<<endl;
    fout<<"/rho_0 "<<physicsParam.rho_0<<" def"<<endl;
    fout<<"/t_max "<<max_iter*dt<<" def"<<endl;
    fout<<"/H "<< physicsParam.canalHeight<<" def"<<endl;
    fout<<"/L "<< physicsParam.canalLength<<" def"<<endl;
    fout<<"/caption-text "<< text<<" def"<<endl;

    fout<<endl;
    
    fout<<"startReadData"<<endl;
    for(int i=0; i<nx; i++){
      rho_i=0.;   ux_i=0.;   uy_i=0.;
      for(int j=0;j<ny;j++){
	rho=0.; u[0]=0.;  u[1]=0.;
	for(int k=0; k<9;k++)rho+=lattice[i][j].fin[k];
    
	for(int k=1; k<9;k++){
	  u[0]+=latticeParam.c[k][0]*lattice[i][j].fin[k];
	  u[1]+=latticeParam.c[k][1]*lattice[i][j].fin[k];
	}
	u[0]=u[0]/rho;        u[1]=u[1]/rho;

	rho_i += rho;	ux_i += u[0];    uy_i+= u[1];

      } // end j loop

      fout<<i*dx<<"  "<<rho_i/ny<<"  "<<v*ux_i/ny<<"   "<<v*uy_i/ny<<"  addData"<<endl;

    } // end i loop
    fout<<"endReadData\n"<<endl;
  } // end PSprintVelocityPressureAlongCanal()


  void printVelocityProfile(int i, int max_iter, ofstream& fout){

    auto u=vector{0.,0.};
    double rho;
    double dx=physicsParam.dx;
    double dt=physicsParam.dt;
    double v=dx/dt;

    fout<<"nx="<<nx<<endl;
    fout<<"ny="<<ny<<endl;

    fout<<"dx="<<physicsParam.dx<<endl;
    fout<<"dt="<<physicsParam.dt<<endl;
    fout<<"Fx="<<physicsParam.F[0]<<endl;
    fout<<"Fy="<<physicsParam.F[1]<<endl;
    fout<<"tau="<<physicsParam.tau<<endl;
    fout<<"rho_0="<<physicsParam.rho_0<<endl;
    fout<<"t_max "<<max_iter*dt<<endl;
    fout<<"canalHeight= "<< physicsParam.canalHeight<<endl;

    fout<<endl;
    
    for(int j=0;j<ny;j++){
      rho=0.; u[0]=0.;  u[1]=0.;
      for(int k=0; k<9;k++)rho+=lattice[i][j].fin[k];
    
      for(int k=1; k<9;k++){
	u[0]+=latticeParam.c[k][0]*lattice[i][j].fin[k];
	u[1]+=latticeParam.c[k][1]*lattice[i][j].fin[k];
      }
      u[0]=u[0]/rho;        u[1]=u[1]/rho;

      fout<<j*dx<<"  "<<rho<<"  "<<v*u[0]<<"   "<<v*u[1]<<endl;
    }
  }  // end printVelocityProfile()

  
  void showFin(int k){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++)cout<<lattice[i][j].fin[k];
      cout<<endl;
    }
  } // end showFin
  
}; // end class D2Q9Fluid


int main()
{
  // in case we need a random number generator (defined as a global function)
  generator = mt19937(time(0));
  rng = uniform_real_distribution<double>(0.0, 1.0);

  auto nx=100;   auto ny=51;  // number of grid points
  auto dt=1.;                // time step
  auto rho_0=2.5;            // initial density
  auto Fx=1.e-4;             // body force along x-axis: rho*Fx
  auto Kx=2e-3;              // body force in the porpus region
  auto Fy=0.;                // same for y axis
  auto visco=5.e-3;          // viscosity in the chosen physical units
  auto canal_height=10;      // height of the canal in physical units. Grid spacing dx is determined from ny
                             // canal length L=nx*dx is determined by the choice of nx
  //  double tau;            // tau is determined according to visco, dx and dt
  int i,j,k;
  


  // instantiate a D2Q9 lattice with the given parameters
  D2Q9Fluid myFluid(nx,ny,dt);


  myFluid.setBounceBackCanal(canal_height,visco,Fx,Fy); // this creates bounceback cells on the upper and lower wall j=0 and j=ny-1

  myFluid.setPorousRegion(Kx, nx/2-10, nx/2+10);

  //  double dx=1;
  //  myFluid.setPeriodicSystem(dx, visco);  // this creates a periodic BGK lattce with given grid spacing and physical visco

  myFluid.init(rho_0);       // initialize the fin[k] as w[k]*rho_0 in each cell [i][j]
  myFluid.printTotalMassAndMomentum();  // to have the initial mass and momentum in the lattice

  // get fin[k] at lattice point (i,j). Note this is valid after initialization or after propagation
  i=nx/2;    j=ny/2;  k=2;
  cout<<"fin[i][j][k]="<<myFluid.getFin(i,j,k)<<endl;
  
  // let us run the fluid simulation for max_iter steps
  int max_iter=25000;  // number of time steps
  for(int iter=0; iter<max_iter;iter++){
    myFluid.collision(); 
    myFluid.propagation();
  }

  // an example to get the velocity in physical units and density at lattice point (i,j)
  i=nx/2;    j=ny/2;
  auto u=myFluid.getVelocity(i,j);
  cout<<"physical velocity at lattice point ("<<i<<","<<j<<")=("<<u[0]<<","<<u[1]<<")"<<endl;

  auto rho=myFluid.getDensity(i,j);
  cout<<"Density at lattice point ("<<i<<","<<j<<")="<<rho<<endl;
    
  // let us print the simulation data (position density and velocity in physical units) and flow parameters in a file
  ofstream fout("velocityPressureAlongCanal.dat");
  fout<<"/Kx "<<Kx<<" def"<<endl;
  myFluid.PSprintVelocityPressureAlongCanal(max_iter, fout);

  //  ofstream fout("velocity-profile.dat");
  //  myFluid.PSprintVelocityProfile(nx/2,max_iter,fout);   // here we write a file for postscript output 
  //  myFluid.printVelocityProfile(nx/2,max_iter,fout);   // here we write a text file
  fout.close();

  // we can check mass consrvation and momentum balance (not concerved with no-slip boundaries 
  myFluid.printTotalMassAndMomentum();
  
  return 0;
}
