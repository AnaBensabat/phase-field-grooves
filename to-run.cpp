#include<bits/stdc++.h>
#include<random>
#include<omp.h>
#include"auxFunctions.h"
#include"grids.h"

void run_simulation(double kappa, double eta, string pwd){
  
  //OPENMP
  omp_set_num_threads(16);
  
  //SET DIMENSIONS OF THE SYSTEM
  double scale = 0.5;
  
  int dimX = 150*scale;
  int dimY = 150*scale;
  int dimZ = 100*scale;
  
  //MAKE GROOVES
  //geometry parameters
  int depth = 0; //to test make it plane!
  int spacing = 0;//10*scale;
  int width = 2*10*scale;

  //dynamic parameters
  double alpha = 0.001;
  double epsilon = 1;
  //construct
  cout<<"----CREATING GROOVES----"<<endl;  
  Groove mygrooves(dimX, dimY, dimZ, depth, spacing, width, epsilon, alpha);
  cout<<"----GROOVES DONE----"<<endl;  
  double dt = 0.0001;
  double time = 5;
  cout<<"----GROOVE STABILIZATION----"<<endl;
  mygrooves.stabilize(dt, time, 5);
  saveGridToVTI(pwd+"environment_0.vti",mygrooves.grid);
  cout<<"----GROOVE STABILIZATION----"<<endl;

  //MAKE CELL
  //geometry parameters
  vector<double> radius = {20*scale,10*scale}; //for the right scale radius = {15,15}  radius_nucleus = {9}
  vector<double> radius_nucleus = {0*scale,0*scale};
  depth = width;
  vector<double> center_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth+(double)radius[1]*2};
  vector<double> nucleus_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth+(double)radius[1]*2};
  
  cout<<"----CREATING CELL----"<<endl;  
  Cell mycell (dimX, dimY, dimZ, center_coordinates, radius, nucleus_coordinates, radius_nucleus, true);
  saveGridToVTI(pwd+"cell_0.vti",mycell.grid);
  saveGridToVTI(pwd+"nucleus_0.vti",mycell.grid_nucleus);
  cout<<"----CELL DONE----"<<endl;  

  
  //evolve cell
  mycell.epsilon = 0.5;         //surf->membrane thickness coefficient
  mycell.alpha = 0.001/12;      //volume coefficient
  mycell.alpha_s = 0.001;       //surface area coefficient
  mycell.gamma = 1;             //repulsion term
  mycell.eta = 2;         //adhesion coefficient 
  mycell.kappa = kappa;         //bending rigidity fraction
  double velocity = 0.;         //drag velocity
  double velocity_nuc = 0;

  string s1 = "cell_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  string s2 = "nucleus_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  string s3 = "environment_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);

  cout<<"----CELL STABILIZATION----"<<endl;  
  time = 3;
  evolve(mycell, mygrooves.grid, 0, 0, dt, time, pwd, true);
  cout<<"----CELL STABILIZATION----"<<endl;

  s1 = "cell_01_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s2 = "nucleus_01_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s3 = "environment_01_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);

  double v0_cell = vol(mycell.grid,1e-6);
  double v0_nucleus = vol(mycell.grid_nucleus,1e-6);

  double s0_cell = area(mycell.grid,mycell.epsilon);
  double s0_nucleus = area(mycell.grid_nucleus,mycell.epsilon);
  
  cout<<"Initial Cell Volume: "<<v0_cell<<'\n';
  cout<<"Initial Nucleus Volume: "<<v0_nucleus<<'\n';
  cout<<"Initial Cell Surface Area: "<<s0_cell<<'\n';
  cout<<"Initial Nucleus Surface Area: "<<s0_nucleus<<'\n';
  
  cout<<"----STARTING TIME EVOLUTION----"<<endl;
  time = 1;
  evolve(mycell, mygrooves.grid, velocity, velocity_nuc,  dt, time, pwd, 100000000);
  cout<<"----TIME EVOLUTION DONE----"<<endl;
  
  s1 = "cell_F_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s2 = "nucleus_F_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s3 = "environment_F_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);

  return ;
}


int main(){

  string pwd = filesystem::current_path();
  filesystem::create_directory("output");
  pwd+="/output/";
  clearDirectory(pwd);
  
  run_simulation(1,1,pwd);

  return 0;
}
