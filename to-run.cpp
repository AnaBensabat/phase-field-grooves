#include<bits/stdc++.h>
#include<random>
#include<omp.h>
#include"auxFunctions.h"
#include"grids.h"

void run_simulation(double kappa, double eta, string pwd){
  
  //OPENMP
  omp_set_num_threads(32);
  
  //SET DIMENSIONS OF THE SYSTEM
  double scale = 0.75;
  
  int dimX = 150*scale;
  int dimY = 150*scale;
  int dimZ = 70*scale;
  
  //MAKE GROOVES
  //geometry parameters
  int depth = 11*scale;   //to test make it plane!
  int spacing = 13*scale;
  int width = 13*scale;

  //dynamic parameters
  double alpha = 0.001;
  double epsilon = 1;
  //construct
  cout<<"----CREATING GROOVES----"<<endl;  
  Groove mygrooves(dimX, dimY, dimZ, depth, spacing, width, epsilon, alpha);
  cout<<"----GROOVES DONE----"<<endl;  
  double dt = 0.0001;
  double time = 0.1;
  cout<<"----GROOVE STABILIZATION----"<<endl;
  mygrooves.stabilize(dt, time, 5);
  saveGridToVTI(pwd+"environment_0.vti",mygrooves.grid);
  cout<<"----GROOVE STABILIZATION----"<<endl;

  //MAKE CELL
  //geometry parameters
  vector<double> radius = {27*scale,17*scale}; //for the right scale radius = {15,15}  radius_nucleus = {9}
  vector<double> radius_nucleus = {0,0}; //{10*scale,7*scale};
  depth = width;
  vector<double> center_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth*2-2};
  vector<double> nucleus_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth*2};
  
  cout<<"----CREATING CELL----"<<endl;  
  Cell mycell (dimX, dimY, dimZ, center_coordinates, radius, nucleus_coordinates, radius_nucleus, true, true);
  saveGridToVTI(pwd+"cell_0.vti",mycell.grid);
  saveGridToVTI(pwd+"nucleus_0.vti",mycell.grid_nucleus);
  cout<<"----CELL DONE----"<<endl;
  
  //evolve cell
  mycell.epsilon = 0.5;         //surf->membrane thickness coefficient
  mycell.alpha = 0.005/12;      //volume coefficient
  mycell.alpha_nuc = 0.005;
  mycell.alpha_s = 0.002;       //surface area coefficient
  mycell.gamma = 1;             //repulsion term
  mycell.eta = eta;             //adhesion coefficient 
  mycell.kappa = 1;             //bending rigidity fraction
  double velocity = 0.;         //drag velocity
  double velocity_nuc = 0;

  string s1 = "cell_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  string s2 = "nucleus_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  string s3 = "environment_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);
  
  cout<<"BEFORE STABILIZATION"<<'\n';
  double v0_cell = vol(mycell.grid,1e-6);
  double v0_nucleus = vol(mycell.grid_nucleus,1e-6);
  
  double s0_cell = area(mycell.grid,mycell.epsilon);
  double s0_nucleus = area(mycell.grid_nucleus,mycell.epsilon);
  
  cout<<"Initial Cell Volume: "<<v0_cell<<'\n';
  cout<<"Initial Nucleus Volume: "<<v0_nucleus<<'\n';
  cout<<"Initial Cell Surface Area: "<<s0_cell<<'\n';
  cout<<"Initial Nucleus Surface Area: "<<s0_nucleus<<'\n';
  
  cout<<"----CELL STABILIZATION----"<<endl;  
  time = 0.01;
  evolve(mycell, mygrooves.grid, 0, 0, dt, time, pwd+"eta="+to_string(eta)+"_k="+to_string(kappa)+"_", true);
  cout<<"----CELL STABILIZATION----"<<endl;

  /*
  s1 = "cell_01_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s2 = "nucleus_01_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s3 = "environment_01_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);
  */
  
  //mycell.grid = loadGridFromVTI("cell_S.vti");
  //mycell.grid_nucleus = loadGridFromVTI("nucleus_S.vti");
  
  v0_cell = vol(mycell.grid,1e-6);
  v0_nucleus = vol(mycell.grid_nucleus,1e-6);

  s0_cell = area(mycell.grid,mycell.epsilon);
  s0_nucleus = area(mycell.grid_nucleus,mycell.epsilon);
  
  cout<<"Initial Cell Volume: "<<v0_cell<<'\n';
  cout<<"Initial Nucleus Volume: "<<v0_nucleus<<'\n';
  cout<<"Initial Cell Surface Area: "<<s0_cell<<'\n';
  cout<<"Initial Nucleus Surface Area: "<<s0_nucleus<<'\n';

  cout<<"----STARTING TIME EVOLUTION----"<<endl;
  //time = 0.01;
  //evolve(mycell, mygrooves.grid, velocity, velocity_nuc,  dt, time, pwd, false, true, 1000);
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

  vector<double> etas={0.5,1};
  for(auto eta:etas){
    run_simulation(1,eta,pwd);
  }
  

  return 0;
}
