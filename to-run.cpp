#include<bits/stdc++.h>
#include<random>
#include<omp.h>
#include"auxFunctions.h"
#include"grids.h"

vector<vector<double>> run_simulation(double kappa, double eta, string pwd){
  
  // threads
  int thread_num;
  omp_set_num_threads(16);
  thread_num = omp_get_max_threads ();
  
  //SET DIMENSIONS OF THE SYSTEM
  double scale = 1;
  
  int dimX = 150*scale;
  int dimY = 150*scale;
  int dimZ = 100*scale;

  // 5 microns for the grooves -> 10 points
  // cell radius of 25 -> 50 points
  // nucleus radius of 10 -> 20 points >.<
  
  //MAKE GROOVES
  //geometry parameters
  int depth = 0;//10*scale;
  int spacing = 10*scale;
  int width = 2*10*scale;

  //dynamic parameters
  double alpha = 0.001;
  double epsilon = 1;
  //construct
  cout<<"----CREATING GROOVES----"<<endl;  
  Groove mygrooves(dimX, dimY, dimZ, depth, spacing, width, epsilon, alpha);
  cout<<"----GROOVES DONE----"<<endl;  
  double dt = 0.007;
  double time = 5;
  cout<<"----GROOVE STABILIZATION----"<<endl;
  mygrooves.stabilize(dt, time, 5);
  //saveGridToVTI(pwd+"environment_0.vti",mygrooves.grid);
  cout<<"----GROOVE STABILIZATION----"<<endl;

  //MAKE CELL
  //geometry parameters
  vector<double> radius = {15*scale,15*scale};
  double radius_nucleus = 9*scale;
  depth = width;
  vector<double> center_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth+(double)radius[1]+1};
  vector<double> nucleus_coordinates = {(double)dimX/2,(double)dimY/2,(double)depth+(double)radius[1]+1};
  
  cout<<"----CREATING CELL----"<<endl;  
  Cell mycell (dimX, dimY, dimZ, center_coordinates, radius, nucleus_coordinates, radius_nucleus);
  //saveGridToVTI(pwd+"cell_0.vti",mycell.grid);
  //saveGridToVTI(pwd+"nucleus_0.vti",mycell.grid_nucleus);
  cout<<"----CELL DONE----"<<endl;  

  //evolve cell  
  time = 50;

  // repetir este com t=50 e eps = 0.5
  mycell.epsilon = 0.5;         //surf->membrane thickness coefficient
  mycell.alpha = 0.001/12;    //volume coefficient
  mycell.gamma = 2;         //repulsion term
  mycell.eta = eta*6/6;          //adhesion coefficient 
  mycell.kappa = kappa;        //surface tension fraction // este é o rácio
  double velocity = 0.;       //drag velocity
  double velocity_nuc = 0;
  
  cout<<"----CELL STABILIZATION----"<<endl;  
  evolve(mycell, mygrooves.grid, 0, 0, dt, 10, pwd);
  cout<<"----CELL STABILIZATION----"<<endl;

  string s1 = "cell_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  string s2 = "nucleus_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  string s3 = "environment_0_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);
  
  double v0_cell = vol(mycell.grid,1e-6);
  double v0_nucleus = vol(mycell.grid_nucleus,1e-6);

  double blows = 0;
  if (v0_cell<1e-3) {
    blows = 1;
  }
  else{
    cout<<"before "<<vol(mycell.grid,1e-6)<<'\n';
    cout<<"----STARTING TIME EVOLUTION----"<<endl;
    evolve(mycell, mygrooves.grid, velocity, velocity_nuc,  dt, time, pwd, 100000000);
    cout<<"----TIME EVOLUTION DONE----"<<endl;
  }
  
  s1 = "cell_F_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s2 = "nucleus_F_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  s3 = "environment_F_eta="+to_string(eta)+"_k="+to_string(kappa)+".vti";
  saveGridToVTI(pwd+s1,mycell.grid);
  saveGridToVTI(pwd+s2,mycell.grid_nucleus);
  saveGridToVTI(pwd+s3,mygrooves.grid);
  
  cout<<"\n";
  cout<<"----POST-PROCESSING----\n";
  vector<double> volumes = DFS(mycell.grid, 2*depth, 0.4);
  vector<double> volumes_nuc = DFS(mycell.grid_nucleus, 2*depth, 0.4);
  double vf_cell = vol(mycell.grid,1e-6);
  double vf_nucleus = vol(mycell.grid_nucleus,1e-6);
  double vol_groove_cell = vol_in_grooves(mycell.grid, depth*2, 1e-6)/vf_cell;
  double vol_groove_nucleus = vol_in_grooves(mycell.grid_nucleus, depth*2, 1e-6)/vf_nucleus;
  double vol_loss_cell = 1-vf_cell/v0_cell;
  double vol_loss_nucleus = 1-vf_nucleus/v0_nucleus;
  double split = isSplit(mycell.grid, 0.4);
  double split_nucleus = isSplit(mycell.grid_nucleus, 0.4);
  
  double blows_nucleus = 0;
  if (vf_cell<1e-3) blows = 1;
  if (vf_nucleus<1e-3) blows_nucleus = 1;
  
  cout<<"----POST-PROCESSING DONE----\n";
  cout<<'\n';  

  vector<double> data_cell = {vol_loss_cell, vol_groove_cell, (double)volumes.size(), split, blows};
  vector<double> data_nucleus = {vol_loss_nucleus, vol_groove_nucleus, (double)volumes_nuc.size(), split_nucleus, blows_nucleus};

  return {data_cell,data_nucleus};
}

int main(){

  string pwd = filesystem::current_path();
  filesystem::create_directory("output");
  pwd+="/output/";
  
  vector<double> etas = {2};
  vector<double> kappas = {2};

  vector<vector<double>> volume_losses_cell(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> volume_grooves_cell(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> N_grooves_cell(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> Split_cell(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> Blows_cell(kappas.size(), vector<double>(etas.size()));

  
  vector<vector<double>> volume_losses_nucleus(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> volume_grooves_nucleus(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> N_grooves_nucleus(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> Split_nucleus(kappas.size(), vector<double>(etas.size()));
  vector<vector<double>> Blows_nucleus(kappas.size(), vector<double>(etas.size()));

  for(int i=0;i<kappas.size();i++){
    for(int j=0;j<etas.size();j++){
      auto x = run_simulation(kappas[i],etas[j],pwd);

      volume_losses_cell[i][j] = x[0][0];
      volume_grooves_cell[i][j] = x[0][1];
      N_grooves_cell[i][j] = x[0][2];
      Split_cell[i][j] = x[0][3];
      Blows_cell[i][j] = x[0][4];
      
      volume_losses_nucleus[i][j] = x[1][0];
      volume_grooves_nucleus[i][j] = x[1][1];
      N_grooves_nucleus[i][j] = x[1][2];
      Split_nucleus[i][j] = x[1][3];
      Blows_nucleus[i][j] = x[1][4];
      
      ofstream volume_losses_info;
      volume_losses_info.open (pwd+"volume_losses.txt");
      for(auto x:volume_losses_cell){
	for(auto y:x) volume_losses_info << y <<'\t';
	volume_losses_info << '\n';
      }
      volume_losses_info << '\n';
      for(auto x:volume_losses_nucleus){
	for(auto y:x) volume_losses_info << y <<'\t';
	volume_losses_info << '\n';
      }
      volume_losses_info.close();

      ofstream volume_grooves_info;
      volume_grooves_info.open (pwd+"volume_grooves.txt");
      for(auto x:volume_grooves_cell){
	for(auto y:x) volume_grooves_info << y <<'\t';
	volume_grooves_info << '\n';
      }
      volume_grooves_info << '\n';
      for(auto x:volume_grooves_nucleus){
	for(auto y:x) volume_grooves_info << y <<'\t';
	volume_grooves_info << '\n';
      }
      volume_grooves_info.close();

      ofstream N_grooves_info;
      N_grooves_info.open (pwd+"N_grooves.txt");
      for(auto x:N_grooves_cell){
	for(auto y:x) N_grooves_info << y <<'\t';
	N_grooves_info << '\n';
      }
      N_grooves_info << '\n';
      for(auto x:N_grooves_nucleus){
	for(auto y:x) N_grooves_info << y <<'\t';
	N_grooves_info << '\n';
      }
      N_grooves_info.close();
      
      ofstream Split_info;
      Split_info.open (pwd+"Split.txt");
      for(auto x:Split_cell){
	for(auto y:x) Split_info << y <<'\t';
	Split_info << '\n';
      }
      Split_info << '\n';
      for(auto x:Split_nucleus){
	for(auto y:x) Split_info << y <<'\t';
	Split_info << '\n';
      }
      Split_info.close();

      ofstream Blows_info;
      Blows_info.open (pwd+"Blows.txt");
      for(auto x:Blows_cell){
	for(auto y:x) Blows_info << y <<'\t';
	Blows_info << '\n';
      }
      Blows_info << '\n';
      for(auto x:Blows_nucleus){
	for(auto y:x) Blows_info << y <<'\t';
	Blows_info << '\n';
      }
      Blows_info.close();
    }
  }
}
