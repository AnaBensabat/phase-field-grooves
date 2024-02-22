#include<bits/stdc++.h>
#include<random>
#include<omp.h>
#include"auxFunctions.h"
#include"grids.h"

int main(){

  string s_cell = "";
  string s_nucleus = "";
  string s_env = "";
  vector<vector<vector<double>>> cell_grid = loadGridFromVTI(s_cell);
  vector<vector<vector<double>>> nucleus_grid = loadGridFromVTI(s_nucleus);
  vector<vector<vector<double>>> env_grid = loadGridFromVTI(s_env);

  cout<<"----POST-PROCESSING----\n";
  double depth = 0;
  double v0_cell = 0;
  double v0_nucleus = 0;

  vector<double> volumes = DFS(cell_grid, 2*depth, 0.4);
  vector<double> volumes_nuc = DFS(nucleus_grid, 2*depth, 0.4);
  double vf_cell = vol(cell_grid,1e-6);
  double vf_nucleus = vol(nucleus_grid,1e-6);
  double vol_groove_cell = vol_in_grooves(cell_grid, depth*2, 1e-6)/vf_cell;
  double vol_groove_nucleus = vol_in_grooves(nucleus_grid, depth*2, 1e-6)/vf_nucleus;
  double vol_loss_cell = 1-vf_cell/v0_cell;
  double vol_loss_nucleus = 1-vf_nucleus/v0_nucleus;
  double split = isSplit(cell_grid, 0.4);
  double split_nucleus = isSplit(nucleus_grid, 0.4);
  
  double blows_nucleus = 0;
  double blows = 0;
  
  if (vf_cell<1e-3) blows = 1;
  if (vf_nucleus<1e-3) blows_nucleus = 1;
  
  cout<<"----POST-PROCESSING DONE----\n";
  cout<<'\n';  

  vector<double> data_cell = {vol_loss_cell, vol_groove_cell, (double)volumes.size(), split, blows};
  vector<double> data_nucleus = {vol_loss_nucleus, vol_groove_nucleus, (double)volumes_nuc.size(), split_nucleus, blows_nucleus};

  return 0;
}
