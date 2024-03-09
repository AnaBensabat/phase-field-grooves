using namespace std;
using matrix = vector<vector<vector<double>>>; 

// AUX
inline double h_function(double x){
  return x*x*(3-2*x);
}

inline double laplacian_h(matrix &grid, int i, int j, int k){
  int dimX = grid.size();
  int dimY = grid[0].size();
  int dimZ = grid[0][0].size();

  double up = grid[i][j][((k+1)%dimZ+dimZ)%dimZ];
  double down = grid[i][j][((k-1)%dimZ+dimZ)%dimZ];
  double left = grid[i][((j+1)%dimY+dimY)%dimY][k];
  double right = grid[i][((j-1)%dimY+dimY)%dimY][k];
  double X1 = grid[((i+1)%dimX+dimX)%dimX][j][k];
  double X2 = grid[((i-1)%dimX+dimX)%dimX][j][k];
  
  double lap = 0;
  lap += h_function(up) + h_function(down) + h_function(left) + h_function(right) + h_function(X1) + h_function(X2) - 6*h_function(grid[i][j][k]);
  
  return lap;	  
}

// CLASSES
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class Cell{
public:
  Cell(int dimX, int dimY, int dimZ, vector<double> center, vector<double> radius, vector<double> center_nucleus, vector<double> radius_nucleus, bool tan = false, bool bomboca = false);
  matrix grid;
  matrix grid_nucleus;
  double epsilon; 
  double alpha;     //volume conservation
  double alpha_nuc; 
  double alpha_s;   //surface conservation
  double eta;       //adhesion
  double gamma;     //repulsion
  double kappa;     //ration between bending rigidities
};

Cell::Cell (int dimX, int dimY, int dimZ, vector<double> center, vector<double> radius, vector<double> center_nucleus, vector<double> radius_nucleus, bool tan, bool bomboca) {
  grid = matrix(dimX, vector<vector<double>>(dimY, vector<double>(dimZ,0)));
  grid_nucleus = matrix(dimX, vector<vector<double>>(dimY, vector<double>(dimZ,0)));

  if (tan){
    for(int i=0;i<dimX;i++){
      for(int j=0;j<dimY;j++){
	for(int k=0;k<dimZ;k++){
	  
	  double dist = sqrt(
			     (i-center[0])*(i-center[0]) +
			     (j-center[1])*(j-center[1]) +
			     (k-center[2])*(k-center[2])
			     );
	  double theta = atan((k-center[2])/sqrt((i-center[0])*(i-center[0])+(j-center[1])*(j-center[1])));
	  double R =  sqrt(1/
			   (sin(theta)*sin(theta)/(radius[1]*radius[1]) + cos(theta)*cos(theta)/(radius[0]*radius[0]))
			   );
	  grid[i][j][k] = (-tanh(( dist-R )
				/epsilon
				)+1)/2;

	  
	  //nucleus	  
	  double dist_n = sqrt(
			     (i-center_nucleus[0])*(i-center_nucleus[0]) +
			     (j-center_nucleus[1])*(j-center_nucleus[1]) +
			     (k-center_nucleus[2])*(k-center_nucleus[2])
			     );
	  double theta_n = atan((k-center_nucleus[2])/sqrt((i-center_nucleus[0])*(i-center_nucleus[0])+(j-center_nucleus[1])*(j-center_nucleus[1])));
	  double R_n =  sqrt(1/
			   (sin(theta)*sin(theta)/(radius_nucleus[1]*radius_nucleus[1]) + cos(theta)*cos(theta)/(radius_nucleus[0]*radius_nucleus[0]))
			   );
	  grid_nucleus[i][j][k] = (radius_nucleus[0]==0)?0:(-tanh(( dist_n-R_n )
				/epsilon
				)+1)/2;

	}
      }
    }

    for(int i=center[0]-1;i<=center[0]+1;i++){
      for(int j=center[1]-1;j<=center[1]+1;j++){
	for(int k=center[2]-1;k<=center[2]+1;k++){
	  grid[i][j][k] = 1;
	}
      }
    }
    
    if (bomboca){
      for(int i=0;i<dimX;i++){
	for(int j=0;j<dimY;j++){
	  for(int k=0;k<center[2];k++){
	    grid[i][j][k] = 0;
	  }
	}
      } 
    }
    
  }
  else{
    # pragma omp parallel for
    for(int i=0;i<dimX;i++){
      for(int j=0;j<dimY;j++){
	for(int k=0;k<dimZ;k++){
	  if (((i-center[0])*(i-center[0])/(radius[0]*radius[0]) + (j-center[1])*(j-center[1])/(radius[0]*radius[0]) + (k-center[2])*(k-center[2])/(radius[1]*radius[1]) )<1)
	    grid[i][j][k]=1;
	  if (((i-center_nucleus[0])*(i-center_nucleus[0])/(radius_nucleus[0]*radius_nucleus[0]) + (j-center_nucleus[1])*(j-center_nucleus[1])/(radius_nucleus[0]*radius_nucleus[0]) + (k-center_nucleus[2])*(k-center_nucleus[2])/(radius_nucleus[1]*radius_nucleus[1]) )<1)
	    grid_nucleus[i][j][k]=1;
	}
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class Groove{
public:
  Groove(int dimX, int dimY, int dimZ, int depth, int spacing, int width, double epsilon, double alpha);
  matrix grid;
  double epsilon;
  double alpha;  
  
  void stabilize(double dt, double time, int Nframes, string s="stabilize") {

    int t0 = 0;
    int N = time/dt;
    int dimX = grid.size();
    int dimY = grid[0].size();
    int dimZ = grid[0][0].size();

    matrix grid_groove2 = grid;
    double eps = 1e-6;
    double vtarget_groove = vol(grid, eps);
    double cur_vol_groove = vtarget_groove;
    
    string id = s+to_string(t0)+".vti";  
    //saveGridToVTI("fiber_"+id,grid);
    
    for(int step=0;step<N;step++){
      # pragma omp parallel for
      for(int i=0;i<dimX;i++){
	for(int j=0;j<dimY;j++){
	  for(int k=0;k<dimZ;k++){

	    grid_groove2[i][j][k] = grid[i][j][k] + dt*(
						       0.5*laplacian(grid,i,j,k) +
						       grid[i][j][k]*(1-grid[i][j][k])*(grid[i][j][k]-0.5 + alpha*(vtarget_groove-cur_vol_groove))
						       );        
	  }
	}
      }
    
      id = "stabilize"+to_string(t0+step+1)+".vti";
      if (step%(int)(N/Nframes)==0) {
	//saveGridToVTI("fiber_"+id,grid);
      }
      swap(grid, grid_groove2);
    }
  }
};

//constructor
Groove::Groove (int dimX, int dimY, int dimZ, int depth, int spacing, int width, double epsilon, double alpha) {

  grid = matrix(dimX, vector<vector<double>>(dimY, vector<double>(dimZ)));

  //construct base
  # pragma omp parallel for
  for(int i=0;i<dimX;i++){
    for(int j=0;j<dimY;j++){
      for(int k=0;k<depth;k++){
	grid[i][j][k] = 1;
      }
    }
  }

  //construct grooves
  int i = width;
  if (depth!=0){
    while(i<dimX-width){
      if (i%(width+spacing)==0) i += width;
      # pragma omp parallel for
      for(int j=0;j<dimY;j++){
	for(int k=depth;k<2*depth;k++){
	  grid[i][j][k] = 1;
	}
      }
      i++;
    }
  }
  else{
    # pragma omp parallel for
    for(int i=0;i<dimX;i++){
      for(int j=0;j<dimY;j++){
	for(int k=0;k<width;k++){
	  grid[i][j][k] = 1;
	}
      }
    }
  }  
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void evolve(Cell &cell, matrix environment, double velocity, double velocity_nuc, double dt, double time, string dir="", bool stab = false, bool control_vol = false, int Nframes=100){
  
  int t0 = 1;
  int N = time/dt;
  int dimX = cell.grid.size();
  int dimY = cell.grid[0].size();
  int dimZ = cell.grid[0][0].size();

  matrix grid_cell2 = cell.grid;
  matrix grid_nucleus2 = cell.grid_nucleus;

  matrix grid_cell_aux = cell.grid;
  matrix grid_nucleus_aux = cell.grid_nucleus;

  double eps = 1e-3;
  double vtarget_cell = vol(cell.grid, eps);
  double cur_vol_cell = vtarget_cell;
  double o_vtarget_cell = vtarget_cell;
  double vtarget_nucleus = vol(cell.grid_nucleus, eps);
  double cur_vol_nucleus = vtarget_nucleus;

  double atarget_cell = area(cell.grid, cell.epsilon);
  double cur_area_cell = atarget_cell;
  double atarget_nucleus = area(cell.grid_nucleus, cell.epsilon);
  double cur_area_nucleus = atarget_nucleus;  
  
  string id = to_string(t0)+".vti";
  //saveGridToVTI(dir+"cell_"+id,cell.grid);
  //saveGridToVTI(dir+"environment_"+id,environment);
  //saveGridToVTI(dir+"nucleus_"+id,cell.grid_nucleus);

  double Mk = 1;

  ofstream volume_vs_time;
  ofstream area_vs_time;
  string l1 = "volume_vs_time";
  string l2 = "area_vs_time";

  ofstream volume_nucleus_vs_time;
  ofstream area_nucleus_vs_time;
  string l3 = "volume_nucleus_vs_time";
  string l4 = "area_nucleus_vs_time";
  if (stab){
    l1+="_stab.txt";
    l2+="_stab.txt";
    l3+="_stab.txt";
    l4+="_stab.txt";
  }
  else{
    l1+=".txt";
    l2+=".txt";
    l3+=".txt";
    l4+=".txt";
  }
  volume_vs_time.open(l1);
  area_vs_time.open(l2);
  volume_nucleus_vs_time.open(l3);
  area_nucleus_vs_time.open(l4);
  
  for(int step=0;step<N;step++){  
    
    cur_vol_cell = vol(cell.grid, eps);
    // squeeze cell
    if (control_vol && step%100==0) vtarget_cell = o_vtarget_cell/(1+double(step)/double(N));
    cur_vol_nucleus = vol(cell.grid_nucleus, eps);

    cur_area_cell = area(cell.grid, cell.epsilon);
    cur_area_nucleus = area(cell.grid_nucleus, cell.epsilon);

    volume_vs_time << cur_vol_cell << endl;
    area_vs_time << cur_area_cell << endl;
    
    volume_nucleus_vs_time << cur_vol_nucleus << endl;
    area_nucleus_vs_time << cur_area_nucleus << endl;
      
    if (step%100==0) {
       cout<<"Step "<<step<<"/"<<N<<endl;
       cout<<"Volume cell "<<cur_vol_cell<<"/"<<vtarget_cell<<" with error "<<abs(vtarget_cell-cur_vol_cell)/vtarget_cell<<endl;
       cout<<"Volume nucleus "<<cur_vol_nucleus<<"/"<<vtarget_nucleus<<" with error of "<<abs(vtarget_nucleus-cur_vol_nucleus)/vtarget_nucleus<<endl;
    }

    # pragma omp parallel for
    for(int i=0;i<dimX;i++){
      for(int j=0;j<dimY;j++){
	for(int k=0;k<dimZ;k++){
	  grid_cell_aux[i][j][k] = cell.grid[i][j][k]*(cell.grid[i][j][k]-1)*(2*cell.grid[i][j][k]-1)-cell.epsilon*cell.epsilon*laplacian(cell.grid,i,j,k);
	  grid_nucleus_aux[i][j][k] = cell.grid_nucleus[i][j][k]*(cell.grid_nucleus[i][j][k]-1)*(2*cell.grid_nucleus[i][j][k]-1)-cell.epsilon*cell.epsilon*laplacian(cell.grid_nucleus,i,j,k);
	}
      }
    } 
    
    # pragma omp parallel for
    for(int i=0;i<dimX;i++){
      for(int j=0;j<dimY;j++){
	for(int k=0;k<dimZ;k++){ 
	  
	  double h_environment = environment[i][j][k]*environment[i][j][k]*(3-2*environment[i][j][k]); //interaction cell-environment
	  double h_cell = (1-cell.grid[i][j][k])*(1-cell.grid[i][j][k])*(3-2*(1-cell.grid[i][j][k])); //repulsion of the nucleus at the border
	  double h_nucleus = cell.grid_nucleus[i][j][k]*cell.grid_nucleus[i][j][k]*(3-2*cell.grid_nucleus[i][j][k]);

	  grid_cell2[i][j][k] = cell.grid[i][j][k] + dt*(
							 -velocity*(cell.grid[i][j][((k+1)%dimZ+dimZ)%dimZ] - cell.grid[i][j][((k-1)%dimZ+dimZ)%dimZ])*0.5 - //advective term
							 Mk*(
							    2*(1+6*cell.grid[i][j][k]*(cell.grid[i][j][k]-1))*grid_cell_aux[i][j][k] -
							    2*cell.epsilon*cell.epsilon*laplacian(grid_cell_aux, i, j, k) -
							    6*cell.eta*cell.grid[i][j][k]*(1-cell.grid[i][j][k])*laplacian_h(environment,i,j,k) +
							    6*cell.gamma*cell.grid[i][j][k]*(1-cell.grid[i][j][k])*(h_environment-h_nucleus)
							     )
        );
	  
	  grid_nucleus2[i][j][k] = cell.grid_nucleus[i][j][k] + dt*(
							 -velocity_nuc*(cell.grid_nucleus[i][j][((k+1)%dimZ+dimZ)%dimZ] - cell.grid_nucleus[i][j][((k-1)%dimZ+dimZ)%dimZ])*0.5 - //advective term
							 Mk*(
							     2*cell.kappa*(1+6*cell.grid_nucleus[i][j][k]*(cell.grid_nucleus[i][j][k]-1))*grid_nucleus_aux[i][j][k] -
									   2*cell.kappa*cell.epsilon*cell.epsilon*laplacian(grid_nucleus_aux, i, j, k)
							     + 6*cell.gamma*cell.grid_nucleus[i][j][k]*(1-cell.grid_nucleus[i][j][k])*h_cell
							     )
        );
	  
	  if (!stab)
	    {
	    grid_cell2[i][j][k] -= Mk*dt*(
				       24*cell.alpha_s*cell.epsilon*(atarget_cell - cur_area_cell)*laplacian(cell.grid, i, j, k) -
				       12*cell.alpha*cell.grid[i][j][k]*(1-cell.grid[i][j][k])*(vtarget_cell-cur_vol_cell)
				       );
	    
	    grid_nucleus2[i][j][k] -= Mk*dt*(
					  24*cell.alpha_s*cell.epsilon*(atarget_nucleus - cur_area_nucleus)*laplacian(cell.grid_nucleus, i, j, k) -
					  12*cell.alpha_nuc*cell.grid_nucleus[i][j][k]*(1-cell.grid_nucleus[i][j][k])*(vtarget_nucleus-cur_vol_nucleus)
					  );
	    }
	}
      }
    }
  
    id = to_string(t0+step+1)+".vti";
    //printf("\r[%3d%%]",double(step)/N*100);
    if (step%Nframes==0) {
      saveGridToVTI(dir+"cell_"+id,cell.grid);
      //saveGridToVTI(dir+"environment_"+id,environment);
      saveGridToVTI(dir+"nucleus_"+id,cell.grid_nucleus);
    }
    swap(cell.grid, grid_cell2);
    swap(cell.grid_nucleus, grid_nucleus2);
  }
  volume_vs_time.close();
  area_vs_time.close();
  volume_nucleus_vs_time.close();
  area_nucleus_vs_time.close();
  t0 +=N;
}



