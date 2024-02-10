using namespace std;
using matrix = vector<vector<vector<double>>>; 

inline double laplacian(matrix &grid, int i, int j, int k){
  int dimX = grid.size();
  int dimY = grid[0].size();
  int dimZ = grid[0][0].size();

  double lap = 0;
  lap += grid[((i+1)%dimX+dimX)%dimX][j][k] + grid[((i-1)%dimX+dimX)%dimX][j][k] +
         grid[i][((j+1)%dimY+dimY)%dimY][k] + grid[i][((j-1)%dimY+dimY)%dimY][k] +
         grid[i][j][((k+1)%dimZ+dimZ)%dimZ] + grid[i][j][((k-1)%dimZ+dimZ)%dimZ]
         - 6*grid[i][j][k];
  
  return lap;	  
}

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

class Cell{
public:
  Cell(int dimX, int dimY, int dimZ, vector<double> center, vector<double> radius, vector<double> center_nucleus, double radius_nucleus);
  matrix grid;
  matrix grid_nucleus;
  double epsilon; //diffusion coefficient
  double alpha;   //volume coefficient
  double eta; //adhesion coefficient
  double gamma; //repulsion term
  double kappa;
};

Cell::Cell (int dimX, int dimY, int dimZ, vector<double> center, vector<double> radius, vector<double> center_nucleus, double radius_nucleus) {
  grid = matrix(dimX, vector<vector<double>>(dimY, vector<double>(dimZ)));
  grid_nucleus = matrix(dimX, vector<vector<double>>(dimY, vector<double>(dimZ)));
  # pragma omp parallel for
  for(int i=0;i<dimX;i++){
    for(int j=0;j<dimY;j++){
      for(int k=0;k<dimZ;k++){
	if (((i-center[0])*(i-center[0])/(radius[0]*radius[0]) + (j-center[1])*(j-center[1])/(radius[0]*radius[0]) + (k-center[2])*(k-center[2])/(radius[1]*radius[1]) )<1)
	  grid[i][j][k]=1;
	if (((i-center_nucleus[0])*(i-center_nucleus[0]) + (j-center_nucleus[1])*(j-center_nucleus[1]) + (k-center_nucleus[2])*(k-center_nucleus[2]) )<radius_nucleus*radius_nucleus)
	  grid_nucleus[i][j][k]=1;
      }
    }
  }
}


//-----------

void make_random_fiber(matrix &grid, int dimX, int dimY, int dimZ, double mu_theta, double mu_phi){
  
  //random points  
  random_device rd; //random number from uniform distribuition as a seed
  mt19937 gen(rd()); //mersenne_twister_engine
  uniform_real_distribution<> disx(0, dimX-1);
  uniform_real_distribution<> disy(0, dimY-1);
  uniform_real_distribution<> disz(0, dimZ-1);
  vector<double> line_point = {disx(gen), disy(gen), disz(gen)};

  //random theta
  random_device rd2{};
  mt19937 gen2{rd2()};
  double sigma_theta=0.3;//acos(-1)/6;
  normal_distribution d_theta{mu_theta, sigma_theta};
  double theta = d_theta(gen2);

  //random phi
  random_device rd4{};
  mt19937 gen4{rd4()};
  double sigma_phi=0.6;//acos(-1)/6;
  normal_distribution d_phi{mu_phi, sigma_phi};
  double phi = d_phi(gen4);
  
  //random point
  random_device rd3{};
  mt19937 gen3{rd3()};
  double mu_radius=3;
  double sigma_radius=0.5;
  normal_distribution d_radius{mu_radius, sigma_radius};
  double radius = d_radius(gen3);
    
  //create
  vector<double> versor = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
  # pragma omp parallel for
  for(int i=0;i<dimX;i++){
    for(int j=0;j<dimY;j++){
      for(int k=0;k<dimZ;k++){
	vector<double> vec = {i-line_point[0], j-line_point[1], k-line_point[2]};
	vector<double> vec2 = cross(vec, versor);
	double dist = norm(vec2)/norm(versor);
	if (dist<radius) grid[i][j][k]=1;
      }
    }
  }
}


class Fiber{
public:
  Fiber(int dimX, int dimY, int dimZ, double mu_theta, double mu_phi, double rho, double epsilon, double alpha);
  matrix grid;
  double epsilon; //diffusion coefficient
  double alpha;   //volume coefficient 
  
  void stabilize(double dt, double time, int Nframes, string s="stabilize") {

    int t0 = 0;
    int N = time/dt;
    int dimX = grid.size();
    int dimY = grid[0].size();
    int dimZ = grid[0][0].size();

    matrix grid_fiber2 = grid;
    double eps = 1e-7;
    double vtarget_fiber = vol(grid, eps);
    double cur_vol_fiber = vtarget_fiber;
    
    string id = s+to_string(t0)+".vti";  
    //saveGridToVTI("fiber_"+id,grid);
    
    for(int step=0;step<N;step++){
      # pragma omp parallel for
      for(int i=0;i<dimX;i++){
	for(int j=0;j<dimY;j++){
	  for(int k=0;k<dimZ;k++){

	    grid_fiber2[i][j][k] = grid[i][j][k] + dt*(
						       laplacian(grid,i,j,k) +
						       grid[i][j][k]*(1-grid[i][j][k])*(grid[i][j][k]-0.5 + alpha*(vtarget_fiber-cur_vol_fiber))
						       );        
	  }
	}
      }
    
      id = "stabilize"+to_string(t0+step+1)+".vti";
      if (step%(int)(N/Nframes)==0) {
	//saveGridToVTI("fiber_"+id,grid);
      }

      swap(grid, grid_fiber2);
    }
  }

  void clear(matrix &grid_cell){
    double eps=1e-7;
    # pragma omp parallel for
    for(int i=0;i<grid_cell.size();i++){
      for(int j=0;j<grid_cell[0].size();j++){
	for(int k=0;k<grid_cell[0][0].size();k++){
	  if (eps<grid[i][j][k] && eps<grid_cell[i][j][k]) grid[i][j][k]=0;
	}
      }
    }
  }
};

//constructor
Fiber::Fiber (int dimX, int dimY, int dimZ, double mu_theta, double mu_phi, double rho, double epsilon, double alpha) {

  grid = matrix(dimX, vector<vector<double>>(dimY, vector<double>(dimZ)));

  double curRho=0;
  while(curRho<rho){ //creates fibers until target density is met
    make_random_fiber(grid, dimX, dimY, dimZ, mu_theta, mu_phi);
    curRho=density(grid);
  }
}


class Groove{
public:
  Groove(int dimX, int dimY, int dimZ, int depth, int spacing, int width, double epsilon, double alpha);
  matrix grid;
  double epsilon; //diffusion coefficient
  double alpha;   //volume coefficient 
  
  void stabilize(double dt, double time, int Nframes, string s="stabilize") {

    int t0 = 0;
    int N = time/dt;
    int dimX = grid.size();
    int dimY = grid[0].size();
    int dimZ = grid[0][0].size();

    matrix grid_groove2 = grid;
    double eps = 1e-7;
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


//-----------


void evolve(Cell &cell, matrix environment, double velocity, double velocity_nuc, double dt, int time, string dir="", int Nframes=100){
  
  int t0 = 1;
  int N = time/dt;
  int dimX = cell.grid.size();
  int dimY = cell.grid[0].size();
  int dimZ = cell.grid[0][0].size();

  matrix grid_cell2 = cell.grid;
  matrix grid_nucleus2 = cell.grid_nucleus;

  double eps = 1e-3;
  double vtarget_cell = vol(cell.grid, eps);
  //double vtarget_cell = vol_citoplasm(cell.grid, cell.grid_nucleus, eps);

  double cur_vol_cell = vtarget_cell;
  double vtarget_nucleus = vol(cell.grid_nucleus, eps);
  double cur_vol_nucleus = vtarget_nucleus;
  
  string id = to_string(t0)+".vti";
  //saveGridToVTI(dir+"cell_"+id,cell.grid);
  //saveGridToVTI(dir+"environment_"+id,environment);
  //saveGridToVTI(dir+"nucleus_"+id,cell.grid_nucleus);
    
  for(int step=0;step<N;step++){

    cur_vol_cell = vol(cell.grid, eps);
    cur_vol_nucleus = vol(cell.grid_nucleus, eps);
    
    // if (step%100) {
    //   cout<<"Step "<<step<<"/"<<N<<endl;
    //   cout<<"Volume cell "<<cur_vol_cell<<"/"<<vtarget_cell<<endl;
    //   cout<<"Volume nucleus "<<cur_vol_nucleus<<"/"<<vtarget_nucleus<<endl;
    // }
    
    # pragma omp parallel for
    for(int i=0;i<dimX;i++){
      for(int j=0;j<dimY;j++){
	for(int k=0;k<dimZ;k++){ 
	  
	  double h_environment = environment[i][j][k]*environment[i][j][k]*(3-2*environment[i][j][k]); //interaction cell-environment
	  double h_cell = (1-cell.grid[i][j][k])*(1-cell.grid[i][j][k])*(3-2*(1-cell.grid[i][j][k])); //repulsion of the nucleus at the border
	  double h_nucleus = cell.grid_nucleus[i][j][k]*cell.grid_nucleus[i][j][k]*(3-2*cell.grid_nucleus[i][j][k]);

	  grid_cell2[i][j][k] = cell.grid[i][j][k] + dt*(
							 -velocity*(cell.grid[i][j][((k+1)%dimZ+dimZ)%dimZ] - cell.grid[i][j][((k-1)%dimZ+dimZ)%dimZ])*0.5 + //advective term
							 cell.epsilon*cell.epsilon*laplacian(cell.grid,i,j,k) + //surface 
							 cell.grid[i][j][k]*(1-cell.grid[i][j][k])*(cell.grid[i][j][k]-0.5
												    -6*cell.gamma*h_environment //repulsion environment
												    +6*cell.gamma*h_nucleus //repulsion nucleus 
												    +6*cell.eta*laplacian_h(environment,i,j,k) //adhesion
												    +12*cell.alpha*(vtarget_cell-cur_vol_cell) //volume conservation
												    )
							 );
	  
	  grid_nucleus2[i][j][k] = cell.grid_nucleus[i][j][k] + dt*(
								    -velocity_nuc*(cell.grid[i][j][((k+1)%dimZ+dimZ)%dimZ] - cell.grid[i][j][((k-1)%dimZ+dimZ)%dimZ])*0.5 + //advective term
								    cell.kappa*cell.epsilon*cell.epsilon*laplacian(cell.grid_nucleus,i,j,k) + //surface
								    cell.grid_nucleus[i][j][k]*(1-cell.grid_nucleus[i][j][k])*( cell.kappa*(cell.grid_nucleus[i][j][k]-0.5)
																-6*cell.gamma*h_cell //repulsion from outside of cell
																+12*cell.alpha*(vtarget_nucleus-cur_vol_nucleus) //volume conservation
																)
								    );
	}
      }
    }
  
    id = to_string(t0+step+1)+".vti";
    if (step%Nframes==0) {
      cout<<"On step "<<step<<"/"<<N<<endl;
      //saveGridToVTI(dir+"cell_"+id,cell.grid);
      //saveGridToVTI(dir+"environment_"+id,environment);
      //saveGridToVTI(dir+"nucleus_"+id,cell.grid_nucleus);
    }
    swap(cell.grid, grid_cell2);
    swap(cell.grid_nucleus, grid_nucleus2);
  }
  t0 +=N;
}



