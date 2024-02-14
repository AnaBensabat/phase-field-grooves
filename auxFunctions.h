using namespace std;
using matrix = vector<vector<vector<double>>>; 

// MATH
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

double norm(vector<double> &vec){
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

vector<double> cross(vector<double> &u, vector<double> &v){
  return {u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-v[0]*u[1]};
}

// PHYSICAL 
double density(matrix &fiber_grid){
  double rho=0;
  for(int i=0;i<fiber_grid.size();i++){
    for(int j=0;j<fiber_grid[0].size();j++){
      for(int k=0;k<fiber_grid[0][0].size();k++){
	rho+=fiber_grid[i][j][k];
      }
    }
  }
  return rho/(fiber_grid.size()*fiber_grid[0].size()*fiber_grid[0][0].size());
}

double vol_in_grooves(matrix &grid, int height, float eps){
  double volume=0;
  int dimx = grid.size();
  int dimy = grid[0].size();
  int dimz = grid[0][0].size();
  for(int i=0;i<dimx;i++){
    for(int j=0;j<dimy;j++){
      for(int k=0;k<=height;k++){
	double h = grid[i][j][k]*grid[i][j][k]*(3-2*grid[i][j][k]);
	volume += (h > eps) ? h : 0;
      }
    }
  }
  return volume;
}

double vol(matrix &grid, float eps){

  double volume=0;
  int dimx = grid.size();
  int dimy = grid[0].size();
  int dimz = grid[0][0].size();
  for(int i=0;i<dimx;i++){
    for(int j=0;j<dimy;j++){
      for(int k=0;k<dimz;k++){
	double h_phi =  grid[i][j][k]*grid[i][j][k]*(3-2*grid[i][j][k]);
	volume += (h_phi > eps) ? h_phi : 0;
      }
    }
  }
  return volume;
}

double area(matrix &grid, float epsilon){
  
  int dimX = grid.size();
  int dimY = grid[0].size();
  int dimZ = grid[0][0].size();

  double A = 0;
  for(int i=0;i<dimX;i++){
    for(int j=0;j<dimY;j++){
      for(int k=0;k<dimZ;k++){
        double dx = (grid[((i+1)%dimX+dimX)%dimX][j][k] - grid[((i-1)%dimX+dimX)%dimX][j][k])*0.5;
        double dy = (grid[i][((j+1)%dimY+dimY)%dimY][k] - grid[i][((j-1)%dimY+dimY)%dimY][k])*0.5;
        double dz = (grid[i][j][((k+1)%dimZ+dimZ)%dimZ] - grid[i][j][((k-1)%dimZ+dimZ)%dimZ])*0.5;
	A += dx*dx + dy*dy + dz*dz;
      }
    }
  }
  
  return 6*epsilon*A;
}

vector<double> CM(matrix &cell_grid){
  double mass=0;
  vector<double> cm(3);
  for(int i=0;i<cell_grid.size();i++){
    for(int j=0;j<cell_grid[0].size();j++){
      for(int k=0;k<cell_grid[0][0].size();k++){
        mass+=cell_grid[i][j][k];
	cm[0]+=i*cell_grid[i][j][k];
	cm[1]+=j*cell_grid[i][j][k];
	cm[2]+=k*cell_grid[i][j][k];
      }
    }
  }
  cm[0]/=mass;
  cm[1]/=mass;
  cm[2]/=mass;
  
  return cm;
}

double lowest(matrix &grid, int depth){
  double bottom=grid[0][0].size();
  for(int i=0;i<grid.size();i++){
    for(int j=0;j<grid[0].size();j++){
      for(double k=0;k<depth*2;k++){
	if (grid[i][j][int(k)]>1e-3) bottom=min(bottom,k);
      }
    }
  }
  return bottom;
}


vector<double> DFS(matrix &original_grid, int height, double eps){
  matrix grid = original_grid;
  int dimX = grid.size();
  int dimY = grid[0].size();
  int dimZ = grid[0][0].size();
  # pragma omp parallel for
  for(int i=0;i<dimX;i++){
    for(int j=0;j<dimY;j++){
      for(int k=height;k<dimZ;k++){
	grid[i][j][k] = 0;
      }
    }
  }

  stack<array<int,3>> stack;
  vector<double> components;
  const vector<array<int,3>> deltas {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
  for(int i=0;i<dimX;i++){
    for(int j=0;j<dimY;j++){
      for(int k=0;k<height;k++){
	if (grid[i][j][k] > eps){
	  stack.push({i,j,k});
	  double volume = 0;
	  while(stack.size()>0){
	    //cout<<stack.size()<<'\n';
	    auto [x,y,z]  = stack.top();
	    stack.pop();
	    volume += grid[x][y][z];
	    grid[x][y][z] = 0; //mark as visited
	    for(auto[dx,dy,dz]:deltas){
	      if (x+dx<0 || x+dx>=dimX || y+dy<0 || y+dy>=dimY || z+dz<0 || z+dz>=dimZ) continue;
	      if(grid[x+dx][y+dy][z+dz]>eps) stack.push({x+dx,y+dy,z+dz});
	    }
	  }
	  components.push_back(volume);
	}
      }
    }
  }
  return components;
}

bool isSplit(matrix &original_grid, double eps){
  
  matrix grid = original_grid;
  int dimX = grid.size();
  int dimY = grid[0].size();
  int dimZ = grid[0][0].size();

  stack<array<int,3>> stack;
  vector<double> components;
  const vector<array<int,3>> deltas {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
  for(int i=0;i<dimX;i++){
    for(int j=0;j<dimY;j++){
      for(int k=0;k<dimZ;k++){
	if (grid[i][j][k] > eps){
	  stack.push({i,j,k});
	  double volume = 1;
	  while(stack.size()>0){
	    //cout<<stack.size()<<'\n';
	    auto [x,y,z]  = stack.top();
	    stack.pop();
	    grid[x][y][z] = 0; //mark as visited}
	    for(auto[dx,dy,dz]:deltas){
	      if (x+dx<0 || x+dx>=dimX || y+dy<0 || y+dy>=dimY || z+dz<0 || z+dz>=dimZ) continue;
	      if(grid[x+dx][y+dy][z+dz]>eps) stack.push({x+dx,y+dy,z+dz});
	    }
	  }
	  components.push_back(volume);
	}
      }
    }
  }
  if (components.size()<=1) return false;
  else return true;
}



// MISC
void clearDirectory(const std::filesystem::path& dir)
{
    for (const auto& entry : std::filesystem::directory_iterator(dir)) 
        std::filesystem::remove_all(entry.path());
}

void saveGridToVTI(string filename, matrix &grid) {

  int dimX = grid.size();
  int dimY = grid[0].size();
  int dimZ = grid[0][0].size();
  FILE* file = fopen(filename.c_str(), "w");
  if (file == NULL) {
    printf("Error opening file %s\n", filename);
    return;
  }
 
  // Write VTI header
  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
  fprintf(file, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\">\n", dimX - 1, dimY - 1, dimZ - 1);
  fprintf(file, "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n", dimX - 1, dimY - 1, dimZ - 1);
  fprintf(file, "      <PointData>\n");
  fprintf(file, "        <DataArray type=\"Float64\" Name=\"Grid\" NumberOfComponents=\"1\" format=\"ascii\">\n");

  // Write grid values
  for (int k = 0; k < dimZ; k++) {
    for (int j = 0; j < dimY; j++) {
      for (int i = 0; i < dimX; i++) {
	fprintf(file, "%.4f ", grid[i][j][k]);
      }
      //fprintf(file, "\n");
    }
    //fprintf(file, "\n");
  }
  
  // Write VTI footer
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </PointData>\n");
  fprintf(file, "    </Piece>\n");
  fprintf(file, "  </ImageData>\n");
  fprintf(file, "</VTKFile>\n");
 
  fclose(file);
}

matrix loadGridFromVTI(string filename){

  ifstream file(filename, ios::binary);
  if (!file.is_open()) {
    cout<<"Failed to open the VTI file."<<endl;
    exit(1);
  }

  string line;
  int dimX, dimY, dimZ;
  while (getline(file, line)) {
    if (line.find("<ImageData WholeExtent") != std::string::npos){
      sscanf(line.c_str(), "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\">\n", &dimX, &dimY, &dimZ);
      break;
    }
  }
  dimX++;
  dimY++;
  dimZ++;

  matrix grid(dimX, vector<vector<double>>(dimY, vector<double>(dimZ)));
  
  // Read the binary data
  bool data=false;
  while (std::getline(file, line)) {    
    if (data){      
      istringstream lineStream(line);
      for (int z = 0; z < dimZ; z++) {
	for (int y = 0; y < dimY; y++) {
	  for (int x = 0; x < dimX; x++) {
	    double value;
	    lineStream >> value;
	    grid[x][y][z] = value;
	  }
	}
      }
      data=false;
    }
    
    if (line.find("<DataArray") != std::string::npos) {
      data=true;
    }
  }
  
  return grid;
}


