#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <pthread.h>
#include "Constants.cpp"

typedef struct {
  float phi[2][N_x][N_y][N_z];
  float vel[2][N_x][N_y][N_z]; /* Phi is complex and we track the 
			    real and imaginary parts separately. */
  float E[dim][N_x][N_y][N_z];
  float A[dim][N_x][N_y][N_z]; /* The link variables are complex, but 
				 there is only one degree of freedom. */
} Coords; 

typedef struct {
  int thread_no;
} ThreadData;

unsigned long int step = 0;
float Energy_vals[frames_total][E_params];
float temp_Energy[E_params];
Coords* const T1 = new Coords;
Coords* const T2 = new Coords;
pthread_barrier_t barr;

inline void periodic_boundaries(int ix, int iy, int iz, int& corr_ix, int& corr_iy, int& corr_iz){
  (ix < 0) ? (corr_ix = N_x+ix):((ix > N_x-1) ? (corr_ix = ix - N_x):corr_ix = ix);
  (iy < 0) ? (corr_iy = N_y+iy):((iy > N_y-1) ? (corr_iy = iy - N_y):corr_iy = iy);
  (iz < 0) ? (corr_iz = N_z+iz):((iz > N_z-1) ? (corr_iz = iz - N_z):corr_iz = iz);
}

void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " "){
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos){
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
}

inline void read_data(Coords& temp_pos){
  int ix, iy, iz, counter;
  std::string infile("initial_data.dat");
  std::string c_line;
  std::vector<std::string> line_elem;
  std::ifstream fin;
  fin.open(infile.c_str());
  if (fin.is_open()){
    counter = 0;
    while(fin){
      getline(fin, c_line);
      Tokenize(c_line, line_elem);
      if(dim + scal_dim <= line_elem.size()){
	ix = atoi(line_elem[0].c_str());
	iy = atoi(line_elem[1].c_str());
	iz = atoi(line_elem[2].c_str());
	for(int pos = dim; pos < dim + scal_dim; ++pos){
	  temp_pos.phi[pos - dim][ix][iy][iz] = scale_factor*atof(line_elem[pos].c_str());
	}
	++counter;
      }
      line_elem.clear();
    }
  }
  else{
    std::cout << " Error: could not open " << infile << std::endl;
    exit(-1);
  }
  fin.close();
  if(counter < pow(N_x, 3)){
    std::cerr << "Not enough data points in: " << infile << std::endl;
    exit(-1);
  }
}

void syncT1T2MT(void* thread_data){
  ThreadData* local_thread_data =  static_cast<ThreadData*>(thread_data);
  int low_x_Lim   = static_cast<int>((N_x*local_thread_data->thread_no)/NoOfThreads);
  int high_x_Lim  = static_cast<int>((N_x*(local_thread_data->thread_no + 1))/NoOfThreads);
  for(int i = low_x_Lim; i < high_x_Lim; ++i){
    for(int j = 0; j < N_y; ++j){
      for(int k = 0; k < N_z; ++k){
	for(int p = 0; p < 2; ++p){
	  T1->phi[p][i][j][k] = T2->phi[p][i][j][k];
	  T1->vel[p][i][j][k] = T2->vel[p][i][j][k];
	}
	for(int p = 0; p < dim; ++p){
	  T1->A[p][i][j][k] = T2->A[p][i][j][k];
	  T1->E[p][i][j][k] = T2->E[p][i][j][k];
	} 
      }
    }
  }
#ifdef DEBUG
  std::cout << "syncT1T2MT: low_x_Lim = " << low_x_Lim << " high_x_Lim " << high_x_Lim << std::endl;
#endif
}

inline float PlaqX1(Coords& Curr_pos, int ix, int iy, int iz){

  int c_ix_pxy_1, c_ix_pxy_2;
  int c_iy_pxy_1, c_iy_pxy_2; 
  int c_iz_pxy_1, c_iz_pxy_2; /* Coordinated of 2 corners 
				 of the plaquette in the xy plane. */
  float tmp_plaq_1x;
  /* Link along x direction with origin at [ix][iy][iz]: A[0][ix][iy][iz]. 
     The 4 plaquettes that contribute are 
     [ix][iy  ][iz] -> [ix+1][iy][iz] -> [ix+1][iy+1][iz] -> [ix  ][iy+1][iz] -> [ix][iy  ][iz] 
     [ix][iy-1][iz] -> [ix  ][iy][iz] -> [ix+1][iy  ][iz] -> [ix+1][iy-1][iz] -> [ix][iy-1][iz] 
     [ix][iy][iz  ] -> [ix+1][iy][iz] -> [ix+1][iy][iz+1] -> [ix  ][iy][iz+1] -> [ix][iy][iz  ]
     [ix][iy][iz-1] -> [ix  ][iy][iz] -> [ix+1][iy][iz  ] -> [ix+1][iy][iz-1] -> [ix][iy][iz-1]
  */
  
  /* First Plaquette xy plane. */

  /* [ix+1][iy  ][iz] -> [ix+1][iy+1][iz] */
  periodic_boundaries(ix+1, iy  , iz, c_ix_pxy_1, c_iy_pxy_1, c_iz_pxy_1); 
  /* [ix+1][iy+1][iz] -> [ix  ][iy+1][iz] */
  periodic_boundaries(ix  , iy+1, iz, c_ix_pxy_2, c_iy_pxy_2, c_iz_pxy_2);     
  
  tmp_plaq_1x = sin(a*(Curr_pos.A[0][ix][iy][iz]+
		       Curr_pos.A[1][c_ix_pxy_1][c_iy_pxy_1][c_iz_pxy_1]-
		       Curr_pos.A[0][c_ix_pxy_2][c_iy_pxy_2][c_iz_pxy_2]-
		       Curr_pos.A[1][ix][iy][iz]))/pow(a, 3);
  
  return(tmp_plaq_1x);
}
/***************************************************************************************/
inline float PlaqX2(Coords& Curr_pos, int ix, int iy, int iz){
  int c_ix_pxy_1, c_ix_pxy_2, c_ix_pxy_3;
  int c_iy_pxy_1, c_iy_pxy_2, c_iy_pxy_3; 
  int c_iz_pxy_1, c_iz_pxy_2, c_iz_pxy_3; /* Coordinated of 2 corners 
					     of the plaquette in the xy plane. */
  float tmp_plaq_2x;
  /* Second Plaquette xy plane*/
  
  /* [ix][iy - 1][iz] -> [ix][iy][iz] */
  periodic_boundaries(ix  , iy-1, iz, c_ix_pxy_1, c_iy_pxy_1, c_iz_pxy_1); 
  /* [ix + 1][iy][iz] -> [ix + 1][iy - 1][iz] */
  periodic_boundaries(ix+1, iy-1, iz, c_ix_pxy_2, c_iy_pxy_2, c_iz_pxy_2);  
  /* [ix + 1][iy - 1][iz] -> [ix][iy - 1][iz] */
  periodic_boundaries(ix  , iy-1, iz, c_ix_pxy_3, c_iy_pxy_3, c_iz_pxy_3);  
  
  tmp_plaq_2x = sin(a*(Curr_pos.A[1][c_ix_pxy_1][c_iy_pxy_1][c_iz_pxy_1]+
		       Curr_pos.A[0][ix][iy][iz]-
		       Curr_pos.A[1][c_ix_pxy_2][c_iy_pxy_2][c_iz_pxy_2]-
		       Curr_pos.A[0][c_ix_pxy_3][c_iy_pxy_3][c_iz_pxy_3]))/pow(a, 3);
  
  return(tmp_plaq_2x);
}
/***************************************************************************************/
inline float PlaqX3(Coords& Curr_pos, int ix, int iy, int iz){

  int c_ix_pxz_1, c_ix_pxz_2;
  int c_iy_pxz_1, c_iy_pxz_2;
  int c_iz_pxz_1, c_iz_pxz_2; /* Coordinated of 2 corners 
				 of the plaquette in the xz plane. */
  float tmp_plaq_3x;
  /* Third Plaquette xz plane  */
  
  /* [ix+1][iy][iz  ] -> [ix+1][iy][iz+1] */
  periodic_boundaries(ix+1, iy  , iz  , c_ix_pxz_1, c_iy_pxz_1, c_iz_pxz_1); 
  /* [ix+1][iy][iz+1] -> [ix  ][iy][iz+1] */
  periodic_boundaries(ix  , iy  , iz+1, c_ix_pxz_2, c_iy_pxz_2, c_iz_pxz_2);     
  
  tmp_plaq_3x = sin(a*(Curr_pos.A[0][ix][iy][iz]+
		       Curr_pos.A[2][c_ix_pxz_1][c_iy_pxz_1][c_iz_pxz_1]-
		       Curr_pos.A[0][c_ix_pxz_2][c_iy_pxz_2][c_iz_pxz_2]-
		       Curr_pos.A[2][ix][iy][iz]))/pow(a, 3);
  
  return(tmp_plaq_3x);
}
/***************************************************************************************/
inline float PlaqX4(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_pxz_1, c_ix_pxz_2, c_ix_pxz_3;
  int c_iy_pxz_1, c_iy_pxz_2, c_iy_pxz_3;
  int c_iz_pxz_1, c_iz_pxz_2, c_iz_pxz_3; /* Coordinated of 2 corners 
					     of the plaquette in the xz plane. */
  float tmp_plaq_4x;
  /* Fourth Plaquette xz plane */
  
  /* [ix][iy][iz-1] -> [ix  ][iy][iz] */
  periodic_boundaries(ix  , iy  , iz-1, c_ix_pxz_1, c_iy_pxz_1, c_iz_pxz_1); 
  /* [ix+1][iy][iz  ] -> [ix+1][iy][iz-1] */
  periodic_boundaries(ix+1, iy  , iz-1, c_ix_pxz_2, c_iy_pxz_2, c_iz_pxz_2);  
  /* [ix+1][iy][iz-1] -> [ix][iy][iz-1] */
  periodic_boundaries(ix  , iy  , iz-1, c_ix_pxz_3, c_iy_pxz_3, c_iz_pxz_3);  
  
  tmp_plaq_4x = sin(a*(Curr_pos.A[2][c_ix_pxz_1][c_iy_pxz_1][c_iz_pxz_1]+
		       Curr_pos.A[0][ix][iy][iz]-
		       Curr_pos.A[2][c_ix_pxz_2][c_iy_pxz_2][c_iz_pxz_2]-
		       Curr_pos.A[0][c_ix_pxz_3][c_iy_pxz_3][c_iz_pxz_3]))/pow(a, 3);
  
  return(tmp_plaq_4x);
}
/***************************************************************************************/


inline float PlaqY1(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_pxy_1, c_ix_pxy_2;
  int c_iy_pxy_1, c_iy_pxy_2; 
  int c_iz_pxy_1, c_iz_pxy_2; /* Coordinated of 2 corners 
				 of the plaquette in the xy plane. */
  float tmp_plaq_1y;
  /* Link along y direction whit origin at [ix][iy][iz]: A[1][ix][iy][iz]. 
     The four plaquettes that contribute are 
     [ix  ][iy][iz] -> [ix][iy+1][iz] -> [ix+1][iy+1][iz] -> [ix+1][iy  ][iz] -> [ix  ][iy][iz] 
     [ix-1][iy][iz] -> [ix][iy  ][iz] -> [ix  ][iy+1][iz] -> [ix-1][iy+1][iz] -> [ix-1][iy][iz]  
     [ix][iy][iz  ] -> [ix][iy+1][iz] -> [ix][iy+1][iz+1] -> [ix][iy  ][iz+1] -> [ix][iy][iz  ] 
     [ix][iy][iz-1] -> [ix][iy  ][iz] -> [ix][iy+1][iz  ] -> [ix][iy+1][iz-1] -> [ix][iy][iz-1]  
  */
  
  /* First Plaquette xy plane */
	
  /* [ix  ][iy+1][iz] -> [ix+1][iy+1][iz] */
  periodic_boundaries(ix  , iy+1 , iz, c_ix_pxy_1, c_iy_pxy_1, c_iz_pxy_1); 
  /* [ix+1][iy+1][iz] -> [ix+1][iy  ][iz] */
  periodic_boundaries(ix+1, iy   , iz, c_ix_pxy_2, c_iy_pxy_2, c_iz_pxy_2);  
  
  tmp_plaq_1y = sin(a*(Curr_pos.A[1][ix][iy][iz]+
		       Curr_pos.A[0][c_ix_pxy_1][c_iy_pxy_1][c_iz_pxy_1]-
		       Curr_pos.A[1][c_ix_pxy_2][c_iy_pxy_2][c_iz_pxy_2]-
		       Curr_pos.A[0][ix][iy][iz]))/pow(a, 3);
  
  return(tmp_plaq_1y);
}
/***************************************************************************************/
inline float PlaqY2(Coords& Curr_pos, int ix, int iy, int iz){
  int c_ix_pxy_1, c_ix_pxy_2, c_ix_pxy_3;
  int c_iy_pxy_1, c_iy_pxy_2, c_iy_pxy_3; 
  int c_iz_pxy_1, c_iz_pxy_2, c_iz_pxy_3; /* Coordinated of 2 corners 
					     of the plaquette in the xy plane. */
  float tmp_plaq_2y;

  /* Second Plaquette xy plane*/
  
  /* [ix-1][iy  ][iz] -> [ix  ][iy  ][iz] */
  periodic_boundaries(ix-1, iy  , iz  , c_ix_pxy_1, c_iy_pxy_1, c_iz_pxy_1); 
  /* [ix  ][iy+1][iz] -> [ix-1][iy+1][iz] */
  periodic_boundaries(ix-1, iy+1, iz  , c_ix_pxy_2, c_iy_pxy_2, c_iz_pxy_2); 
  /* [ix-1][iy+1][iz] -> [ix-1][iy  ][iz] */
  periodic_boundaries(ix-1, iy  , iz  , c_ix_pxy_3, c_iy_pxy_3, c_iz_pxy_3); 
  
  tmp_plaq_2y = sin(a*(Curr_pos.A[0][c_ix_pxy_1][c_iy_pxy_1][c_iz_pxy_1]+
		       Curr_pos.A[1][ix][iy][iz]-
		       Curr_pos.A[0][c_ix_pxy_2][c_iy_pxy_2][c_iz_pxy_2]-
		       Curr_pos.A[1][c_ix_pxy_3][c_iy_pxy_3][c_iz_pxy_3]))/pow(a, 3);
  
  return(tmp_plaq_2y);
}
/***************************************************************************************/
inline float PlaqY3(Coords& Curr_pos, int ix, int iy, int iz){

  int c_ix_pyz_1, c_ix_pyz_2; 
  int c_iy_pyz_1, c_iy_pyz_2; 
  int c_iz_pyz_1, c_iz_pyz_2;/* Coordinated of 2 corners 
				of the plaquette in the yz plane. */
  float tmp_plaq_3y;
  /* Third Plaquette yz plane */
  
  /* [ix][iy+1][iz  ] -> [ix][iy+1][iz+1] */
  periodic_boundaries(ix  , iy+1, iz  , c_ix_pyz_1, c_iy_pyz_1, c_iz_pyz_1); 
  /* [ix][iy+1][iz+1] -> [ix][iy  ][iz+1] */
  periodic_boundaries(ix  , iy  , iz+1, c_ix_pyz_2, c_iy_pyz_2, c_iz_pyz_2);  
  
  tmp_plaq_3y = sin(a*(Curr_pos.A[1][ix][iy][iz]+
		       Curr_pos.A[2][c_ix_pyz_1][c_iy_pyz_1][c_iz_pyz_1]-
		       Curr_pos.A[1][c_ix_pyz_2][c_iy_pyz_2][c_iz_pyz_2]-
		       Curr_pos.A[2][ix][iy][iz]))/pow(a, 3);
  
  return(tmp_plaq_3y);
}
/***************************************************************************************/
inline float PlaqY4(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_pyz_1, c_ix_pyz_2, c_ix_pyz_3; 
  int c_iy_pyz_1, c_iy_pyz_2, c_iy_pyz_3; 
  int c_iz_pyz_1, c_iz_pyz_2, c_iz_pyz_3;/* Coordinated of 2 corners 
					    of the plaquette in the yz plane. */
  float tmp_plaq_4y;
  /* Fourth Plaquette yz plane*/
  
  /* [ix][iy  ][iz-1] -> [ix][iy  ][iz  ] */
  periodic_boundaries(ix  , iy  , iz-1, c_ix_pyz_1, c_iy_pyz_1, c_iz_pyz_1); 
  /* [ix][iy+1][iz  ] -> [ix][iy+1][iz-1] */
  periodic_boundaries(ix  , iy+1, iz-1, c_ix_pyz_2, c_iy_pyz_2, c_iz_pyz_2); 
  /* [ix][iy+1][iz-1] -> [ix][iy  ][iz-1] */
  periodic_boundaries(ix  , iy  , iz-1, c_ix_pyz_3, c_iy_pyz_3, c_iz_pyz_3); 

  tmp_plaq_4y = sin(a*(Curr_pos.A[2][c_ix_pyz_1][c_iy_pyz_1][c_iz_pyz_1]+
		       Curr_pos.A[1][ix][iy][iz]-
		       Curr_pos.A[2][c_ix_pyz_2][c_iy_pyz_2][c_iz_pyz_2]-
		       Curr_pos.A[1][c_ix_pyz_3][c_iy_pyz_3][c_iz_pyz_3]))/pow(a, 3);
  
  return(tmp_plaq_4y);
}
/***************************************************************************************/


inline float PlaqZ1(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_pxz_1, c_ix_pxz_2;
  int c_iy_pxz_1, c_iy_pxz_2;
  int c_iz_pxz_1, c_iz_pxz_2; /* Coordinated of 2 corners 
				 of the plaquette in the xz plane. */
  float tmp_plaq_1z;
  /* Link along y direction whit origin at [ix][iy][iz]: A[2][ix][iy][iz]. 
     The four plaquettes that contribute are 
     [ix  ][iy][iz] -> [ix][iy][iz+1] -> [ix+1][iy][iz+1] -> [ix+1][iy][iz  ] -> [ix  ][iy][iz] 
     [ix-1][iy][iz] -> [ix][iy][iz  ] -> [ix][iy  ][iz+1] -> [ix-1][iy][iz+1] -> [ix-1][iy][iz] 
     [ix][iy  ][iz] -> [ix][iy][iz+1] -> [ix][iy+1][iz+1] -> [ix][iy+1][iz  ] -> [ix][iy  ][iz] 
     [ix][iy-1][iz] -> [ix][iy][iz  ] -> [ix][iy  ][iz+1] -> [ix][iy-1][iz+1] -> [ix][iy-1][iz] 
  */
  
  /* First Plaquette xz plane */
  
  /* [ix  ][iy][iz+1] -> [ix+1][iy][iz+1] */
  periodic_boundaries(ix  , iy  , iz+1, c_ix_pxz_1, c_iy_pxz_1, c_iz_pxz_1); 
  /* [ix+1][iy][iz+1] -> [ix+1][iy][iz  ] */
  periodic_boundaries(ix+1, iy  , iz  , c_ix_pxz_2, c_iy_pxz_2, c_iz_pxz_2);  
  
  tmp_plaq_1z = sin(a*(Curr_pos.A[2][ix][iy][iz]+
		       Curr_pos.A[0][c_ix_pxz_1][c_iy_pxz_1][c_iz_pxz_1]-
		       Curr_pos.A[2][c_ix_pxz_2][c_iy_pxz_2][c_iz_pxz_2]-
		       Curr_pos.A[0][ix][iy][iz]))/pow(a, 3);
  return(tmp_plaq_1z);
}
/***************************************************************************************/
inline float PlaqZ2(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_pxz_1, c_ix_pxz_2, c_ix_pxz_3;
  int c_iy_pxz_1, c_iy_pxz_2, c_iy_pxz_3;
  int c_iz_pxz_1, c_iz_pxz_2, c_iz_pxz_3; /* Coordinated of 2 corners 
					     of the plaquette in the xz plane. */
  float tmp_plaq_2z;
  /* Second Plaquette xz plane*/
  
  /* [ix-1][iy][iz  ] -> [ix  ][iy][iz  ] */
  periodic_boundaries(ix-1, iy  , iz  , c_ix_pxz_1, c_iy_pxz_1, c_iz_pxz_1); 
  /* [ix  ][iy][iz+1] -> [ix-1][iy][iz+1] */
  periodic_boundaries(ix-1, iy  , iz+1, c_ix_pxz_2, c_iy_pxz_2, c_iz_pxz_2); 
  /* [ix-1][iy][iz+1] -> [ix-1][iy][iz  ] */
  periodic_boundaries(ix-1, iy  , iz  , c_ix_pxz_3, c_iy_pxz_3, c_iz_pxz_3); 
  
  tmp_plaq_2z = sin(a*(Curr_pos.A[0][c_ix_pxz_1][c_iy_pxz_1][c_iz_pxz_1]+
		       Curr_pos.A[2][ix][iy][iz]-
		       Curr_pos.A[0][c_ix_pxz_2][c_iy_pxz_2][c_iz_pxz_2]-
		       Curr_pos.A[2][c_ix_pxz_3][c_iy_pxz_3][c_iz_pxz_3]))/pow(a, 3);
  return(tmp_plaq_2z);
}
/***************************************************************************************/
inline float PlaqZ3(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_pyz_1, c_ix_pyz_2; 
  int c_iy_pyz_1, c_iy_pyz_2; 
  int c_iz_pyz_1, c_iz_pyz_2;/* Coordinated of 2 corners 
				of the plaquette in the yz plane. */
  float tmp_plaq_3z;
  /* Third Plaquette yz plane */
  
  /* [ix][iy  ][iz+1] -> [ix][iy+1][iz+1] */
  periodic_boundaries(ix  , iy  , iz+1, c_ix_pyz_1, c_iy_pyz_1, c_iz_pyz_1); 
  /* [ix][iy+1][iz+1] -> [ix][iy+1][iz  ] */
  periodic_boundaries(ix  , iy+1, iz  , c_ix_pyz_2, c_iy_pyz_2, c_iz_pyz_2);  
  
  tmp_plaq_3z = sin(a*(Curr_pos.A[2][ix][iy][iz]+
		       Curr_pos.A[1][c_ix_pyz_1][c_iy_pyz_1][c_iz_pyz_1]-
		       Curr_pos.A[2][c_ix_pyz_2][c_iy_pyz_2][c_iz_pyz_2]-
		       Curr_pos.A[1][ix][iy][iz]))/pow(a, 3);
  return(tmp_plaq_3z);
}
/***************************************************************************************/
inline float PlaqZ4(Coords& Curr_pos, int ix, int iy, int iz){

  int c_ix_pyz_1, c_ix_pyz_2, c_ix_pyz_3; 
  int c_iy_pyz_1, c_iy_pyz_2, c_iy_pyz_3; 
  int c_iz_pyz_1, c_iz_pyz_2, c_iz_pyz_3;/* Coordinated of 2 corners 
					    of the plaquette in the yz plane. */
  float tmp_plaq_4z;
  /* Fourth Plaquette yz plane*/
  
  /*  */
  periodic_boundaries(ix  , iy-1, iz  , c_ix_pyz_1, c_iy_pyz_1, c_iz_pyz_1); 
  /*  */
  periodic_boundaries(ix  , iy-1, iz+1, c_ix_pyz_2, c_iy_pyz_2, c_iz_pyz_2); 
  /*  */
  periodic_boundaries(ix  , iy-1, iz  , c_ix_pyz_3, c_iy_pyz_3, c_iz_pyz_3); 
  
  tmp_plaq_4z = sin(a*(Curr_pos.A[1][c_ix_pyz_1][c_iy_pyz_1][c_iz_pyz_1]+
		       Curr_pos.A[2][ix][iy][iz]-
		       Curr_pos.A[1][c_ix_pyz_2][c_iy_pyz_2][c_iz_pyz_2]-
		       Curr_pos.A[2][c_ix_pyz_3][c_iy_pyz_3][c_iz_pyz_3]))/pow(a, 3);
  return(tmp_plaq_4z);
}
/***************************************************************************************/

/***************************************************************************************/
inline float ScalX(Coords& Curr_pos, int ix, int iy, int iz){
  int c_ix_scal, c_iy_scal, c_iz_scal; /* Coordinates for the scalar fields 
					  at the other end of the link. 
					  (The first end being [ix][iy][iz])*/
  float tmp_plaq_scalx;
  
  /* Scalar field interaction term. The scalars that 
     contribute are the ones at the ends of 
     the A[0][ix][iy][iz] link field, 
     namely the ones at [ix][iy][iz] and [ix+1][iy][iz]. */
  periodic_boundaries(ix+1, iy, iz, c_ix_scal, c_iy_scal, c_iz_scal);
	
  tmp_plaq_scalx = (cos(a*Curr_pos.A[0][ix][iy][iz])*
		    (Curr_pos.phi[0][ix][iy][iz]*Curr_pos.phi[1][c_ix_scal][c_iy_scal][c_iz_scal] - 
		     Curr_pos.phi[1][ix][iy][iz]*Curr_pos.phi[0][c_ix_scal][c_iy_scal][c_iz_scal])+
		    sin(a*Curr_pos.A[0][ix][iy][iz])*
		    (Curr_pos.phi[0][ix][iy][iz]*Curr_pos.phi[0][c_ix_scal][c_iy_scal][c_iz_scal] + 
		     Curr_pos.phi[1][ix][iy][iz]*Curr_pos.phi[1][c_ix_scal][c_iy_scal][c_iz_scal]))*
    (pow(g, 2)/(2*a));
  
  return(tmp_plaq_scalx);
}
/***************************************************************************************/
inline float ScalY(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_scal, c_iy_scal, c_iz_scal; /* Coordinates for the scalar fields 
					  at the other end of the link. 
					  (The first end being [ix][iy][iz])*/
  float tmp_plaq_scaly;
  /* Scalar field interaction term. The scalars that 
     contribure are the ones at the ends of 
     the A[1][ix][iy][iz] link field, 
     namely the ones at [ix][iy][iz] and [ix][iy+1][iz]. */
  periodic_boundaries(ix, iy+1, iz, c_ix_scal, c_iy_scal, c_iz_scal);
  
  tmp_plaq_scaly = (cos(a*Curr_pos.A[1][ix][iy][iz])*
		    (Curr_pos.phi[0][ix][iy][iz]*Curr_pos.phi[1][c_ix_scal][c_iy_scal][c_iz_scal] - 
		     Curr_pos.phi[1][ix][iy][iz]*Curr_pos.phi[0][c_ix_scal][c_iy_scal][c_iz_scal])+
		    sin(a*Curr_pos.A[1][ix][iy][iz])*
		    (Curr_pos.phi[0][ix][iy][iz]*Curr_pos.phi[0][c_ix_scal][c_iy_scal][c_iz_scal] + 
		     Curr_pos.phi[1][ix][iy][iz]*Curr_pos.phi[1][c_ix_scal][c_iy_scal][c_iz_scal]))
    *(pow(g, 2)/(2*a));
  
  return(tmp_plaq_scaly);
}
/***************************************************************************************/
inline float ScalZ(Coords& Curr_pos, int ix, int iy, int iz){
  
  int c_ix_scal, c_iy_scal, c_iz_scal; /* Coordinates for the scalar fields 
					  at the other end of the link. 
					  (The first end being [ix][iy][iz])*/
  float tmp_plaq_scalz;
  /* Scalar field interaction term. The scalars that 
     contribure are the ones at the ends of 
     the A[2][ix][iy][iz] link field, 
     namely the ones at [ix][iy][iz] and [ix][iy][iz+1]. */
  periodic_boundaries(ix, iy, iz+1, c_ix_scal, c_iy_scal, c_iz_scal);
  
  tmp_plaq_scalz = (cos(a*Curr_pos.A[2][ix][iy][iz])*
		    (Curr_pos.phi[0][ix][iy][iz]*Curr_pos.phi[1][c_ix_scal][c_iy_scal][c_iz_scal] - 
		     Curr_pos.phi[1][ix][iy][iz]*Curr_pos.phi[0][c_ix_scal][c_iy_scal][c_iz_scal])+
		    sin(a*Curr_pos.A[2][ix][iy][iz])*
		    (Curr_pos.phi[0][ix][iy][iz]*Curr_pos.phi[0][c_ix_scal][c_iy_scal][c_iz_scal] + 
		     Curr_pos.phi[1][ix][iy][iz]*Curr_pos.phi[1][c_ix_scal][c_iy_scal][c_iz_scal]))
    *(pow(g, 2)/(2*a));
  
  return(tmp_plaq_scalz);
}
/***************************************************************************************/
inline float ScalarRealPart(Coords& Curr_pos, int ix, int iy, int iz){

  int c_ix_P, c_iy_P, c_iz_P; /* m + j "m Plus  j"*/
  int c_ix_M, c_iy_M, c_iz_M; /* m - j "m Minus j"*/
  float tmp0; /* For the real and imaginary 
		       parts of the scalar field. */

  /* tmp0 corresponds to the acceleration of the real part of the scalar field. */
  tmp0 = -lambda*Curr_pos.phi[0][ix][iy][iz]*
    (pow(Curr_pos.phi[0][ix][iy][iz], 2) + pow(Curr_pos.phi[1][ix][iy][iz], 2) - pow(sigma, 2));
  
  periodic_boundaries(ix+1, iy, iz, c_ix_P, c_iy_P, c_iz_P);      
  periodic_boundaries(ix-1, iy, iz, c_ix_M, c_iy_M, c_iz_M);
  
  /* Contribution to Re\phi[ix][iy][iz] from the gauge fields along the x-direction. */
  tmp0 = tmp0 + 
    pow(1/a, 2)*(cos(a*Curr_pos.A[0][ix][iy][iz])*Curr_pos.phi[0][c_ix_P][c_iy_P][c_iz_P] - 
		 sin(a*Curr_pos.A[0][ix][iy][iz])*Curr_pos.phi[1][c_ix_P][c_iy_P][c_iz_P] - 
		 2*Curr_pos.phi[0][ix][iy][iz] + 
		 cos(a*Curr_pos.A[0][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[0][c_ix_M][c_iy_M][c_iz_M]+
		 sin(a*Curr_pos.A[0][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[1][c_ix_M][c_iy_M][c_iz_M]);
  
  periodic_boundaries(ix, iy+1, iz, c_ix_P, c_iy_P, c_iz_P);      
  periodic_boundaries(ix, iy-1, iz, c_ix_M, c_iy_M, c_iz_M);
  
  /* Contribution to Re\phi[ix][iy][iz] from the gauge fields along the y-direction. */
  tmp0 = tmp0 + 
    pow(1/a, 2)*(cos(a*Curr_pos.A[1][ix][iy][iz])*Curr_pos.phi[0][c_ix_P][c_iy_P][c_iz_P] - 
		 sin(a*Curr_pos.A[1][ix][iy][iz])*Curr_pos.phi[1][c_ix_P][c_iy_P][c_iz_P] - 
		 2*Curr_pos.phi[0][ix][iy][iz] + 
		 cos(a*Curr_pos.A[1][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[0][c_ix_M][c_iy_M][c_iz_M]+
		 sin(a*Curr_pos.A[1][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[1][c_ix_M][c_iy_M][c_iz_M]);

  periodic_boundaries(ix, iy, iz+1, c_ix_P, c_iy_P, c_iz_P);      
  periodic_boundaries(ix, iy, iz-1, c_ix_M, c_iy_M, c_iz_M);
  
  /* Contribution to Re\phi[ix][iy][iz] from the gauge fields along the y-direction. */
  tmp0 = tmp0 + 
    pow(1/a, 2)*(cos(a*Curr_pos.A[2][ix][iy][iz])*Curr_pos.phi[0][c_ix_P][c_iy_P][c_iz_P] - 
		 sin(a*Curr_pos.A[2][ix][iy][iz])*Curr_pos.phi[1][c_ix_P][c_iy_P][c_iz_P] - 
		 2*Curr_pos.phi[0][ix][iy][iz] + 
		 cos(a*Curr_pos.A[2][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[0][c_ix_M][c_iy_M][c_iz_M]+
		 sin(a*Curr_pos.A[2][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[1][c_ix_M][c_iy_M][c_iz_M]);
  
  return(tmp0);
}
/***************************************************************************************/
inline float ScalarImagPart(Coords& Curr_pos, int ix, int iy, int iz){

  int c_ix_P, c_iy_P, c_iz_P; /* m + j "m Plus  j"*/
  int c_ix_M, c_iy_M, c_iz_M; /* m - j "m Minus j"*/
  float tmp1; /* For the real and imaginary 
		       parts of the scalar field. */

  /* tmp0 corresponds to the acceleration of the real part of the scalar field. */
  tmp1 = -lambda*Curr_pos.phi[1][ix][iy][iz]*
    (pow(Curr_pos.phi[0][ix][iy][iz], 2) + pow(Curr_pos.phi[1][ix][iy][iz], 2) - pow(sigma, 2));
  /* Time evolve the real part of the scalar field. */
  
  periodic_boundaries(ix+1, iy, iz, c_ix_P, c_iy_P, c_iz_P);      
  periodic_boundaries(ix-1, iy, iz, c_ix_M, c_iy_M, c_iz_M);
  
  /* Contribution to Im\phi[ix][iy][iz] from the gauge fields along the x-direction. */
  
  tmp1 = tmp1 + 
    pow(1/a, 2)*(cos(a*Curr_pos.A[0][ix][iy][iz])*Curr_pos.phi[1][c_ix_P][c_iy_P][c_iz_P] +
		 sin(a*Curr_pos.A[0][ix][iy][iz])*Curr_pos.phi[0][c_ix_P][c_iy_P][c_iz_P] -
		 2*Curr_pos.phi[1][ix][iy][iz] +
		 cos(a*Curr_pos.A[0][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[1][c_ix_M][c_iy_M][c_iz_M] -
		 sin(a*Curr_pos.A[0][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[0][c_ix_M][c_iy_M][c_iz_M]);
  
  periodic_boundaries(ix, iy+1, iz, c_ix_P, c_iy_P, c_iz_P);      
  periodic_boundaries(ix, iy-1, iz, c_ix_M, c_iy_M, c_iz_M);

  /* Contribution to Im\phi[ix][iy][iz] from the gauge fields along the y-direction. */
  
  tmp1 = tmp1 + 
    pow(1/a, 2)*(cos(a*Curr_pos.A[1][ix][iy][iz])*Curr_pos.phi[1][c_ix_P][c_iy_P][c_iz_P] +
		 sin(a*Curr_pos.A[1][ix][iy][iz])*Curr_pos.phi[0][c_ix_P][c_iy_P][c_iz_P] -
		 2*Curr_pos.phi[1][ix][iy][iz] +
		 cos(a*Curr_pos.A[1][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[1][c_ix_M][c_iy_M][c_iz_M] -
		 sin(a*Curr_pos.A[1][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[0][c_ix_M][c_iy_M][c_iz_M]);
  
  periodic_boundaries(ix, iy, iz+1, c_ix_P, c_iy_P, c_iz_P);      
  periodic_boundaries(ix, iy, iz-1, c_ix_M, c_iy_M, c_iz_M);
  
  /* Contribution to Im\phi[ix][iy][iz] from the gauge fields along the y-direction. */
  
  tmp1 = tmp1 + 
    pow(1/a, 2)*(cos(a*Curr_pos.A[2][ix][iy][iz])*Curr_pos.phi[1][c_ix_P][c_iy_P][c_iz_P] +
		 sin(a*Curr_pos.A[2][ix][iy][iz])*Curr_pos.phi[0][c_ix_P][c_iy_P][c_iz_P] -
		 2*Curr_pos.phi[1][ix][iy][iz] +
		 cos(a*Curr_pos.A[2][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[1][c_ix_M][c_iy_M][c_iz_M] -
		 sin(a*Curr_pos.A[2][c_ix_M][c_iy_M][c_iz_M])*Curr_pos.phi[0][c_ix_M][c_iy_M][c_iz_M]);
  
  return(tmp1);
}
/***************************************************************************************/
inline void UpdateScalarVelReal(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Update the momenta and the field: Real part .*/
  Future_pos.vel[0][ix][iy][iz] = Curr_pos.vel[0][ix][iy][iz] + 
    dt*(ScalarRealPart(Curr_pos, ix, iy, iz) - H*Curr_pos.vel[0][ix][iy][iz]);
}
/***************************************************************************************/
inline void UpdateScalarVelImag(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Update the momenta and the field: Imaginary part .*/
  Future_pos.vel[1][ix][iy][iz] = Curr_pos.vel[1][ix][iy][iz] + 
    dt*(ScalarImagPart(Curr_pos, ix, iy, iz) - H*Curr_pos.vel[1][ix][iy][iz] );
}
/***************************************************************************************/
inline void UpdateScalarReal(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Update the momenta and the field: Real part .*/
  Future_pos.phi[0][ix][iy][iz] = Curr_pos.phi[0][ix][iy][iz] + dt*Future_pos.vel[0][ix][iy][iz]; 
}
/***************************************************************************************/
inline void UpdateScalarImag(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Update the momenta and the field: Imaginary part .*/
  Future_pos.phi[1][ix][iy][iz] = Curr_pos.phi[1][ix][iy][iz] + dt*Future_pos.vel[1][ix][iy][iz];
}
/***************************************************************************************/
inline void UpdateEFieldX(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Finally update the electric field and the gauge potential. 
     The minus sign must be the correct sign. */
  if(Lock_Gauge_Field){
    Future_pos.A[0][ix][iy][iz] = 0;
    Future_pos.E[0][ix][iy][iz] = 0;
  }
  else{
    Future_pos.E[0][ix][iy][iz] = Curr_pos.E[0][ix][iy][iz] - 
      dt*(PlaqX1(Curr_pos, ix, iy, iz) + PlaqX2(Curr_pos, ix, iy, iz) + 
	  PlaqX3(Curr_pos, ix, iy, iz) + PlaqX4(Curr_pos, ix, iy, iz) + 
	  ScalX(Curr_pos, ix, iy, iz)  + H*Curr_pos.E[0][ix][iy][iz]);
  }
}
/***************************************************************************************/
inline void UpdateEFieldY(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Finally update the electric field and the gauge potential.
     The minus sign must be the correct sign. */ 
  if(Lock_Gauge_Field){
    Future_pos.A[1][ix][iy][iz] = 0;
    Future_pos.E[1][ix][iy][iz] = 0;
  }
  else{
    Future_pos.E[1][ix][iy][iz] = Curr_pos.E[1][ix][iy][iz] - 
      dt*(PlaqY1(Curr_pos, ix, iy, iz) + PlaqY2(Curr_pos, ix, iy, iz) + 
	  PlaqY3(Curr_pos, ix, iy, iz) + PlaqY4(Curr_pos, ix, iy, iz) + 
	  ScalY(Curr_pos, ix, iy, iz)  + H*Curr_pos.E[1][ix][iy][iz]);
  }
}
/***************************************************************************************/
inline void UpdateEFieldZ(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Finally update the electric field and the gauge potential.
     The minus sign must be the correct sign. */
  if(Lock_Gauge_Field){
    Future_pos.A[2][ix][iy][iz] = 0;
    Future_pos.E[2][ix][iy][iz] = 0;
  }
  else{
    Future_pos.E[2][ix][iy][iz] = Curr_pos.E[2][ix][iy][iz] - 
      dt*(PlaqZ1(Curr_pos, ix, iy, iz) + PlaqZ2(Curr_pos, ix, iy, iz) + 
	  PlaqZ3(Curr_pos, ix, iy, iz) + PlaqZ4(Curr_pos, ix, iy, iz) + 
	  ScalZ(Curr_pos, ix, iy, iz)  + H*Curr_pos.E[2][ix][iy][iz]); 
  }
}

/***************************************************************************************/
inline void UpdateLinkX(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Finally update the electric field and the gauge potential. 
     The minus sign must be the correct sign. */
  Future_pos.A[0][ix][iy][iz] = Curr_pos.A[0][ix][iy][iz] + dt*Future_pos.E[0][ix][iy][iz];
}
/***************************************************************************************/
inline void UpdateLinkY(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Finally update the electric field and the gauge potential.
     The minus sign must be the correct sign. */ 
  Future_pos.A[1][ix][iy][iz] = Curr_pos.A[1][ix][iy][iz] + dt*Future_pos.E[1][ix][iy][iz]; 
}
/***************************************************************************************/
inline void UpdateLinkZ(Coords& Future_pos, Coords& Curr_pos, int ix, int iy, int iz){
  /* Finally update the electric field and the gauge potential.
     The minus sign must be the correct sign. */
  Future_pos.A[2][ix][iy][iz] = Curr_pos.A[2][ix][iy][iz] + dt*Future_pos.E[2][ix][iy][iz];
}

void* EvolveMT(void* thread_data){
  ThreadData* local_thread_data =  static_cast<ThreadData*>(thread_data);
  int low_x_Lim   = static_cast<int>((N_x*local_thread_data->thread_no)/NoOfThreads);
  int high_x_Lim  = static_cast<int>((N_x*(local_thread_data->thread_no + 1))/NoOfThreads);
  for(int ix = low_x_Lim; ix < high_x_Lim; ++ix){
    for(int iy = 0; iy < N_y; ++iy){
      for(int iz = 0; iz < N_z; ++iz){//Start looping 
		
	UpdateScalarVelReal(*T2, *T1, ix, iy, iz);
	UpdateScalarVelImag(*T2, *T1, ix, iy, iz);
	
	UpdateEFieldX(*T2, *T1, ix, iy, iz);
	UpdateEFieldY(*T2, *T1, ix, iy, iz);
	UpdateEFieldZ(*T2, *T1, ix, iy, iz);
      }
    }
  }
  for(int ix = low_x_Lim; ix < high_x_Lim; ++ix){
    for(int iy = 0; iy < N_y; ++iy){
      for(int iz = 0; iz < N_z; ++iz){//Start looping 
	
	UpdateScalarReal(*T2, *T1, ix, iy, iz);
	UpdateScalarImag(*T2, *T1, ix, iy, iz);

	UpdateLinkX(*T2, *T1, ix, iy, iz);
	UpdateLinkY(*T2, *T1, ix, iy, iz);
	UpdateLinkZ(*T2, *T1, ix, iy, iz);
      }
    }
  }
#ifdef DEBUG
  std::cout << "EvolveMT : low_x_Lim = " << low_x_Lim << " high_x_Lim = " << high_x_Lim << std::endl;
#endif
}

inline void Check_Gauss(Coords& Curr_pos, unsigned int& local_frame){
  int c_ix, c_iy, c_iz;
  float temp_x, temp_y, temp_z, temp_scal, temp_gauss;
  char buffer[16];
  sprintf(buffer, "%03d", local_frame);
  std::string curr_gauss_out("gauss" + std::string(buffer) + ".txt");
  std::ofstream gauss_out(curr_gauss_out.c_str());
  for(int ix = 0; ix < N_x; ++ix){
    for(int iy = 0; iy < N_y; ++iy){
      for(int iz = 0; iz < N_z; ++iz){
	if((ix%grid_reduction == 0)&&
	   (iy%grid_reduction == 0)&&
	   (iz%grid_reduction == 0)){
	  /* Contribution from the electric field along the x direction. 
	     The links that contribute  are the following ones: 
	     [ix-1][iy][iz] -> [ix][iy][iz] and [ix][iy][iz] -> [ix+1][iy][iz]. 
	     They correspond to E[0][ix-1][iy][iz] and E[0][ix][iy][iz] */
	  
	  periodic_boundaries(ix-1, iy, iz, c_ix, c_iy, c_iz);
	  temp_x = Curr_pos.E[0][ix][iy][iz] - Curr_pos.E[0][c_ix][c_iy][c_iz];
	  
	  /* Contribution from the electric field along the y direction. 
	     The links that contribute  are the following ones: 
	     [ix][iy-1][iz] -> [ix][iy][iz] and [ix][iy][iz] -> [ix][iy+1][iz]. 
	     They correspond to E[1][ix][iy-1][iz] and E[1][ix][iy][iz] */
	  
	  periodic_boundaries(ix, iy-1, iz, c_ix, c_iy, c_iz);
	  temp_y = Curr_pos.E[1][ix][iy][iz] - Curr_pos.E[1][c_ix][c_iy][c_iz];

	   /* Contribution from the electric field along the y direction. 
	     The links that contribute  are the following ones: 
	     [ix][iy][iz-1] -> [ix][iy][iz] and [ix][iy][iz] -> [ix][iy][iz+1]. 
	     They correspond to E[1][ix][iy][iz-1] and E[1][ix][iy][iz] */
	  
	  periodic_boundaries(ix, iy, iz-1, c_ix, c_iy, c_iz);
	  temp_z = Curr_pos.E[2][ix][iy][iz] - Curr_pos.E[2][c_ix][c_iy][c_iz];
	  
	  temp_scal = (Curr_pos.phi[0][ix][iy][iz]*Curr_pos.vel[1][ix][iy][iz] - 
		       Curr_pos.phi[1][ix][iy][iz]*Curr_pos.vel[0][ix][iy][iz])*(a*pow(g, 2)/2);
	  
	  temp_gauss = temp_x + temp_y + temp_z + temp_scal;
	  gauss_out << ix << " " << iy << " " << iz << " " <<  temp_gauss << std::endl;
	}
      }
    }
  }
  gauss_out.close();
}

inline void Energy(Coords& Future_pos, Coords& Curr_pos, float (&tmp_Energy)[6]){

  int c_ix_pxy_1, c_iy_pxy_1, c_iz_pxy_1;
  int c_ix_pxy_2, c_iy_pxy_2, c_iz_pxy_2;
  
  int c_ix_pyz_1, c_iy_pyz_1, c_iz_pyz_1;
  int c_ix_pyz_2, c_iy_pyz_2, c_iz_pyz_2;

  int c_ix_pxz_1, c_iy_pxz_1, c_iz_pxz_1;
  int c_ix_pxz_2, c_iy_pxz_2, c_iz_pxz_2; /* Coordinated of 2 corners 
					     of the 3 plaquettes contributing 
					     to the B-field energy.*/
  
  int c_ix_1, c_ix_2, c_ix_3; 
  int c_iy_1, c_iy_2, c_iy_3;
  int c_iz_1, c_iz_2, c_iz_3;

  float Mag_Energy, El_Energy, Grad_phi_Energy, Kin_phi_Energy, Pot_Energy, Tot_Energy;

  Mag_Energy      = 0; 
  El_Energy       = 0; 
  Grad_phi_Energy = 0; 
  Kin_phi_Energy  = 0; 
  Pot_Energy      = 0;

  for(int ix = 0; ix < N_x; ++ix){
    for(int iy = 0; iy < N_y; ++iy){
      for(int iz = 0; iz < N_z; ++iz){
	
	/* Magnetic energy, plaquette term. */
	
	/* First Plaquette - xy plane. */
	
	/* [ix+1][iy  ][iz] -> [ix+1][iy+1][iz] */
	periodic_boundaries(ix+1, iy  , iz, c_ix_pxy_1, c_iy_pxy_1, c_iz_pxy_1); 
	/* [ix+1][iy+1][iz] -> [ix  ][iy+1][iz] */
	periodic_boundaries(ix  , iy+1, iz, c_ix_pxy_2, c_iy_pxy_2, c_iz_pxy_2);     
	
	Mag_Energy += -pow(sin(a*(Curr_pos.A[0][ix][iy][iz]+
				  Curr_pos.A[1][c_ix_pxy_1][c_iy_pxy_1][c_iz_pxy_1]-
				  Curr_pos.A[0][c_ix_pxy_2][c_iy_pxy_2][c_iz_pxy_2]-
				  Curr_pos.A[1][ix][iy][iz])), 2)/(pow(a*a*g, 2));
	
	/* Second Plaquette - xz plane. */

	/* [ix+1][iy][iz  ] -> [ix+1][iy][iz+1] */
	periodic_boundaries(ix+1, iy  , iz  , c_ix_pxz_1, c_iy_pxz_1, c_iz_pxz_1); 
	/* [ix+1][iy][iz+1] -> [ix  ][iy][iz+1] */
	periodic_boundaries(ix  , iy  , iz+1, c_ix_pxz_2, c_iy_pxz_2, c_iz_pxz_2);     
	
	Mag_Energy += -pow(sin(a*(Curr_pos.A[0][ix][iy][iz]+
				  Curr_pos.A[2][c_ix_pxz_1][c_iy_pxz_1][c_iz_pxz_1]-
				  Curr_pos.A[0][c_ix_pxz_2][c_iy_pxz_2][c_iz_pxz_2]-
				  Curr_pos.A[2][ix][iy][iz])), 2)/(pow(a*a*g, 2));
	
	/* Third Plaquette - yz plane. */

	/* [ix][iy+1][iz  ] -> [ix][iy+1][iz+1] */
	periodic_boundaries(ix  , iy+1, iz  , c_ix_pyz_1, c_iy_pyz_1, c_iz_pyz_1); 
	/* [ix][iy+1][iz+1] -> [ix][iy  ][iz+1] */
	periodic_boundaries(ix  , iy  , iz+1, c_ix_pyz_2, c_iy_pyz_2, c_iz_pyz_2);  
	
	Mag_Energy += -pow(sin(a*(Curr_pos.A[1][ix][iy][iz]+
				  Curr_pos.A[2][c_ix_pyz_1][c_iy_pyz_1][c_iz_pyz_1]-
				  Curr_pos.A[1][c_ix_pyz_2][c_iy_pyz_2][c_iz_pyz_2]-
				  Curr_pos.A[2][ix][iy][iz])), 2)/(pow(a*a*g, 2));
	
	// 	Mag_Energy += -pow((Curr_brane.A[0][i][j]+
// 			    Curr_brane.A[1][c_i_px_1][c_j_px_1]-
// 			    Curr_brane.A[0][c_i_px_2][c_j_px_2]-
// 			    Curr_brane.A[1][i][j]), 2)/(pow(a*g, 2));
	
      /* Energy of the electric field. 
	 It takes the average of the 
	 electric field values at
	 t + \delta t/2 and  t - \delta t/2, 
	 for each field component. */
	
	El_Energy += -(pow((Future_pos.E[0][ix][iy][iz] + 
			    Curr_pos.E[0][ix][iy][iz])/2, 2)+
		       pow((Future_pos.E[1][ix][iy][iz] + 
			    Curr_pos.E[1][ix][iy][iz])/2, 2)+
		       pow((Future_pos.E[2][ix][iy][iz] + 
			    Curr_pos.E[2][ix][iy][iz])/2, 2))/(pow(g, 2));
	
	Kin_phi_Energy += (2*Curr_pos.phi[0][ix][iy][iz]*Future_pos.phi[0][ix][iy][iz] + 
			   2*Curr_pos.phi[1][ix][iy][iz]*Future_pos.phi[1][ix][iy][iz] - 
			   (pow(Curr_pos.phi[0][ix][iy][iz], 2) + 
			    pow(Curr_pos.phi[1][ix][iy][iz], 2)) - 
			   (pow(Future_pos.phi[0][ix][iy][iz], 2) + 
			    pow(Future_pos.phi[1][ix][iy][iz], 2)))/(2*pow(dt, 2));
	
	periodic_boundaries(ix+1, iy  , iz  , c_ix_1, c_iy_1, c_iz_1);
	periodic_boundaries(ix  , iy+1, iz  , c_ix_2, c_iy_2, c_iz_2);
	periodic_boundaries(ix  , iy  , iz+1, c_ix_3, c_iy_3, c_iz_3);
	
	Grad_phi_Energy += (2*cos(a*Curr_pos.A[0][ix][iy][iz])*
			    (Curr_pos.phi[0][c_ix_1][c_iy_1][c_iz_1]*Curr_pos.phi[0][ix][iy][iz] + 
			     Curr_pos.phi[1][c_ix_1][c_iy_1][c_iz_1]*Curr_pos.phi[1][ix][iy][iz]) - 
			    2*sin(a*Curr_pos.A[0][ix][iy][iz])*
			    (Curr_pos.phi[1][c_ix_1][c_iy_1][c_iz_1]*Curr_pos.phi[0][ix][iy][iz] - 
			     Curr_pos.phi[0][c_ix_1][c_iy_1][c_iz_1]*Curr_pos.phi[1][ix][iy][iz]) +
			    2*cos(a*Curr_pos.A[1][ix][iy][iz])*
			    (Curr_pos.phi[0][c_ix_2][c_iy_2][c_iz_2]*Curr_pos.phi[0][ix][iy][iz] + 
			     Curr_pos.phi[1][c_ix_2][c_iy_2][c_iz_2]*Curr_pos.phi[1][ix][iy][iz]) -
			    2*sin(a*Curr_pos.A[1][ix][iy][iz])*
			    (Curr_pos.phi[1][c_ix_2][c_iy_2][c_iz_2]*Curr_pos.phi[0][ix][iy][iz] - 
			     Curr_pos.phi[0][c_ix_2][c_iy_2][c_iz_2]*Curr_pos.phi[1][ix][iy][iz]) +
			    2*cos(a*Curr_pos.A[2][ix][iy][iz])*
			    (Curr_pos.phi[0][c_ix_3][c_iy_3][c_iz_3]*Curr_pos.phi[0][ix][iy][iz] + 
			     Curr_pos.phi[1][c_ix_3][c_iy_3][c_iz_3]*Curr_pos.phi[1][ix][iy][iz]) -
			    2*sin(a*Curr_pos.A[2][ix][iy][iz])*
			    (Curr_pos.phi[1][c_ix_3][c_iy_3][c_iz_3]*Curr_pos.phi[0][ix][iy][iz] - 
			     Curr_pos.phi[0][c_ix_3][c_iy_3][c_iz_3]*Curr_pos.phi[1][ix][iy][iz])
			    -3*(pow(Curr_pos.phi[0][ix][iy][iz], 2) +  
				pow(Curr_pos.phi[1][ix][iy][iz], 2))- 
			    (pow(Curr_pos.phi[0][c_ix_1][c_iy_1][c_iz_1], 2) +  
			     pow(Curr_pos.phi[1][c_ix_1][c_iy_1][c_iz_1], 2))- 
			    (pow(Curr_pos.phi[0][c_ix_2][c_iy_2][c_iz_2], 2) +  
			     pow(Curr_pos.phi[1][c_ix_2][c_iy_2][c_iz_2], 2))-
			    (pow(Curr_pos.phi[0][c_ix_3][c_iy_3][c_iz_3], 2) +  
			     pow(Curr_pos.phi[1][c_ix_3][c_iy_3][c_iz_3], 2))
			    )/(2*pow(a, 2));
	
	Pot_Energy += -(lambda/4.0)*pow(pow(Curr_pos.phi[0][ix][iy][iz], 2) + 
					pow(Curr_pos.phi[1][ix][iy][iz], 2) - 
					pow(sigma, 2), 2);
	
	
      }
    }
  }
  
  Tot_Energy = Mag_Energy + El_Energy + Grad_phi_Energy + Kin_phi_Energy + Pot_Energy;

  tmp_Energy[0] = Mag_Energy;
  tmp_Energy[1] = El_Energy;
  tmp_Energy[2] = Grad_phi_Energy;
  tmp_Energy[3] = Kin_phi_Energy;
  tmp_Energy[4] = Pot_Energy;
  tmp_Energy[5] = Tot_Energy;
}

inline void reset_coords(Coords& curr_pos){
  for(int i = 0; i < N_x; ++i){
    for(int j = 0; j < N_y; ++j){
      for(int k = 0; k < N_z; ++k){
	for(int p = 0; p < 2; ++p){
	  curr_pos.phi[p][i][j][k] = 0;
	  curr_pos.vel[p][i][j][k] = 0;
	}
	for(int p = 0; p < dim; ++p){
	  curr_pos.A[p][i][j][k] = 0; 
	  curr_pos.E[p][i][j][k] = 0; 
	}
      }
    }
  }
}

inline void sync_coords(Coords& new_pos, Coords& old_pos){ 
  for(int i = 0; i < N_x; ++i){
    for(int j = 0; j < N_y; ++j){
      for(int k = 0; k < N_z; ++k){
	for(int p = 0; p < 2; ++p){
	  new_pos.phi[p][i][j][k] = old_pos.phi[p][i][j][k];
	  new_pos.vel[p][i][j][k] = old_pos.vel[p][i][j][k];
	}
	for(int p = 0; p < dim; ++p){
	  new_pos.A[p][i][j][k] = old_pos.A[p][i][j][k];
	  new_pos.E[p][i][j][k] = old_pos.E[p][i][j][k];
	} 
      }
    }
  }
}

void* F(void* thread_data){
  int rc = 0;
  for(unsigned int k = 0; k < reduction; ++k){
    EvolveMT(thread_data);
#ifdef DEBUG
    ThreadData* local_thread_data =  static_cast<ThreadData*>(thread_data);
    std::cout << "F; EvolveMT finished, thread = " << local_thread_data->thread_no << std::endl;
#endif
    rc = pthread_barrier_wait(&barr);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
      std::cout << "Could not wait on barrier" << std::endl;
      exit(-1);
    }
    syncT1T2MT(thread_data);
    rc = pthread_barrier_wait(&barr);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
      std::cout << "Could not wait on barrier" << std::endl;
      exit(-1);
    }
  }
  pthread_exit(NULL);
}

void writeScalars(unsigned int& local_frame){
  char buffer[16];
  sprintf(buffer, "%03d", local_frame);
  std::string curr_scalar_out("scalar" + std::string(buffer) + ".txt");
  std::ofstream scalar_out(curr_scalar_out.c_str());
  scalar_out.precision(4);
  for(int ix = 0; ix < N_x; ++ix){
    for(int iy = 0; iy < N_y; ++iy){
      for(int iz = 0; iz < N_z; ++iz){
	if((ix%grid_reduction == 0)&&(iy%grid_reduction == 0)&&(iz%grid_reduction == 0)){
	  scalar_out << ix << " " << iy << " " << iz << " "
		     <<  sigma - sqrt(pow(T1->phi[0][ix][iy][iz], 2) + 
				      pow(T1->phi[1][ix][iy][iz], 2)) << " "  << std::endl;  
	  
	}
      }
    }
  }
  scalar_out.close();
}

void writeLinks(unsigned int& local_frame){
  char buffer[16];
  sprintf(buffer, "%03d", local_frame);
  std::string curr_link_out("link" + std::string(buffer) + ".txt");	  
  std::ofstream link_out(curr_link_out.c_str());
  for(int ix = 0; ix < N_x; ++ix){
    for(int iy = 0; iy < N_y; ++iy){
      for(int iz = 0; iz < N_z; ++iz){
	if((ix%grid_reduction == 0)&&(iy%grid_reduction == 0)&&(iz%grid_reduction == 0)){
	  link_out << ix << " " << iy << " " << iz << " " ;
	  for(int p = 0; p < dim; ++p){
	    link_out << T1->A[p][iz][iy][iz] << " ";
	  }
	  for(int p = 0; p < dim; ++p){
	    link_out << T1->E[p][ix][iy][iz] << " ";
	  }
	  link_out << std::endl; 
	}
      }
    }
  }
  link_out.close();
}

int main(int argc, char* argv[]){  
  unsigned int local_frame;
  pthread_attr_t attr;
  
  unsigned int thr_id[NoOfThreads];
  ThreadData* TT[NoOfThreads];
  pthread_t  p_thread[NoOfThreads];

#ifdef LINUX
  cpu_set_t cpuset;
  size_t cpu_set_size = sizeof(cpuset);
#endif

  for(unsigned int i = 0; i < NoOfThreads; ++i){
    TT[i] = new ThreadData;
    TT[i]->thread_no = i;
  }

  reset_coords(*T1);
  reset_coords(*T2); 

  read_data(*T1);
  sync_coords(*T2, *T1);

  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if(pthread_barrier_init(&barr, NULL, NoOfThreads)){
    std::cout << "Could not create a barrier" << std::endl;
    exit(-1);
  }

#ifdef LINUX
  CPU_ZERO(&cpuset);
  for (unsigned int j = 0; j < NoOfThreads; j++){
    CPU_SET(j, &cpuset);
  }
#endif
  while(step < frames_total){
    for(unsigned int i = 0; i < NoOfThreads; ++i){
      /* Force each thread to run on only one CPU. */
#ifdef LINUX
      CPU_ZERO(&cpuset);
      CPU_SET(i, &cpuset);
      int s = pthread_attr_setaffinity_np(&attr, cpu_set_size, &cpuset);
      if (s != 0){ std::cerr << "pthread_attr_setaffinity_np() , s = " << s << std::endl; }
#endif
      thr_id[i] = pthread_create(&(p_thread)[i], &attr, F, static_cast<void*>(TT[i]));

#ifdef LINUX
#ifdef DEBUG
      int aff = pthread_attr_getaffinity_np(&attr,  cpu_set_size, &cpuset);
      if(aff == 0){
	for (int j = 0; j < CPU_SETSIZE; j++){
	  if (CPU_ISSET(j, &cpuset)){  printf(" attr: CPU %d\n", j);  }
	}
      }
#endif
      /* Uncomment this section if you want to set the affinity _after_ the thread is created. */
      //int s = pthread_setaffinity_np(p_thread[i], sizeof(cpu_set_t), &cpuset);
      //if (s != 0){ std::cerr << "pthread_setaffinity_np() , s = " << s << std::endl; }
#ifdef DEBUG
      //s = pthread_getaffinity_np(p_thread[i], sizeof(cpu_set_t), &cpuset);
      //printf("Set returned by pthread_getaffinity_np() contained:\n");

      for (int j = 0; j < CPU_SETSIZE; j++){
	//if (CPU_ISSET(j, &cpuset)){  printf(" CPU %d\n", j);  }
      }
#endif
#endif
    }
    
    for(unsigned int i = 0; i < NoOfThreads; ++i){
      pthread_join(p_thread[i], NULL);
    }
    
    if(frames_yes){
      local_frame = step;
      if(energy_yes){
	if(local_frame < frames_total){
	  Energy(*T2, *T1, temp_Energy);
	  for(int uu = 0; uu < E_params; ++uu){
	    Energy_vals[local_frame][uu] = temp_Energy[uu];
	  }
	}
      }
      if(gauss_yes){//Check if the constraint on the initial conditions is satisfied.
	Check_Gauss(*T2, local_frame);
      }
      if(link_yes){// Writing out the frames for the links. 
	writeLinks(local_frame);
      }
      if(scalars_yes){// Writing out the frames for the scalar field 
	writeScalars(local_frame);
      }
    }
    ++step;
    std::cout << " done " << step << "/" << frames_total << std::endl;
  }

  /* Delete the pointers to the heap. */
  delete T1;
  delete T2;
  pthread_attr_destroy(&attr);
  pthread_barrier_destroy(&barr);
  for(unsigned int i = 0; i < NoOfThreads; ++i){
    delete TT[i];
  }
  
  std::ofstream energy_out("Energy.txt");
  for(unsigned int i = 1; i < frames_total; ++i){
    energy_out << i << " ";
    for(int uu = 0; uu < E_params; ++uu){
      energy_out << pow(a, 3)*Energy_vals[i][uu] << " ";
    }
    energy_out << std::endl;
  }
  energy_out.close();
  return(0);
}

  
