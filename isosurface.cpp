#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <GL/glut.h>
#include <GL/glu.h>

const GLint dim = 3;
const GLint N_x = 32;
const GLint N_y = N_x;
const GLint N_z = N_x;

GLfloat eye_x = 0;
GLfloat eye_y = 0;
const GLfloat darken = 0.5;

typedef struct {
  GLfloat X[N_x][N_y][N_z];
  GLfloat Y[N_x][N_y][N_z];
  GLfloat Z[N_x][N_y][N_z];
  GLfloat F[N_x][N_y][N_z];
} Coords;

typedef struct {
  int LEnds[dim][dim]; /* The indices for the 3 other links of the triad. */
  GLint cross[dim]; /* 1 if the link is crossed, 0 if it is not */
} Triad;

typedef struct {
  GLfloat T_Corners[6][4][dim];
  int T_CornerIndices[6][4][dim];
} SixTetrahedra;

typedef struct{
  GLfloat V[dim];
  int End1Index[dim];
  int End2Index[dim];
} Vec3D;

typedef std::vector<Vec3D> VertexCoords;

typedef struct{
  bool BulkLink;
  unsigned int Corner1Index;
  unsigned int Corner2Index;
} VolumeLink;

GLint isosurface_landscape = 1;
GLint volume_edges   = 2;

/* This one works for smooth motions of the mouse */
void mauz(GLint x, GLint y){
  eye_x = x;
  eye_y = y;
  glutPostRedisplay();
}

void init_lighting(void)
{
  GLfloat twilight = 0.1;
  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[] = { 2.0 };
  GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
  GLfloat ambient_light[] = {twilight, twilight, twilight, 1.0};

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient_light);
  
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
  glColor3f(1.0, 0.0, 0.0);
  
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_COLOR_MATERIAL);
  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
}

void read_data(Coords& TT, char* infile, GLfloat scale_length, GLfloat treshold, GLint red){
  std::ifstream fin;
  fin.open(infile);
  if (fin.is_open()){
    for(int ix = 0; ix < N_x; ++ix){
      for(int iy = 0; iy < N_y; ++iy){
	for(int iz = 0; iz < N_z; ++iz){
	  fin >> TT.X[ix][iy][iz];
	  fin >> TT.Y[ix][iy][iz];
	  fin >> TT.Z[ix][iy][iz];
	  fin >> TT.F[ix][iy][iz];

	  TT.X[ix][iy][iz] = scale_length*(TT.X[ix][iy][iz] - red*(N_x - 1.0)/2);
	  TT.Y[ix][iy][iz] = scale_length*(TT.Y[ix][iy][iz] - red*(N_y - 1.0)/2);
	  TT.Z[ix][iy][iz] = scale_length*(TT.Z[ix][iy][iz] - red*(N_z - 1.0)/2);
	  
	  TT.F[ix][iy][iz] = TT.F[ix][iy][iz] - treshold;
	  if(TT.F[ix][iy][iz] > 0){
	  }
	}
      }
    }
  }
  else{
    std::cout << " Error: could not open " << infile << std::endl;
    exit(-1);
  }
  fin.close();
}

void inter_point(Coords& TT, 
		 int ix_1, int iy_1, int iz_1, 
		 int ix_2, int iy_2, int iz_2, GLfloat& X, GLfloat& Y, GLfloat&Z){
  GLfloat Dx, Dy, Dz;
  Dx = TT.X[ix_2][iy_2][iz_2] - TT.X[ix_1][iy_1][iz_1];
  Dy = TT.Y[ix_2][iy_2][iz_2] - TT.Y[ix_1][iy_1][iz_1];
  Dz = TT.Z[ix_2][iy_2][iz_2] - TT.Z[ix_1][iy_1][iz_1];
  X = TT.X[ix_1][iy_1][iz_1] + Dx*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1]) + fabs(TT.F[ix_2][iy_2][iz_2]));
  Y = TT.Y[ix_1][iy_1][iz_1] + Dy*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1]) + fabs(TT.F[ix_2][iy_2][iz_2]));
  Z = TT.Z[ix_1][iy_1][iz_1] + Dz*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1]) + fabs(TT.F[ix_2][iy_2][iz_2]));
}

void vector_product(GLfloat (&v1)[3], GLfloat (&v2)[3], GLfloat (&result)[3]){
  GLfloat norm = 0;
  
  result[0] = v1[1]*v2[2] - v1[2]*v2[1];
  result[1] = v1[2]*v2[0] - v1[0]*v2[2];
  result[2] = v1[0]*v2[1] - v1[1]*v2[0];

  norm = sqrt(pow(result[0], 2) + pow(result[1], 2) + pow(result[2], 2));

  result[0] = darken*result[0]/norm;
  result[1] = darken*result[1]/norm;
  result[2] = darken*result[2]/norm;
}

void find_crossings(GLfloat ixx, GLfloat iyy, GLfloat izz, 
		    int ix, int iy, int iz, 
		    SixTetrahedra& T6, GLfloat scale_length){

  GLfloat CubeCorners[2][2][2][dim];
  int CubeCornerIndices[2][2][2][dim];
  int ShiftCorners[2][2][2][dim];
  GLfloat LocalPos[dim];
  int LocalIndices[dim];
  
  LocalPos[0] = ixx;
  LocalPos[1] = iyy;
  LocalPos[2] = izz;
  LocalIndices[0] = ix;
  LocalIndices[1] = iy;
  LocalIndices[2] = iz;

  for(int i0 = 0; i0 < 2; ++i0){
    for(int i1 = 0; i1 < 2; ++i1){
      for(int i2 = 0; i2 < 2; ++i2){
	ShiftCorners[i0][i1][i2][0] = i0;
	ShiftCorners[i0][i1][i2][1] = i1;
	ShiftCorners[i0][i1][i2][2] = i2;
	for(int k = 0; k < dim; ++k){
	  CubeCorners[i0][i1][i2][k] = LocalPos[k] + scale_length*ShiftCorners[i0][i1][i2][k];
	  CubeCornerIndices[i0][i1][i2][k] = LocalIndices[k] + ShiftCorners[i0][i1][i2][k];
	}
      }
    }
  }
  
  for(int k = 0; k < dim; ++k){    
    /***********    Tetrahedron # 1    ************/
    T6.T_Corners[0][0][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[0][1][k] = CubeCorners[1][1][0][k];
    T6.T_Corners[0][2][k] = CubeCorners[1][1][1][k];
    T6.T_Corners[0][3][k] = CubeCorners[0][1][1][k];

    T6.T_CornerIndices[0][0][k] = CubeCornerIndices[0][0][0][k];
    T6.T_CornerIndices[0][1][k] = CubeCornerIndices[1][1][0][k];
    T6.T_CornerIndices[0][2][k] = CubeCornerIndices[1][1][1][k];
    T6.T_CornerIndices[0][3][k] = CubeCornerIndices[0][1][1][k];

    /***********    Tetrahedron # 2    ************/
    T6.T_Corners[1][0][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[1][1][k] = CubeCorners[0][1][0][k];
    T6.T_Corners[1][2][k] = CubeCorners[1][1][0][k];
    T6.T_Corners[1][3][k] = CubeCorners[0][1][1][k];

    T6.T_CornerIndices[1][0][k] = CubeCornerIndices[0][0][0][k];
    T6.T_CornerIndices[1][1][k] = CubeCornerIndices[0][1][0][k];
    T6.T_CornerIndices[1][2][k] = CubeCornerIndices[1][1][0][k];
    T6.T_CornerIndices[1][3][k] = CubeCornerIndices[0][1][1][k];

    /***********    Tetrahedron # 3    ************/
    T6.T_Corners[2][0][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[2][1][k] = CubeCorners[0][0][1][k];
    T6.T_Corners[2][2][k] = CubeCorners[1][1][1][k];
    T6.T_Corners[2][3][k] = CubeCorners[0][1][1][k];

    T6.T_CornerIndices[2][0][k] = CubeCornerIndices[0][0][0][k];
    T6.T_CornerIndices[2][1][k] = CubeCornerIndices[0][0][1][k];
    T6.T_CornerIndices[2][2][k] = CubeCornerIndices[1][1][1][k];
    T6.T_CornerIndices[2][3][k] = CubeCornerIndices[0][1][1][k];

    /***********    Tetrahedron # 4    ************/
    T6.T_Corners[3][0][k] = CubeCorners[1][0][0][k];
    T6.T_Corners[3][1][k] = CubeCorners[0][0][1][k];
    T6.T_Corners[3][2][k] = CubeCorners[1][0][1][k];
    T6.T_Corners[3][3][k] = CubeCorners[1][1][1][k];

    T6.T_CornerIndices[3][0][k] = CubeCornerIndices[1][0][0][k];
    T6.T_CornerIndices[3][1][k] = CubeCornerIndices[0][0][1][k];
    T6.T_CornerIndices[3][2][k] = CubeCornerIndices[1][0][1][k];
    T6.T_CornerIndices[3][3][k] = CubeCornerIndices[1][1][1][k];

    /***********    Tetrahedron # 5    ************/
    T6.T_Corners[4][0][k] = CubeCorners[1][0][0][k];
    T6.T_Corners[4][1][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[4][2][k] = CubeCorners[1][1][0][k];
    T6.T_Corners[4][3][k] = CubeCorners[1][1][1][k];

    T6.T_CornerIndices[4][0][k] = CubeCornerIndices[1][0][0][k];
    T6.T_CornerIndices[4][1][k] = CubeCornerIndices[0][0][0][k];
    T6.T_CornerIndices[4][2][k] = CubeCornerIndices[1][1][0][k];
    T6.T_CornerIndices[4][3][k] = CubeCornerIndices[1][1][1][k];

    /***********    Tetrahedron # 6    ************/
    T6.T_Corners[5][0][k] = CubeCorners[1][0][0][k];
    T6.T_Corners[5][1][k] = CubeCorners[0][0][1][k];
    T6.T_Corners[5][2][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[5][3][k] = CubeCorners[1][1][1][k];

    T6.T_CornerIndices[5][0][k] = CubeCornerIndices[1][0][0][k];
    T6.T_CornerIndices[5][1][k] = CubeCornerIndices[0][0][1][k];
    T6.T_CornerIndices[5][2][k] = CubeCornerIndices[0][0][0][k];
    T6.T_CornerIndices[5][3][k] = CubeCornerIndices[1][1][1][k];
  }
}

void draw_isosurface_triangles(Coords& TT, SixTetrahedra& T6){
  GLfloat v1[dim], v2[dim], n1[dim];
  Vec3D tmp_Vertex;
  VertexCoords VTT;
  VolumeLink VLink[6];
  bool CoincidentCorners[4];

  for(int t_n = 0; t_n < 6; ++t_n){
    for(int c1 = 0; c1 < 4; ++c1){
      for(int c2 = c1; c2 < 4; ++c2){

	if(TT.F[T6.T_CornerIndices[t_n][c1][0]][T6.T_CornerIndices[t_n][c1][1]][T6.T_CornerIndices[t_n][c1][2]]*
	   TT.F[T6.T_CornerIndices[t_n][c2][0]][T6.T_CornerIndices[t_n][c2][1]][T6.T_CornerIndices[t_n][c2][2]] 
	   <= 0){//Find out if the isosurface crosses the line between the two vertices. 
	  inter_point(TT,
		      T6.T_CornerIndices[t_n][c1][0], 
		      T6.T_CornerIndices[t_n][c1][1],
		      T6.T_CornerIndices[t_n][c1][2],
		      T6.T_CornerIndices[t_n][c2][0], 
		      T6.T_CornerIndices[t_n][c2][1],
		      T6.T_CornerIndices[t_n][c2][2], 
		      tmp_Vertex.V[0], 
		      tmp_Vertex.V[1], 
		      tmp_Vertex.V[2]);

	  for(int kki1 = 0; kki1 < dim; ++kki1){
	    tmp_Vertex.End1Index[kki1] = T6.T_CornerIndices[t_n][c1][kki1];
	    tmp_Vertex.End2Index[kki1] = T6.T_CornerIndices[t_n][c2][kki1];
	  }
	  VTT.push_back(tmp_Vertex);
	}
      }
    }

    if(VTT.size() == 3){//In this case the isosurface separates one node of the tetrahedron from the other three.
      for(unsigned int j0 = 0; j0 < VTT.size(); ++j0){
	for(unsigned int j1 = j0+1; j1 < VTT.size(); ++j1){
	  for(unsigned int j2 = j1+1; j2 < VTT.size(); ++j2){
	    for(int i = 0; i < dim; ++i){
	      v1[i] =  VTT[j1].V[i] - VTT[j0].V[i];
	      v2[i] =  VTT[j2].V[i] - VTT[j0].V[i];
	    }
	    vector_product(v1, v2, n1);
	    
	    glNormal3fv(n1);
	    glVertex3f(VTT[j0].V[0], VTT[j0].V[1], VTT[j0].V[2]);
	    glNormal3fv(n1);
	    glVertex3f(VTT[j1].V[0], VTT[j1].V[1], VTT[j1].V[2]);
	    glNormal3fv(n1);
	    glVertex3f(VTT[j2].V[0], VTT[j2].V[1], VTT[j2].V[2]);
	  }
	}
      }
    }
    if(VTT.size() == 4){//In this case two nodes of the tetrahedron are on one side and two on the other side of the isosurface. 
      int LinkCounter = 0;
      for(unsigned int j0 = 0; j0 < VTT.size(); ++j0){
	for(unsigned int j1 = j0+1; j1 < VTT.size(); ++j1){

	  VLink[LinkCounter].Corner1Index = j0;
	  VLink[LinkCounter].Corner2Index = j1;
	  
	  CoincidentCorners[0] = ((VTT[j0].End1Index[0] == VTT[j1].End1Index[0])&&
				  (VTT[j0].End1Index[1] == VTT[j1].End1Index[1])&&
				  (VTT[j0].End1Index[2] == VTT[j1].End1Index[2]));
	  
	  CoincidentCorners[1] = ((VTT[j0].End1Index[0] == VTT[j1].End2Index[0])&&
				  (VTT[j0].End1Index[1] == VTT[j1].End2Index[1])&&
				  (VTT[j0].End1Index[2] == VTT[j1].End2Index[2]));

	  CoincidentCorners[2] = ((VTT[j0].End2Index[0] == VTT[j1].End1Index[0])&&
				  (VTT[j0].End2Index[1] == VTT[j1].End1Index[1])&&
				  (VTT[j0].End2Index[2] == VTT[j1].End1Index[2]));
	  
	  CoincidentCorners[3] = ((VTT[j0].End2Index[0] == VTT[j1].End2Index[0])&&
				  (VTT[j0].End2Index[1] == VTT[j1].End2Index[1])&&
				  (VTT[j0].End2Index[2] == VTT[j1].End2Index[2]));

	  VLink[LinkCounter].BulkLink = !(CoincidentCorners[0]||
					  CoincidentCorners[1]||
					  CoincidentCorners[2]||
					  CoincidentCorners[3]);
	  ++LinkCounter;
	}
      }
      
      LinkCounter = 0;
      
      for(unsigned int j0 = 0; j0 < VTT.size(); ++j0){
	for(unsigned int j1 = j0+1; j1 < VTT.size(); ++j1){
	  
	  if(VLink[LinkCounter].BulkLink){
	    for(unsigned int j2 = 0; j2 < VTT.size(); ++j2){

	      if((j2 != VLink[LinkCounter].Corner1Index)&&(j2 != VLink[LinkCounter].Corner2Index)){
		for(int i = 0; i < dim; ++i){
		  v1[i] =  -VTT[j2].V[i] + VTT[VLink[LinkCounter].Corner2Index].V[i];
		  v2[i] =   VTT[j2].V[i] - VTT[VLink[LinkCounter].Corner1Index].V[i];
		}
		vector_product(v1, v2, n1);
		
		glNormal3fv(n1);
		glVertex3f(VTT[VLink[LinkCounter].Corner1Index].V[0], 
			   VTT[VLink[LinkCounter].Corner1Index].V[1], 
			   VTT[VLink[LinkCounter].Corner1Index].V[2]);
		glNormal3fv(n1);
		glVertex3f(VTT[VLink[LinkCounter].Corner2Index].V[0], 
			   VTT[VLink[LinkCounter].Corner2Index].V[1], 
			   VTT[VLink[LinkCounter].Corner2Index].V[2]);
		glNormal3fv(n1);
		glVertex3f(VTT[j2].V[0], 
			   VTT[j2].V[1], 
			   VTT[j2].V[2]);
	      }
	    }
	    j0 = VTT.size();
	    j1 = VTT.size();
	  }
	  ++LinkCounter;
	}
      }
    }
    VTT.clear();
  }
}

void create_isosurface_landscape_list(Coords& TT, GLfloat scale_length){
  
  SixTetrahedra T6;
  GLfloat mat_specular[] = { 1.0, 1.0, 0.0, 1.0 };
  
  glNewList(isosurface_landscape, GL_COMPILE);
  glColor3f(1.0, 1.0, 0.0);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glBegin(GL_TRIANGLES);
  
  for(int ix = 0; ix < N_x-1; ++ix){
    for(int iy = 0; iy < N_y-1; ++iy){
      for(int iz = 0; iz < N_z-1; ++iz){
		
	find_crossings(TT.X[ix][iy][iz], TT.Y[ix][iy][iz], TT.Z[ix][iy][iz], ix, iy, iz, T6, scale_length);
	draw_isosurface_triangles(TT, T6);
      }
    }
  }
  glEnd();
  glColor3f(1.0, 1.0, 1.0);
  glEndList();
}

void draw_volume_edges(GLfloat scale_length, GLint reduction){
  glNewList(volume_edges, GL_COMPILE);
  glColor3f(1.0, 0.0, 1.0);
  glBegin(GL_LINES);
  
  /************************************/
  glVertex3f(-scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(-scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  
  glVertex3f(scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);

  glVertex3f(-scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(-scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  
  glVertex3f(scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
    
  /************************************/
  glVertex3f(-scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(-scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  
  glVertex3f(scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);

  glVertex3f(-scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  glVertex3f(-scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  
  glVertex3f(scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);

  /************************************/
  glVertex3f(-scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  
  glVertex3f(-scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, -scale_length*reduction*N_z/2);

  glVertex3f(-scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, -scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  
  glVertex3f(-scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);
  glVertex3f(scale_length*reduction*N_x/2, scale_length*reduction*N_y/2, scale_length*reduction*N_z/2);

  glEnd();
  glEndList();
}

void landscape(void){
    glCallList(isosurface_landscape);
    glCallList(volume_edges);
}

void display(void){
    /* Called when OpenGL needs to update the display */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);/* Clears the window */
    glLoadIdentity();
    gluLookAt(10, 10, 25.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);    
    glPushMatrix();
    glRotatef(eye_x, 0, 1, 0);
    glRotatef(-eye_y, 0, 0, 1);
    landscape();
    glPopMatrix();
    glutSwapBuffers();
    glFlush();
}

void keyboard(unsigned char key, int x, int y){
  /* Called when a key is pressed */
  if(key == 27){
    exit(0); /* 27 is the excape key */
  }
  else{
    printf("You pressed %c\n", key);
  }
}

void reshape(int w, int h){
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION); /* Select the projection matrix */
  glLoadIdentity(); /* Initialize it */
  gluPerspective(30, (GLfloat)w/(GLfloat)h, 10.0, 50.0);
  glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv){ 
  Coords* const T = new Coords;
  if(argc < 5){
    std::cout << "Usage ./isosurface <datafile> <length_scale> <treshold> <reduction>" << std::endl;
    return(0);
  }
  read_data(*T, argv[1], atof(argv[2]), atof(argv[3]), atoi(argv[4]));

  glutInit(&argc, argv); /* Initializes OpenGL */
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_DEPTH| GLUT_RGBA);
  glutInitWindowSize(800, 600); /* Set the window size */
  glutInitWindowPosition(100, 100); /* Set the window position */
  glutCreateWindow("OpenGL_3D_data"); /* Create a window */
  
  init_lighting();  
  create_isosurface_landscape_list(*T, atof(argv[2]));
  draw_volume_edges(atof(argv[2]), atoi(argv[4]));
    
  glutDisplayFunc(display); /* Register the "display" function */
  glutKeyboardFunc(keyboard); /* Register the "keyboard" function */
  glutReshapeFunc(reshape); /* Register the "reshape" function */
  glutMotionFunc(mauz); /* Register the "mauz" function */
  glutMainLoop(); /* Enter the OpenGL main loop */
  return(0);
}

