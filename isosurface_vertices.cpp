#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <math.h>
#include <GL/glut.h>
#include <GL/glu.h>

const GLint dim = 3;
const GLint N_x = 32;
const GLint N_y = N_x;
const GLint N_z = N_x;

GLfloat eye_x = 0;
GLfloat eye_y = 0;

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

GLint wire_landscape = 1;
GLint volume_edges   = 2;

void mouse_motion(GLint x, GLint y){
  eye_x = x;
  eye_y = y;
  glutPostRedisplay();
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
		 int ix_2, int iy_2, int iz_2, 
		 GLfloat& X, GLfloat& Y, GLfloat&Z){
  GLfloat Dx, Dy, Dz;
  Dx = TT.X[ix_2][iy_2][iz_2] - TT.X[ix_1][iy_1][iz_1];
  Dy = TT.Y[ix_2][iy_2][iz_2] - TT.Y[ix_1][iy_1][iz_1];
  Dz = TT.Z[ix_2][iy_2][iz_2] - TT.Z[ix_1][iy_1][iz_1];
  X = TT.X[ix_1][iy_1][iz_1] + Dx*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1]) + fabs(TT.F[ix_2][iy_2][iz_2]));
  Y = TT.Y[ix_1][iy_1][iz_1] + Dy*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1]) + fabs(TT.F[ix_2][iy_2][iz_2]));
  Z = TT.Z[ix_1][iy_1][iz_1] + Dz*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1]) + fabs(TT.F[ix_2][iy_2][iz_2]));
}

void split_cube(GLfloat ixx, GLfloat iyy, GLfloat izz, SixTetrahedra& T6, GLfloat scale_length){
  GLfloat CubeCorners[2][2][2][dim];
  int     ShiftCorners[2][2][2][dim];
  GLfloat LocalPos[dim];

  LocalPos[0] = ixx;
  LocalPos[1] = iyy;
  LocalPos[2] = izz;

  for(int i0 = 0; i0 < 2; ++i0){
    for(int i1 = 0; i1 < 2; ++i1){
      for(int i2 = 0; i2 < 2; ++i2){
	ShiftCorners[i0][i1][i2][0] = i0;
	ShiftCorners[i0][i1][i2][1] = i1;
	ShiftCorners[i0][i1][i2][2] = i2;
	for(int k = 0; k < dim; ++k){
	  CubeCorners[i0][i1][i2][k] = LocalPos[k] + scale_length*ShiftCorners[i0][i1][i2][k];
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

    /***********    Tetrahedron # 2    ************/
    T6.T_Corners[1][0][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[1][1][k] = CubeCorners[0][1][0][k];
    T6.T_Corners[1][2][k] = CubeCorners[1][1][0][k];
    T6.T_Corners[1][3][k] = CubeCorners[0][1][1][k];

    /***********    Tetrahedron # 3    ************/
    T6.T_Corners[2][0][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[2][1][k] = CubeCorners[0][0][1][k];
    T6.T_Corners[2][2][k] = CubeCorners[1][1][1][k];
    T6.T_Corners[2][3][k] = CubeCorners[0][1][1][k];

    /***********    Tetrahedron # 4    ************/
    T6.T_Corners[3][0][k] = CubeCorners[1][0][0][k];
    T6.T_Corners[3][1][k] = CubeCorners[0][0][1][k];
    T6.T_Corners[3][2][k] = CubeCorners[1][0][1][k];
    T6.T_Corners[3][3][k] = CubeCorners[1][1][1][k];

    /***********    Tetrahedron # 5    ************/
    T6.T_Corners[4][0][k] = CubeCorners[1][0][0][k];
    T6.T_Corners[4][1][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[4][2][k] = CubeCorners[1][1][0][k];
    T6.T_Corners[4][3][k] = CubeCorners[1][1][1][k];

    /***********    Tetrahedron # 5    ************/    
    T6.T_Corners[5][0][k] = CubeCorners[1][0][0][k];
    T6.T_Corners[5][1][k] = CubeCorners[0][0][1][k];
    T6.T_Corners[5][2][k] = CubeCorners[0][0][0][k];
    T6.T_Corners[5][3][k] = CubeCorners[1][1][1][k];
  }
}

void draw_tetrahedron_edges(SixTetrahedra& T6){
  for(int t_n = 0; t_n < 6; ++t_n){
    for(int c1 = 0; c1 < 4; ++c1){
      for(int c2 = c1; c2 < 4; ++c2){
	glVertex3f(T6.T_Corners[t_n][c1][0], T6.T_Corners[t_n][c1][1], T6.T_Corners[t_n][c1][2]);
	glVertex3f(T6.T_Corners[t_n][c2][0], T6.T_Corners[t_n][c2][1], T6.T_Corners[t_n][c2][2]);
      }
    }
  } 
}

void find_crossings(GLfloat ixx, GLfloat iyy, GLfloat izz, 
		    int ix, int iy, int iz, 
		    SixTetrahedra& T6, GLfloat scale_length){

  GLfloat CubeCorners[2][2][2][dim];
  int     CubeCornerIndices[2][2][2][dim];
  int     ShiftCorners[2][2][2][dim];
  GLfloat LocalPos[dim];
  int     LocalIndices[dim];
  
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

void draw_tetrahedron_triangles(Coords& TT, SixTetrahedra& T6){
  
  GLfloat L_X, L_Y, L_Z;
  
  for(int t_n = 0; t_n < 6; ++t_n){
    for(int c1 = 0; c1 < 4; ++c1){
      for(int c2 = c1; c2 < 4; ++c2){

	if(TT.F[T6.T_CornerIndices[t_n][c1][0]][T6.T_CornerIndices[t_n][c1][1]][T6.T_CornerIndices[t_n][c1][2]]*
	   TT.F[T6.T_CornerIndices[t_n][c2][0]][T6.T_CornerIndices[t_n][c2][1]][T6.T_CornerIndices[t_n][c2][2]] 
	   <= 0){//Isosurface crosses the line connecting the two points.
	  inter_point(TT,
		      T6.T_CornerIndices[t_n][c1][0], 
		      T6.T_CornerIndices[t_n][c1][1],
		      T6.T_CornerIndices[t_n][c1][2],
		      T6.T_CornerIndices[t_n][c2][0], 
		      T6.T_CornerIndices[t_n][c2][1],
		      T6.T_CornerIndices[t_n][c2][2], L_X, L_Y, L_Z);
	  
	  glVertex3f(L_X, L_Y, L_Z);
	}
      }
    }
  } 
}

void create_vertices_landscape_list(Coords& TT, GLfloat scale_length){

  GLfloat currLCoor[dim][dim];
  Triad K;
  SixTetrahedra T6;
  
  glNewList(wire_landscape, GL_COMPILE);
  glColor3f(1.0, 1.0, 0.0);
  glBegin(GL_POINTS);
  
  for(int ix = 0; ix < N_x-1; ++ix){
    for(int iy = 0; iy < N_y-1; ++iy){
      for(int iz = 0; iz < N_z-1; ++iz){
	
	for(int kk2 = 0; kk2 < dim; ++kk2){
	  if(kk2 == 0){
	    K.LEnds[kk2][0] = ix + 1;
	    K.LEnds[kk2][1] = iy;
	    K.LEnds[kk2][2] = iz;
	  }
	  else{
	    if(kk2 == 1){
	      K.LEnds[kk2][0] = ix;
	      K.LEnds[kk2][1] = iy + 1;
	      K.LEnds[kk2][2] = iz;
	    }
	    else{
	      K.LEnds[kk2][0] = ix;
	      K.LEnds[kk2][1] = iy;
	      K.LEnds[kk2][2] = iz + 1;
	    }
	  }
	}
	/********* Determine which links are crossed by the surface **********/
	for(int kk2 = 0; kk2 < dim; ++kk2){
	  if(TT.F[ix][iy][iz]*
	     TT.F[K.LEnds[kk2][0]][K.LEnds[kk2][1]][K.LEnds[kk2][2]] <= 0){
	    inter_point(TT, ix, iy, iz, K.LEnds[kk2][0], K.LEnds[kk2][1], K.LEnds[kk2][2], 
			currLCoor[kk2][0], currLCoor[kk2][1], currLCoor[kk2][2]);
	    K.cross[kk2] = 1;
	  }
	  else{
	    K.cross[kk2] = 0;
	  }
	}
	find_crossings(TT.X[ix][iy][iz], TT.Y[ix][iy][iz], TT.Z[ix][iy][iz],ix, iy, iz, T6, scale_length);
	draw_tetrahedron_triangles(TT, T6);
	
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
  glColor3f(1.0, 1.0, 1.0);
  glEndList();
}

void landscape(void){
    glCallList(wire_landscape);
    glCallList(volume_edges);
}

void display(void){
    /* Called when OpenGL needs to update the display */
    glClear(GL_COLOR_BUFFER_BIT);/* Clears the window */
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
    std::cout << "Usage ./isosurface_vertices.srl <datafile> <length_scale> <treshold> <reduction>" << std::endl;
    return(-1);
  }
  read_data(*T, argv[1], atof(argv[2]), atof(argv[3]), atoi(argv[4]));

  glutInit(&argc, argv); /* Initializes OpenGL */
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_DEPTH| GLUT_RGBA);
  glutInitWindowSize(800, 600); /* Set the window size */
  glutInitWindowPosition(100, 100); /* Set the window position */
  glutCreateWindow("OpenGL_3D_data"); /* Create a window */

  create_vertices_landscape_list(*T, atof(argv[2]));
  draw_volume_edges(atof(argv[2]), atoi(argv[4]));
    
  glutDisplayFunc(display); /* Register the "display" function */
  glutKeyboardFunc(keyboard); /* Register the "keyboard" function */
  glutReshapeFunc(reshape); /* Register the "reshape" function */
  glutMotionFunc(mouse_motion); /* Register the "mouse_motion" function */
  glutMainLoop(); /* Enter the OpenGL main loop */
  return(0);
}

