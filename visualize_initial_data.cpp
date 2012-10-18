#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <GL/glut.h>
#include <GL/glu.h>

//For the surface
const GLint dim = 3;

const GLint N_x = 64;
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

GLint once = 1;
GLint point_landscape = 1;
GLint volume_edges   = 2;
bool halt_program = false;

/* This one works for smooth motions of the mouse */
void mouse(GLint x, GLint y){
  eye_x = x;
  eye_y = y;
  glutPostRedisplay();
}


void read_data(Coords& TT, char* infile, GLfloat scale_length, GLfloat treshold){
  std::ifstream fin;
  GLfloat T1, T2;
  int ix, iy, iz;
  fin.open(infile);
  if (fin.is_open()){
    while(!fin.eof()){
      
      fin >> ix;
      fin >> iy;
      fin >> iz;
      fin >> T1;
      fin >> T2;

      TT.X[ix][iy][iz] = ix;
      TT.Y[ix][iy][iz] = iy;
      TT.Z[ix][iy][iz] = iz;

      TT.F[ix][iy][iz] = sqrt(pow(T1, 2) + pow(T2, 2));
      
      TT.X[ix][iy][iz] = scale_length*(TT.X[ix][iy][iz] - N_x/2);
      TT.Y[ix][iy][iz] = scale_length*(TT.Y[ix][iy][iz] - N_y/2);
      TT.Z[ix][iy][iz] = scale_length*(TT.Z[ix][iy][iz] - N_x/2);
      
      TT.F[ix][iy][iz] = TT.F[ix][iy][iz] - treshold;
    
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
		 int ix_2, int iy_2, int iz_2, float& X, float& Y, float&Z){
  float Dx, Dy, Dz;
  Dx = TT.X[ix_2][iy_2][iz_2] - TT.X[ix_1][iy_1][iz_1];
  Dy = TT.Y[ix_2][iy_2][iz_2] - TT.Y[ix_1][iy_1][iz_1];
  Dz = TT.Z[ix_2][iy_2][iz_2] - TT.Z[ix_1][iy_1][iz_1];
  X = TT.X[ix_1][iy_1][iz_1] + 
    Dx*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1])+
				     fabs(TT.F[ix_2][iy_2][iz_2]));
  Y = TT.Y[ix_1][iy_1][iz_1] + 
    Dy*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1])+
				     fabs(TT.F[ix_2][iy_2][iz_2]));
  Z = TT.Z[ix_1][iy_1][iz_1] + 
    Dz*fabs(TT.F[ix_1][iy_1][iz_1])/(fabs(TT.F[ix_1][iy_1][iz_1])+
				     fabs(TT.F[ix_2][iy_2][iz_2]));
}

void create_point_landscape_list(Coords& TT){
 
  float currLCoor[dim][dim];
  Triad K;
  
  glNewList(point_landscape, GL_COMPILE);
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
	/********* Link indiscriminately all crossing points inside a cube *********/
	
	for(int kk3 = 0; kk3 < dim; ++kk3){
	  if(K.cross[kk3] == 1){
	    glVertex3f(currLCoor[kk3][0], currLCoor[kk3][1], currLCoor[kk3][2]);
	  }
	}
	/****************************************************************/
      }
    }
  }
  glEnd();
  glColor3f(1.0, 1.0, 1.0);
  glEndList();
}

void draw_volume_edges(GLfloat scale_length){
  glNewList(volume_edges, GL_COMPILE);
  glColor3f(1.0, 0.0, 1.0);
  glBegin(GL_LINES);
  
  /************************************/
  glVertex3f(-scale_length*N_x/2, -scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(-scale_length*N_x/2, -scale_length*N_y/2, scale_length*N_z/2);
  
  glVertex3f(scale_length*N_x/2, -scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, -scale_length*N_y/2, scale_length*N_z/2);

  glVertex3f(-scale_length*N_x/2, scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(-scale_length*N_x/2, scale_length*N_y/2, scale_length*N_z/2);
  
  glVertex3f(scale_length*N_x/2, scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, scale_length*N_y/2, scale_length*N_z/2);
    
  /************************************/
  glVertex3f(-scale_length*N_x/2, -scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(-scale_length*N_x/2, scale_length*N_y/2, -scale_length*N_z/2);
  
  glVertex3f(scale_length*N_x/2, -scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, scale_length*N_y/2, -scale_length*N_z/2);

  glVertex3f(-scale_length*N_x/2, -scale_length*N_y/2, scale_length*N_z/2);
  glVertex3f(-scale_length*N_x/2, scale_length*N_y/2, scale_length*N_z/2);
  
  glVertex3f(scale_length*N_x/2, -scale_length*N_y/2, scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, scale_length*N_y/2, scale_length*N_z/2);

  /************************************/
  glVertex3f(-scale_length*N_x/2, -scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, -scale_length*N_y/2, -scale_length*N_z/2);
  
  glVertex3f(-scale_length*N_x/2, scale_length*N_y/2, -scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, scale_length*N_y/2, -scale_length*N_z/2);

  glVertex3f(-scale_length*N_x/2, -scale_length*N_y/2, scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, -scale_length*N_y/2, scale_length*N_z/2);
  
  glVertex3f(-scale_length*N_x/2, scale_length*N_y/2, scale_length*N_z/2);
  glVertex3f(scale_length*N_x/2, scale_length*N_y/2, scale_length*N_z/2);

  glEnd();
  glColor3f(1.0, 1.0, 1.0);
  glEndList();
  
}

void landscape(void){
    glCallList(point_landscape);
    glCallList(volume_edges);
}

void display(void){
    /* Called when OpenGL needs to update the display */
    //glClearColor(1.0, 1.0, 1.0, 0.0);
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
  if(argc < 4){
    std::cout << "Too few arguments: usage " 
	      << "./view_initial_data.srl <datafile> <length_scale> <treshold>" << std::endl;
    return(0);
  }
  read_data(*T, argv[1], atof(argv[2]), atof(argv[3]));
  if(halt_program){
    return(0);
  }
  glutInit(&argc, argv); /* Initializes OpenGL */
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_DEPTH| GLUT_RGBA);
  glutInitWindowSize(800, 600); /* Set the window size */
  glutInitWindowPosition(100, 100); /* Set the window position */
  glutCreateWindow("OpenGL_3D_initial_data"); /* Create a window */
  create_point_landscape_list(*T);
  draw_volume_edges(atof(argv[2]));
  glutDisplayFunc(display); /* Register the "display" function */
  glutKeyboardFunc(keyboard); /* Register the "keyboard" function */
  glutReshapeFunc(reshape); /* Register the "reshape" function */
  glutMotionFunc(mouse); /* Register the "mauz" function */
  glutMainLoop(); /* Enter the OpenGL main loop */
  return(0);
}

