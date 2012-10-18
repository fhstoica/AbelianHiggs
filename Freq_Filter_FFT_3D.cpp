#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

/*
This program performs a 1/f^p filtering of the frequency 
followed by a 3D Fast Fourier Transform. The initial Fourier 
Coefficients must be given in an external file
*/

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

const int n = 64;

typedef struct{
  double all_data[n][n][n][2];
  double data[2*n];
} FourierStruct;

FourierStruct* const tmp_fourier = new FourierStruct;


int main (int argc, char* argv[]){
  int ix, iy, iz; 
  float decay_pow;
  float damping_factor;
  int central_freq;
  
  if(argc < 4){
    std::cout << "Usage: ./Freq_Filter_FFT_3D.srl ";
    std::cout << "infile outfile central_freq power " << std::endl;
    return(0);
  }

  central_freq = atoi(argv[3]);
  decay_pow    = atof(argv[4]);

  gsl_fft_complex_wavetable* wavetable; /* Pointers to structures. */
  gsl_fft_complex_workspace* workspace;

  std::ifstream fin;
  fin.open(argv[1]);
  if (fin.is_open()){
    while(!fin.eof()){
      fin >> ix;
      fin >> iy;
      fin >> iz;
	  
      fin >> tmp_fourier->all_data[ix][iy][iz][0];
      fin >> tmp_fourier->all_data[ix][iy][iz][1];
      
      damping_factor = pow(sqrt(pow(ix-central_freq-0.5, 2) + 
				pow(iy-central_freq-0.5, 2) + 
				pow(iz-central_freq-0.5, 2)), decay_pow);
      
      tmp_fourier->all_data[ix][iy][iz][0] = tmp_fourier->all_data[ix][iy][iz][0]/damping_factor;
      tmp_fourier->all_data[ix][iy][iz][1] = tmp_fourier->all_data[ix][iy][iz][1]/damping_factor;
    }
  }
  else{
    std::cout << " Could not open " << argv[1] << std::endl;
    return(0);
  }
  fin.close();

  wavetable = gsl_fft_complex_wavetable_alloc (n);
  workspace = gsl_fft_complex_workspace_alloc (n);

  /* Pass 1; Perform N^2 one-dimensional FFTs along the Z direction, one for each point in the X-Y plane. */
  for(int i = 0; i < n; ++i){
    for (int j = 0; j < n; j++){
      for (int k = 0; k < n; k++){
	REAL(tmp_fourier->data,k) = tmp_fourier->all_data[i][j][k][0];
	IMAG(tmp_fourier->data,k) = tmp_fourier->all_data[i][j][k][1];
      }
      
      gsl_fft_complex_forward (tmp_fourier->data, 1, n, wavetable, workspace);
      
      for (int k = 0; k < n; k++){
	tmp_fourier->all_data[i][j][k][0] = REAL(tmp_fourier->data,k);
	tmp_fourier->all_data[i][j][k][1] = IMAG(tmp_fourier->data,k);
      }
    }
  }
  /* Pass 2; Perform N^2 one-dimensional FFTs along the Y direction, one for each point in the X-Z plane. */
  for(int i = 0; i < n; ++i){
    for (int k = 0; k < n; k++){
      for (int j = 0; j < n; j++){
	REAL(tmp_fourier->data,j) = tmp_fourier->all_data[i][j][k][0];
	IMAG(tmp_fourier->data,j) = tmp_fourier->all_data[i][j][k][1];
      }
      
      gsl_fft_complex_forward (tmp_fourier->data, 1, n, wavetable, workspace);
      
      for (int j = 0; j < n; j++){
	tmp_fourier->all_data[i][j][k][0] = REAL(tmp_fourier->data,j);
	tmp_fourier->all_data[i][j][k][1] = IMAG(tmp_fourier->data,j);
      }
    }
  }
  /* Pass 3; Perform N^2 one-dimensional FFTs along the X direction, one for each point in the Y-Z plane. */
  
  for (int j = 0; j < n; j++){
    for (int k = 0; k < n; k++){
      for(int i = 0; i < n; ++i){
	REAL(tmp_fourier->data,i) = tmp_fourier->all_data[i][j][k][0];
	IMAG(tmp_fourier->data,i) = tmp_fourier->all_data[i][j][k][1];
      }
      
      gsl_fft_complex_forward (tmp_fourier->data, 1, n, wavetable, workspace);
      
      for (int i = 0; i < n; i++){
	tmp_fourier->all_data[i][j][k][0] = REAL(tmp_fourier->data,i);
	tmp_fourier->all_data[i][j][k][1] = IMAG(tmp_fourier->data,i);
      }
    }
  }
  
  std::ofstream f_out;
  f_out.open(argv[2]);
  if (f_out.is_open()){
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
	for (int k = 0; k < n; k++){
	  f_out << i << " " << j << " " << k << " "
		<< tmp_fourier->all_data[i][j][k][0] << " " 
		<< tmp_fourier->all_data[i][j][k][1] << std::endl; 
	  
	}
      }
    }
  }
  else{
    std::cout << " Could not open " << argv[2] << std::endl;
    return(0);
  }

  f_out.close();
  
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
  return 0;
}
