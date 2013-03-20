const int dim      = 3; //Number of space dimensions
const int scal_dim = 2; /* Number of components of the 
			   complex scalar field. */
const unsigned int NoOfThreads = 4;//Make sure it is not larger than the number of CPUs when setting the affinity!

const int N_x = 64;  //Number of intervals X-Axis. Easier if it is a power of 2
const int N_y = N_x; //Number of intervals Y-Axis. Easier if it is a power of 2
const int N_z = N_x; //Number of intervals Z-Axis. Easier if it is a power of 2

unsigned long int reduction = 100;/* Write out the data only once every "reduction" timesteps. */
unsigned long int grid_reduction = 1; /*For very large lattices write out only one out 
					of "grid_reduction" points along each direction */

const unsigned long int frames_total = 20;
unsigned long int Steps = reduction*frames_total; //Total number of steps

const float dt = 1.0e-2; //Time  Interval

const float a = 0.45; /* Lattice Spacing  */ 

float scale_factor  =  0.2;//Overall scaling of the initial fluctuations. 
const int E_params = 6;



/* Parameters of \phi^4 lagrangian. */
const float g = 1.0; //Gauge coupling constant.
float sigma   = 1.0; //V.E.V of the scalar field
float lambda  = 1.0; //Strength of the potential
float H       = 0.4; //Damping factor.

bool Lock_Gauge_Field = false;//true; /* Turn off the gauge fields to get global strings with long-range interactions. */
bool frames_yes       = true;
bool scalars_yes      = true;
bool link_yes         = false;//true;
bool gauss_yes        = false;//true;
bool energy_yes       = true;

bool Debug_ON         = false;//true;
