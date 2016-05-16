//this is a lite version of a GPU accelerated N-body simulation. Has to be run on an NVIDIA machine with CUDA enabled.
//the interaction is just gravitation
//the simulation trajectory is to be visualized in VMD
#include <cstdio>
#include <cstdlib>
#include <cmath>

#define N 9999  // number of bodies
#define MASS 0  // row in array for mass
#define X_POS 1 // row in array for x position
#define Y_POS 2 // row in array for y position
#define Z_POS 3 // row in array for z position
#define X_VEL 4 // row in array for x velocity
#define Y_VEL 5 // row in array for y velocity
#define Z_VEL 6 // row in array for z velocity
#define G 10    // "gravitational constant" (not really)

float dt = 0.05; // time interval

// each thread computes new position of one body
__global__ void nbody(float *dev_body, float dt) {
int i=threadIdx.x + blockIdx.x* blockDim.x;
int j;
if(i<N) 
{    // force calculation
	float Fx_dir;
 	float Fy_dir;
	float Fz_dir;
    // initialize forces to zero
      Fx_dir = 0.0;
      Fy_dir = 0.0;
      Fz_dir = 0.0; 

 	for(j=0;j<N&&j!=i;j++)
 	{ 
// force on body x due to all other bodies 
	float x_diff, y_diff, z_diff;

	  x_diff = dev_body[i*7+X_POS] - dev_body [j*7+X_POS];  // difference in x direction
	  y_diff = dev_body[i*7+Y_POS] - dev_body [j*7+Y_POS];  // difference in y direction
	  z_diff = dev_body[i*7+Z_POS] - dev_body [j*7+Z_POS];  // difference in z direction

	  // calculate distance (r)
	  float rr = (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
	  float r = sqrt(rr);

	  // force between bodies i and x
	  float F = G * dev_body[i*7+MASS] * dev_body[j*7+MASS] / r;

	  // if sufficiently far away, gravitation force
	  if (r > 10.0) {
	    Fx_dir += -F * x_diff / r;  // resolve forces in x and y directions
	    Fy_dir += -F * y_diff / r;  // and accumulate forces
	    Fz_dir += -F * z_diff / r;  // 
	  } 
	  else if(r<=10.0&&r>0.5){//avoid extremely large acceleration due to long time interval
	    // if too close, anti-gravitational force
	    Fx_dir -= -F * x_diff / r;  // resolve forces in x and y directions
	    Fy_dir -= -F * y_diff / r;  // and accumulate forces
	    Fz_dir -= -F * z_diff / r;  // 
	  }
	  }


    // update postions and velocity in array
        
	// update velocities
	dev_body[i*7+X_VEL] += Fx_dir * dt / dev_body[i*7+MASS];
	dev_body[i*7+Y_VEL] += Fy_dir * dt / dev_body[i*7+MASS];
	dev_body[i*7+Z_VEL] += Fz_dir * dt / dev_body[i*7+MASS];
	// update positions
	dev_body[i*7+X_POS] += dev_body[i*7+X_VEL] * dt;
	dev_body[i*7+Y_POS] += dev_body[i*7+Y_VEL] * dt;
	dev_body[i*7+Z_POS] += dev_body[i*7+Z_VEL] * dt;
	}
 
}


int main(int argc, char **argv) {
  float *body; // host data array of bodies
  float *dev_body; // device data array of bodies

  int tmax = 0;

  if (argc != 2) {
    fprintf(stderr, "Format: %s { number of timesteps }\n", argv[0]);
    exit (-1);
  }

  tmax = atoi(argv[1]);

  // allocate memory size for the body
  int bodysize = N * 7 * sizeof(float);    
  body = (float *)malloc(bodysize);
  cudaMalloc((void**) &dev_body, bodysize);

  // assign each body a random position
  for (int i = 0; i < N; i++) {
    body[i * 7 + MASS] =  i%1001?1:1000;//create several heavy regions
    body[i * 7 + X_POS] = (i%2?0.0:200.0)+drand48() * 50.0;//define two galaxis
    body[i * 7 + Y_POS] = drand48() * 50.0;
    body[i * 7 + Z_POS] = drand48() * 20.0;//a plate-like distribution
    body[i * 7 + X_VEL] = drand48() * 0.1/body[i * 7 + MASS];
    body[i * 7 + Y_VEL] = (i%2?-10:10)+drand48() * 1.0/body[i * 7 + MASS];//angular momentum
    body[i * 7 + Z_VEL] = drand48() * 0.1/body[i * 7 + MASS];
  }

  // print out initial positions in PDB format

  printf("MODEL %8d\n", 0);
  for (int i = 0; i < N; i++) {
	printf("%s%7d  %s %s %s%4d     %7.0f %7.0f %7.0f  %4.2f  %4.3f\n",
           "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i * 7 + X_POS], body[i * 7 + Y_POS], body[i * 7 + Z_POS], 1.00, 0.00);
  }
  printf("TER\nENDMDL\n");
 // copy nbody info over to GPU
    cudaMemcpy(dev_body, body, bodysize, cudaMemcpyHostToDevice);

  // step through each time step
  for (int t = 0; t < tmax; t++) {
 
    dim3 blockDim(1024);
    dim3 gridDim((int)ceil(N*1.0 / blockDim.x));

    // run nbody calculation
    nbody<<<gridDim, blockDim>>>(dev_body, dt);
    cudaThreadSynchronize();
if(!(t%1))//change output frequency by the mod factor. help to determine the robusticity quickly
{
    // copy nbody info back to CPU
    cudaMemcpy(body, dev_body, bodysize, cudaMemcpyDeviceToHost);

    // print out positions in PDB format
    printf("MODEL %8d\n", t+1);
    for (int i = 0; i < N; i++) {
	printf("%s%7d  %s %s %s%4d     %7.0f %7.0f %7.0f  %4.2f  %4.3f\n",
               "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i * 7 + X_POS], body[i * 7 + Y_POS], body[i * 7 + Z_POS], 1.00, 0.00);
    }
    printf("TER\nENDMDL\n");
   }

  }  // end of time period loop
  free(body);
  cudaFree(body);
}
