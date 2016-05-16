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

float body[100000][7]; // data array of bodies

int main(int argc, char **argv) {

  int tmax = 0;
  float Fx_dir[N];
  float Fy_dir[N];
  float Fz_dir[N];

  if (argc != 2) {
    fprintf(stderr, "Format: %s { number of timesteps }\n", argv[0]);
    exit (-1);
  }

  tmax = atoi(argv[1]);

  // assign each body a random position
  for (int i = 0; i < N; i++) {
    body[i][MASS] =i%500?1.0:100.0;
    body[i][X_POS] = drand48() * 100.0;
    body[i][Y_POS] = drand48() * 100.0;
    body[i][Z_POS] = drand48() * 100.0;
    body[i][X_VEL] = drand48() * 100.0/body[i][MASS];
    body[i][Y_VEL] = drand48() * 100.0/body[i][MASS];
    body[i][Z_VEL] = drand48() * 100.0/body[i][MASS];
  }

  // print out initial positions in PDB format
  printf("MODEL %8d\n", 0);
  for (int i = 0; i < N; i++) {
    printf("%s%7d  %s %s %s%4d    %8.3f%8.3f%8.3f  %4.2f  %4.3f\n",
           "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i][X_POS], body[i][Y_POS], body[i][Z_POS], 1.00, 0.00);
  }
  printf("TER\nENDMDL\n");

  // step through each time step
  for (int t = 0; t < tmax; t++) {
    // force calculation

    // initialize forces to zero
    for (int i = 0; i < N; i++) {
      Fx_dir[i] = 0.0;
      Fy_dir[i] = 0.0;
      Fz_dir[i] = 0.0;
    }

    for (int x = 0; x < N; x++) {  // force on body x due to
      for (int i = 0; i < N; i++) {   // all other bodies 
	float x_diff, y_diff, z_diff;

	if (i != x) {
	  x_diff = body[i][X_POS] - body [x][X_POS];  // difference in x direction
	  y_diff = body[i][Y_POS] - body [x][Y_POS];  // difference in y direction
	  z_diff = body[i][Z_POS] - body [x][Z_POS];  // difference in z direction

	  // calculate distance (r)
	  float rr = (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
	  float r = sqrt(rr);

	  // force between bodies i and x
	  float F = G * body[i][MASS] * body[x][MASS] / r;

	  // if sufficiently far away, gravitation force
	  if (r > 10.0) {
	    Fx_dir[x] += F * x_diff / r;  // resolve forces in x and y directions
	    Fy_dir[x] += F * y_diff / r;  // and accumulate forces
	    Fz_dir[x] += F * z_diff / r;  // 
	  } else {
	    // if too close, anti-gravitational force
	    Fx_dir[x] -= F * x_diff / r;  // resolve forces in x and y directions
	    Fy_dir[x] -= F * y_diff / r;  // and accumulate forces
	    Fz_dir[x] -= F * z_diff / r;  // 
	  }
	}
      }
    }

    // update postions and velocity in array
    for (int i = 0; i < N; i++) {

        // update velocities
	body[i][X_VEL] += Fx_dir[i] * dt / body[i][MASS];
	body[i][Y_VEL] += Fy_dir[i] * dt / body[i][MASS];
	body[i][Z_VEL] += Fz_dir[i] * dt / body[i][MASS];
	// update positions
	body[i][X_POS] += body[i][X_VEL] * dt;
	body[i][Y_POS] += body[i][Y_VEL] * dt;
	body[i][Z_POS] += body[i][Z_VEL] * dt;
    }

    // print out positions in PDB format
    printf("MODEL %8d\n", t+1);
    for (int i = 0; i < N; i++) {
	printf("%s%7d  %s %s %s%4d    %8.3f%8.3f%8.3f  %4.2f  %4.3f\n",
               "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i][X_POS], body[i][Y_POS], body[i][Z_POS], 1.00, 0.00);
    }
    printf("TER\nENDMDL\n");
  }  // end of time period loop
}
