/**
 * MaxFile name: MovingAverageSimple
 * Using fourth order runge kutta method to calculate next point of y value 
 * Calcualte next to next for 16 y points given equation 3y^2 + 2
 * Compare actual value and the expected value of output
 */

#include "Maxfiles.h"
#include <MaxSLiCInterface.h>

// initial values
const int size = 16;
float h = 2.0;
float x = 1.0;
const float k[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const float y[16] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
float y_output[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
float expected_y_output[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// implement derivative method to calculate the dervative function for every y value
// dy/dx = 3y^2 + 2
static void calc_derivation(float x, float *y, float *k, int size){
  for(int j = 0 ; j < size ; j++){
    k[j] = 3.0 * y[j] * y[j] + 2.0;
  }  
}

// implement 4th order runge kutta method to calculate the actual expected output values
static void calc_order4_runge_kutta(const float y[], const float k[], int size, 
                                    float x, float h, float y_output[]){
  
  float h_h, h_six, x_h, *y_t, *dy_t, *dy_m;
  
  h_h = h * 0.5;
  h_six = h / 6.0;
  x_h = x + h_h;
  
  // add vector 
  y_t = malloc(sizeof(void*) * size);
  dy_t = malloc(sizeof(void*) * size);
  dy_m = malloc(sizeof(void*) * size);
  
  for(int j = 0 ; j < size ; j++){
    y_t[j] = y[j] + h_h * k[j];             //yt[i]=y[i]+hh*dydx[i]
  }
  calc_derivation(x_h, y_t, dy_t, size);

  for(int j = 0 ; j < size ; j++){
    y_t[j] = y[j] + h_h * dy_t[j];          //yt[i]=y[i]+hh*dyt[i] 
  }
  calc_derivation(x_h, y_t, dy_m, size);
    
  for(int j = 0 ; j < size ; j++){
    y_t[j]  = y[j] + h * dy_m[j];           //yt[i]=y[i]+h*dym[i]
    dy_m[j] = dy_m[j] + dy_t[j];            //dym[i] += dyt[i]
  }
  calc_derivation(x + h, y_t, dy_t, size);

  for(int j = 0 ; j < size ; j++){
    // yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i])
    y_output[j] = y[j] + h_six * (k[j] + dy_t[j] + 2.0 * dy_m[j]);
  }
}

int main()
{
    printf("Running DFE Fourth Order Runge Kutta Method.\n");
    
    // get actual y values using fourth order runge kutta method
	MovingAverageSimple(size, h, x, k, y, y_output);
	
	// get expected y values using fourth order runge kutta method
	calc_order4_runge_kutta(y, k, size, x, h, expected_y_output);
	
	// compare expected and actual output of y values
	for(int j=0; j < size; j++)
	{
		if(y_output[j] != expected_y_output[j]) {
			fprintf(stderr, "@ point %d: output y is %1.8g; expected output y is %1.8g. small 7th decimal error. \n",
				    j, y_output[j], expected_y_output[j]);
		}else {
		    fprintf(stderr, "@ point %d: output y is %1.8g; expected output y is %1.8g. no error. \n",
				    j, y_output[j], expected_y_output[j]); 
		}
	}
	return 0;
}

