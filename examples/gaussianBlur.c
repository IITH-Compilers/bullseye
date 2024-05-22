/**
Gaussian Blur Kernel
Taken from
Pankaj Kukreja
github.com/proton0001
Indian Institute of Technology Hyderabad
*/

#include <math.h> // For M_PI and exp
#define HEIGHT 3840
#define WIDTH 2160
#define OFFSET 54
#define BOX_SIZE 59
#define SIGMA 59

int main() {
  int height = HEIGHT;
  int width = WIDTH;
  int inputImage[HEIGHT][WIDTH];
  int outputImage[HEIGHT][WIDTH];
  float sigma = 9;
  float s = 2.0 * SIGMA * SIGMA;
  int offset = (BOX_SIZE - 1) / 2;

  float sum = 0;
  float gaussianFilter[BOX_SIZE][BOX_SIZE] = {0};

#pragma scop
  for (int x = -1 * OFFSET; x <= OFFSET; x++) {
    for (int y = -1 * OFFSET; y <= OFFSET; y++) {
      gaussianFilter[x + OFFSET][y + OFFSET] =
          (exp(-(x * x + y * y) / s)) / (M_PI * s);
      sum += gaussianFilter[x + OFFSET][y + OFFSET];
    }
  }

  float sum_in_current_frame = 0;
  for (int i = OFFSET; i < HEIGHT - OFFSET; i++) {
    for (int j = OFFSET; j < WIDTH - OFFSET; j++) {
      /* Computing sum of (elements * corresponding gaussianFilter) in window
       * centered at i,j */
      sum_in_current_frame = 0;
      for (int k = -1 * OFFSET; k <= OFFSET; k++) {
        for (int l = -1 * OFFSET; l <= OFFSET; l++) {
          sum_in_current_frame += (inputImage[i + k][j + l] *
                                   gaussianFilter[k + OFFSET][l + OFFSET]);
        }
      }
      outputImage[i][j] = (sum_in_current_frame) / sum;
    }
  }
#pragma endscop
}
