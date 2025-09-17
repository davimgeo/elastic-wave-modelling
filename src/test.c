#include <stdlib.h>

#include "bin.h"

#define PATH "data/input/geometry.txt"

int main(void)
{
  int f_lines = 10;
  int f_cols = 4;

  read_geometry(PATH, f_lines, f_cols);

}
