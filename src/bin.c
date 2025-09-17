#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bin.h"

#define BUFFER_SIZE 64

typedef struct
{
  float *sIdx, sIdz;
  float *rIdx, rIdz;
} geom_t;

void read_geometry(const char *path, int f_lines, int f_cols)
{
  FILE* fptr = fopen(path, "r");
  if (fptr == NULL) 
  {
      printf("Could not open the file.\n");
      exit(-1);
  }

  float* result = (float *)malloc(f_lines * f_cols * sizeof(float));

  char* buff[BUFFER_SIZE] = "";
  char ch;
  int len = 0;

  while ((ch = fgetc(fptr)) != EOF && ch != '#' && ch != ',')
  {
    buff[len++] = ch;

    buff[len] = "\0";
    printf("%s", buff);
    //result[i * f_lines + j] = atof(buff);  

  }

  // for (int i = 0; i < f_lines; i++) 
  // {
  //   for (int j = 0; j < f_cols; j++) 
  //   {
  //     printf("%f ", result[i * f_cols + j]);
  //   }
  //   printf("\n");
  // }

  fclose(fptr);
}

void read2D(const char* PATH, void* arr, size_t type, int row, int column) 
{
  FILE* bin_data = fopen(PATH, "rb"); 
  if (bin_data == NULL) 
  {
      printf("Could not open the file.\n");
      exit(-1);
  }

  fread(arr, type, row * column, bin_data); 

  fclose(bin_data);   
}

void write2D(const char* PATH, void* arr, size_t type, int row, int column) 
{
  FILE* bin_data = fopen(PATH, "wb"); 
  if (bin_data == NULL) 
  {
      printf("Could not write the file.\n");
      exit(-1);
  }

  fwrite(arr, type, row * column, bin_data); 

  fclose(bin_data);   
}

