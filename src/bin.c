#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bin.h"

#define BUFFER_SIZE 64

void read_geometry(const char *path, int f_lines, int f_cols, geom_t *geom)
{
  FILE* fptr = fopen(path, "r");
  if (fptr == NULL) 
  {
      printf("Could not open the file.\n");
      exit(-1);
  }

  int total_size = f_lines * f_cols;

  float* result = (float *)malloc(total_size * sizeof(float));

  char buff[BUFFER_SIZE];
  int len = 0, idx = 0;

  char ch;
  while ((ch = fgetc(fptr)) != EOF && idx < total_size) 
  {
    if (ch == '#') 
    {
      while (ch != '\n' && ch != EOF) 
        ch = fgetc(fptr);
      continue;
    }

    if (ch == ' ' || ch == '\n' || ch == ',') 
      continue;

    len = 0;
    while (ch != ',' && ch != '\n' && 
        ch != EOF && ch != ' ') 
    {
        buff[len++] = ch;
        ch = fgetc(fptr);
    }
    buff[len] = '\0';

    result[idx++] = atof(buff);
  }

  // rec pos
  int col1 = 0;
  for (int i = 0; i < f_lines; i++) 
  {
    geom->rIdx[i] = result[i * f_lines];
  }

  // rec depth
  int col2 = 1;
  for (int i = 0; i < f_lines; i++) 
  {
    geom->rIdz[i] = result[i * f_lines + col2];
  }

  // src pos
  int col3 = 2;
  for (int i = 0; i < f_lines; i++) 
  {
    geom->sIdx[i] = result[i * f_lines + col3];
  }

  // src pos
  int col4 = 3;
  for (int i = 0; i < f_lines; i++) 
  {
    geom->sIdz[i] = result[i * f_lines + col4];
  }

  for (int i = 0; i < f_lines; i++) 
  {
    for (int j = 0; j < f_cols; j++) 
    {
      printf("%f ", result[i * f_cols + j]);
    }
    printf("\n");
  }

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

