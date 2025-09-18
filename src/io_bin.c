#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_bin.h"

#define BUFFER_SIZE 64

static void read_file_separed_by_comma(FILE* fptr, int total_size, float* result)
{
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
}

void read_receivers(const char *PATH, int f_lines, int f_cols, float *rIdx, float *rIdz)
{
  FILE* fptr = fopen(PATH, "r");
  if (fptr == NULL) 
  {
      printf("Could not open the file.\n");
      exit(-1);
  }

  int total_size = f_lines * f_cols;

  float* result = (float *)malloc(total_size * sizeof(float));

  read_file_separed_by_comma(fptr, total_size, result);

  int col1 = 0, col2 = 1;
  for (int i = 0; i < f_lines; i++) 
  {
    rIdx[i] = result[i * f_cols + col1];
    rIdz[i] = result[i * f_cols + col2];
  }

  fclose(fptr);
}

void read_sources(const char *PATH, int f_lines, int f_cols, float *sIdx, float *sIdz)
{
  FILE* fptr = fopen(PATH, "r");
  if (fptr == NULL) 
  {
      printf("Could not open the file.\n");
      exit(-1);
  }

  int total_size = f_lines * f_cols;

  float* result = (float *)malloc(total_size * sizeof(float));

  read_file_separed_by_comma(fptr, total_size, result);

  int col1 = 0, col2 = 1;
  for (int i = 0; i < f_lines; i++) 
  {
    sIdx[i] = result[i * f_cols + col1];
    sIdz[i] = result[i * f_cols + col2];
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
