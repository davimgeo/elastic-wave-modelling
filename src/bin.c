#include <stdio.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>

#include "bin.h"

#define GEOM_COL 4
#define BUFFER_SIZE 1024

#define PATTERN_COMMA "^[0-9]+(\\.[0-9]+)?(,[0-9]+(\\.[0-9]+)?)*$"

void read_comma_separed_file(const char *path, int nx)
{
  FILE* fptr = fopen(path, "r");

  regex_t regex;
  int value;

  value = regcomp(&regex, PATTERN_COMMA, REG_EXTENDED);
  if (value) 
  {
    fprintf(stderr, "Could not compile regex\n");
    exit(1);
  }

  char msgbuf[BUFFER_SIZE];

  while (fgets(msgbuf, BUFFER_SIZE, fptr) != NULL) 
  {
    if (msgbuf[0] != '#') 
    {
      msgbuf[strcspn(msgbuf, "\n")] = 0;

      value = regexec(&regex, msgbuf, 0, NULL, 0);
      if (!value) 
      {
        puts("Match");
      } 
      else 
      {
        char errbuf[BUFFER_SIZE];
        regerror(value, &regex, errbuf, sizeof(errbuf));
        fprintf(stderr, "Regex match failed: %s\n", errbuf);
        exit(1);
      }
    }
  }

  regfree(&regex);
  fclose(fptr);
}

float* read_f32_bin_model(const char *path, int nx, int nz)
{
  FILE *model_file = fopen(path, "rb"); 
  if (model_file == NULL) 
  {
    printf("Model file could not be opened\n");
    return NULL;
  }

  size_t n = nx * nz;
  float *model = (float *)calloc(n, sizeof(float));

  for (int i = 0; i < nx; i++) 
  {
    for (int j = 0; j < nz; j++) 
    {
      if (fread(&model[i * nz + j], sizeof(float), 1, model_file) != 1) 
      {
        printf("Error reading from file at (%d,%d)\n", i, j);
        free(model);
        fclose(model_file);
        return NULL;
      }
    }
  }

  fclose(model_file);  
  return model;
}

void write_f32_bin_model(const char *path, float *model, int nx, int nz)
{
  FILE *model_file = fopen(path, "wb");
  if (model_file == NULL) 
  {
      printf("Model file could not be created\n");
      return;
  }

  for (int i = 0; i < nx; i++) 
  {
      for (int j = 0; j < nz; j++) 
      {
        if (fwrite(&model[i * nz + j], sizeof(float), 1, model_file) != 1) 
        {
          printf("Error writing to file at (%d,%d)\n", i, j);
          fclose(model_file);
          return;
        }
      }
  }

  fclose(model_file);
}

