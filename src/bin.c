#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bin.h"

#define GEOM_COL 4
#define BUFFER_SIZE 1024

void read2D(const char* PATH, void* arr, size_t type, int row, int column) 
{
    FILE* bin_data = fopen(PATH, "rb"); 
    if (bin_data == NULL) 
    {
        printf("Could not read binary file.\n");
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
        printf("Could not write binary file.\n");
        exit(-1);
    }

    fwrite(arr, type, row * column, bin_data); 

    fclose(bin_data);   
}

