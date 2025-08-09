#include <stdio.h>

void print_f32_2d_arr(float *arr, int row, int col)               
{                                              
    for (int i = 0; i < row; i++) {               
        for (int j = 0; j < col; j++) {           
            printf("%f ", arr[i + j*col]);        
        }                                         
        printf("\n");                             
    }                                             
} 

void print_f32_1d_arr(float *arr, int size)
{
    for (int i = 0; i < size; i++) {
        printf("%f ", arr[i]);
    }
}

void f32arr_to_txt(const char* output_path, float* arr, int arr_size) 
{
    FILE *fptr = fopen(output_path, "w");
    if (fptr == NULL) {
        return;
    }

    for (int i = 0; i < arr_size; i++) {
        fprintf(fptr, "%f\n", arr[i]);  
    }

    fclose(fptr);
}