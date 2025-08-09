#ifndef DEBUG_H
#define DEBUG_H

void print_f32_2d_arr(float *arr, int row, int col);

void print_f32_1d_arr(float *arr, int size);

void f32arr_to_txt(const char* output_path, float* arr, int arr_size);

#endif // DEBUG_H