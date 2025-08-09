#include <stdlib.h>

#include "debug/debug.h"

#include "src/par.h"
#include "src/read_bin.h"
#include "src/FD.h"

int main(int argc, char** argv)
{
    modelPar mpar = {
       "data/model_vp_2d_1150x648.bin",
        1150,
        648
    };

    waveletPar wpar = {
        1001,
        0.004f,
        10
    };

    float *ricker_wavelet = ricker(wpar.nt, wpar.dt, wpar.fmax);

    f32arr_to_txt("data/ricker_nt1001_dt4e-3.txt", ricker_wavelet, wpar.nt);

    float *model = read_f32_bin_model(mpar.path, mpar.nx, mpar.nz);

    // print_f32_2d_arr(model, msize.nx, msize.nz);

    free(model);
    return 0;
}

