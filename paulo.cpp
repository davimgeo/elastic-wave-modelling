# include <cmath>
# include <omp.h>
# include <chrono>
# include <fstream>
# include <iostream>

# define FDM8E1 6.97545e-4f 
# define FDM8E2 9.57031e-3f 
# define FDM8E3 7.97526e-2f 
# define FDM8E4 1.19628906f 

# define FDM6E1 4.6875e-3f 
# define FDM6E2 6.5104e-2f 
# define FDM6E3 1.171875f 

# define FDM4E1 4.166666e-2f 
# define FDM4E2 1.125f 

# define FDM2E1 1.0f 

void export_binary_float(std::string path, float * array, int n)
{
    std::ofstream file(path, std::ios::out);
    file.write((char *) array, n * sizeof(float));
    file.close();
}

int main()
{
    int nx = 501;
    int nz = 501;
    int nt = 501;

    float dh = 10.0f;
    float dt = 1e-3f;

    float fmax = 30.0f;

    int sIdx = 250;
    int sIdz = 250;

    float vp = 1500.0f;
    float vs = 0.0f;
    float ro = 1000.0f;

    float M = vs*vs*ro;
    float L = vp*vp*ro - 2.0f*M; 
    float B = 1.0f / ro;
    
    float t0 = 2.0f * sqrtf(M_PI) / fmax;
    float fc = fmax / (3.0f * sqrtf(M_PI));
    
    float * Vx = new float[nx*nz]();
    float * Vz = new float[nx*nz]();
    float * Txx = new float[nx*nz]();
    float * Tzz = new float[nx*nz]();
    float * Txz = new float[nx*nz]();

    auto ti = std::chrono::system_clock::now();

    for (int n = 0; n < nt; n++)
    {
        float arg = M_PI*M_PI*M_PI*fc*fc*(n*dt - t0)*(n*dt - t0);
        float wavelet = (1.0f - 2.0f*arg) * expf(-arg);

        Txx[sIdz + sIdx*nz] += wavelet / (dh*dh);
        Tzz[sIdz + sIdx*nz] += wavelet / (dh*dh);

        # pragma omp parallel for
        for (int index = 0; index < nx*nz; index++)
        {
            int i = (int)(index % nz);
            int j = (int)(index / nz);

            if((i >= 0) && (i < nz-1) && (j > 0) && (j < nx)) 
            {
               
                float dTxx_dx = (FDM2E1*(Txx[i + j*nz]     - Txx[i + (j-1)*nz])) / dh;

                float dTxz_dz = (FDM2E1*(Txz[(i+1) + j*nz] - Txz[i + j*nz])) / dh;

                Vx[index] += dt*B*(dTxx_dx + dTxz_dz); 
            }

            if((i > 0) && (i < nz) && (j >= 0) && (j < nx-1)) 
            {
                
                float dTxz_dx = (FDM2E1*(Txz[i + (j+1)*nz] - Txz[i + j*nz])) / dh;

                float dTzz_dz = (FDM2E1*(Tzz[i + j*nz]     - Tzz[(i-1) + j*nz])) / dh;

                Vz[index] += dt*B*(dTxz_dx + dTzz_dz); 
            }            
        }
        
        # pragma omp parallel for
        for (int index = 0; index < nx*nz; index++)
        {
            int i = (int)(index % nz);
            int j = (int)(index / nz);

            if((i >= 0) && (i < nz-1) && (j >= 0) && (j < nx-1)) 
            {    
                
                float dVx_dx = (FDM2E1*(Vx[i + (j+1)*nz] - Vx[i + j*nz])) / dh;

                float dVz_dz = (FDM2E1*(Vz[(i+1) + j*nz] - Vz[i + j*nz])) / dh;
                                
                Txx[index] += dt*((L + 2*M)*dVx_dx + L*dVz_dz);
                Tzz[index] += dt*((L + 2*M)*dVz_dz + L*dVx_dx);                    
            }

            if((i > 0) && (i < nz) && (j > 0) && (j < nx)) 
            {
                float dVx_dz = (FDM2E1*(Vx[i + j*nz]     - Vx[(i-1) + j*nz])) / dh;

                float dVz_dx = (FDM2E1*(Vz[i + j*nz]     - Vz[i + (j-1)*nz])) / dh;

                Txz[index] += dt*M*(dVx_dz + dVz_dx);
            }
        }
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;

    export_binary_float("Vx.bin", Vx, nx*nz);
    export_binary_float("Vz.bin", Vz, nx*nz);
    export_binary_float("Txx.bin", Txx, nx*nz);
    export_binary_float("Tzz.bin", Tzz, nx*nz);
    export_binary_float("Txz.bin", Txz, nx*nz);

    return 0;
}
