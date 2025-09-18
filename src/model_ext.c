#include <stdio.h>
#include <stdlib.h>

#include "cfg_par.h"

void set_boundary(const config_t *cfg, model_t *m)
{
  int nb  = cfg->nb;
  int nx  = cfg->nx;
  int nz  = cfg->nz;
  int nxx = cfg->nxx;
  int nzz = cfg->nzz;

  size_t n = (size_t)nxx * (size_t)nzz;

  float *vp_ext  = (float *)calloc(n, sizeof(float));
  float *vs_ext  = (float *)calloc(n, sizeof(float));
  float *rho_ext = (float *)calloc(n, sizeof(float));

  if (!vp_ext || !vs_ext || !rho_ext) {
    fprintf(stderr, "set_boundary: calloc failed\n");
    exit(EXIT_FAILURE);
  }

  /* copy original arr into ext */
  for (int j = 0; j < nx; j++)
  {
    for (int i = 0; i < nz; i++)
    {
      vp_ext[(i + nb) + (j + nb) * nzz]  = m->vp[i + j * nz];
      vs_ext[(i + nb) + (j + nb) * nzz]  = m->vs[i + j * nz];
      rho_ext[(i + nb) + (j + nb) * nzz] = m->rho[i + j * nz];
    }
  }

  /* pad top/bottom (rows) */
  for (int j = nb; j < nx + nb; j++)
  {
    for (int i = 0; i < nb; i++)
    {
      vp_ext[i + j * nzz]  = vp_ext[nb + j * nzz];
      vs_ext[i + j * nzz]  = vs_ext[nb + j * nzz];
      rho_ext[i + j * nzz] = rho_ext[nb + j * nzz];

      vp_ext[(nz + nb + i) + j * nzz]  = vp_ext[(nz + nb - 1) + j * nzz];
      vs_ext[(nz + nb + i) + j * nzz]  = vs_ext[(nz + nb - 1) + j * nzz];
      rho_ext[(nz + nb + i) + j * nzz] = rho_ext[(nz + nb - 1) + j * nzz];
    }
  }

  /* pad left/right (columns) */
  for (int i = 0; i < nzz; i++)
  {
    for (int j = 0; j < nb; j++)
    {
      vp_ext[i + j * nzz]  = vp_ext[i + nb * nzz];
      vs_ext[i + j * nzz]  = vs_ext[i + nb * nzz];
      rho_ext[i + j * nzz] = rho_ext[i + nb * nzz];

      vp_ext[i + (nx + nb + j) * nzz]  = vp_ext[i + (nx + nb - 1) * nzz];
      vs_ext[i + (nx + nb + j) * nzz]  = vs_ext[i + (nx + nb - 1) * nzz];
      rho_ext[i + (nx + nb + j) * nzz] = rho_ext[i + (nx + nb - 1) * nzz];
    }
  }

  /* swap pointers to new arr */
  free(m->vp);  m->vp  = vp_ext;
  free(m->vs);  m->vs  = vs_ext;
  free(m->rho); m->rho = rho_ext;
}
