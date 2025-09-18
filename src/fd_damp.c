#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cfg_par.h"
#include "fd_core.h"

static void get_damp_x(damping_t *damp, const config_t *cfg)
{
  for (int j = 0; j < cfg->nxx; j++)
  {
    if ((unsigned)(j - cfg->nb) < (cfg->nx))
    {
      damp->x[j] = 1.0f;
    }
    else if (j < cfg->nb)
    {
      int d = cfg->nb - j;
      damp->x[j] = expf(- (cfg->factor * d) * (cfg->factor * d));
    }
    else
    {
      int d = j - (cfg->nb + cfg->nx - 1);
      damp->x[j] = expf(- (cfg->factor * d) * (cfg->factor * d));
    }
  }
}

static void get_damp_z(damping_t *damp, const config_t *cfg)
{
  for (int i = 0; i < cfg->nzz; i++)
  {
    if ((unsigned)(i - cfg->nb) < (cfg->nz))
    {
      damp->z[i] = 1.0f;
    }
    else if (i < cfg->nb)
    {
      int d = cfg->nb - i;
      damp->z[i] = expf(- (cfg->factor * d) * (cfg->factor * d));
    }
    else
    {
      int d = i - (cfg->nb + cfg->nz - 1);
      damp->z[i] = expf(- (cfg->factor * d) * (cfg->factor * d));
    }
  }
}

damping_t get_damp(const config_t *cfg)
{
  damping_t damp;

  damp.x = (float *)malloc(cfg->nxx * sizeof(float));
  damp.z = (float *)malloc(cfg->nzz * sizeof(float));

  if (!damp.x || !damp.z)
  {
    perror("Could not allocate damping\n");
    exit(EXIT_FAILURE);
  }

  get_damp_x(&damp, cfg);
  get_damp_z(&damp, cfg);

  return damp;
}
