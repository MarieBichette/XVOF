#include <stdlib.h>

typedef struct MieGruneisenParameters MieGruneisenParameters_t;

struct MieGruneisenParameters {
        double c_zero;
        double s1;
        double s2;
        double s3;
        double rho_zero;
        double gamma_zero;
        double coeff_b;
        double e_zero;
        void (*solve)(MieGruneisenParameters_t*, const double, const double, double *, double *, double *);
};

void solveVolumeEnergy(MieGruneisenParameters_t*, const double specific_volume, const double internal_energy, double *pressure, double *gamma_per_vol, double *c_son);

 
