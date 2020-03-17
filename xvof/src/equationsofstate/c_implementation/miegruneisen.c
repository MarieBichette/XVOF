#include <math.h>
#include <stdio.h>
#include <signal.h>
#include "miegruneisen.h"

void solveVolumeEnergy(MieGruneisenParameters_t *params, const double specific_volume, const double internal_energy,
                      double* pressure, double* gamma_per_vol, double* c_son)
{
        // phi : pression sur l hugoniot
        // einth : energie interne specifique sur l hugoniot
        // dpdv : dp/dv
        double c_zero_2 = 0., dgam = 0., epsv = 0., epsv2 = 0., redond_a = 0., phi = 0., einth = 0., dpdv = 0.;
        double denom = 0., dphi = 0., deinth = 0., vson_2 = 0.;
        // Carré de la vitesse du son initiale
        c_zero_2 = params->c_zero * params->c_zero;
        // Dérivee du coefficient de gruneisen
        dgam = params->rho_zero * (params->gamma_zero - params->coeff_b);
        epsv = 1.0 - params->rho_zero * specific_volume;
        epsv2 = epsv * epsv;
        // Coefficient de gruneisen
        *gamma_per_vol = (params->gamma_zero * (1.0 - epsv) + params->coeff_b * epsv) / specific_volume;
        redond_a = params->s1 + 2. * params->s2 * epsv + 3. * params->s3 * epsv2;
        if (epsv > 0) {
            // ============================================================
            // si epsv > 0, la pression depend de einth et phi.
            // denom : racine du denominateur de phi
            // dphi : derivee de ph par rapport a
            // deinth : derivee de einth par rapport a v
            // ============================================================
            denom = (1. - (params->s1 + params->s2 * epsv + params->s3 * epsv2) * epsv);
            phi = params->rho_zero * c_zero_2 * epsv / (denom * denom); 
            einth = params->e_zero + phi * epsv / (2. * params->rho_zero);
            //
            dphi = phi * params->rho_zero * (-1. / epsv - 2. * redond_a / denom);
            //
            deinth = phi * (-1. - epsv * redond_a / denom);
            //
            dpdv = dphi + (dgam - *gamma_per_vol) *  (internal_energy - einth) / specific_volume - *gamma_per_vol * deinth;
        } else {
            // ============================================================
            // traitement en tension : epsv < 0
            // la pression depend d une fonction de v : phi
            // et
            // de e0 appelee einth
            // ============================================================
            phi = params->rho_zero * c_zero_2 * epsv / (1. - epsv);
            // einth ---> e0
            einth = params->e_zero;
            //
            dphi = -c_zero_2 / (specific_volume * specific_volume);
            //
            dpdv = dphi + (dgam - *gamma_per_vol) * (internal_energy - einth) / specific_volume;
        }
        // ****************************
        // Pression :
        // ****************************
        *pressure = phi + (*gamma_per_vol) * (internal_energy - einth);
        // ======================================
        // Carre de la vitesse du son :
        // ======================================
        vson_2 = specific_volume * specific_volume * (*pressure * *gamma_per_vol - dpdv);
        if (vson_2 < 0) {
            printf("Carré de la vitesse du son < 0\n");
            printf("specific_volume = %15.9g\n", specific_volume);
            printf("pressure = %15.9g\n", *pressure);
            printf("dpsurde = %15.9g\n", *gamma_per_vol);
            printf("dpdv = %15.9g\n", dpdv);
            raise(SIGABRT);
        }
        *c_son = sqrt(vson_2);
        //
        if ((*c_son <= 0.) && (*c_son < 10000.)) {
            *c_son = 0.;
        }
}

//int main() 
//{
//        double czero = 3980.0, S1 = 1.58, S2 = 0., S3 = 0., rhozero = 8129.0, grunzero = 1.6, b = 0.5, ezero = 0.;
//        MieGruneisenParameters_t params= {czero, S1, S2, S3, rhozero, grunzero, b, ezero, solveVolumeEnergy};
//        double rho = 9.000000000001e+03, e_int = 1.0e+04;
//        double pression = 0., gruneisen = 0., vitesse_du_son = 0.;
//        params.solve(&params, 1./rho, e_int, &pression, &gruneisen, &vitesse_du_son);
//        printf("Pression = %15.9g \nGruneisen = %15.9g \nVitesse du son = %15.9g\n", pression, gruneisen, vitesse_du_son);
//}
