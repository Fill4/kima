#include <cmath>

#include "elliptic.h"

// TODO: Get this out of here
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

template <typename T>
inline T quadratic_ld(const T& b, const T& p, const T& z, const T& c1, const T& c2) {

    const T omega = 1.0 - c1/3.0 - c2/6.0;
    const T tol = 1.0e-14; //< double precision equality tolerance for corner case issues

    T kap0 = 0.0, kap1 = 0.0;
    T lambdad=0.0, lambdae=0.0, etad=0.0;

    // allow for negative impact parameters. use d instead of b to keep the same
    double d = fabs(b);

    // check the corner cases
    if(fabs(p - d) < tol) {
        d = p;
    }
    if(fabs(p - 1.0 - d) < tol) {
        d = p - 1.0;
    }
    if(fabs(1.0 - p - d) < tol) {
        d = 1.0 - p;
    }
    if(d < tol) {
        d = 0.0;
    }

    double x1 = pow((p - d), 2.0);
    double x2 = pow((p + d), 2.0);
    double x3 = p*p - d*d;

    //< source is unocculted
    if (d >= 1.0 + p) {
        return 1.0;
    }
    //< source is completely occulted
    if (p >= 1.0 && d <= p - 1.0) {
        lambdad = 0.0;
        etad = 0.5; // error in Fortran code corrected here, following Jason Eastman's python code
        lambdae = 1.0;
        return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*(lambdad + 2.0/3.0) + c2*etad)/omega;
    }
    //< source is partly occulted and occulting object crosses the limb
    if (d >= fabs(1.0 - p) && d <= 1.0 + p) {
        kap1 = acos(MIN((1.0 - p*p + d*d)/2.0/d, 1.0));
        kap0 = acos(MIN((p*p + d*d - 1.0)/2.0/p/d, 1.0));
        lambdae = p*p*kap0 + kap1;
        lambdae = (lambdae - 0.50*sqrt(MAX(4.0*d*d - pow((1.0 + d*d - p*p), 2.0), 0.0)))/M_PI;
    }

    //< edge of the occulting star lies at the origin
    if (d == p) {
        // zone 5
        if (d < 0.5) {
            // zone 5.2
            T q = 2.0 * p;
            T Kk = K(q);
            T Ek = E(q);
            lambdad = 1.0/3.0 + 2.0/9.0/M_PI*(4.0*(2.0*p*p - 1.0)*Ek + (1.0 - 4.0*p*p)*Kk);
            etad = p*p/2.0*(p*p + 2.0*d*d);
            return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
        }
        else if (d > 0.5) {
            // zone 5.1
            T q = 0.5 / p;
            T Kk = K(q);
            T Ek = E(q);
            lambdad = 1.0/3.0 + 16.0*p/9.0/M_PI*(2.0*p*p - 1.0)*Ek -  \
                        (32.0*pow(p, 4.0) - 20.0*p*p + 3.0)/9.0/M_PI/p*Kk;
            etad = 1.0/2.0/M_PI*(kap1 + p*p*(p*p + 2.0*d*d)*kap0 -  \
                                (1.0 + 5.0*p*p + d*d)/4.0*sqrt((1.0 - x1)*(x2 - 1.0)));
        }
        else {
            // zone 6
            lambdad = 1.0/3.0 - 4.0/M_PI/9.0;
            etad = 3.0/32.0;
            return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
        }

        return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
    }

    // occulting star partly occults the source and crosses the limb
    //if((d > 0.5 + fabs(p  - 0.5) && d < 1.0 + p) || (p > 0.5 && d > fabs(1.0 - p)*1.0001 && d < p))  
    // the factor of 1.0001 is from the Mandel/Agol Fortran routine, but gave bad output for d near fabs(1-p)
    if ((d > 0.5 + fabs(p  - 0.5) && d < 1.0 + p) || (p > 0.5 && d > fabs(1.0 - p) && d < p)) {
        // zone 3.1
        T q = sqrt((1.0 - x1)/4.0/d/p);
        T Kk = K(q);
        T Ek = E(q);
        T n = 1.0/x1 - 1.0;
        T Pk = ell_bulirsch(n, q);
        lambdad = 1.0/9.0/M_PI/sqrt(p*d)*(((1.0 - x2)*(2.0*x2 +  \
                x1 - 3.0) - 3.0*x3*(x2 - 2.0))*Kk + 4.0*p*d*(d*d +  \
                7.0*p*p - 4.0)*Ek - 3.0*x3/x1*Pk);
        if(d < p) lambdad += 2.0/3.0;
        etad = 1.0/2.0/M_PI*(kap1 + p*p*(p*p + 2.0*d*d)*kap0 -  \
            (1.0 + 5.0*p*p + d*d)/4.0*sqrt((1.0 - x1)*(x2 - 1.0)));
        return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
    }

    // occulting star transits the source
    if(p <= 1.0  && d <= (1.0 - p)) {
        etad = p*p/2.0*(p*p + 2.0*d*d);
        lambdae = p*p;

        // zone 4.1
        T q = sqrt((x2 - x1)/(1.0 - x1));
        T Kk = K(q);
        T Ek = E(q);
        T n = x2/x1 - 1.0;
        T Pk = ell_bulirsch(n, q);

        lambdad = 2.0/9.0/M_PI/sqrt(1.0 - x1)*((1.0 - 5.0*d*d + p*p +  \
                    x3*x3)*Kk + (1.0 - x1)*(d*d + 7.0*p*p - 4.0)*Ek - 3.0*x3/x1*Pk);

        // edge of planet hits edge of star
        if (fabs(p + d - 1.0) <= tol) {
            lambdad = 2.0/3.0/M_PI*acos(1.0 - 2.0*p) - 4.0/9.0/M_PI* \
                        sqrt(p*(1.0 - p))*(3.0 + 2.0*p - 8.0*p*p);
        }
        if (d < p) {
            lambdad += 2.0/3.0;
        }
    }
    return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;

}