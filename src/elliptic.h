#pragma once

#include <cmath>

/**
    Computes polynomial approximation for the complete elliptic integral of the second kind (Hasting's approximation)

    Obtained from batman: https://github.com/lkreidberg/batman/blob/master/c_src/_quadratic_ld.c
*/
template <typename T>
inline T E(T k) {
    T m1, a1, a2, a3, a4, b1, b2, b3, b4, ee1, ee2, ellec;
    m1 = 1.0 - k*k;
    a1 = 0.44325141463;
    a2 = 0.06260601220;
    a3 = 0.04757383546;
    a4 = 0.01736506451;
    b1 = 0.24998368310;
    b2 = 0.09200180037;
    b3 = 0.04069697526;
    b4 = 0.00526449639;
    ee1 = 1.0 + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4)));
    ee2 = m1*(b1 + m1*(b2 + m1*(b3 + m1*b4)))*log(1.0/m1);
    ellec = ee1 + ee2;
    return ellec;
}

/**
    Computes polynomial approximation for the complete elliptic integral of the first kind (Hasting's approximation)
    
    Obtained from batman: https://github.com/lkreidberg/batman/blob/master/c_src/_quadratic_ld.c
*/
template <typename T>
inline T K(T k) {
    T a0, a1, a2, a3, a4, b0, b1, b2, b3, b4, ellk,  ek1, ek2, m1;
    m1 = 1.0 - k*k;
    a0 = 1.38629436112;
    a1 = 0.09666344259;
    a2 = 0.03590092383;
    a3 = 0.03742563713;
    a4 = 0.01451196212;
    b0 = 0.5;
    b1 = 0.12498593597;
    b2 = 0.06880248576;
    b3 = 0.03328355346;
    b4 = 0.00441787012;
    ek1 = a0 + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4)));
    ek2 = (b0 + m1*(b1 + m1*(b2 + m1*(b3 + m1*b4))))*log(m1);
    ellk = ek1 - ek2;
    return ellk;
}

/**
    Computes the complete elliptical integral of the third kind using the algorithm of Bulirsch (1965):
    
    Bulirsch 1965, Numerische Mathematik, 7, 78
    Bulirsch 1965, Numerische Mathematik, 7, 353

    Obtained from batman: https://github.com/lkreidberg/batman/blob/master/c_src/_quadratic_ld.c
*/
template <typename T>
inline T ell_bulirsch(T n, T k) {
    T kc = sqrt(1.-k*k);
    T p = sqrt(n + 1.);
    T m0 = 1.;
    T c = 1.;
    T d = 1./p;
    T e = kc;
    T f, g;

    uint32_t nit = 0;
    while(nit < 10000) {
        f = c;
        c = d/p + c;
        g = e/p;
        d = 2.*(f*g + d);
        p = g + p;
        g = m0;
        m0 = kc + m0;
        if (fabs(1.-kc/g) > 1.0e-8) {
            kc = 2.*sqrt(e);
            e = kc*m0;
        }
        else {
            return 0.5*M_PI*(c*m0+d)/(m0*(m0+p));
        }
        nit++;
    }
    return 0;
}