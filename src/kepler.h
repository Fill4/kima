#pragma once

#include <cmath>

/**
    Calculates the mean anomaly at time t by solving Kepler's equation.

    \param t time at current timestep
    \param tp time of periastron passage
    \param n mean motion
    \return mean anomaly
*/
template <typename T>
inline T mean_anomaly(T t, T tp, T n) {
    return n*(t - tp);
}

/**
    Provides a starting value to solve Kepler's equation.
    See "A Practical Method for Solving the Kepler Equation", Marc A. Murison, 2006

    @param M mean anomaly (in radians)
    @param e eccentricity of the orbit
    @return starting value for the eccentric anomaly
*/
template <typename T>
inline T keplerstart3(T M, T ecc) {
    T t34 = ecc*ecc;
    T t35 = ecc*t34;
    T t33 = cos(M);
    return M + (-0.5*t35 + ecc + (t34 + 1.5*t33*t35)*t33)*sin(M);
}

/**
    An iteration (correction) method to solve Kepler's equation.
    See "A Practical Method for Solving the Kepler Equation", Marc A. Murison, 2006

    @param M mean anomaly (in radians)
    @param e the eccentricity of the orbit
    @param x starting value for the eccentric anomaly
    @return corrected value for the eccentric anomaly
*/
template <typename T>
inline T eps3(T M, T ecc, T x) {
    T t1 = cos(x);
    T t2 = -1 + ecc*t1;
    T t3 = sin(x);
    T t4 = ecc*t3;
    T t5 = -x + t4 + M;
    T t6 = t5/(0.5*t5*t4/t2+t2);

    return t5/((0.5*t3 - 1/6*t1*t6)*ecc*t6+t2);
}

/**
    Calculates the eccentric anomaly at time t by solving Kepler's equation.
    See "A Practical Method for Solving the Kepler Equation", Marc A. Murison, 2006

    @param M mean anomaly
    @param ecc eccentricity of the orbit
    @return eccentric anomaly
*/
template <typename T>
inline T ecc_anomaly(T M, T ecc) {
    if (ecc <= 1.0e-10) {
        return M;
    }
    
    T tol;
    if (ecc < 0.8) {
        tol = 1.0e-14;
    }
    else {
        tol = 1.0e-10;
    }    

    T Mnorm = fmod(M, 2.*M_PI);
    T E0 = keplerstart3(Mnorm, ecc);
    T dE = tol + 1;
    T E;
    uint32_t count = 0;
    while (dE > tol) {
        E = E0 - eps3(Mnorm, ecc, E0);
        dE = abs(E-E0);
        E0 = E;
        count++;
        // failed to converge, this only happens for nearly parabolic orbits
        if (count == 100) break;
    }
    return E;
}

/**
    Calculates the true anomaly at time t.
    See Eq. 2.6 of The Exoplanet Handbook, Perryman 2010

    @param E eccentric anomaly
    @param ecc eccentricity of the orbit
    @return true anomaly
*/
template <typename T>
inline T true_anomaly(T E, T ecc) {
    if (ecc <= 1.0e-10) {
        return E;
    }
    T f = 2.0*atan(sqrt((1.0+ecc)/(1.0-ecc))*tan(E/2.0));
    return f;
}

