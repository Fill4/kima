#include <cmath>

#include "DNest4.h"
#include "spdlog/spdlog.h"

#include "planet.h"

Planet::Planet(bool has_rv, bool has_transit) {
    // For now all the priors of the parameters are setup when the instances are initialized 
    // before the System constructor is called
    // TODO: Setup the priors before the constructor so that we don't have to rewrite on top

    // Common priors
    // P.prior = std::make_shared<DNest4::LogUniform>(1., 5.);
    // ecc.prior = std::make_shared<DNest4::Uniform>(0., 1);
    // w.prior = std::make_shared<DNest4::Uniform>(0., 2.*M_PI);
    // phi.prior = std::make_shared<DNest4::Uniform>(0., 2.*M_PI);

    // RV priors
    // K.prior = std::make_shared<DNest4::ModifiedLogUniform>(1., 1e3);

    // Transit priors  //! Right now the values are absurd and should be ignored
    // t0.prior = std::make_shared<DNest4::Uniform>(0, 1);
    // RpRs.prior = std::make_shared<DNest4::LogUniform>(0.001, 0.1);
    // aR.prior = std::make_shared<DNest4::Uniform>(1., 5.);
    // inc.prior = std::make_shared<DNest4::Uniform>(0., M_PI/2.);

    // This should never trigger as we have the check before constructor in system.cpp
    if ((!has_rv) || (!has_transit)) {
        spdlog::error("Cannot create planet with no RV and transit data");
        exit(4);
    }
    //! No identification of the parameters right now. Maybe add a name to the Parameter class
    // Common parameters always get added
    // Period P
    parameters.push_back(Parameter(std::make_shared<DNest4::LogUniform>(1., 5.)));
    // Eccentricity ecc
    parameters.push_back(Parameter(std::make_shared<DNest4::Uniform>(0., 1)));
    // Argument of periastron w
    parameters.push_back(Parameter(std::make_shared<DNest4::Uniform>(0., 2.*M_PI)));
    // Phi
    parameters.push_back(Parameter(std::make_shared<DNest4::Uniform>(0., 2.*M_PI)));

    // Parameters in case of RV
    if (has_rv) {
        // Semi-amplitude K
        parameters.push_back(Parameter(std::make_shared<DNest4::ModifiedLogUniform>(1., 1e3)));
    }
    if (has_transit) {
        // Ratio of radii RpRs
        parameters.push_back(Parameter(std::make_shared<DNest4::LogUniform>(0.001, 0.1)));
        // Semi-major axis a
        parameters.push_back(Parameter(std::make_shared<DNest4::Uniform>(1., 5.)));
        // Inclination inc
        parameters.push_back(Parameter(std::make_shared<DNest4::Uniform>(0., M_PI/2.)));
    }
}

// TODO: The class has no hyperparameters right now so the function is empty
void Planet::from_prior(DNest4::RNG& rng) {
    return;
}

// TODO: The class has no hyperparameters right now so the function is empty
double Planet::perturb_hyperparameters(DNest4::RNG& rng) {
    double logH = 0.0;
    return logH;
}

// TODO: No hyperpriors for now. Also only working for RV parameters, ignoring the transit ones
// TODO: There are no checks for boundaries for now. Should not be necessary for uniform priors but
// TODO: might not be true for the remaining distributions
double Planet::log_pdf(const std::vector<double>& vec) const {
    double logpdf = 0.;
    uint32_t i = 0;
    for (const Parameter& parameter : parameters) {
        logpdf += parameter.prior->log_pdf(vec[i]);
        ++i;
    }
    return logpdf;
    // return P.prior->log_pdf(vec[0]) + 
    //        ecc.prior->log_pdf(vec[1]) + 
    //        w.prior->log_pdf(vec[2]) + 
    //        phi.prior->log_pdf(vec[3]) +  
    //        RpRs.prior->log_pdf(vec[4]) + 
    //        aR.prior->log_pdf(vec[5]) + 
    //        inc.prior->log_pdf(vec[6]);
}

// TODO: No hyperpriors for now. Also only working for RV parameters, ignoring the transit ones
void Planet::from_uniform(std::vector<double>& vec) const {
    uint32_t i = 0;
    for (const Parameter& parameter : parameters) {
        vec[i] = parameter.prior->cdf_inverse(vec[i]);
        ++i;
    }
    // vec[0] = P.prior->cdf_inverse(vec[0]);
    // vec[1] = ecc.prior->cdf_inverse(vec[1]);
    // vec[2] = w.prior->cdf_inverse(vec[2]);
    // vec[3] = phi.prior->cdf_inverse(vec[3]);
    // vec[4] = RpRs.prior->cdf_inverse(vec[4]);
    // vec[5] = aR.prior->cdf_inverse(vec[5]);
    // vec[6] = inc.prior->cdf_inverse(vec[6]);
}

// TODO: No hyperpriors for now. Also only working for RV parameters, ignoring the transit ones
void Planet::to_uniform(std::vector<double>& vec) const {
    uint32_t i = 0;
    for (const Parameter& parameter : parameters) {
        vec[i] = parameter.prior->cdf(vec[i]);
        ++i;
    }
    // vec[0] = P.prior->cdf(vec[0]);
    // vec[1] = ecc.prior->cdf(vec[1]);
    // vec[2] = w.prior->cdf(vec[2]);
    // vec[3] = phi.prior->cdf(vec[3]);
    // vec[4] = RpRs.prior->cdf(vec[4]);
    // vec[5] = aR.prior->cdf(vec[5]);
    // vec[6] = inc.prior->cdf(vec[6]);
}

// TODO: The class has no hyperparameters right now so the function is empty
void Planet::print(std::ostream& out) const {
    return;
}