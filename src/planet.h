#pragma once

#include <vector>
#include <ostream>

#include "RNG.h"
#include "RJObject/ConditionalPriors/ConditionalPrior.h"

#include "parameter.h"

struct Planet : public DNest4::ConditionalPrior
{
    Planet() {}
    Planet(bool has_rv, bool has_transit);

    double perturb_hyperparameters(DNest4::RNG& rng);

    /// Generate a point from the prior
    void from_prior(DNest4::RNG& rng);

    double log_pdf(const std::vector<double>& vec) const;
    /// Get parameter sample from a uniform sample (CDF)
    void from_uniform(std::vector<double>& vec) const;
    /// Get uniform sample from a parameter sample (inverse CDF)
    void to_uniform(std::vector<double>& vec) const;

    void print(std::ostream& out) const;

    static const int weight_parameter = 1;

    std::vector<Parameter> parameters;

    //--------------------------------------------------------
    // Priors for all planet parameters

    // Parameters common for both RVs and transits
    /// Orbital period
    Parameter P;
    /// Eccentricity
    Parameter ecc;
    /// Argument of the periastron
    Parameter w;
    /// Phase - Can be substituted by epoch. Used to determine periastron time.
    Parameter phi;

    // Parameters for RVs
    /// Semi-amplitude (in m/s)
    // Parameter K;

    // Parameters for transits
    /// Epoch
    // Parameter t0;
    /// Ratio between stellar and planetary radius
    Parameter RpRs;
    /// Semi-major axis (in stellar radii units)
    Parameter aR;
    /// Inclination 
    Parameter inc;
    //-------------------------------------------------------
};