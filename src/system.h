# pragma once

#include <cstdint>
#include <ostream>
#include <vector>
#include <string>
#include <functional>

#include "RJObject/RJObject.h"
#include "RNG.h"

#include "parameter.h"
#include "planet.h"

struct System
{
    System();

    // Methods for DNest4
    /// @brief Generate a point from the prior.
    void from_prior(DNest4::RNG& rng);

    /// @brief Do Metropolis-Hastings proposals
    double perturb(DNest4::RNG& rng);

    // Likelihood function
    double log_likelihood();

    // Print parameters to stream
    void print(std::ostream& out) const;

    // Return string with column information
    std::string description() const;


    // Flags and general definitions
    /// Fix the number of planets? (by default, yes)
    bool fix {true};
    /// Maximum number of planets
    uint32_t npmax {1};


    // -----------------------------------------------------------
    // Data specific flags
    /// Whether the model includes a linear trend in the RVs
    bool trend {false};
    /// Whether the data comes from multiple instruments
    bool multi_instrument {false};
    /// Whether the data comes from after the change in HARPS fibers
    bool obs_after_HARPS_fibers {false};

    // TODO: No moving average (MA) of Gaussian process (GP) is included yet

    // Planets
    // DNest4::RJObject<Planet> planets = DNest4::RJObject<Planet>(7, npmax, fix, Planet());
    DNest4::RJObject<Planet> planets;

    // Non-planet parameters - RV
    /// Systemic velocity
    Parameter background;
    /// RV Jitter (extra white noise)
    Parameter rv_jitter;
    /// Slope for the trend (used if 'trend == true')
    Parameter slope;
    /// HARPS fiber RV offset
    Parameter fiber_offset;
    /// Multi instrument offsets
    Parameter offset;

    // Non-planet parameters - Transit
    // Quadratic limb-darkening coeficcients
    Parameter u1, u2;
    /// Transit Jitter (extra white noise)
    Parameter transit_jitter;

    // Helper arrays
    // Offsets between instruments
    std::vector<double> offsets;
    // Jitter for each instrument
    std::vector<double> rv_jitters;
    // "slopes" for each indicator
    std::vector<double> betas;

    void calculate_rv(const double& n, const double& tp, const double& ecc, const double& w, const double& K);
    void calculate_transit(const double& n, const double& tp, const double& ecc, const double& w, 
                        const double& RpRs, const double& aR, const double& inc);

    double likelihood_rv();
    double likelihood_transit();

    // RV model array for all planets
    std::vector<double> rv;
    // Transit model array for all planets
    std::vector<double> transit;
    // -----------------------------------------------------------
};
