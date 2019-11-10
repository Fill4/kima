#pragma once

#include <memory>

#include "Distributions/ContinuousDistribution.h"

struct Parameter {
    /// Parameter value
    double value {0};
    /// Pointer to parameter prior distribution
    std::shared_ptr<DNest4::ContinuousDistribution> prior = nullptr;
    
    // Contructors
    Parameter() {};
    Parameter(std::shared_ptr<DNest4::ContinuousDistribution> prior) : prior(prior) {}
};
