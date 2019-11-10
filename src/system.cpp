#include <cstdint>
#include <ostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <iomanip>

#include "DNest4.h"
#include "spdlog/spdlog.h"

#include "system.h"
#include "data.h"
#include "kepler.h"
#include "limbdark.h"

const double halflog2pi = 0.5*log(2.*M_PI);

System::System() {
    Data& data = Data::instance();
    // TODO: Flag for which data we init and use should start here
    
    // ---------------------------- RV -----------------------------
    // Init the priors for all the common parameters
    background.prior = std::make_shared<DNest4::Uniform>(data.rv_y_min(), data.rv_y_max());
    rv_jitter.prior = std::make_shared<DNest4::ModifiedLogUniform>(1.0, 100.);
    slope.prior = std::make_shared<DNest4::Uniform>(-data.rv_top_slope(), data.rv_top_slope());
    fiber_offset.prior = std::make_shared<DNest4::Uniform>(0, 50);
    offset.prior = std::make_shared<DNest4::Uniform>(-data.rv_y_span(), data.rv_y_span());

    // Init all arrays based on the processed data
    offsets = std::vector<double>(data.rv.number_instruments - 1);
    rv_jitters = std::vector<double>(data.rv.number_instruments);
    betas = std::vector<double>(data.rv.number_indicators);

    rv = std::vector<double>(data.rv.t.size());
    
    // ----------------------------- Transit -------------------------
    u1.prior = std::make_shared<DNest4::Uniform>(0.4, 0.6);
    u2.prior = std::make_shared<DNest4::Uniform>(0.05, 0.15);
    transit_jitter.prior = std::make_shared<DNest4::Uniform>(1.0, 400.0);

    transit = std::vector<double>(data.transit.t.size());

    // TODO: Init all the planets priors
    // Get the number of parameters in planet according to the data we have
    uint32_t npar = data.rv.has_rv ? 5 : 4;
    npar += data.transit.has_transit ? 3 : 0;
    if (npar == 4) {
        spdlog::error("Need to have either RV or transit data");
        exit(3);
    }
    planets = DNest4::RJObject<Planet>(npar, npmax, fix, Planet(data.rv.has_rv, data.transit.has_transit));
}

// TODO: No active indicators, moving average or GPs included yet
void System::from_prior(DNest4::RNG& rng) {

    // Sample the priors of all the planet parameters
    planets.from_prior(rng);
    planets.consolidate_diff();

    // Sample RV global parameters
    background.value = background.prior->generate(rng);

    if (multi_instrument) {
        for(uint32_t i=0; i<offsets.size(); i++) {
            offsets[i] = offset.prior->generate(rng);
        }
        for(uint32_t i=0; i<rv_jitters.size(); i++) {
            rv_jitters[i] = rv_jitter.prior->generate(rng);
        }
    } 
    else {
        rv_jitter.value = rv_jitter.prior->generate(rng);
    }

    if (obs_after_HARPS_fibers) {
        fiber_offset.value = fiber_offset.prior->generate(rng);
    }
    if (trend) {
        slope.value = slope.prior->generate(rng);
    }

    // Sample transit global parameters
    u1.value = u1.prior->generate(rng);
    u2.value = u2.prior->generate(rng);
    transit_jitter.value = transit_jitter.prior->generate(rng);
}

// TODO: No active indicators, moving average or GPs included yet
double System::perturb(DNest4::RNG& rng) {
    double logH = 0.0;

    if (rng.rand() <= 0.75) {
        logH += planets.perturb(rng);
        planets.consolidate_diff();
        // Planets parameters have been changed so orbits need to be updated
        // update_orbits = true; // TODO: Not implemented yet. Need to save planet arrays
    } 
    else if (rng.rand() <= 0.5) {
        if (multi_instrument) {
            for(int i=0; i<rv_jitters.size(); i++) {
                rv_jitter.prior->perturb(rv_jitters[i], rng);
            }
        }
        else {
            rv_jitter.prior->perturb(rv_jitter.value, rng);
            transit_jitter.prior->perturb(transit_jitter.value, rng);
        }
    }
    else {
        background.prior->perturb(background.value, rng);

        // Propose new instrument offsets
        if (multi_instrument) {
            for(uint32_t j=0; j<offsets.size(); j++) {
                offset.prior->perturb(offsets[j], rng);
            }
        }

        // Propose new fiber offset
        if (obs_after_HARPS_fibers) {
            fiber_offset.prior->perturb(fiber_offset.value, rng);
        }

        // propose new slope
        if (trend) {
            slope.prior->perturb(slope.value, rng);
        }

        u1.prior->perturb(u1.value, rng);
        u2.prior->perturb(u2.value, rng);
    }

    return logH;
}

void System::print(std::ostream& out) const {
    Data& data = Data::instance();
    // output precision
    out.setf(std::ios::fixed,std::ios::floatfield);
    out.precision(8);

    if (multi_instrument) {
        for(int j=0; j<rv_jitters.size(); j++) {
            out << rv_jitters[j] << ' ';
        }
    }
    else {
        out << rv_jitter.value << ' ';
    }

    if (trend) {
        out << slope.value << ' ';
    }

    if (obs_after_HARPS_fibers) {
        out << fiber_offset.value << ' ';
    }

    if (multi_instrument) {
        for(int j=0; j<offsets.size(); j++) {
            out << offsets[j] << ' ';
        }
    }

    if (data.rv.indicator_correlations){
        for(size_t j=0; j<data.rv.number_indicators; j++){
            out << betas[j] << ' ';
        }
    }

    out << transit_jitter.value << ' ';
    out << u1.value << ' ';
    out << u2.value << ' ';

    planets.print(out);

    // out << ' ' << staleness << ' ';
    out << background.value;
}

std::string System::description() const {
    Data& data = Data::instance();
    std::string desc;

    if (multi_instrument) {
        for(int j=0; j<rv_jitters.size(); j++) {
           desc += "rv_jitter" + std::to_string(j+1) + "  ";
        }
    }
    else {
        desc += "extra_sigma  ";
    }

    if (trend) {
        desc += "slope  ";
    }

    if (obs_after_HARPS_fibers) {
        desc += "fiber_offset  ";
    }

    if (multi_instrument) {
        for(unsigned j=0; j<offsets.size(); j++) {
            desc += "offset" + std::to_string(j+1) + "  ";
        }
    }

    if (data.rv.indicator_correlations) {
        for(size_t j=0; j<data.rv.number_indicators; j++) {
            desc += "beta" + std::to_string(j+1) + "  ";
        }
    }

    desc += "tr_jitter  ";
    desc += "u1  ";
    desc += "u2  ";

    desc += "ndim  maxNp  ";
    desc += "Np  ";

    if (planets.get_max_num_components()>0) {
        desc += "P  ecc  w  phi  K  RpRs  aR  inc  u1  u2  ";
    }

    desc += "staleness  vsys";

    return desc;
}

double System::log_likelihood() {
    // Get the data
    Data& data = Data::instance();

    // Get planet components  // TODO: Right now no staleness calculation is done
    const std::vector<std::vector<double>>& components = planets.get_components();

    // For each of the planets in the list. Define planet params and other necessary variables
    double P, ecc, w, phi, K, t0, RpRs, aR, inc;
    double n, fc, Ec, Mc, tp, tc;
    for(size_t j=0; j<components.size(); j++) {
        P = components[j][0];
        ecc = components[j][1];
        w = components[j][2];
        phi = components[j][3];
        K = components[j][4];
        // t0 = components[j][5];
        RpRs = components[j][5];
        aR = components[j][6];
        inc = components[j][7];
        
        // Mean motion
        n = 2*M_PI/P;

        // Time of periastron passage for RVs from phi (which seems to be the mean anomaly at t[0])
        //! Right now requires phi to define the time of periastron, which can't work when there is only a transit
        tp = data.rv.t[0] - phi/n;

        //! This method is only valid for non-extreme eccentricities. Should be checked beforehand
        // True anomaly at inferior conjunction
        fc = M_PI/2.0 - w;
        // Eccentric anomaly at inferior conjuction
        Ec = 2.0*atan(sqrt((1.0 - ecc)/(1.0 + ecc))*tan(fc/2.0));
        // Mean anomaly at inferior conjuction
        Mc = Ec - ecc*sin(Ec);
        // Time of inferior conjuction from time of periastron passage for transits
        tc = tp + Mc/n;

        // Calculate the RV and transit signals for this planet and add them to the arrays
        //? We are right now assuming that no transiting planets cross in front of each other
        // calculate_rv(n, tp, ecc, w, K);
        calculate_transit(n, tp, ecc, w, RpRs, aR, inc);
    }
    
    double loglike_rv = 0.0;
    double loglike_tr = 0.0;
    // Then get the likelihoods from RV and transit
    loglike_rv = likelihood_rv();
    loglike_tr = likelihood_transit();

    double loglike = 0;
    if (std::isnan(loglike_rv) || std::isinf(loglike_rv) || std::isnan(loglike_tr) || std::isinf(loglike_tr)) {
        loglike = std::numeric_limits<double>::infinity();
    } else {
        loglike = loglike_rv + loglike_tr;
    }
    return loglike;
}

// Calculate the RVs for all times
void System::calculate_rv(const double& n, const double& tp, const double& ecc, const double& w, const double& K) {
    Data& data = Data::instance();
    double M, E, f;
    for(size_t i=0; i<data.rv.t.size(); i++) {
        // Mean anomaly //! Using tp_rv
        M = mean_anomaly(data.rv.t[i], tp, n);
        // Eccentric anomaly
        E = ecc_anomaly(M, ecc);
        // True anomaly
        f = true_anomaly(E, ecc);
        // Radial velocity // TODO: Need to init this rvs array to 0
        rv[i] += K*(cos(w + f) + ecc*(cos(w)));

        // Add other contributions to rv signal
        //! Right now there are no trends, multi_instrument or anything other than the rvs
    }
}

// Calculate transit depth for all times
void System::calculate_transit(const double& n, const double& tp, const double& ecc, const double& w, 
                            const double& RpRs, const double& aR, const double& inc) {
    Data& data = Data::instance();
    double M, E, f, r, b, x, y, z;
    for(size_t i=0; i<data.transit.t.size(); i++) {
        // Mean anomaly //! Using tp_tr
        M = mean_anomaly(data.transit.t[i], tp, n);
        // Eccentric anomaly
        E = ecc_anomaly(M, ecc);
        // True anomaly
        f = true_anomaly(E, ecc);
        // Distance between centers (normalized to Rstar as aR is normalized already)
        r = aR*(1 - ecc*ecc) / (1 + ecc*cos(f));
        // Instantaneous impact parameter
        b = r * sqrt(1 - sin(w + f) * sin(w + f) * sin(inc) * sin(inc));
        // Cartesian coordinates x, y and z //! Only z is being calculated
        z = r * sin(w + f) * sin(inc);
        //! Right now we are following batman and pysyzygy and will ignore secondaries
        transit[i] = quadratic_ld(b, RpRs, z, u1.value, u2.value);
    }
}

double System::likelihood_rv() {
    Data& data = Data::instance();
    double loglike = 0.0;

    double jit, var;
    for(size_t i=0; i<data.rv.t.size(); i++) {
        if(multi_instrument) {
            jit = rv_jitters[data.rv.instrument_ids[i]-1];
            var = data.rv.y_err[i]*data.rv.y_err[i] + jit*jit;
        }
        else {
            var = data.rv.y_err[i]*data.rv.y_err[i] + rv_jitter.value*rv_jitter.value;
        }
        loglike += -halflog2pi - 0.5*log(var) - 0.5*(pow(data.rv.y[i] - rv[i], 2)/var);
    }

    return loglike;
}

double System::likelihood_transit() {
    Data& data = Data::instance();
    double loglike = 0.0;

    for(size_t i=0; i<data.transit.t.size(); i++) {
        double var = data.transit.y_err[i]*data.transit.y_err[i] + transit_jitter.value*transit_jitter.value;
        loglike += -halflog2pi - 0.5*log(var) - 0.5*(pow(data.transit.y[i] - transit[i], 2)/var);
    }

    return loglike;
}
