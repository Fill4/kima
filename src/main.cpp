#include "kima.h"

// const bool obs_after_HARPS_fibers = false;
// const bool GP = false;
// const bool MA = false;
// const bool hyperpriors = false;
// const bool trend = false;
// const bool multi_instrument = false;

int main(int argc, char** argv)
{
    Data& data = Data::instance();
    // Set the RV and transit data files
    rv_files = {"../examples/BL2009/BL2009_dataset1.kms.rv"};
    transit_files = {"../examples/sim-1560290/sim-1560290.lc"};

    // Load the files (RVs are in km/s). Don't skip any lines in the header.
    data.load_rv(rv_files, "kms", 0);
    // Load the transit files (Fractional flux)
    data.load_transit(transit_files, "frac", 0);

    // Check the dataset was read correctly. DEBUG only
    // std::cout << std::setprecision(12);
    // std::cout << "RV Dataset: " << data.rv.t.size() << " points" << std::endl;
    // for (size_t i = 0; i < data.rv.t.size(); i++) {
    //     std::cout << "t: " << data.rv.t[i] << " y: " << data.rv.y[i] << " y_err: " << data.rv.y_err[i] << std::endl;
    // }
    // std::cout << "Transit Dataset: " << data.transit.t.size() << " points" << std::endl;
    // for (size_t i = 0; i < 10; i++) {
    //     std::cout << "t: " << data.transit.t[i] << " y: " << data.transit.y[i] << " y_err: " << data.transit.y_err[i] << std::endl;
    // }
    // exit(0);

    // set the sampler and run it!
    // Sampler<RVmodel> sampler = setup<RVmodel>(argc, argv);
    DNest4::Sampler<System> sampler = DNest4::setup<System>(argc, argv);
    sampler.run(50);

    return 0;
}
