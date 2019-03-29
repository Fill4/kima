#include "DNest4.h"
#include "Data.h"
#include "RVmodel.h"

using namespace DNest4;

#include "default_priors.h"

const bool obs_after_HARPS_fibers = false;
const bool GP = true;
const bool hyperpriors = false;
const bool trend = false;
const bool multi_instrument = false;

// options for the model
RVmodel::RVmodel():fix(false),npmax(5)
{}


int main(int argc, char** argv)
{
    /* set the RV data file */
    // kima skips the first 2 lines in the header
    // and reads the first 3 columns into time, vrad and svrad
    char* datafile = "corot7.txt";

    Data::get_instance().load(datafile, "ms");
    
    // set the sampler and run it!
    Sampler<RVmodel> sampler = setup<RVmodel>(argc, argv);
    sampler.run();

    return 0;
}
