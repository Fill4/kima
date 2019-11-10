#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <vector>
#include <string>
#include <set>
#include <cmath>

#include "data.h"

// Reads the contents from a csv file into a std::vector<std::vector<double>>
// TODO: Move this into the header to template it. Maybe
void Data::read_csv(const std::string& filename, std::vector<std::vector<double>>& data) {
    std::ifstream file(filename);

    std::string line;
    std::vector<double> row;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double f;
        while (ss >> f) {
            row.push_back(f);
        }
        data.push_back(row);
        row.clear();
    }
    
    // Complain if something went wrong.
    if (!file.eof()) {
        std::cout <<"Could not read data file (" << filename << ")!" << std::endl;
        exit(1);
    }
    file.close();
}

// std::istream& operator >> ( std::istream& ins, record_t& record ) {
//     // Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
//     //-----------------------------------------------------------------------------

//     // make sure that the returned record contains only the stuff we read now
//     record.clear();

//     // read the entire line into a string (a CSV record is terminated by a newline)
//     std::string line;
//     getline(ins, line);

//     // now we'll use a stringstream to separate the fields out of the line
//     std::stringstream ss( line );

//     // convert each field to a double and 
//     // add the newly-converted field to the end of the record
//     double f;
//     while (ss >> f) {
//         record.push_back(f);
//     }

//     // Now we have read a single line, converted into a list of fields, converted the fields
//     // from strings to doubles, and stored the results in the argument record, so
//     // we just return the argument stream as required for this kind of input overload function.
//     return ins;
// }


// std::istream& operator >> ( std::istream& ins, data_t& data ) {
//     // Let's likewise overload the stream input operator to read a list of CSV records.
//     // This time it is a little easier, just because we only need to worry about reading
//     // records, and not fields.
//     //-----------------------------------------------------------------------------
    
//     // make sure that the returned data only contains the CSV data we read here
//     // data.clear();

//     // For every record we can read from the file, append it to our resulting data
//     // except if it's in the header
//     int i = 0;
//     record_t record;
//     while (ins >> record) {
//         data.push_back(record);
//     }

//     // Again, return the argument stream
//     return ins;  
// }


// /**
//  * @brief Load RV data from a file.
//  *
//  * Read a tab/space separated file with columns  
//  * ```
//  *   time  vrad  error  
//  *   ...   ...   ...
//  * ```
//  * 
//  * @param filename   the name of the file
//  * @param units      units of the RVs and errors, either "kms" or "ms"
//  * @param skip       number of lines to skip in the beginning of the file (default = 2)
//  * @param indicators
// */
// void Data::load_rv(const std::string filename, const std::string units, uint32_t skip, const std::vector<std::string>& indicators) {
//     data_t data;

//     // Empty the vectors
//     rv.t.clear();
//     rv.y.clear();
//     rv.y_err.clear();
    
//     // Check for indicator correlations and store stuff
//     size_t nempty = count(indicators.begin(), indicators.end(), "");
//     rv.number_indicators = static_cast<unsigned int>(indicators.size() - nempty);
//     rv.indicator_correlations = rv.number_indicators > 0;
//     rv.indicator_names = indicators;
//     rv.indicator_names.erase(
//         std::remove( rv.indicator_names.begin(), rv.indicator_names.end(), "" ), 
//         rv.indicator_names.end() );

//     // Empty the indicator vectors as well
//     rv.activity_indicators.clear();
//     rv.activity_indicators.resize(rv.number_indicators);
//     for (size_t n = 0; n < rv.number_indicators; n++) {
//         rv.activity_indicators[n].clear();
//     }

//     // Read the file into the data container
//     std::ifstream infile(filename);
//     infile >> data;

//     // Complain if something went wrong.
//     if (!infile.eof()) {
//         std::cout <<"Could not read data file (" << filename << ")!" << std::endl;
//         exit(1);
//     }
//     infile.close();

//     rv.filenames.push_back(filename);
//     rv.units = units;
//     rv.skip = skip;
//     rv.number_instruments = 1;


//     double factor = 1.0;
//     if (units == "kms") factor = 1E3;
//     int j;

//     for (size_t n = 0; n < data.size(); n++) {
//         if (n<skip) continue;
//         rv.t.push_back(data[n][0]);
//         rv.y.push_back(data[n][1] * factor);
//         rv.y_err.push_back(data[n][2] * factor);
    
//         if (rv.indicator_correlations) {
//             j = 0;
//             for (size_t i = 0; i < rv.number_indicators + nempty; i++) {
//                 if (indicators[i] == "") {
//                     continue; // skip column
//                 }
//                 else {
//                     rv.activity_indicators[j].push_back(data[n][3+i] * factor);
//                     j++;
//                 }
//             }
//         }
//     }
    

//     // subtract means from activity indicators
//     if (rv.indicator_correlations) {
//         double mean;
//         for (auto& i : rv.activity_indicators) { // use auto& instead of auto to modify i
//             mean = std::accumulate(i.begin(), i.end(), 0.0) / rv.y.size();
//             std::for_each(i.begin(), i.end(), [mean](double& d) { d -= mean;});
//         }
//     }

//     // How many points did we read?
//     std::cout << "# Loaded " << rv.t.size() << " data points from file" << filename << std::endl;
//     // Did we read activity indicators? how many?
//     if (rv.indicator_correlations) {
//         std::cout << "# Loaded " << rv.activity_indicators[0].size() << " observations of "
//                   << rv.activity_indicators.size() << " activity indicators: " << std::endl;
//         for (const auto i: indicators) {
//             if (i != "") {
//                 std::cout << "'" << i << "'" << std::endl;
//                 (i != indicators.back()) ? std::cout << ", " : std::cout << " ";
//             }
//         }
//         std::cout << std::endl;
//     }
//     // What are the units?
//     if (units == "kms") {
//         std::cout << "# Multiplied all RVs by 1000; units are now m/s." << std::endl;
//     }

//     for(size_t i=0; i<data.size(); i++) {
//         if (rv.t[i] > 57170.) {
//             rv.index_fibers = static_cast<uint32_t>(i);
//             break;
//         }
//     }
// }


// /**
//  * @brief Load RV data from a multi-instrument file.
//  * 
//  * Read a tab/space separated file with columns  
//  * ```
//  *   time  vrad  error  obs
//  *   ...   ...   ...    ...
//  * ```
//  * The `obs` column should be an integer identifying the instrument.
//  * 
//  * @param filename   the name of the file
//  * @param units      units of the RVs and errors, either "kms" or "ms"
//  * @param skip       number of lines to skip in the beginning of the file (default = 2)
// */
// void Data::load_rv_multi(const std::string filename, const std::string units, uint32_t skip) {
//     data_t data;

//     // Empty the vectors
//     rv.t.clear();
//     rv.y.clear();
//     rv.y_err.clear();
//     rv.instrument_ids.clear();

//     // Read the file into the data container
//     std::ifstream infile(filename);
//     infile >> data;

//     // Complain if something went wrong.
//     if (!infile.eof()) {
//         std::cout <<"Could not read data file (" << filename << ")!" << std::endl;
//         exit(1);
//     }
//     infile.close();

//     rv.filenames.push_back(filename);
//     rv.units = units;
//     rv.skip = skip;

//     double factor = 1.;
//     if (units == "kms") factor = 1E3;

//     for (size_t n = 0; n < data.size(); n++) {
//         if (n<skip) continue;
//         rv.t.push_back(data[n][0]);
//         rv.y.push_back(data[n][1] * factor);
//         rv.y_err.push_back(data[n][2] * factor);
//         rv.instrument_ids.push_back(static_cast<int32_t>(data[n][3]));
//     }

//     // How many points did we read?
//     std::cout << "# Loaded " << rv.t.size() << " data points from file" << filename << std::endl;

//     // Of how many instruments?
//     std::set<int> s( rv.instrument_ids.begin(), rv.instrument_ids.end() );
//     std::cout << "# RVs come from " << s.size() << " different instruments." << std::endl;
//     rv.number_instruments = static_cast<uint32_t>(s.size());
    
//     if (units == "kms") {
//         std::cout << "# Multiplied all RVs by 1000; units are now m/s." << std::endl;
//     }

//     for(size_t i=0; i<data.size(); i++) {
//         if (rv.t[i] > 57170.) {
//             rv.index_fibers = static_cast<uint32_t>(i);
//             break;
//         }
//     }
// }



// /**
//  * @brief Load RV data from a multiple files.
//  * 
//  * Read a tab/space separated files, each with columns  
//  * ```
//  *   time  vrad  error
//  *   ...   ...   ...
//  * ```
//  * All files should have the same structure and values in the same units.
//  * 
//  * @param filenames  the names of the files
//  * @param units      units of the RVs and errors, either "kms" or "ms"
//  * @param skip       number of lines to skip in the beginning of the file (default = 2)
//  * @param indicators
// */
// void Data::load_rv_multi(std::vector<std::string> filenames, const std::string units, uint32_t skip, const std::vector<std::string>& indicators) {
//     data_t data;

//     // Empty the vectors
//     rv.t.clear();
//     rv.y.clear();
//     rv.y_err.clear();
//     rv.instrument_ids.clear();

//     std::string dump; // to dump the first skip lines of each file
//     uint32_t filecount = 1;
//     size_t last_file_size = 0;

//     // Read the files into the data container
//     for (auto &filename : filenames) {
//         std::ifstream infile(filename);
//         for (int i=0; i<skip; i++) { // skip the first `skip` lines of each file
//             getline(infile, dump);
//         }
//         infile >> data;

//         // Complain if something went wrong.
//         if (!infile.eof()) {
//             std::cout <<"Could not read data file (" << filename << ")!" << std::endl;
//             exit(1);
//         }
//         infile.close();

//         // Assign instrument int identifier to obsi
//         for(size_t i=last_file_size; i<data.size(); i++) {
//             rv.instrument_ids.push_back(filecount);
//         }

//         last_file_size = data.size();
//         filecount++;
//     }

//     for (auto filename : filenames) {
//         rv.filenames.push_back(filename);
//     }
//     rv.units = units;
//     rv.skip = skip;

//     double factor = 1.;
//     if (units == "kms") factor = 1E3;

//     for (size_t n=0; n<data.size(); n++) {
//         // if (n<skip) continue;
//         rv.t.push_back(data[n][0]);
//         rv.y.push_back(data[n][1] * factor);
//         rv.y_err.push_back(data[n][2] * factor);
//     }

//     // How many points did we read?
//     std::cout << "# Loaded " << rv.t.size() << " data points from all files" << std::endl;
//     std::cout << "# ";
//     for (auto f: filenames) {
//         std::cout << f << " ; ";
//     }
//     std::cout << std::endl;

//     // Of how many instruments?
//     std::set<int> s( rv.instrument_ids.begin(), rv.instrument_ids.end() );
//     // set<int>::iterator iter;
//     // for(iter=s.begin(); iter!=s.end();++iter) {  cout << (*iter) << endl;}
//     std::cout << "# RVs come from " << s.size() << " different instruments." << std::endl;
//     rv.number_instruments = static_cast<uint32_t>(s.size());

//     if(units == "kms") {
//         std::cout << "# Multiplied all RVs by 1000; units are now m/s." << std::endl;
//     }

//     if (rv.number_instruments > 1) {
//         // We need to sort t because it comes from different instruments
//         size_t N = rv.t.size();
//         std::vector<double> tt(N), yy(N);
//         std::vector<double> sigsig(N);
//         std::vector<int32_t> obsiobsi(N);
//         std::vector<uint32_t> order(N);

//         // order = argsort(t)
//         uint32_t x=0;
//         std::iota(order.begin(), order.end(), x++);
//         sort( order.begin(),order.end(), [&](uint32_t i,uint32_t j){return rv.t[i] < rv.t[j];} );

//         for(size_t i=0; i<N; i++) {
//             tt[i] = rv.t[order[i]];
//             yy[i] = rv.y[order[i]];
//             sigsig[i] = rv.y_err[order[i]];
//             obsiobsi[i] = rv.instrument_ids[order[i]];
//         }

//         for(size_t i=0; i<N; i++) {
//             rv.t[i] = tt[i];
//             rv.y[i] = yy[i];
//             rv.y_err[i] = sigsig[i];
//             rv.instrument_ids[i] = obsiobsi[i];
//         }

//         // debug
//         // for(std::vector<uint32_t>::size_type i = 0; i != t.size(); i++) {
//         //     std::cout << rv.t[i] << "\t" << rv.y[i] << "\t" << rv.y_err[i] << "\t" << rv.instrument_ids[i] <<  std::endl;
//         // }
//     }

//     for(size_t i=0; i<data.size(); i++) {
//         if (rv.t[i] > 57170.) {
//             rv.index_fibers = static_cast<uint32_t>(i);
//             break;
//         }
//     }
// }

// TODO: Right now doesnt handle any of the extra data (multi_instrument, slope, etc.)
void Data::load_rv(const std::vector<std::string>& filenames, const std::string& units, uint32_t skip) {
    // Get data from files
    std::vector<std::vector<double>> data;
    for (std::string filename : filenames) {
        read_csv(filename, data);
        rv.filenames.push_back(filename);
    }

    rv.units = units;
    rv.skip = skip;
    rv.number_instruments = 1;

    // Prepare the vectors for the data
    size_t num_points = data.size() - skip;
    rv.t.reserve(num_points);
    rv.y.reserve(num_points);
    rv.y_err.reserve(num_points);

    double factor = 1.0;
    if (units == "kms") factor = 1E3;
    // Fill the arrays
    for (size_t n = skip; n < data.size(); n++) {
        // if (n<skip) continue;
        rv.t.push_back(data[n][0]);
        rv.y.push_back(data[n][1]);
        rv.y_err.push_back(data[n][2]);
    }

    rv.has_rv = true;
}

void Data::load_transit(const std::vector<std::string>& filenames, const std::string& units, uint32_t skip) {
    // Get data from files
    std::vector<std::vector<double>> data;
    for (std::string filename : filenames) {
        read_csv(filename, data);
        transit.filenames.push_back(filename);
    }
    transit.units = "frac";
    transit.skip = skip;

    // Prepare the vectors for the data
    size_t num_points = data.size() - skip;
    transit.t.reserve(num_points);
    transit.y.reserve(num_points);
    transit.y_err.reserve(num_points);

    // Fill the arrays
    for (size_t n = skip; n < data.size(); n++) {
        // if (n<skip) continue;
        transit.t.push_back(data[n][0]);
        transit.y.push_back(data[n][1]);
        transit.y_err.push_back(data[n][2]);
    }

    transit.has_transit = true;
}


double Data::rv_t_min() const { 
    return *std::min_element(rv.t.begin(), rv.t.end());
}
double Data::rv_t_max() const {
    return *std::max_element(rv.t.begin(), rv.t.end());
}
double Data::rv_t_middle() const {
    return rv_t_min() + 0.5*(rv_t_max() - rv_t_min());
}
double Data::rv_timespan() const {
    return rv_t_max() - rv_t_min();
}
double Data::rv_y_min() const {
    return *std::min_element(rv.y.begin(), rv.y.end());
}
double Data::rv_y_max() const {
    return *std::max_element(rv.y.begin(), rv.y.end());
}
double Data::rv_y_span() const {
    return rv_y_max() - rv_y_min();
}
double Data::rv_y_var() const {
    double sum = std::accumulate(std::begin(rv.y), std::end(rv.y), 0.0);
    double mean =  sum / rv.y.size();

    double accum = 0.0;
    std::for_each (std::begin(rv.y), std::end(rv.y), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });
    return accum / (rv.y.size()-1);
}
double Data::rv_y_std() const {
    return std::sqrt(rv_y_var());
}


double Data::rv_top_slope() const {
    return std::abs(rv_y_max() - rv_y_min()) / (rv.t.back() - rv.t.front());
}