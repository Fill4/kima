// (c) 2019 João Faria
// This file is part of kima, which is licensed under the MIT license (see LICENSE for details)

#pragma once

#include <vector>
#include <string>

struct Data
{
	struct {
		bool has_rv = false;

		std::vector<std::string> filenames;
		std::string units;
		uint32_t skip;

		std::vector<double> t, y, y_err;
		std::vector<int32_t> instrument_ids;
		std::vector<std::vector<double>> activity_indicators;

		uint32_t index_fibers;
		
		bool indicator_correlations;
		uint32_t number_indicators;
		std::vector<std::string> indicator_names;

		uint32_t number_instruments;
	} rv;

	struct {
		bool has_transit = false;

		std::vector<std::string> filenames;
		std::string units;
		uint32_t skip;

		std::vector<double> t, y, y_err;
	} transit;

	void read_csv(const std::string& filename, std::vector<std::vector<double>>& data);
	// Load up the RV data from a 3 column file. Can handle multiple instruments from multiple files
	void load_rv(const std::vector<std::string>& filenames, const std::string& units, uint32_t skip=2);
	// Load up the transit data from a 3 column file
	void load_transit(const std::vector<std::string>& filenames, const std::string& units, uint32_t skip=2);

	/// @brief Get the mininum (starting) RV time @return double
	double rv_t_min() const;
	/// @brief Get the maximum (ending) RV time @return double
	double rv_t_max() const;
	/// @brief Get the middle (average) RV time @return double
	double rv_t_middle() const;
	/// @brief Get the RV timespan @return double
	double rv_timespan() const;

	/// @brief Get the mininum RV @return double
	double rv_y_min() const;
	/// @brief Get the maximum RV @return double
	double rv_y_max() const;
	/// @brief Get the RV span @return double
	double rv_y_span() const;
	/// @brief Get the variance of the RVs @return double
	double rv_y_var() const;
	/// @brief Get the standard deviation of the RVs @return double
	double rv_y_std() const;
	
	/// @brief Get the RV span, adjusted for multiple instruments @return double
	double rv_adjusted_y_span() const;
	/// @brief Get the RV variance, adjusted for multiple instruments @return double
	double rv_adjusted_y_var() const;
	/// @brief Get the RV standard deviation, adjusted for multiple instruments @return double
	double rv_adjusted_y_std() const;
	
	/// @brief Get the maximum slope allowed by the data. @return double
	double rv_top_slope() const;

	// Singleton
	static Data& instance() { static Data theData; return theData; }
	// Destroy all ways to create another Data instance
	private: 
		Data() {} // Private constructor
	public:
		Data(Data const&) = delete;
        void operator=(Data const&) = delete;
};
