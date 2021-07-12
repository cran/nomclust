#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> lin_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w) {
	std::vector<double> lin(r * r, 0);
	double agreement_sum = 0;
	double lin_weight_sum = 0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
				int d = static_cast<int>(dat_vec[j + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
					agreement_sum = agreement_sum + (w[k]*(2*log(freq_rel[c + (freq_rel_r * k)])));
				}
				else {
					agreement_sum = agreement_sum + (w[k] * (2 * log(freq_rel[c + (freq_rel_r * k)] + freq_rel[d + (freq_rel_r * k)])));
				}
				lin_weight_sum = lin_weight_sum + ((log(freq_rel[c + (freq_rel_r * k)]) + log(freq_rel[d + (freq_rel_r * k)]))*w[k]);
			}
			if (i == j) {
				lin[j + (i * r)] = 0;
			} 
			else {
				lin[i + (j * r)] = (1/(1 / lin_weight_sum * agreement_sum )) - 1;
				lin[j + (i * r)] = lin[i + (j * r)];
			}
			agreement_sum = 0;
			lin_weight_sum = 0;
		}
	}
	return lin;
}
