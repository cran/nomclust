#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> of_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_abs, const double freq_abs_r, const double n, const std::vector<double>& w, const double sum_w) {
	std::vector<double> of(r * r, 0);
	double agreement_sum = 0.0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
				int d = static_cast<int>(dat_vec[j + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
					agreement_sum = agreement_sum + w[k];
				} 
				else {
					agreement_sum = agreement_sum + (w[k]/(1.0 + log(n/freq_abs[c + (freq_abs_r * k)]) * log(n/freq_abs[d + (freq_abs_r * k)])));
				}
				
			}
			of[i + (j * r)] = (1.0/(1.0 / sum_w * agreement_sum)) - 1.0;
			of[j + (i * r)] = of[i + (j * r)];
			agreement_sum = 0.0;
		}
	}
	return of;
}
