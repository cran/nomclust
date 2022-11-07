#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> ve_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& norm_entropy, const std::vector<double>& w, const double sum_w) {
	std::vector<double> vm(r * r, 0);
	double agreement_sum = 0.0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
					agreement_sum = agreement_sum + (norm_entropy[k] * w[k]);
				}
			}
			if (i == j) {
				vm[i + (j * r)] = 0.0;
			}
			else {
				vm[i + (j * r)] = (1.0 - (1.0 / sum_w * agreement_sum));
				vm[j + (i * r)] = vm[i + (j * r)];
			}
			agreement_sum = 0.0;
		}
	}
	return vm;
}
