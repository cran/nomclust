#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> gambaryan_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const double num_cat_sum, const std::vector<double>& w, const double sum_w) {
	std::vector<double> gambaryan(r * r, 0);
	double agreement_sum = 0.0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
				  agreement_sum = agreement_sum - (w[k]*(freq_rel[c + (freq_rel_r * k)] * log2(freq_rel[c + (freq_rel_r * k)]) + (1.0- freq_rel[c + (freq_rel_r * k)]) * log2(1.0-freq_rel[c + (freq_rel_r * k)])));
				} 
				
			}
			if(i==j) {
			  gambaryan[i + (j * r)] = 0.0;
			}
			else {
			  gambaryan[i + (j * r)] = 1.0- (agreement_sum / num_cat_sum);
			  gambaryan[j + (i * r)] = gambaryan[i + (j * r)];
			}
			agreement_sum = 0.0;
		}
	}
	return gambaryan;
}