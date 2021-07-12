#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> good3_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w) {
	std::vector<double> good3(r * r, 0);
	std::vector<double> logic(s, 0);
	double agreement_sum = 0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
					agreement_sum = agreement_sum + (w[k]*(1 - (freq_rel[c + (freq_rel_r * k)]*freq_rel[c + (freq_rel_r * k)])));
				}
			}
			if (i == j) {
				good3[j + (i * r)] = 0;
			}
			else {
				good3[i + (j * r)] = (1 - (1 / sum_w * agreement_sum));
				good3[j + (i * r)] = good3[i + (j * r)];
			}
			agreement_sum = 0;
		}
	}
	return good3;
}
