#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> lin1_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w) {
	std::vector<double> lin1(r * r, 0);
	std::vector<double> logic(s, 0);
	double agreement_sum = 0.0;
	double lin1_weight_sum = 0.0;
	double logic_sum = 0.0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
				int d = static_cast<int>(dat_vec[j + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
					for (int ctg = 0; ctg < freq_rel_r; ++ctg) {
						if (freq_rel[ctg + (freq_rel_r * k)] == freq_rel[c + (freq_rel_r * k)]) {
							logic_sum = logic_sum + freq_rel[ctg + (freq_rel_r * k)];
						}
					}
					agreement_sum = agreement_sum + (w[k]*(log(logic_sum)));
				}
				else {
					if (freq_rel[c + (freq_rel_r * k)] >= freq_rel[d + (freq_rel_r * k)]) { 
						for (int ctg = 0; ctg < freq_rel_r; ++ctg) {
							if ((freq_rel[ctg + (freq_rel_r * k)] >= freq_rel[d + (freq_rel_r * k)]) && (freq_rel[ctg + (freq_rel_r * k)] <= freq_rel[c + (freq_rel_r * k)])) {
								logic_sum = logic_sum + freq_rel[ctg + (freq_rel_r * k)];
							}
						}
						agreement_sum = agreement_sum + (w[k]*(2*log(logic_sum)));
					}
					else {
						for (int ctg = 0; ctg < freq_rel_r; ++ctg) {
							if ((freq_rel[ctg + (freq_rel_r * k)] >= freq_rel[c + (freq_rel_r * k)]) && (freq_rel[ctg + (freq_rel_r * k)] <= freq_rel[d + (freq_rel_r * k)])) {
								logic_sum = logic_sum + freq_rel[ctg + (freq_rel_r * k)];
							}
						}
						agreement_sum = agreement_sum + (w[k]*(2.0*log(logic_sum)));
					}
				}
				lin1_weight_sum = lin1_weight_sum + ((log(freq_rel[c + (freq_rel_r * k)]) + log(freq_rel[d + (freq_rel_r * k)]))*w[k]);
				logic_sum = 0.0;
			}
			//std::cout << agreement_sum << std::endl;
			lin1[i + (j * r)] = (1.0/(1.0 / lin1_weight_sum * agreement_sum )) - 1.0;
			lin1[j + (i * r)] = lin1[i + (j * r)];
			agreement_sum = 0.0;
			lin1_weight_sum = 0.0;
		}
	}
	return lin1;
}
