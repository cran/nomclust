#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> anderberg_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& num_cat, const std::vector<double>& w, const double sum_w) {
	std::vector<double> anderberg(r * r, 0);
	double nominator_sum = 0.0;
	double denominator_sum = 0.0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
			  int d = static_cast<int>(dat_vec[j + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
				  nominator_sum = nominator_sum + ((1.0/freq_rel[c + (freq_rel_r * k)])*(1.0/freq_rel[c + (freq_rel_r * k)]) * 2.0 / num_cat[k] / (num_cat[k]+1.0));
				} 
				else {
				  denominator_sum = denominator_sum + (1.0/2.0/freq_rel[c + (freq_rel_r * k)]/freq_rel[d + (freq_rel_r * k)] * 2.0  / num_cat[k] / (num_cat[k]+1.0));
				}
				
			}
			
			
			
			if(i==j) {
			  anderberg[i + (j * r)] = 0.0;
			}
			else {
			  anderberg[i + (j * r)] = 1.0 - (nominator_sum / (nominator_sum+denominator_sum));
			  anderberg[j + (i * r)] = anderberg[i + (j * r)];
			}
			nominator_sum = 0.0;
			denominator_sum = 0.0;
		}
	}
	return anderberg;
}
