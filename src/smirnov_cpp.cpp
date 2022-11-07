#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> smirnov_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_abs, const double freq_abs_r, const double num_cat_sum, const std::vector<double>& w, const double sum_w) {
	std::vector<double> smirnov(r * r, 0);
	double agreement_sum = 0.0;
	double dominator_sum = 0.0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
			  int d = static_cast<int>(dat_vec[j + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
				  for(int cc = 0; cc < freq_abs_r; ++cc) {
				    if(cc!=c) {
				      dominator_sum = dominator_sum + freq_abs[cc + (freq_abs_r * k)] / (r - freq_abs[cc + (freq_abs_r * k)]);
				    }
				  }
				  agreement_sum = agreement_sum + (w[k]*(2.0+(r-freq_abs[c + (freq_abs_r * k)])/freq_abs[c + (freq_abs_r * k)] + dominator_sum));
				  dominator_sum = 0.0;
				} else {
				  for(int cc = 0; cc < freq_abs_r; ++cc) {
				    if(cc!=c && cc!=d) {
				      dominator_sum = dominator_sum + freq_abs[cc + (freq_abs_r * k)] / (r - freq_abs[cc + (freq_abs_r * k)]);
				    }
				  }
				  agreement_sum = agreement_sum + dominator_sum*w[k];
				  dominator_sum = 0.0;
				}
				
			}
			smirnov[i + (j * r)] = 1.0/(1.0/num_cat_sum*(agreement_sum)+1.0);
			smirnov[j + (i * r)] = smirnov[i + (j * r)];
			agreement_sum = 0.0;
		}
	}
	return smirnov;
}
