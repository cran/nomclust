#include <Rcpp.h>
#include <iostream>
#include <vector> 
#include <fstream>

// [[Rcpp::export]]
std::vector<double> burnaby_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w) {
	std::vector<double> burnaby(r * r, 0);
	double agreement_sum = 0.0;
	double nominator = 0.0;
	for (int i = 0; i < (r - 1); ++i) {
		for (int j = (1 + i); j < r; ++j) {
			for (int k = 0; k < s; ++k) {
				int c = static_cast<int>(dat_vec[i + (r * k)] - 1);
				int d = static_cast<int>(dat_vec[j + (r * k)] - 1);
				if (dat_vec[i + (r * k)] == dat_vec[j + (r * k)]) {
					agreement_sum = agreement_sum + w[k];
				} 
				else {
				  for (int ctg = 0; ctg < freq_rel_r; ++ctg) {
				    nominator = nominator + (2.0*log(1.0-freq_rel[ctg + (freq_rel_r * k)]));
				  }
				  
					agreement_sum = agreement_sum + (w[k]*(nominator/(log( freq_rel[c + (freq_rel_r * k)] * freq_rel[d + (freq_rel_r * k)] / (1.0-freq_rel[c + (freq_rel_r * k)]) / (1.0-freq_rel[d + (freq_rel_r * k)])) + nominator)));
				  nominator = 0;
				}
				
			}
			if(i==j) {
			  burnaby[i + (j * r)] = 0.0;
			}
			else {
			  burnaby[i + (j * r)] = 1.0- (agreement_sum / sum_w);
		    burnaby[j + (i * r)] = burnaby[i + (j * r)];
			}
			agreement_sum = 0.0;
		}
	}
	return burnaby;
}
