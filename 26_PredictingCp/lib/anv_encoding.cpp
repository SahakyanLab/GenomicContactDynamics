#include <Rcpp.h>
#include <map>

using namespace Rcpp;
using namespace std;


float base_covariances(Rcpp::NumericVector aif_1, Rcpp::NumericVector aif_2) {
    // Calculate covariances between two AIF rows repesenting nucleotide distribution
    float this_seq_len = aif_1.size();
    float theta_1 = accumulate(aif_1.begin(), aif_1.end(), 0.0) / this_seq_len;
    float theta_2 = accumulate(aif_2.begin(), aif_2.end(), 0.0) / this_seq_len;

    // Denominator: Product of base counts
    float denom = aif_1[this_seq_len-1] * aif_2[this_seq_len-1];
    
    // Numerator: Sum of products of (aif - theta) for both bases
    float numerator = 0;
    for (int i = 0; i < this_seq_len; i++) {
        numerator += (aif_1[i] - theta_1) * (aif_2[i] - theta_2);
    }

    float covariance = numerator / denom;
    return covariance;
}


// [[Rcpp::export]]
Rcpp::NumericVector rcpp_anv_encode( Rcpp::CharacterVector char_string ) {

    // Create mapping between base-pairs and corners to move towards
    std::map<char, Rcpp::NumericVector> base_onehot = {
        { 'A', {1., 0., 0., 0.} },
        { 'C', {0., 1., 0., 0.} },
        { 'G', {0., 0., 1., 0.} },
        { 'T', {0., 0., 0., 1.} }
    };

    // Convert R character vector to string
    std::string dna_string = Rcpp::as<std::string>(char_string);

    // Extract length
    float sequence_length = dna_string.size();
    float base_count = 4;

    // Generate the indicator functions (IF) and accumulated IF (AIF)
    Rcpp::NumericMatrix IF(4, sequence_length);
    Rcpp::NumericMatrix padded_AIF(4, sequence_length+1);
    for (int ind=0; ind<sequence_length; ++ind) {
        // Get current character
        char current_base = dna_string[ind];
        // Get corresponding one-hot vector
        Rcpp::NumericVector onehot = base_onehot[current_base];
        // Update IF and padded_AIF
        IF(_,ind) = onehot;
        padded_AIF(_,ind+1) = padded_AIF(_,ind) + onehot;
    }
    Rcpp::NumericMatrix AIF = padded_AIF(_,Range(1,sequence_length));

    // Extract base counts
    Rcpp::NumericVector all_base_counts = AIF(_,sequence_length-1);
		
    // Extract average distances from end -> sum(AIF) / base count
	Rcpp::NumericVector average_distance(4);
	for (int base=0; base<base_count; ++base)  {
		Rcpp::NumericVector aif_base = AIF(base,_);
		float base_ave_dist = accumulate(aif_base.begin(), aif_base.end(), 0.0) / all_base_counts[base];
	    average_distance[base] = base_ave_dist;
	}
	
    // Return combined Accumulated Natural Vector
    Rcpp::NumericVector anv = NumericVector::create(
        all_base_counts[0],
        all_base_counts[1],
        all_base_counts[2],
        all_base_counts[3],
		
        average_distance[0],  							
        average_distance[1],
        average_distance[2],
        average_distance[3],
		
        base_covariances(AIF(0,_), AIF(0,_)), // A divergences[0],
        base_covariances(AIF(1,_), AIF(1,_)), // C divergences[1],
        base_covariances(AIF(2,_), AIF(2,_)), // G divergences[2],
        base_covariances(AIF(3,_), AIF(3,_)), // T divergences[3],
		
        base_covariances(AIF(0,_), AIF(1,_)), // AC
        base_covariances(AIF(0,_), AIF(2,_)), // AG
        base_covariances(AIF(0,_), AIF(3,_)), // AT
        base_covariances(AIF(1,_), AIF(2,_)), // CG
        base_covariances(AIF(1,_), AIF(3,_)), // CT
        base_covariances(AIF(2,_), AIF(3,_))  // GT
    );
    return anv;
	
}
