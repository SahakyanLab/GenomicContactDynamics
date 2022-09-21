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
Rcpp::NumericVector rcpp_anv_encode_J( Rcpp::CharacterVector char_string ) {

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
        char current_base = dna_string[ind];
        Rcpp::NumericVector onehot = base_onehot[current_base];
        IF(_,ind) = onehot;
        padded_AIF(_,ind+1) = padded_AIF(_,ind) + onehot;
    }
    Rcpp::NumericMatrix AIF = padded_AIF(_,Range(1,sequence_length));

    // Extract base counts
    Rcpp::NumericVector all_base_counts = AIF(_,sequence_length-1);

    // Extract average distances from end 
    // JAMES MODIFIED (COPY MATLAB)
    Rcpp::NumericVector omega(4);
    for (int ind=0; ind<sequence_length; ind++) {
        for(int inner_ind=0; inner_ind<=ind; inner_ind++) {
            for (int base=0; base<base_count; ++base)  {
                omega[base] += IF(base, inner_ind);
            }
        }
    }
    Rcpp::NumericVector average_distance = omega / all_base_counts;

    // LIEZELS
    // Rcpp::NumericVector average_distance(4);
	// for (int base=0; base<base_count; ++base)  {
	// 	Rcpp::NumericVector aif_base = AIF(base,_);
	// 	float base_ave_dist = accumulate(aif_base.begin(), aif_base.end(), 0.0) / all_base_counts[base];
	//     average_distance[base] = base_ave_dist;
	// }

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
        
        base_covariances(AIF(0,_), AIF(0,_)), // A
        base_covariances(AIF(1,_), AIF(1,_)), // C
        base_covariances(AIF(2,_), AIF(2,_)), // G
        base_covariances(AIF(3,_), AIF(3,_)), // T

        base_covariances(AIF(0,_), AIF(1,_)), // AC
        base_covariances(AIF(0,_), AIF(2,_)), // AG
        base_covariances(AIF(0,_), AIF(3,_)), // AT
        base_covariances(AIF(1,_), AIF(2,_)), // CG
        base_covariances(AIF(1,_), AIF(3,_)), // CG
        base_covariances(AIF(2,_), AIF(3,_))  // GT
    );
    return anv;
}



// [[Rcpp::export]]
Rcpp::NumericVector rcpp_anv_oneloop_encode( Rcpp::CharacterVector char_string ) {

    // Create mapping between base-pairs important vectors
    std::map<char, int> base_index_map = {
        { 'A', 0 },
        { 'C', 1 },
        { 'G', 2 },
        { 'T', 3 }
    };

    // Convert R character vector to string
    std::string dna_string = Rcpp::as<std::string>(char_string);

    // Extract length
    float sequence_length = dna_string.size();
    float num_bases = 4;

    // Create storage for all the outputs
    Rcpp::NumericVector all_base_counts(num_bases);
    Rcpp::NumericVector sum_base_distances(num_bases);
    Rcpp::NumericVector aif_sums(num_bases);

    // Iterate along the sequence
    Rcpp::NumericMatrix AIF(num_bases, sequence_length);
    for (int ind=0; ind<sequence_length; ++ind) {
        // Get current character
        int this_base_index = base_index_map[dna_string[ind]];

        // Update running total, and storage
        all_base_counts[this_base_index] += 1;
        AIF(_,ind) = all_base_counts;

        // Update distance and sum storage
        sum_base_distances[this_base_index] += ind;
        aif_sums += all_base_counts;
    }

    // Average the positions from start
    Rcpp::NumericVector average_position = sum_base_distances / all_base_counts;

    // Calculate the position-wise divergences
    Rcpp::NumericVector divergences(num_bases);
    for (int base=0; base<num_bases; ++base)  {
        Rcpp::NumericVector position_divergences =  AIF(base,_) - (aif_sums[base]/sequence_length);
        Rcpp::NumericVector divergence = pow(position_divergences / all_base_counts[base], 2);
        divergences[base] = accumulate(divergence.begin(), divergence.end(), 0.0);
    }

    // Return combined Accumulated Natural Vector
    Rcpp::NumericVector anv = NumericVector::create(
        all_base_counts[0],
        all_base_counts[1],
        all_base_counts[2],
        all_base_counts[3],
        average_position[0],
        average_position[1],
        average_position[2],
        average_position[3],
        divergences[0],
        divergences[1],
        divergences[2],
        divergences[3],
        base_covariances(AIF(0,_), AIF(1,_)), // AC
        base_covariances(AIF(0,_), AIF(2,_)), // AG
        base_covariances(AIF(0,_), AIF(3,_)), // AT
        base_covariances(AIF(1,_), AIF(2,_)), // CG
        base_covariances(AIF(1,_), AIF(3,_)), // CT
        base_covariances(AIF(2,_), AIF(3,_))  // GT
    );
    return anv;
}
