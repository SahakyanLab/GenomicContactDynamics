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
Rcpp::NumericVector rcpp_anv_oneloop_encode( Rcpp::CharacterVector char_string ) {

    // Create mapping between base-pairs important vectors
    std::map<char, int> base_index_map = {
        { 'A', 0 },
        { 'C', 1 },
        { 'G', 2 },
        { 'T', 3 }
    };
    float num_bases = 4;

    // Convert R character vector to string
    std::string dna_string = Rcpp::as<std::string>(char_string);
    float sequence_length = dna_string.size();

    // Create storage for all the outputs
    Rcpp::NumericVector all_base_counts(num_bases);
    Rcpp::NumericVector aif_sums(num_bases);

    // Iterate along the sequence
    Rcpp::NumericMatrix AIF(num_bases, sequence_length);
    for (int ind=0; ind<sequence_length; ++ind) {
        int this_base_index = base_index_map[dna_string[ind]];
        all_base_counts[this_base_index] += 1;
        AIF(_,ind) = all_base_counts;
        aif_sums += all_base_counts;
    }

    // Average the positions from start
    Rcpp::NumericVector average_position = aif_sums / all_base_counts;

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
        base_covariances(AIF(0,_), AIF(0,_)), // A
        base_covariances(AIF(1,_), AIF(1,_)), // G
        base_covariances(AIF(2,_), AIF(2,_)), // C
        base_covariances(AIF(3,_), AIF(3,_)), // T
        base_covariances(AIF(0,_), AIF(1,_)), // AC
        base_covariances(AIF(0,_), AIF(2,_)), // AG
        base_covariances(AIF(0,_), AIF(3,_)), // AT
        base_covariances(AIF(1,_), AIF(2,_)), // CG
        base_covariances(AIF(1,_), AIF(3,_)), // CT
        base_covariances(AIF(2,_), AIF(3,_))  // GT
    );
    return anv;
}
