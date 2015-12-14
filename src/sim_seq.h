#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;

class twosubtrees {
public:
    double left_length, right_length;
    // lengths of the branches going to the left and the right subtree

    string left, right;
    // the left and the right subtree in newick format

    bool left_is_leaf, right_is_leaf;
    // are the subtrees just leaves?
};

class F84 {
    // parameters of the F84 substitution model, with a branch length and transition probabilities
public:
    vector<double> nucfreq; // nucleotide freqency A,C,G,T
    double ti_rate, tv_rate; // transition and transversion rate
    double t; // branch length
    F84(double freqA, double freqC, double freqG, double ti_rate, double tv_rate);
    F84(double freqA, double freqC, double freqG, double ti_rate, double tv_rate, double t);
    inline void set_t(double nt) {t=nt;};
};

twosubtrees  split_tree(const string & tree);
// tree must be a binary tree in newick format

void continue_simulation(string tree, const vector<unsigned> & rootseq, F84 f84,
			 vector<vector<unsigned> > * seq, vector<string> * tax);
// simulates sequences along the tree, starting with rootseq. Will return
// simulated sequences bei putting them coded as vectors over {0,1,2,3} in
// &seq. Taxa will go into tax.

vector<unsigned> simulate_branch(const vector<unsigned> & rootseq,F84 f84);

pair<vector<string>,vector<vector<unsigned> > > simulate_sequences_F84(string tree, unsigned length, double freqA, double freqC, double freqG,
								     double tirate, double tvrate, bool only_segregating);

inline unsigned int sample_int(const unsigned int max) {
  return(static_cast<int>(R::runif(0, 1) * max));
}
