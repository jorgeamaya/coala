#include "sim_seq.h"

F84::F84(double freqA, double freqC, double freqG,
         double nti_rate,  double ntv_rate) :
  ti_rate(nti_rate), tv_rate(ntv_rate), t(0) {
  nucfreq = std::vector<double>(4);
  nucfreq[0] = freqA;
  nucfreq[1] = freqC;
  nucfreq[2] = freqG;
  nucfreq[3] = 1 - freqA - freqC - freqG;
}

F84::F84(double freqA, double freqC, double freqG,
         double nti_rate,  double ntv_rate,  double nt) :
  ti_rate(nti_rate), tv_rate(ntv_rate) {
  nucfreq = std::vector<double>(4);
  nucfreq[0] = freqA;
  nucfreq[1] = freqC;
  nucfreq[2] = freqG;
  nucfreq[3] = 1 - freqA - freqC - freqG;
  set_t(nt);
}

twosubtrees  split_tree(const string & tree) {
  // tree must be a binary tree in newick format

  twosubtrees ts;
  int par_counter=0; // keeps track of level of parentheses, assuming
  // tree[0]=='(', not counting for the very outside
  ts.left_is_leaf=true;
  ts.right_is_leaf=true;
  unsigned int treebegin=0, treeend=0;
  bool left_side=true; // are we still in the left tree?
  for(unsigned i=1; i<tree.length(); ++i) {
    switch(tree[i]) {
    case '(' :
           ++par_counter;
      if(left_side)
        ts.left_is_leaf=false;
      else {
        ts.right_is_leaf=false;
      }
      break;
    case ')' :
      --par_counter;
      if(par_counter<0) {
        // we are at the very end
        ts.right_length=atof(tree.substr(treeend+1,i-treeend).c_str());
        ts.right=tree.substr(treebegin,treeend-treebegin);
        return ts;
      }
      break;
    default:
      if(par_counter==0) {
        switch(tree[i]) {
        case ':' :
          treeend=i;
          break;
        case ',' :
          ts.left_length=atof(tree.substr(treeend+1,i-treeend).c_str());
          ts.left=tree.substr(1,treeend-1);
          treebegin=i+1;
          left_side=false;
          break;
        default:
          break;
        }
      }
      break;
    }
  }
  return ts; // should actually never be reached
};

vector<unsigned> simulate_branch(const vector<unsigned> & rootseq,F84 f84) {
  vector<unsigned> newseq(rootseq);
  const unsigned n=rootseq.size();
  //poisson_distribution<> randtv(n*f84.t*f84.tv_rate), randti(n*f84.t*f84.ti_rate);
  //uniform_int_distribution<> r_unif_n(0,n-1);
  //uniform_real_distribution<> r_unif(0,1);

  const unsigned rtv = R::rpois(n*f84.t*f84.tv_rate),
                 rti = R::rpois(n*f84.t*f84.ti_rate);

  if(f84.t*(f84.tv_rate+f84.ti_rate)>0.4) {
    vector<bool> mutated(rootseq.size(),false);
    for(unsigned i=0; i<rtv; ++i) {
      // pepper potential transversions into the sequence
      unsigned site = sample_int(n);
      if(!mutated[site]) {
        double re = R::runif(0, 1) - f84.nucfreq[0];
        unsigned b=0;
        while(re>0) {
          re-=f84.nucfreq[++b];
        }
        newseq[site]=b;
        mutated[site]=true;
      }
    }
    for(unsigned i=0; i<rti; ++i) {
      // pepper potential transitions into the sequence
      unsigned site= sample_int(n);
      if(!mutated[site]) {
        switch(newseq[site]) {
        case 0 :
          if(R::runif(0, 1) <f84.nucfreq[2]/(f84.nucfreq[0]+f84.nucfreq[2])) {
            newseq[site]=2;
          }
          break;
        case 1:
          if(R::runif(0, 1) <f84.nucfreq[3]/(f84.nucfreq[1]+f84.nucfreq[3])) {
            newseq[site]=3;
          }
          break;
        case 2:
          if(R::runif(0, 1) <f84.nucfreq[0]/(f84.nucfreq[0]+f84.nucfreq[2])) {
            newseq[site]=0;
          }
          break;
        case 3:
          if(R::runif(0, 1) <f84.nucfreq[3]/(f84.nucfreq[1]+f84.nucfreq[3])) {
            newseq[site]=3;
          }
          break;
        default:
          Rcpp::stop("unknown nucleotide in vector<unsigned> simulate_branch(const vector<unsigned>&,F84)\n");
        break;
        }
        mutated[site]=true;
      }
    }
  } else {
    for(unsigned i=0; i<rtv; ++i) {
      // pepper potential transversions into the sequence
      unsigned site=sample_int(n);
      double re=R::runif(0, 1) -f84.nucfreq[0];
      unsigned b=0;
      while(re>0) {
        re-=f84.nucfreq[++b];
      }
      newseq[site]=b;
    }
    for(unsigned i=0; i<rti; ++i) {
      // pepper potential transitions into the sequence
      unsigned site=sample_int(n);
      switch(newseq[site]) {
      case 0 :
        if(R::runif(0, 1) <f84.nucfreq[2]/(f84.nucfreq[0]+f84.nucfreq[2])) {
          newseq[site]=2;
        }
        break;
      case 1:
        if(R::runif(0, 1) <f84.nucfreq[3]/(f84.nucfreq[1]+f84.nucfreq[3])) {
          newseq[site]=3;
        }
        break;
      case 2:
        if(R::runif(0, 1) <f84.nucfreq[0]/(f84.nucfreq[0]+f84.nucfreq[2])) {
          newseq[site]=0;
        }
        break;
      case 3:
        if(R::runif(0, 1) <f84.nucfreq[3]/(f84.nucfreq[1]+f84.nucfreq[3])) {
          newseq[site]=3;
        }
        break;
      default:
        Rcpp::stop("unknown nucleotide in vector<unsigned> simulate_branch(const vector<unsigned>&,F84)\n");
      break;
      }
    }
  }
  return newseq;
}

void continue_simulation(string tree, const vector<unsigned> & rootseq, F84 f84,
                         vector<vector<unsigned> > * seq, vector<string> * tax) {
  // simulates sequences along the tree, starting with rootseq. Will return
  // simulated sequences bei putting them coded as vectors over {0,1,2,3} in
  // &seq. Taxa will go into tax.

  twosubtrees ts=split_tree(tree);
  f84.set_t(ts.left_length);
  vector<unsigned> newseq=simulate_branch(rootseq,f84);
  if(ts.left_is_leaf) {
    seq -> push_back(newseq);
    tax -> push_back(ts.left);
  } else {
    continue_simulation(ts.left,newseq,f84,seq,tax);
  }
  f84.set_t(ts.right_length);
  newseq=simulate_branch(rootseq,f84);
  if(ts.right_is_leaf) {
    seq -> push_back(newseq);
    tax -> push_back(ts.right);
  } else {
    continue_simulation(ts.right,newseq,f84,seq,tax);
  }
}

pair<vector<string>,vector<vector<unsigned> > > simulate_sequences_F84(string tree, unsigned length,
                                                                     double freqA, double freqC, double freqG,
                                                                     double tirate, double tvrate,
                                                                     bool only_segregating) {
  F84 f84(freqA, freqC, freqG, tirate, tvrate);
  vector<vector<unsigned> > seq(0);
  vector<string> tax(0);
  vector<unsigned> rootseq(length);
  for(unsigned i=0; i<length; ++i) {
    double re=R::runif(0, 1) -freqA;
    if(re>0) {
      re-=freqC;
      if(re>0) {
        re-=freqG;
        if(re>0) rootseq[i]=3;
        else rootseq[i]=2;
      } else rootseq[i]=1;
    } else rootseq[i]=0;
  }
  continue_simulation(tree, rootseq, f84, &seq, &tax);
  if(only_segregating) {
    vector<vector<unsigned> > s(seq.size());
    for(unsigned i=0; i<length; ++i) {
      unsigned b=seq[0][i];
      for(unsigned j=1; j<seq.size(); ++j) {
        if(b!=seq[j][i]) {
          for(unsigned k=0; k<seq.size(); ++k)
            s[k].push_back(seq[k][i]);
          break;
        }
      }
    }
    return pair<vector<string>,vector<vector<unsigned> > >(tax, s);
  }
  return pair<vector<string>,vector<vector<unsigned> > >(tax, seq);
}


// [[Rcpp::export]]
Rcpp::NumericVector sim_seq(const Rcpp::CharacterVector trees) {
  Rcpp::RNGScope scope;

  std::string line, tree;
  size_t length, digits;
  pair<vector<string>, vector<vector<unsigned> > > p;

  for (int i = 0; i < trees.size(); ++i) {
    line = Rcpp::as<std::string>(trees.at(i));
    digits = line.find("]")-1;
    length = std::atoi(line.substr(1, digits).c_str());
    tree = line.substr(digits+2, std::string::npos);

    p = simulate_sequences_F84(tree, length, 0.2, 0.3, 0.2, 0.01, 0.01, false);
    for(unsigned i=0; i<p.first.size(); ++i) {
      Rcpp::Rcout << p.first[i] << "  ";
      for(unsigned j=0; j<p.second[i].size(); ++j) {
        Rcpp::Rcout << p.second[i][j];
      }
      Rcpp::Rcout << endl;
    }
  }

  return(Rcpp::NumericVector(0));
}
