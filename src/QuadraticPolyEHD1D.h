#ifndef QuadraticPoly_EHD1D_H_
#define QuadraticPoly_EHD1D_H_

#include "./QuadraticPoly1D.h"

namespace PCMBaseCpp {

typedef SPLITT::OrderedTree<SPLITT::uint, LengthAndRegime> EHD1DTreeType;

template<class TreeType, class DataType>
struct CondGaussianEHD1D: public CondGaussianOmegaPhiV1D {
  
  TreeType const& ref_tree_;
  
  // number of regimes;
  uint R_;
  
  bool transpose_Sigma_x = false;
  
  // EHD parameters for 1D case
  std::vector<double> X0;
  std::vector<double> rho;      // Early burst parameter
  std::vector<double> H;        // Selection strength (scalar)
  std::vector<double> Theta;    // Optimum value (scalar)
  std::vector<double> Sigma_x;  // Variance parameter
  std::vector<double> Sigmae_x; // Measurement error variance
  
  CondGaussianEHD1D(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->R_ = R;
    this->transpose_Sigma_x = ref_data.transpose_Sigma_x;
  }
  
  CondGaussianEHD1D(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->R_ = ref_data.R_;
    this->transpose_Sigma_x = ref_data.transpose_Sigma_x;
  }
  
  arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) {
    using namespace arma;
    
    uint npar = R_*(5);  // X0 + rho + H + Theta + Sigma_x + Sigmae_x
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os<<"QuadraticPolyEHD1D.h:CondEHD1D.SetParameter:: The length of the parameter vector minus offset ("<<par.size() - offset<<
        ") should be at least of R*5, where "+
          " R="<<R_<<" is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    X0.assign(R_, 0.0);
    rho.assign(R_, 0.0);
    H.assign(R_, 0.0);
    Theta.assign(R_, 0.0);
    Sigma_x.assign(R_, 0.0);
    Sigmae_x.assign(R_, 0.0);
    
    for(uword r = 0; r < R_; r++) {
      X0[r] = par[offset + r];
      rho[r] = par[offset + R_ + r];
      H[r] = par[offset + 2*R_ + r];
      Theta[r] = par[offset + 3*R_ + r];
      Sigma_x[r] = par[offset + 4*R_ + r];
      Sigmae_x[r] = par[offset + 5*R_ + r];
    }
    
    return npar;
  }
  
  void CalculateOmegaPhiV(uint i, arma::uword ri, double& omega, double& phi, double& V) {
    double t = this->ref_tree_.LengthOfBranch(i).length_;
    
    // Get node height for EB component (simplified for now)
    double node_height = 0.0; // TODO: Implement proper node height calculation
    
    // Early Burst scaling factor
    double eb_factor = (rho[ri] == 0.0) ? 1.0 : (1.0 - exp(-rho[ri] * node_height)) / rho[ri];
    double sigma2 = Sigma_x[ri] * Sigma_x[ri] * eb_factor;
    
    if(std::abs(H[ri]) < 1e-8) {
      // BM case with EB scaling
      phi = 1.0;
      omega = Theta[ri] * t;
      V = sigma2 * t;
    } else {
      // OU case with EB scaling
      double exp_Ht = exp(-H[ri] * t);
      phi = exp_Ht;
      omega = (1.0 - exp_Ht) * Theta[ri];
      V = sigma2 * (1.0 - exp(-2.0 * H[ri] * t)) / (2.0 * H[ri]);
    }
    
    // Add measurement error for tips
    if(i < this->ref_tree_.num_tips()) {
      V += Sigmae_x[ri] * Sigmae_x[ri];
    }
  }
};

class EHD1D: public QuadraticPoly1D<EHD1DTreeType> {
public:
  typedef EHD1DTreeType TreeType;
  typedef QuadraticPoly1D<TreeType> BaseType;
  typedef EHD1D MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData1D<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;
  
  CondGaussianEHD1D<TreeType, DataType> cond_dist_;
  
  EHD1D(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }
  
  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};

typedef TraversalTaskWrapper<EHD1D> QuadraticPolyEHD1D;
}

#endif // QuadraticPoly_EHD1D_H_