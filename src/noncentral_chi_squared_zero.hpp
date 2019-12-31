#ifndef __NONCENTRAL_CHI_SQUARED_ZERO_HPP__
#define __NONCENTRAL_CHI_SQUARED_ZERO_HPP__

// Include std::random library.
#include <random>

namespace std {

  template<class RealType=double> class noncentral_chi_squared_zero {
  public:
    noncentral_chi_squared_zero(RealType nc_param) 
      : noncentrality_parameter(nc_param), 
        dist_poisson(std::poisson_distribution(static_cast<RealType>(nc_param/2.))),
        dist_normal(std::normal_distribution<RealType>(0., 1.))
    {};

    template<class URNG>
    RealType operator()(URNG& generator) {
      // Generate the number of variables to use.
      int n_normal_variables = dist_poisson(generator);
      // Generate the value.
      RealType value = 0.;
      for (int i=0; i<n_normal_variables; ++i)
        value += internal_sqr(dist_normal(generator));
      // Return the value.
      return value;
    }

    //! \brief Resets the distribution, so that subsequent uses of the object do not depend on values already produced by it.
    //! 
    //! This is here for compatibility, but doesn't currently do anything.
    void reset() {};

  private:

    //! \brief Private squaring function.
    inline RealType internal_sqr(RealType x) { return x*x; }

    //! \brief The noncentrality parameter.
    RealType noncentrality_parameter;

    //! \brief A poisson distribution.
    std::poisson_distribution<int> dist_poisson;
    //! \brief A normal distribution.
    std::normal_distribution<RealType> dist_normal;
  };


}

#endif // __NONCENTRAL_CHI_SQUARED_ZERO_HPP__