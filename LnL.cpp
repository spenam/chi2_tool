#include <cmath>
//#include <iostream>
#include <pybind11/pybind11.h>
//#include <boost/python.hpp>
//#include <TMath.h>
//
namespace py = pybind11;

double LnL0(double data, double mc){ // BEWARE this returns -2 * ( log P(n_data|l_model) - log P(n_data|n_data) )

  if(data<0) return 1e12;

  double MCError = 1e-12;
  if(mc > MCError){

    if(data > 0) return 2*(mc - data + data*log(data/mc));
    else         return 2*mc;

  }

  else{

    if(data > MCError) return 2  * (mc - data + data*log(data/MCError)) +
                              data * (3 - 4*mc/MCError + std::pow(mc/MCError,2));

    else if(data > 0)  return mc*mc/MCError + MCError -
                              2*data*(mc/MCError - log(data/MCError));

    else               return mc*mc/MCError + MCError;

  }

}


double LnL(double data, double mc, double s){

  if(data<0) return 1e12;

  if(s<=0) return LnL0(data, mc);

  double MCError = 1e-12;

  double ss = s*s;

  if(mc > MCError){

    if(data > 0){

      double dmc = 1 - ss*mc;

      double sb = 0.5 * (dmc + std::sqrt(4*data*ss + dmc*dmc)) - 1;

      double mcs = mc * (1 + sb);

      return 2*(mcs - data + data*log(data/mcs)) + sb*sb / ss;

    }
    else{
      if(mc*ss < 1) return mc * (2 - mc*ss);
      else          return 1 / ss;
    }

  }
  else{

    if(data > 0){

      if(data != MCError){

        double dmc = 1 - ss*MCError;

        double sb = 0.5 * (dmc + std::sqrt(4*data*ss + dmc*dmc)) - 1;

        double mcs = MCError * (1 + sb);

        double fmc = (mc*mc + data*(MCError - 2*mc)) / MCError;

        return fmc * (data - mcs)/(data - MCError) + 2*data*log(data/mcs) - sb*sb / ss;

      }
      else return std::pow(mc - MCError,2) / (MCError*(1 + MCError*ss));

    }
    else{

      if(MCError*ss < 1) return MCError + mc*mc*(1/MCError - ss);
      else               return 1 / ss;

    }
  }
}

PYBIND11_MODULE(lnl,handle){
    handle.def("LnL", &LnL);
}

//BOOST_PYTHON_FUNCTION_OVERLOADS(lnl_overloads, LnL, 3, 3)

//BOOST_PYTHON_MODULE(lnl){
//      //boost::python::def("LnL", LnL, lnl_overloads());
//      boost::python::def("LnL", LnL);
//}

/*
int main()
{
    double data;
    double mc;
    double s;

    std::cin >> data;
    std::cin >> mc;
    std::cin >> s;

    std::cout<<("%.16f \n", LnL(data, mc, s)) << std::endl;

    return 0;

}
*/
