// g++ -I/home/ydzhao/cuth/install/boost/include -L/home/ydzhao/cuth/install/boost/lib GF_para.cc -lboost_program_options
#include <stdlib.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

void init(int argc, char **argv, Int_para &para)
{
  po::options_description desc("GFFA options");
  desc.add_options()("help", "help message")
                    // ("box,b", po::value<double>(&para.box_length)->default_value(1.0), "box length. For Cartesian coordinate")
                    ("p_interval,p", po::value<double>(&para.p_interval)->default_value(1), "integration interval of p's norm. For spherical coordinate")
                    ("eps_rel,e", po::value<double>(&para.eps_rel)->default_value(1e-3), "required relative error")
                    ("term1_eps", po::value<double>(&para.term1_eps)->default_value(1e-2), "epsilon for calculating principal part in term 1")
                    ("algorithm,a", po::value<std::string>(&para.algorithm)->default_value(""), "Cuhre or Divonne")
                    ("term,t", po::value<std::string>(&para.term)->default_value(""), "term1, term2 or term3")
                    ("x_for_p1,x", po::value<int>(&para.x_for_p1)->default_value(0), "x for with_p1")
                    ("M_hadron", po::value<double>(&para.M_hadron)->default_value(0), "Intial state hadron mass")
                    ("M_lepton", po::value<double>(&para.M_lepton)->default_value(0), "Final state lepton mass")
                    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options // command line options have higher priority
  po::store(po::parse_config_file<char>("int.ini", desc), vm);
  po::notify(vm);

  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

	std::cout << std::string(20, '*') << std::endl;
	std::cout << "Initial state hadron mass: " << para.M_hadron << std::endl;
	std::cout << "Final state lepton mass: " << para.M_lepton << std::endl;
	std::cout << std::string(20, '*') << std::endl;
	std::cout << "p_interval: " << para.p_interval << std::endl;
	std::cout << "eps_rel: " << para.eps_rel << std::endl;
	std::cout << "term1_eps: " << para.term1_eps << std::endl;
	std::cout << "algorithm: " << para.algorithm << std::endl;
	std::cout << "term: " << para.term << std::endl;
	if(para.term=="with_p1") std::cout << "x_for_p1: " << para.x_for_p1 << std::endl;
	std::cout << std::string(20, '*') << std::endl;

}
