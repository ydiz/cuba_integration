// std::vector<double> Func::spherical_integral_w(const std::vector<double>& vx) const
// {
//   std::vector<double> ans(1);
// 	double p = vx[0] * p_interval; // p: [0, p_interval]
// 	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
// 	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]
//
// 	double d = p_interval; // box size is [-d/2, d/2]
//
// 	double Mpi = 0.135;
// 	double me = 0.000511;
// 	double pe = sqrt(0.5 * Mpi * 0.5 * Mpi - me*me);
//
//   ans[0] = 4 * PI * p_interval * cos( p * (sqrt(1 - c*c)*( cos(phi)*w[0] + sin(phi)*w[1] )  + c*w[2]) ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));
//
//   return ans;
// }

inline std::vector<double> test_integrand4d(const std::vector<double>& vx)
{
  assert(4 == vx.size());
  std::vector<double> ans(1);
  ans[0] = vx[0] * vx[1] * vx[2] * vx[3] + sin(vx[0] * PI) * sin(vx[1] * PI) + sqrt(vx[3]) * sqrt(vx[2]) * sqrt(vx[1]);
  return ans;
}

inline std::vector<double> my_integral(const std::vector<double>& vx, double box_size)
{
  // TIMER_VERBOSE("test_integrand4d");
  // using namespace qlat;
  // assert(3 == vx.size());
  std::vector<double> ans(1);
	std::vector<double> p(vx.size());
	for(size_t i=0; i<vx.size(); ++i) p[i] = vx[i] - 0.5;

	double d = box_size; // box size is [-d/2, d/2]

	double Mpi = 0.135;
	double me = 0;
	double pMinus = 0.5 * Mpi;
	double Ep = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]); 
	double twoEpPlusMpi = 2. * Ep + Mpi / d;
	// double Epe2 = p[0]*p[0] + p[1]*p[1] + (p[2] - pMinus/d) * (p[2] - pMinus/d) + me*me/d/d;

  // ans[0] = 1.0 / d * cos( (p[0] + p[1] + p[2]) * d ) / ((Ep + 0.001) * twoEpPlusMpi * (twoEpPlusMpi * twoEpPlusMpi - 4. * Epe2 + 0.001));
  ans[0] = 1.0 / d * cos( (p[0] + p[1] + p[2]) * d ) / ((Ep + 0.001) * twoEpPlusMpi * (Ep * Mpi + 2. * p[2] * pMinus  + 0.001));
  // ans[0] = cos( (p[0] + p[1] + p[2]) * d ) / (Ep*Ep*Ep + 0.01);
	// std::cout << "p[0]: "  << p[0] << " p[1]: "<< p[1] << " p[2]: " << p[2] << std::endl;
	// std::cout << " Ep: " << Ep << " Epe2: "<< Epe2 << "twoEpPlusMpi: "  << twoEpPlusMpi  << " Ans[0]: " << ans[0] << std::endl;
  // ans[0] = d * cos( (2 * p[0] + 2 * p[1] + 2 * p[2]) * d ) / (p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + 1. / d / d);
  return ans;
}

// std::vector<double> Func::spherical_integral_term_1(const std::vector<double>& vx) const
// {
//   std::vector<double> ans(1);
// 	double p = vx[0]*(p_interval - (Mpi*0.5 + 0.01)) + (Mpi*0.5 + 0.01); // p: [Mpi/2 + eps, p_interval]
// 	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
// 	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]
//
//   ans[0] = 4 * PI * (p_interval - (Mpi*0.5 + 0.01)) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
//
//   return ans;
// }
