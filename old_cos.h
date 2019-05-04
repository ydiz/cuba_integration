

// combining upper and lower integral into one is faster than calculating two separately
std::vector<double> Func::spherical_integral_term_1(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

  ans[0] = 4 * PI * (Mpi*0.5 - eps) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]

	ans[0] += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

  return ans;
}


std::vector<double> Func::spherical_integral_term_2(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

  ans[0] = 4 * PI * p_interval * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));

  return ans;
}

std::vector<double> Func::spherical_integral_term_3(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  ans[0] = 4 * PI * p_interval * exp(- Epe * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

  return ans;
}





std::vector<double> Func::integral_total(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

	double term1, term2, term3;

	double tmp = sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]; // \hat{p} \dot \vec{w} // performace can improve ~25%
	// term 2
  // term2 = 4 * PI * p_interval * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));
  term2 = 4 * PI * p_interval * exp(- p * v[2]) * cos( p * tmp ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  // term3 = 4 * PI * p_interval * exp(- Epe * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));
  term3 = 4 * PI * p_interval * exp(- Epe * v[2]) * cos( p * tmp ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

	// term 1
	p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]

  // term1 = 4 * PI * (Mpi*0.5 - eps) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
  term1 = 4 * PI * (Mpi*0.5 - eps) * exp(- p * v[2]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]

	// term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
	term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * exp(- p * v[2]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	// final result
	// ans[0] = term1 + term2 + term3;
	ans[0] = (1. / Mpi) * exp(0.5*Mpi*v[2]) * term1 + (1. / Mpi) * exp(-0.5*Mpi*v[2]) * term2 + term3;
// 
  return ans;
}

// can be simplified
std::vector<double> Func::integral_total_with_pe0(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

	double term1, term2, term3;

	double tmp = sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]; // \hat{p} \dot \vec{w} // performace can improve ~25%
	// term 2
  term2 = 4 * PI * p_interval * exp(- p * v[2]) * cos( p * tmp ) / (Mpi * 0.5 + pe * c);

	// term 3
	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  term3 = 4 * PI * p_interval * exp(- Epe * v[2]) * cos( p * tmp ) / ((Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

	// term 1
	term1 = 4 * PI * p_interval * exp(- p * v[2]) * cos( p * tmp ) / (-Mpi * 0.5 + pe * c);

	// final result
	ans[0] = - (1. / Mpi) * exp(0.5*Mpi*v[2]) * term1 + (1. / Mpi) * exp(-0.5*Mpi*v[2]) * term2 + term3; // note the sign of term 1

  return ans;
}

std::vector<double> Func::integral_total_with_p3(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

	double term1, term2, term3;

	double tmp = sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]; // \hat{p} \dot \vec{w} // performace can improve ~25%
	// term 2
  term2 = 4 * PI * p_interval * p * c * exp(- p * v[2]) * cos( p * tmp ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  term3 = 4 * PI * p_interval * p * c * exp(- Epe * v[2]) * cos( p * tmp ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

	// term 1
	p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]

  term1 = 4 * PI * (Mpi*0.5 - eps) * p * c * exp(- p * v[2]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]

	term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * p * c * exp(- p * v[2]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	// final result
	// ans[0] = term1 + term2 + term3;
	ans[0] = (1. / Mpi) * exp(0.5*Mpi*v[2]) * term1 + (1. / Mpi) * exp(-0.5*Mpi*v[2]) * term2 + term3;
// 
  return ans;
}

std::vector<double> Func::integral_total_with_p1(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double sin_theta = sqrt(1 - c*c); // sin_theta
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

	double term1, term2, term3;
	double sin_phi, cos_phi;
	sincos(phi, &sin_phi, &cos_phi);

	double tmp = sin_theta*cos_phi*v[0] + sin_theta*sin_phi*v[1] + c*v[2]; // \hat{p} \dot \vec{w} // performace can improve ~25%
	// term 2
  term2 = 4 * PI * p_interval * p * sin_theta * cos_phi * exp(- p * v[3]) * cos( p * tmp ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  term3 = 4 * PI * p_interval * p * sin_theta * cos_phi * exp(- Epe * v[3]) * cos( p * tmp ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

	// term 1
	p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]

  term1 = 4 * PI * (Mpi*0.5 - eps) * p * sin_theta * cos_phi * exp(- p * v[3]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]

	term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * p * sin_theta * cos_phi * exp(- p * v[3]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	// final result
	// ans[0] = term1 + term2 + term3;
	ans[0] = (1. / Mpi) * exp(0.5*Mpi*v[3]) * term1 + (1. / Mpi) * exp(-0.5*Mpi*v[3]) * term2 + term3;
// 
  return ans;
}



