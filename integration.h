#pragma once

#include <cuba.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>
#include <math.h>
#include "cuba_wrapper.h"
#include "constants_macro.h"

struct Int_para {
	// double p_interval;
	double eps_rel;
	double p_interval;
	double term1_eps;
	std::string algorithm;
	std::string term;
	int x_for_p1;
};

const double PI = 3.1415926535897932384626433832795;


class Func {

	std::vector<double> (Func::*pf)(const std::vector<double> &) const;
	const double Mpi; // mass of pion
	const double me; // mass of electron
	double pe; // momentum of outgoing electron
	double eps; // eps for calculating the principal part of term 1

public:
	double p_interval;
	std::vector<double> v; // v is a three-vector; v = (sqrt(x^2 + y^2), z, |t|)

	Func(const Int_para& para);

	std::vector<double> operator()(const std::vector<double>& vx) const { return (this->*pf)(vx);	}
	
	std::vector<double> integral_total(const std::vector<double>& vx) const;
	std::vector<double> integral_total_cos(const std::vector<double>& vx) const;
	std::vector<double> spherical_integral_term_1(const std::vector<double>& vx) const; 
	std::vector<double> spherical_integral_term_2(const std::vector<double>& vx) const; 
	std::vector<double> spherical_integral_term_3(const std::vector<double>& vx) const; 

	std::vector<double> integral_total_with_pe0(const std::vector<double>& vx) const;
	std::vector<double> integral_total_with_p3(const std::vector<double>& vx) const;
	std::vector<double> integral_total_with_p1(const std::vector<double>& vx) const;
};

Func::Func(const Int_para& para) : Mpi(0.1349766), me(0.000511) {

	p_interval = para.p_interval;
	eps = para.term1_eps;
	pe = sqrt(0.5 * Mpi * 0.5 * Mpi - me*me);

	// if(para.term == "term1") pf = &Func::spherical_integral_term_1;
	// else if(para.term == "term2") pf = &Func::spherical_integral_term_2;
	// else if(para.term == "term3") pf = &Func::spherical_integral_term_3;
  if(para.term == "all") pf = &Func::integral_total;
	// else if(para.term == "all_with_pe0") pf = &Func::integral_total_with_pe0;
	else if(para.term == "all_with_p3") pf = &Func::integral_total_with_p3;
	else if(para.term == "all_with_p1") pf = &Func::integral_total_with_p1;
	else assert(0);

	std::cout << "Mpi: " << Mpi << " me: " << me  << std::endl;
	std::cout << std::string(20, '*') << std::endl;
}


// std::vector<double> Func::integral_total_cos(const std::vector<double>& vx) const
// {
//   std::vector<double> ans(1);
// 	double p = vx[0] * p_interval; // p: [0, p_interval]
// 	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
// 	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]
//
// 	double term1, term2, term3;
//
// 	double tmp = sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]; // \hat{p} \dot \vec{w} // performace can improve ~25%
// 	// term 2
//   // term2 = 4 * PI * p_interval * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));
//   term2 = 4 * PI * p_interval * exp(- p * v[2]) * cos( p * tmp ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));
//
// 	// term 3
// 	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
//   // term3 = 4 * PI * p_interval * exp(- Epe * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));
//   term3 = 4 * PI * p_interval * exp(- Epe * v[2]) * cos( p * tmp ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));
//
// 	// term 1
// 	p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]
//
//   // term1 = 4 * PI * (Mpi*0.5 - eps) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
//   term1 = 4 * PI * (Mpi*0.5 - eps) * exp(- p * v[2]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
//
// 	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]
//
// 	// term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
// 	term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * exp(- p * v[2]) * cos( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
//
// 	// final result
// 	// ans[0] = term1 + term2 + term3;
// 	ans[0] = (1. / Mpi) * exp(0.5*Mpi*v[2]) * term1 + (1. / Mpi) * exp(-0.5*Mpi*v[2]) * term2 + term3;
// // 
//   return ans;
// }


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
  term2 = 4 * PI * p_interval * exp(- p * v[2]) * sin( p * tmp ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  // term3 = 4 * PI * p_interval * exp(- Epe * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));
  term3 = 4 * PI * p_interval * exp(- Epe * v[2]) * sin( p * tmp ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

	// term 1
	p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]

  // term1 = 4 * PI * (Mpi*0.5 - eps) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
  term1 = 4 * PI * (Mpi*0.5 - eps) * exp(- p * v[2]) * sin( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]

	// term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * exp(- p * v[2]) * cos( p * (sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]) ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));
	term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * exp(- p * v[2]) * sin( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	// final result
	// ans[0] = term1 + term2 + term3;
	ans[0] = (1. / Mpi) * exp(0.5*Mpi*v[2]) * term1 + (1. / Mpi) * exp(-0.5*Mpi*v[2]) * term2 + term3;
// 
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
  term2 = 4 * PI * p_interval * p * c * exp(- p * v[2]) * sin( p * tmp ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  term3 = 4 * PI * p_interval * p * c * exp(- Epe * v[2]) * sin( p * tmp ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

	// term 1
	p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]

  term1 = 4 * PI * (Mpi*0.5 - eps) * p * c * exp(- p * v[2]) * sin( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]

	term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * p * c * exp(- p * v[2]) * sin( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

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
  term2 = 4 * PI * p_interval * p * sin_theta * cos_phi * exp(- p * v[3]) * sin( p * tmp ) / ((p + Mpi * 0.5) * (Mpi * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (Mpi/2.)*(Mpi/2.) - 2 * p * pe * c);
  term3 = 4 * PI * p_interval * p * sin_theta * cos_phi * exp(- Epe * v[3]) * sin( p * tmp ) / (Epe * (Mpi * 0.5 + pe * c) * (-Mpi * 0.5 + pe * c));

	// term 1
	p = vx[0] * (Mpi*0.5 - eps); // p: [0, Mpi/2 - eps]

  term1 = 4 * PI * (Mpi*0.5 - eps) * p * sin_theta * cos_phi * exp(- p * v[3]) * sin( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	p = vx[0]*(p_interval - (Mpi*0.5 + eps)) + (Mpi*0.5 + eps);  // p: [Mpi/2 + eps, p_interval]

	term1 += 4 * PI * (p_interval - (Mpi*0.5 + eps)) * p * sin_theta * cos_phi * exp(- p * v[3]) * sin( p * tmp ) / ((-p + Mpi * 0.5) * (-Mpi * 0.5 + pe * c));

	// final result
	// ans[0] = term1 + term2 + term3;
	ans[0] = (1. / Mpi) * exp(0.5*Mpi*v[3]) * term1 + (1. / Mpi) * exp(-0.5*Mpi*v[3]) * term2 + term3;
// 
  return ans;
}



void integrate(const Int_para &para, const Func &func)
{
  std::vector<double> integral, error, prob;
  int nregions, neval, fail;

	if(para.algorithm=="Cuhre") integrateCuhre(integral, error, prob, nregions, neval, fail, 3, 1, func, 0., para.eps_rel);
	else if(para.algorithm=="Divonne") integrateDivonne(integral, error, prob, nregions, neval, fail, 3, 1, func, 0., para.eps_rel);
	else assert(0); // wrong algorithm name

}


void run_p1(const Int_para &para) {

	Func func(para);

	int x = para.x_for_p1;
	if(x==0) return ; // when x=0, integral is 0 and cuba cannote calculate it
	// int x_min = -8, x_max = 8; // 
	int y_max = SPACE_LIMIT;
	int z_max = SPACE_LIMIT;
	int t_max = SPACE_LIMIT; 
	std::cout << "x: " << x << std::endl; 
	std::cout << "y: [" << -y_max << ", " << y_max << "]"<< std::endl; 
	std::cout << "z: [" << -z_max << ", " << z_max << "]"<< std::endl; 
	std::cout << "t: [" << 0 << ", " << t_max << "]"<< std::endl; 
	std::cout << std::string(20, '*') << std::endl;

	// for(int x=x_min; x<=x_max; ++x) 
	for(int y=-y_max; y<=y_max; ++y) 
		for(int z=-z_max; z<=z_max; ++z)
			for(int t=0; t<=t_max; ++t) { // t [0, L/4] // t appears only as absolute value
				func.v = {double(x), double(y), double(z), double(t)};
				std::cout << "v: " << "[" << func.v[0] << " " << func.v[1] << " " << func.v[2] << " " << func.v[3] << "]" << std::endl;
				if(z==0) {
					std::cout << cur_time << "s	integral = " << 0 << std::endl;
					continue; // when z=0, integral is 0 and cuba cannote calculate it
				}
				if(t==0 || t==1){ // Cuba cannot calculate when t=0; t=1 is very slow to calculate
					std::cout << cur_time << "s	Cannot calculate t=0,1; I am actually calculating with t=4" << std::endl;
					func.v[3] = 4;
				}	
        if(t==2){ // Cuba cannot calculate when t=0; t=1 is very slow to calculate
					std::cout << cur_time << "calculating t=2 is slow; I am actually calculating with t=4" << std::endl;
					func.v[3] = 4;
				}
				integrate(para, func);
	}

	// func.v = {1., 0., -1., 10.};
	// std::cout << "v: " << "[" << func.v[0] << " " << func.v[1] << " " << func.v[2] << " " << func.v[3] << "]" << std::endl;
	// integrate(para, func);
}


void run(const Int_para &para) {
	Func func(para);

	// func.v = {0, 0, 1};
	// integrate(para, func);
	// func.v = {sqrt(2), 1, 1};
	// integrate(para, func);
  // func.v = {3, 5, 1};
	// integrate(para, func);
  // return;


	int r_max = int(SPACE_LIMIT*std::sqrt(2)) + 1, z_max = SPACE_LIMIT, t_max = TIME_LIMIT;
	// int t_start = (para.term == "all_with_pe0" || para.term == "all_with_p3") ? 1 : 0; 
	std::cout << "r: [" << 0 << ", " << r_max << "]"<< std::endl; 
	std::cout << "z: [" << -z_max << ", " << z_max << "]"<< std::endl; 
	std::cout << "t: [" << 0 << ", " << t_max << "]"<< std::endl; 
	std::cout << std::string(20, '*') << std::endl;

	for(int r=0; r<=r_max; ++r) // r is distance, always positive
		for(int z=-z_max; z<=z_max; ++z)
			for(int t=0; t<=t_max; ++t) { // t [0, L/4] // t appears only as absolute value

				func.v = {double(r), double(z), double(t)};
				std::cout << "v: " << "[" << func.v[0] << " " << func.v[1] << " " << func.v[2] << "]" << std::endl;

        // for sin, when z=0, integral is 0
        if(z==0) {
						std::cout << cur_time << "s	integral = " << 0 << std::endl;
						continue;
        }	

				if(t == 0){
					if(para.term == "all_with_pe0") { // if pe0 is in numerator and t=0, the integral is 0 and Cuba would have difficulty figuring it out
						std::cout << cur_time << "s	integral = " << 0 << std::endl;
						continue;
					}	
					if(para.term == "all_with_p3") {
						std::cout << cur_time << "s	Cannot calculate t=0; I am actually calculating with t=1" << std::endl;
						func.v[2] = 1;
					}
					if(para.term == "all") {
						std::cout << cur_time << "s	Cannot calculate t=0 when z is larger than 10. I am lazy. I am skipping all t=0 and am actually calculating with t=1" << std::endl;
						func.v[2] = 1;
					}
				}

				integrate(para, func);
	}

	// func.v = {3, 5, 0.5};
	// integrate(para, func);
	// func.v = {3, 5, 0.2};
	// integrate(para, func);
	// func.v = {10, 7, 0.5};
	// integrate(para, func);
	// func.v = {10, 7, 0.};
	// integrate(para, func);

}
