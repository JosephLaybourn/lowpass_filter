// LowPassRes.cpp
// -- high order LP filter with resonance
// Joey Laybourn
// cs246 10/19

#include "LowPassRes.h"
#include <cmath>
#include <complex>

#define PI 3.14159265359

///////////////////////////////////////////////////
// SECOND ORDER IMPLEMENTATION
///////////////////////////////////////////////////
LowPassRes2::LowPassRes2(float f, float r, float R) :
  irate(R),
  resonance(r),
  cutoff(f),
  x1(0),
  x2(0),
  y1(0),
  y2(0)
{
  reset();
}

void LowPassRes2::setFrequency(float f)
{
  cutoff = f;

  reset();
}

void LowPassRes2::setResonance(float r)
{
  resonance = r;

  reset();
}

float LowPassRes2::operator()(float x)
{
  float y = (a0 * x) + (a1 * x1) + (a0 * x2) + (b1 * y1) + (b2 * y2);
  y2 = y1;
  y1 = y;
  x2 = x1;
  x1 = x;

  return y;
}

float LowPassRes2::sGain(float f)
{
  std::complex<double> phasor = std::polar(1.0, (2.0 * PI * f) / irate);
  double tau = tan((PI * cutoff) / irate);
  double xi = (PI / 4.0) * (3.0 - resonance);
  std::complex<double> numerator = (tau * tau) * ((phasor * phasor) + (2.0 * phasor) + 1.0);
  std::complex<double> denom1 = (1 - (2 * tau * cos(xi)) + (tau * tau)) * (phasor * phasor);
  std::complex<double> denom2 = (-2 * (1 - (tau * tau))) * phasor;
  std::complex<double> denom3 = 1 + (2 * tau * cos(xi)) + (tau + tau);
  std::complex<double> denominator = denom1 + denom2 + denom3;

  return std::norm(numerator) / std::norm(denominator);
}

void LowPassRes2::reset(void)
{
  double tau = tan((PI * cutoff) / irate);
  double xi = (PI / 4.0) * (3.0 - resonance);
  double gamma = 1.0 - (2.0 * tau * cos(xi)) + (tau * tau);
  a0 = (tau * tau) / gamma;
  a1 = (2.0 * tau * tau) / gamma;
  b1 = (2.0 * (1.0 - (tau * tau))) / gamma;
  b2 = -1.0 * (1.0 + (2.0 * tau * cos(xi) + (tau * tau))) / gamma;
}

///////////////////////////////////////////////////
// THIRD ORDER IMPLEMENTATION
///////////////////////////////////////////////////
LowPassRes3::LowPassRes3(float f, float r, float R) :
  irate(R),
  resonance(r),
  cutoff(f),
  x1(0),
  x2(0),
  x3(0),
  y1(0),
  y2(0),
  y3(0)
{
  reset();
}

void LowPassRes3::setFrequency(float f)
{
  cutoff = f;
  reset();
}

void LowPassRes3::setResonance(float r)
{
  resonance = r;
  reset();
}

float LowPassRes3::operator()(float x)
{
  float y = (a0 * (x + (3 * x1) + (3 * x2) + (3 * x3))) + (b1 * y1) + (b2 * y2) + (b3 * y3);
  x3 = x2;
  x2 = x1;
  x1 = x;
  y3 = y2;
  y2 = y1;
  y1 = y;

  return y;
}

float LowPassRes3::sGain(float f)
{
  std::complex<double> phasor = std::polar(1.0, (2.0 * PI * f) / irate);
  
  double zeta = (PI / 6.0) * (4.0 - resonance);
  double tau = tan((PI * cutoff) / irate);

  std::complex<double> num1 = a0 * phasor * phasor * phasor;
  std::complex<double> num2 = a0 * 3.0 * phasor * phasor;
  std::complex<double> num3 = a0 * 3.0 * phasor;
  std::complex<double> num4(a0, 0.0);

  std::complex<double> denom1 = phasor * phasor * phasor;
  std::complex<double> denom2 = -b1 * phasor * phasor;
  std::complex<double> denom3 = -b2 * phasor;
  std::complex<double> denom4(-b3, 0.0);

  return std::norm(num1 + num2 + num3 + num4) / std::norm(denom1 + denom2 + denom3 + denom4);
}

void LowPassRes3::reset(void)
{
  double zeta = (PI / 6.0) * (4.0 - resonance);
  double tau = tan((PI * cutoff) / irate);
  double delta = 1.0 - (2.0 * cos(zeta));
  double gamma = 1.0 + (delta * tau * (1.0 + tau)) + (tau * tau * tau);
  
  a0 = (tau * tau * tau) / gamma;
  b1 = (3.0 + (delta * tau * (1.0 - tau)) - (3.0 * (tau * tau * tau))) / gamma;
  b2 = -(3.0 - (delta * tau * (1.0 + tau)) + (3.0 * (tau * tau * tau))) / gamma;
  b3 = (1.0 - (delta * tau * (1.0 - tau)) - (tau * tau * tau)) / gamma;
}