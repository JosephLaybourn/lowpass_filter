// LowPassRes.h
// -- high order LP filter with resonance
// cs246 10/19

#ifndef CS246_LOWPASSRES_H
#define CS246_LOWPASSRES_H


#include "FrequencyFilter.h"


class LowPassRes2 : public FrequencyFilter {
  public:
    LowPassRes2(float f=0, float r=0, float R=44100);
    void setFrequency(float f);
    void setResonance(float r);
    float operator()(float x);
    float sGain(float f);
  private:
    float irate, resonance, cutoff;
    double a0, a1,
           b1, b2,
           x1, x2,
           y1, y2;
    void reset(void);
};


class LowPassRes3 : public FrequencyFilter {
  public:
    LowPassRes3(float f=0, float r=0, float R=44100);
    void setFrequency(float f);
    void setResonance(float r);
    float operator()(float x);
    float sGain(float f);
  private:
    float irate, resonance, cutoff;
    double a0, a1, a2,
           b1, b2, b3,
           x1, x2, x3,
           y1, y2, y3;
    void reset(void);
};


#endif

