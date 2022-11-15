// FrequencyFilter.h
// -- interface for filters with adjustable frequency parameter
// cs246 9/19


#ifndef CS246_FREQUENCYFILTER_H
#define CS246_FREQUENCYFILTER_H


#include "Filter.h"


class FrequencyFilter : public Filter {
  public:
    virtual ~FrequencyFilter(void) { }
    virtual void setFrequency(float f) = 0;
};


#endif

