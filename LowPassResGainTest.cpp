// LowPassResGainTest.cpp
// cs246 10/19

#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include "LowPassRes.h"
using namespace std;


const float PI = 4.0f*atan(1.0f);


float filterGain(Filter& F, float f, float R) {
  int count = int(R),
      min_count = int(0.1f*R);
  double darg = 2*PI*f/R,
         arg = 0;
  float maxg = 0;
  for (int i=0; i < count; ++i) {
    float g = abs(F(float(sin(arg))));
    arg += darg;
    if (i > min_count && g > maxg)
      maxg = g;
  }
  return maxg;
}



int main(void) {
  const float MAX = float((1<<15)-1),
              RATE = 44100.0f,
              MINDB = -48.0f,
              MAXDB = 48.0f;
  const int COUNT = 1000;

  LowPassRes2 lp2(500,0.95f,RATE);
  LowPassRes3 lp3(500,0.95f,RATE);

  for (int k=0; k < 2; ++k) {
    short *samples = new short[2*COUNT];

    // compute log(f) vs log(gain) plot
    for (int i=0; i < COUNT; ++i) {
      float frequency = 20.0f  * pow(10,3*float(i)/float(COUNT-1)),
            gain = (k == 0) ? filterGain(lp2,frequency,RATE)
                            : filterGain(lp3,frequency,RATE),
      dB = min(MAXDB,max(MINDB,dB));
      samples[2*i+0] = short(MAX + 2*MAX*(dB-MAXDB)/(MAXDB-MINDB));
      gain = (k == 0) ? lp2.sGain(frequency) : lp3.sGain(frequency);
      dB = min(MAXDB,max(MINDB,20.0f*log10(gain)));
      samples[2*i+1] = short(MAX + 2*MAX*(dB-MAXDB)/(MAXDB-MINDB));
    }

    // write WAVE file
    struct {
      char riff_chunk[4];
      unsigned chunk_size;
      char wave_fmt[4];
      char fmt_chunk[4];
      unsigned fmt_chunk_size;
      unsigned short audio_format;
      unsigned short number_of_channels;
      unsigned sampling_rate;
      unsigned bytes_per_second;
      unsigned short block_align;
      unsigned short bits_per_sample;
      char data_chunk[4];
      unsigned data_chunk_size;
    }
    header = { {'R','I','F','F'},
               36 + 4*COUNT,
               {'W','A','V','E'},
               {'f','m','t',' '},
               16,1,2,unsigned(RATE),4*unsigned(RATE),2,16,
               {'d','a','t','a'},
               4*COUNT
             };
    string fname = "LowPassResGainTest." + to_string(k+2) + ".wav";
    fstream out(fname.c_str(),ios_base::binary|ios_base::out);
    out.write(reinterpret_cast<char*>(&header),44);
    out.write(reinterpret_cast<char*>(samples),4*COUNT);

    delete[] samples;
  }

  return 0;
}

