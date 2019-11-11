#ifndef __EXPERIMENTS_HPP__
#define __EXPERIMENTS_HPP__

#include "predictive.hpp"

typedef pair<float,float> point;

//! \brief This class encapsulates a number of experiments that can be done using the Scheduler.
class Experiments {
public:

  void generate_events(vector<point>&, vector<point>&, float, float, float);

  // --- Experiments

  void generate_scores(vector<float>&, int, float, float);

  void create_system_grid(vector<vector<float> >&, float, float, float, float, int);

  void score_vary_total_time(vector<point>&, float, float, float, float, int, int);

  void score_vary_bin_granularity(vector<pair<int,float> >&, int, int, int, int, float, float);

  // --- Mutators

  void setMacroBinSize(int b) {
    schedule.set_macro_bin_size(b);
  }

  void setMacroBinTime(float t) {
    schedule.set_macro_bin_time(t);
  }

  void setTotalTime(float t) {
    schedule.set_time_sight(t);
    if (t>0) total_time = t;
  }

private:
  // Create scheduler.
  Schedule schedule;
  
  float total_time = 0;

};

#endif // __EXPERIMENTS_HPP__