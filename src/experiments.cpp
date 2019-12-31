#include "experiments.hpp"

Experiments::Experiments() {
  seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
  generator = std::mt19937(seed);
}

point Experiments::score() {
  // Find max and baseline scores.
  float ms = schedule.max_score();
  float bs = schedule.base_line_score();
  // Return max and baseline scores.
  return make_pair(ms, bs);
}

void Experiments::generate_scores(vector<float>& score_record, int trials, float rhoL, float rhoR, bool only_max) {
  vector<point> left, right;
  score_record.clear();
  // Calculate data.
  setDensities(rhoL, rhoR);
  for (int i=0; i<trials; ++i) {
    // Create new events.
    generate_events(left, right, total_time);
    schedule.set_events(left, right);
    // Find the maximum and baseline score.
    float ms = schedule.max_score(), bs = schedule.base_line_score();
    if (only_max) score_record.push_back(ms);
    else score_record.push_back(ms - bs);
  }
}

void Experiments::create_system_grid(vector<vector<float> >& multidata, float minRL, float maxRL, float minRR, float maxRR, int divisions) {
  // System steps and vectors for events.
  float dl = (maxRL-minRL)/(divisions-1);
  float dr = (maxRR-minRR)/(divisions-1);
  vector<point> left, right;
  multidata.clear();
  float macro_bin_time = schedule.get_macro_bin_time();
  // Calculate data.
  for (int nl=1; nl<divisions; ++nl)
    for (int nr=1; nr<divisions; ++nr) {
      vector<float> entry;
      // Choose some densities.
      float rhoL = nl*dl + minRL;
      float rhoR = nr*dr + minRR;
      // Generate events
      setDensities(rhoL, rhoR);
      generate_events(left, right, total_time);
      schedule.set_events(left, right);
      // Find the maximum and baseline score.
      entry.push_back(rhoL*macro_bin_time);
      entry.push_back(rhoR*macro_bin_time);
      entry.push_back(schedule.max_score());
      entry.push_back(schedule.base_line_score());
      multidata.push_back(entry);
    }
}

void Experiments::score_vary_tau(vector<point>& data, float rhoL, float rhoR, float minT, float maxT, int divisions, int trials) {
  // Clear old data vector.
  data.clear();
  // Make sure no numerical errors occur.
  divisions = std::max(2, divisions);
  trials = std::max(1, trials);
  // Calculate time step.
  float dt = (maxT - minT)/(divisions-1);
  // Generate data.
  for (int iter=0; iter<divisions; ++iter) {
    // Calculate and set tau.
    float tau = minT + iter*dt;
    setTau(tau);
    // Calculate average score.
    float ave_score = 0.f;
    for (int tr=0; tr<trials; ++tr) {
      generate_events(rhoL, rhoR);
      ave_score += schedule.max_score();
    }
    data.push_back(point(tau, ave_score/trials));
  }
}

void Experiments::score_vary_total_time(vector<vector<float> >& growth, float rhoL, float rhoR, float minT, float maxT, int divisions, int trials) {
  // Time step and vectors for events.
  float dt = (maxT-minT)/(divisions-1);
  vector<point> left, right;
  growth.clear();
  setDensities(rhoL, rhoR);
  // Calculate data.
  for (int j=0; j<divisions; ++j) {
    float time = j*dt + minT, normalization = 1./(trials*time);;
    float ave_max_score = 0.f, ave_base_line_score = 0.f;
    // Set total time and time sight (they are the same).
    setTotalTime(time);
    setTimeSight(time);
    // Compute average score.
    for (int i=0; i<trials; ++i) {
      generate_events(rhoL, rhoR);
      // Find the maximum and baseline score.
      ave_max_score += schedule.max_score();
      ave_base_line_score += schedule.base_line_score();
    }
    growth.push_back(vector<float>{time, ave_max_score*normalization, ave_base_line_score*normalization});
  }
}

void Experiments::score_vary_system(vector<vector<float> >& data, float rhoL, float minRhoR, float maxRhoR, int divisions, int trials) {
  data.clear();
  float drho = (maxRhoR - minRhoR)/(divisions-1);
  for (int iter=0; iter<divisions; ++iter) {
    float ave_score = 0.f, baseline_score = 0.f;
    float rhoR = minRhoR + drho*iter;
    setDensities(rhoL, rhoR);
    for (int trial=0; trial<trials; ++trial) {
      generate_events(rhoL, rhoR);
      ave_score += schedule.max_score();
      baseline_score += schedule.base_line_score();
    }
    data.push_back(vector<float>{rhoR, ave_score/(trials*total_time), baseline_score/(trials*total_time)});
  }
}

void Experiments::score_vary_bin_granularity(vector<pair<int,float> >& data, int minB, int maxB, int dB, int trials, float rhoL, float rhoR) {
  // Vectors for events.
  vector<point> left, right;
  data.clear();
  setDensities(rhoL, rhoR);
  // Calculate data.
  for (int b=minB; b<=maxB; b += dB) {
    float average_max = 0.f, average_base_line = 0.f;
    // Set macro bin size.
    schedule.set_granularity(b);
    // Do some number of trials.
    for (int i=0; i<trials; ++i) {
      // Create new events.
      generate_events(left, right, total_time);
      schedule.set_events(left, right);
      // Find the maximum and baseline score.
      float max_score = schedule.max_score();
      float base_line_score = schedule.base_line_score();
      // Calculate ratio.
      average_max += max_score;
      average_base_line += base_line_score;
    }
    // Don't divide by zero.
    if (average_base_line>0) data.push_back(make_pair(b, average_max/average_base_line));
  }
}

void Experiments::score_vary_bin_granularity(vector<pair<int,float> >& data, int minB, int maxB, int dB, int trials, 
  float rhoL, float rhoR, vector<pair<int,float> >& max_scores, vector<pair<int,float> >& base_line_scores) {
  // Vectors for events.
  vector<point> left, right;
  data.clear();
  max_scores.clear();
  base_line_scores.clear();
  float time = schedule.get_time_sight();
  setDensities(rhoL, rhoR);
  // Calculate data.
  for (int b=minB; b<=maxB; b += dB) {
    float average_max = 0.f, average_base_line = 0.f;
    // Set macro bin size.
    schedule.set_granularity(b);
    // Do some number of trials.
    for (int i=0; i<trials; ++i) {
      // Create new events.
      generate_events();
      // Find the maximum and baseline score.
      float max_score = schedule.max_score();
      float base_line_score = schedule.base_line_score();
      // Calculate ratio.
      average_max += max_score;
      average_base_line += base_line_score;
    }
    // Don't divide by zero.
    if (average_base_line>0) data.push_back(make_pair(b, average_max/average_base_line));
    max_scores.push_back(make_pair(b, average_max/time/trials));
    base_line_scores.push_back(make_pair(b, average_base_line/time/trials));
  }
}

void Experiments::demon_vary_time_sight(vector<point>& data, float time, float min_t_sight, float max_t_sight, int bins, int trials, float rhoL, float rhoR) {
  // Make sure the data vector is empty.
  data.clear();
  // Run all time sights [trials] number of times.
  float dsight = (max_t_sight - min_t_sight)/(bins-1);
  for (int b=0; b<bins; ++b) {
    float t_sight = min_t_sight + b*dsight;
    float ave_score = 0.f;
    for (int iter=0; iter<trials; ++iter) {
      generate_events(rhoL, rhoR);
      ave_score += simulate_demon(time, t_sight);
    }
    data.push_back(make_pair(t_sight, ave_score/trials));
  }
}

void Experiments::demon_vary_tau(vector<point>& data, float time, float t_sight, float min_tau, float max_tau, int bins, int trials) {
  // Make sure the data vector is empty.
  data.clear();
  // Run all time sights [trials] number of times.
  float dt = (max_tau - min_tau)/(bins-1);
  for (int iter=0; iter<bins; ++iter) {
    float tau = min_tau + iter*dt;
    setTau(tau);
    float ave_score = 0.f;
    for (int trial=0; trial<trials; ++trial) {
      generate_events();
      ave_score += simulate_demon(time, t_sight);
    }
    data.push_back(make_pair(tau, ave_score/trials));
  }
}

float Experiments::simulate_demon(float time, float t_sight) {
  // Set total time and time sight.
  setTotalTime(time);
  setTimeSight(t_sight);
  // Set up / get some data.
  float current_time = 0.f;
  float micro_bin_time = schedule.get_micro_bin_time();
  int granularity = schedule.get_granularity();
  int total_micro_bins = schedule.get_total_micro_bins();

  // Initial best scheduling.
  float max_score = schedule.max_score();
  vector<pair<int, int>> binning = schedule.optimal_binning();

  vector<pair<int, int>> finite_demon_binning;
  point current_bin;

  bool is_open = (!binning.empty() && binning[0].first==0);
  if (is_open) current_bin.first = 0;

  // How long the door has been in the same state and how many time bins have gone by.
  int state_length = granularity, cumulative_bins = 0; 
  float total_score = 0.f;
  while (current_time + t_sight < time) {
    // Find next point when we can recalculate the optimal binning, and what the score is between now and then.
    int move_bins; // Number of microbins to move forward.

    // --- Determine how far ahead in time we have to move.

    // If the door should open (if binning is empty, door should be closed the entire time).
    if (!binning.empty() && binning[0].first==0) {
      if (is_open) { // If the door was open.
        move_bins = 1;
        ++state_length;
      }
      else { // If the door was closed, door should OPEN.
        move_bins = granularity;
        state_length = granularity;
        // Door opens.
        current_bin.first = cumulative_bins;
      }
      is_open = true;
    }
    // If the door should be closed.
    else {
      if (is_open) { // If the door was open, door should CLOSE.
        move_bins = granularity;
        state_length = granularity;
        // Door closes
        current_bin.second = cumulative_bins-1;
        finite_demon_binning.push_back(current_bin);
      }
      else { // Door was already closed
        move_bins = 1;
        ++state_length;
      }
      is_open = false;
    }

    // Tally up score between now and the new time.
    if (is_open)
      for (int i=0; i<move_bins; ++i)
        total_score += schedule.get_score(i);

    // Advance time.
    current_time += move_bins*micro_bin_time;
    cumulative_bins += move_bins;
    
    // Move schedule time forward.
    schedule.set_current_time(current_time);
    schedule.set_time_sight(std::min(t_sight, time - current_time));
    // Compute new max score/optimal binning. We really only need to recompute this if a new event appeared.
    max_score = schedule.max_score(is_open, state_length);
    binning = schedule.optimal_binning(is_open, state_length);
  }

  // --- Add the rest of the binning to the demon's binning.
  int i=0;
  // If door is open, remember that.
  if (is_open) {
    // If the bin extends onwards.
    if (!binning.empty() && binning[0].first==0) {
      current_bin.second = binning[0].second + cumulative_bins;
      ++i;
    }
    else current_bin.second = cumulative_bins-1;
    finite_demon_binning.push_back(current_bin);
  }
  // Add in the rest of the binning.
  for (; i<binning.size(); ++i)
    finite_demon_binning.push_back(make_pair(binning[i].first + cumulative_bins, binning[i].second + cumulative_bins));

  // Total up the rest of the score.
  total_score += max_score;

  // Return the average score rate for the sight limited demon.
  return total_score / time;
}

void Experiments::generate_events(const float rhoL, const float rhoR) {
  // Generate all events.
  vector<point> left, right;
  setDensities(rhoL, rhoR);
  generate_events(left, right, total_time);
  schedule.set_events(left, right);
}

void Experiments::generate_events(vector<point>& left, vector<point>& right, float time) {
  // Clear vectors.
  left.clear();
  right.clear();
  // Generate some events.
  fill_event_vector(left, rho_left, beta_left, time);
  fill_event_vector(right, rho_right, beta_right, time);
}

void Experiments::generate_events() {
  generate_events(rho_left, rho_right);
}

//! \brief Set the macro bin size.
void Experiments::setGranularity(const int b) {
  schedule.set_granularity(b);
}

//! \brief Set the time length of a macrobin.
void Experiments::setTau(const float t) {
  schedule.set_macro_bin_time(t);
}

//! \brief Set the total time that we generate events for.
void Experiments::setTotalTime(const float t) {
  if (t>0) total_time = t;
}

//! \brief Set the time sight of the schedule.
void Experiments::setTimeSight(const float t) {
  if (t>0) {
    time_sight = t;
    schedule.set_time_sight(t);
    if (time_sight>total_time) 
      setTotalTime(time_sight);
  }
}

void Experiments::setDimensionality(const unsigned int d) {
  dimensionality = d;
  if (d==1) gate_area = 1.f;
}

  //! \brief Set the area of the demon's gate (if dimensionality>1).
void Experiments::setArea(const float A) {
  if (dimensionality>1) gate_area = A;
}

void Experiments::resetTimes() {
  setTimeSight(total_time);
  schedule.set_current_time(0);
}

void Experiments::setSeed(unsigned s) {
  seed = s;
  generator = std::mt19937(seed);
}

void Experiments::setDemonMode(DemonMode mode) {
  demon_mode = mode;
}

void Experiments::setDensities(float rl, float rr) {
  rho_left = rl;
  rho_right = rr;
}

void Experiments::setBetas(float bl, float br) {
  beta_left = bl;
  beta_right = br;
}

inline void Experiments::fill_event_vector(vector<point>& list, const float rho, const float beta, const float time) {
  // Clear list.
  list.clear();
  // Determine the number of events that we should create. The number of events is the same for any demon.
  float point_intensity = rho*gate_area/sqrt(2*PI*particle_mass*beta);
  std::poisson_distribution<int> subsystem_poisson(point_intensity*time);
  int number_of_events = subsystem_poisson(generator);
  // Generate events.
  switch (demon_mode) {
    case (DemonMode::Number): {
      for (int i=0; i<number_of_events; ++i)
        list.push_back(std::make_pair(drand48()*time, 1.f));
      break;
    }
    case (DemonMode::Energy): {
      // E ~ Gamma[(d+1)/2, \beta], and c++ <random> uses 1/\beta for \beta
      energy_distribution = std::gamma_distribution<float>(0.5*(dimensionality+1.), 1./beta); 
      for (int i=0; i<number_of_events; ++i) {
        float t1 = drand48()*time;
        float energy = energy_distribution(generator);
        list.push_back(std::make_pair(t1, energy));
      }
      break;
    }
    default: {
      cout << "Invalid demon mode.\n";
      break;
    }
  }

}
