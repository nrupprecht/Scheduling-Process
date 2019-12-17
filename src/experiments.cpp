#include "experiments.hpp"

Experiments::Experiments() {
  seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
  generator = std::mt19937(seed);
}

void Experiments::generate_events(vector<point>& left, vector<point>& right, float rhoL, float rhoR, float time) {
  // Clear vectors.
  left.clear();
  right.clear();
  // Poisson generators.
  std::poisson_distribution<int> left_poisson(rhoL*time);
  std::poisson_distribution<int> right_poisson(rhoR*time);

  // Generate some events.
  int nL = left_poisson(generator);
  int nR = right_poisson(generator);
  for (int i=0; i<nL; ++i) {
    float t1 = drand48()*time;
    left.push_back(std::make_pair(t1, 1.f));
  }
  for (int i=0; i<nR; ++i) {
    float t2 = drand48()*time;
    right.push_back(std::make_pair(t2, 1.f));
  }
}

point Experiments::score() {
  // Find max and baseline scores.
  float ms = schedule.max_score();
  float bs = schedule.base_line_score();
  // Return max and baseline scores.
  return make_pair(ms, bs);
}

void Experiments::generate_scores(vector<float>& score_record, int trials, float rhoL, float rhoR) {
  vector<point> left, right;
  score_record.clear();
  // Calculate data.
  for (int i=0; i<trials; ++i) {
    // Create new events.
    generate_events(left, right, rhoL, rhoR, total_time);
    schedule.set_events(left, right);
    // Find the maximum and baseline score.
    float ms = schedule.max_score(), bs = schedule.base_line_score();
    score_record.push_back(ms - bs);
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
      generate_events(left, right, rhoL, rhoR, total_time);
      schedule.set_events(left, right);
      // Find the maximum and baseline score.
      entry.push_back(rhoL*macro_bin_time);
      entry.push_back(rhoR*macro_bin_time);
      entry.push_back(schedule.max_score());
      entry.push_back(schedule.base_line_score());
      multidata.push_back(entry);
    }
}

void Experiments::score_vary_total_time(vector<vector<float> >& growth, float rhoL, float rhoR, float minT, float maxT, int divisions, int trials) {
  // Time step and vectors for events.
  float dt = (maxT-minT)/(divisions-1);
  vector<point> left, right;
  growth.clear();
  // Calculate data.
  for (int j=0; j<divisions; ++j) {
    float time = j*dt + minT, normalization = 1./(trials*time);;
    float ave_max_score = 0.f, ave_base_line_score = 0.f;
    // Compute average score.
    for (int i=0; i<trials; ++i) {
    // Create new events.
      generate_events(left, right, rhoL, rhoR, time);
      schedule.set_time_sight(time);
      schedule.set_events(left, right);
      // Find the maximum and baseline score.
      ave_max_score += schedule.max_score();
      ave_base_line_score += schedule.base_line_score();
    }
    growth.push_back(vector<float>{time, ave_max_score*normalization, ave_base_line_score*normalization});
  }
}

void Experiments::score_vary_bin_granularity(vector<pair<int,float> >& data, int minB, int maxB, int dB, int trials, float rhoL, float rhoR) {
  // Vectors for events.
  vector<point> left, right;
  data.clear();
  // Calculate data.
  for (int b=minB; b<=maxB; b += dB) {
    float average_max = 0.f, average_base_line = 0.f;
    // Set macro bin size.
    schedule.set_granularity(b);
    // Do some number of trials.
    for (int i=0; i<trials; ++i) {
      // Create new events.
      generate_events(left, right, rhoL, rhoR, total_time);
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
  // Calculate data.
  for (int b=minB; b<=maxB; b += dB) {
    float average_max = 0.f, average_base_line = 0.f;
    // Set macro bin size.
    schedule.set_granularity(b);
    // Do some number of trials.
    for (int i=0; i<trials; ++i) {
      // Create new events.
      generate_events(left, right, rhoL, rhoR, total_time);
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
    max_scores.push_back(make_pair(b, average_max/time/trials));
    base_line_scores.push_back(make_pair(b, average_base_line/time/trials));
  }
}

void Experiments::number_demon_vary_time_sight(vector<point>& data, float time, float min_t_sight, float max_t_sight, int bins, int trials, float rhoL, float rhoR) {
  // Make sure the data vector is empty.
  data.clear();
  // Run all time sights [trials] number of times.
  float dsight = (max_t_sight - min_t_sight)/(bins-1);
  for (int b=0; b<bins; ++b) {
    float t_sight = min_t_sight + b*dsight;
    float ave_score = 0.f;
    for (int iter=0; iter<trials; ++iter) {
      generate_number_events(rhoL, rhoR);
      ave_score += simulate_number_demon(time, t_sight);
    }
    ave_score /= trials;
    data.push_back(pair(t_sight, ave_score));
  }
}

float Experiments::simulate_number_demon(float time, float t_sight) {
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
        total_score += schedule.get_bin_score(i);

    // Advance time.
    current_time += move_bins*micro_bin_time;
    cumulative_bins += move_bins;
    
    // Move schedule time forward.
    schedule.set_current_time(current_time);
    schedule.set_time_sight(std::min(t_sight, time - current_time));
    // Compute new max score/optimal binning.
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
    finite_demon_binning.push_back(pair(binning[i].first + cumulative_bins, binning[i].second + cumulative_bins));

  // Total up the rest of the score.
  total_score += max_score;



  // cout << "Final bit of demon:\n";
  // for (auto p : binning)
  //   cout << "(" << p.first + cumulative_bins << ", " << p.second + cumulative_bins << ") ";
  // cout << endl << "Score: " << max_score << endl << endl;

  // cout << "Final binning:\n";
  // for (auto p : finite_demon_binning)
  //   cout << "(" << p.first << ", " << p.second << ") ";
  // cout << endl;


  // resetTimes();
  // schedule.compute_micro_bin_scores();
  // cout << "Check: " << schedule.get_score(finite_demon_binning) << endl;
  // cout << "Vs: " << total_score << endl;


  // Return the total score for the sight limited demon.
  return total_score;
}

void Experiments::generate_number_events(float rhoL, float rhoR) {
  // Generate all events.
  vector<point> left, right;
  generate_events(left, right, rhoL, rhoR, total_time);
  schedule.set_events(left, right);
}

//! \brief Set the macro bin size.
void Experiments::setGranularity(int b) {
  schedule.set_granularity(b);
}

//! \brief Set the time length of a macrobin.
void Experiments::setMacroBinTime(float t) {
  schedule.set_macro_bin_time(t);
}

//! \brief Set the total time that we generate events for.
void Experiments::setTotalTime(float t) {
  if (t>0) total_time = t;
}

//! \brief Set the time sight of the schedule.
void Experiments::setTimeSight(float t) {
  if (t>0) {
    time_sight = t;
    schedule.set_time_sight(t);
    if (time_sight>total_time) 
      setTotalTime(time_sight);
  }
}

void Experiments::resetTimes() {
  setTimeSight(total_time);
  schedule.set_current_time(0);
}

void Experiments::setSeed(unsigned s) {
  seed = s;
  generator = std::mt19937(seed);
}
