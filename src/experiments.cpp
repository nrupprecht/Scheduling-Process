#include "experiments.hpp"

Experiments::Experiments() {
  if(seed==0) seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
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

point Experiments::run_single(float rhoL, float rhoR) {
  // Generate some events.
  vector<point> left, right;
  generate_events(left, right, rhoL, rhoR, total_time);
  schedule.set_events(left, right);
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

void Experiments::simulate_number_demon(float time, float t_sight, float rhoL, float rhoR) {

  setTotalTime(time);
  setTimeSight(t_sight);
  float current_time = 0.f;

  // Generate all events.
  vector<point> left, right;
  generate_events(left, right, rhoL, rhoR, total_time);
  schedule.set_events(left, right);

  // Generate first guess at a good schedule.
  cout << schedule.max_score() << endl;
  auto binning = schedule.optimal_binning();

  while (current_time + time_sight < total_time) {
    current_time += 1.f;
    schedule.set_current_time(current_time);

    float max_s = schedule.max_score();

    cout << max_s << endl;
  }



}
