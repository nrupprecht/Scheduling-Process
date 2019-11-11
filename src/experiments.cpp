#include "experiments.hpp"

void Experiments::generate_events(vector<point>& left, vector<point>& right, float rhoL, float rhoR, float time) {
  // Clear vectors.
  left.clear();
  right.clear();
  // Generate some events.
  int nL = rhoL*time, nR = rhoR*time;
  for (int i=0; i<nL; ++i) {
    float t1 = drand48()*time;
    left.push_back(std::make_pair(t1, 1.f));
  }
  for (int i=0; i<nR; ++i) {
    float t2 = drand48()*time;
    right.push_back(std::make_pair(t2, 1.f));
  }
}

void Experiments::generate_scores(vector<float>& score_record, int trials, float rhoL, float rhoR) {
  vector<point> left, right;
  for (int i=0; i<trials; ++i) {
    // Create new events.
    generate_events(left, right, rhoL, rhoR, total_time);
    schedule.set_events(left, right);
    // Find the maximum and baseline score.
    float ms = schedule.test_max_score(), bs = schedule.base_line_score();
    score_record.push_back(ms - bs);
  }
}

void Experiments::create_system_grid(vector<vector<float> >& multidata, float minRL, float maxRL, float minRR, float maxRR, int divisions) {

  float dl = (maxRL-minRL)/(divisions-1);
  float dr = (maxRR-minRR)/(divisions-1);
  vector<point> left, right;

  float macro_bin_time = schedule.get_macro_bin_time();

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

void Experiments::score_vary_total_time(vector<point>& growth, float rhoL, float rhoR, float minT, float maxT, int divisions, int trials) {

  float dt = (maxT-minT)/(divisions-1);
  vector<point> left, right;

  for (int j=0; j<divisions; ++j) {
    float time = j*dt + minT;
    float ave_score = 0;
    // Compute average score.
    for (int i=0; i<trials; ++i) {
    // Create new events.
      generate_events(left, right, rhoL, rhoR, time);
      schedule.set_time_sight(time);
      schedule.set_events(left, right);
      // Find the maximum and baseline score.
      ave_score += (schedule.max_score() - schedule.base_line_score());
    }
    growth.push_back(make_pair(time, ave_score/trials/time));
  }
}

void Experiments::score_vary_bin_granularity(vector<pair<int,float> >& data, int minB, int maxB, int dB, int trials, float rhoL, float rhoR) {

  vector<point> left, right;

  for (int b=minB; b<=maxB; b += dB) {
    float average_ratio = 0.f;
    // Set macro bin size.
    schedule.set_macro_bin_size(b);
    // Do some number of trials.
    for (int i=0; i<trials; ++i) {
      // Create new events.
      generate_events(left, right, rhoL, rhoR, total_time);
      schedule.set_events(left, right);
      // Find the maximum and baseline score.
      float max_score = schedule.max_score();
      float base_line_score = schedule.base_line_score();
      // Calculate ratio.
      average_ratio += max_score/base_line_score;
    }
    average_ratio /= trials;
    data.push_back(make_pair(b, average_ratio));
  }
}
