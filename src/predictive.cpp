#include "predictive.hpp"

Schedule::Schedule(float microtime, int macrobin) : micro_bin_time(microtime), macro_bin_size(macrobin) {
  if (micro_bin_time<=0) {
    cout << "Microbin time must be > 0. Exiting.\n";
    exit(0);
  }
}

void Schedule::set_events(const vector<pair<float, float>>& left, const vector<pair<float, float>>& right) {
  left_events = left;
  right_events = right;
  reset();
}

void Schedule::set_time_sight(float t) {
  if (t>0 && t!=time_sight) {
    time_sight = t;
    reset();
  }
  else if (t<=0) cout << "Time sight should be > 0.\n";
}

void Schedule::set_current_time(float t) {
  if (t!=current_time) {
    current_time = t;
    reset();
  }
}

void Schedule::set_macro_bin_time(float t) {
  if (t>0 && t!=macro_bin_time) {
    macro_bin_time = t;
    reset();
  }
  else if (t<=0) cout << "Macro bin time should be > 0.\n";
}

void Schedule::set_granularity(int n) {
  if (n>0) {
    macro_bin_size = n;
    reset();
  }
  else cout << "Macro bin size should be > 0.\n";
}

void Schedule::set_score_mode(ScoreMode mode) {
  score_mode = mode;
}

float Schedule::max_score() {
  if (optimal_score!=-1) return optimal_score;
  if (total_micro_bins<=0) {
    cout << "total_micro_bins is <= 0. Error, exiting.\n";
    exit(0);
  }
  if (macro_bin_size<=0) {
    cout << "macro_bin_size is <= 0. Error, exiting.\n";
    exit(0);
  }

  // If there is no time (aka microbins), return 0.
  if (total_micro_bins==0) return 0.f;

  // Compute the score structure.
  compute_score_structure();

  // Since we allow the first bin to be 1 - [macro_bin_size-1] in length, we well as >= [macro_bin_size] lengths, we need to check what
  // scores a shorter bin would allow us to obtain.
  float cumulative_score = 0.f;
  optimal_score = 0.f;
  for (int i=macro_bin_size-1; 1<=i; --i) {
    // Update cumulative score.
    cumulative_score += get_score(macro_bin_size-i-1);
    // Open
    float s1 = cumulative_score + score_structure[macro_bin_size-i].first;
    // Closed
    float s2 = score_structure[macro_bin_size-i].second;
    optimal_score = std::max(std::max(s1, s2), optimal_score);
  }
  optimal_score = std::max(std::max(score_structure[0].first, score_structure[0].second), optimal_score);

  // Return the optimal score.
  return optimal_score;
}

float Schedule::max_score(bool is_open, int micro_bins_in_current_state) {
  // If an optimal score is stored, then the score structure is as well.
  compute_score_structure();

  // For how many micro bins must the gate remain open/closed.
  int remaining_bins_in_state = std::max(macro_bin_size - micro_bins_in_current_state, 0);

  // If we can choose freely whether we want the gate to be open or closed from the start, then the conditioning does not matter.
  if (remaining_bins_in_state==0 || micro_bins_in_current_state==-1) {
    if (0<=optimal_score) return optimal_score;
    else return max_score();
  }
  // Reset optimal score.
  optimal_score = 0.f;

  // Compute optimal score.
  if (is_open) {
    for (int i=0; i<remaining_bins_in_state && i<micro_bin_scores.size(); ++i) optimal_score += get_score(i);
    // After having to leave the door open for the required amount of time, choose the best thing to do.
    if (remaining_bins_in_state<score_structure.size())
      optimal_score += std::max(score_structure[remaining_bins_in_state].first, score_structure[remaining_bins_in_state].second);
  }
  else { // Gate has been closed
    if (remaining_bins_in_state<score_structure.size())
      optimal_score = score_structure[remaining_bins_in_state].second;
  }

  // Return the optimal score.
  return optimal_score;
}

float Schedule::max_score(bool is_open) {
  // If an optimal score is stored, then the score structure is as well.
  compute_score_structure();

  if (is_open) optimal_score = score_structure[0].first;
  else optimal_score = score_structure[0].second;

  return optimal_score;
}

float Schedule::base_line_score() {
  // Make sure micro bin scores have been computed.
  compute_micro_bin_scores();

  float baseline_score = 0, test_score = 0;
  // The first bin can be 1 - macro_bin_size, so we try them all.
  for (int start=0; start<macro_bin_size; ++start) {
    test_score = 0;
    // Find total for the initial macro bin.
    float score = 0.f;
    for (int t=0; t<=start; ++t) score += get_score(t);

    // Decide whether macro bin should be open or closed.
    test_score += std::max(0.f, score);
    // Find total for each macro bin.
    for (int t=start+1; t<total_micro_bins; t+=macro_bin_size) {
      // Compute score for this macro bin.
      score = 0.f;
      for (int i=t; i<total_micro_bins && i<t+macro_bin_size; ++i) score += get_score(i);
      // Decide whether macro bin should be open or closed.
      test_score += std::max(0.f, score);
    }
    // Done with this test score. Is it better than any previous test score?
    baseline_score = std::max(test_score, baseline_score);
  }
  // Return.
  return baseline_score;
}

float Schedule::get_binning_score(const vector<pair<int,int>>& binning, ScoreMode mode) {
  // Make sure micro bin scores have been computed.
  compute_micro_bin_scores();
  // Save the old score mode, use the specified score mode.
  ScoreMode save = score_mode;
  if (mode!=ScoreMode::Default) score_mode = mode;
  // Compute the score the binning has.
  // Calculate the score.
  float score = 0;
  for (auto pr : binning) {
    int start = std::max(pr.first, 0), end = std::min(pr.second, static_cast<int>(micro_bin_scores.size()-1));
    for (int i=start; i<=end; ++i) score += get_score(i);
  }
  // Reset the score mode to what it was previously.
  score_mode = save;
  // Return the score.
  return score;
}

float Schedule::get_binning_score(const vector<pair<int,int> >& binning, float min_t, float max_t, ScoreMode mode) {
  // Make sure micro bin scores have been computed.
  compute_micro_bin_scores();
  // Calculate the score.
  float score = 0;
  int min_bin = (min_t - current_time)/micro_bin_time;
  int max_bin = (max_t - current_time)/micro_bin_time;
  // Save the old score mode, use the specified score mode.
  ScoreMode save = score_mode;
  if (mode!=ScoreMode::Default) score_mode = mode;
  // Compute the score the binning has.
  for (auto pr : binning) {
    int start = std::max(pr.first, 0), end = std::min(pr.second, static_cast<int>(micro_bin_scores.size()-1));
    for (int i=std::max(start, min_bin); i<=std::min(end, max_bin); ++i) score += get_score(i);
  }
  // Reset the score mode to what it was previously.
  score_mode = save;
  // Return the score.
  return score;
}

point Schedule::get_both_binning_score(const vector<pair<int,int> >& binning) {
  float score = get_binning_score(binning, ScoreMode::Score);
  float net_events = get_binning_score(binning, ScoreMode::Number);
  return point(score, net_events);
}

float Schedule::get_score(const int index) {
  switch (score_mode) {
    case ScoreMode::Default:
    case ScoreMode::Score:
      return micro_bin_scores[index];
    case ScoreMode::Number:
      return micro_bin_counts[index];
    default: {
      cout << "Error. Undefined score mode.\n";
      return 0;
    }
  }
}

bool Schedule::valid_binning(const vector<pair<int,int> >& binning) {
  // Check if any of the bins are too small.
  for (int i=1; i<binning.size(); ++i) {
    if (binning[i].first - binning[i-1].second < macro_bin_size ||
      binning[i].second - binning[i].first < macro_bin_size) 
      return false;
  }
  // Otherwise, the binning is valid.
  return true;
}

vector<pair<int, int> > Schedule::optimal_binning(bool is_open, int micro_bins_in_current_state) {
  // If the optimal score has not yet been calculated, calculate it.
  max_score(is_open, micro_bins_in_current_state);

  // Keep track of whether the door is open or closed.
  bool open = false;
  // s0 is how long the door has been open before the start of time.
  int s0 = 0, time = 0, open_time = 0;

  // If allowed, find the optimal state to start in (door can be open less than a full macrobin length).
  if (micro_bins_in_current_state<0 || macro_bin_size<=micro_bins_in_current_state) {
    float cumulative_score = 0.f;
    for (int i=0; i<macro_bin_size; ++i) {
      // Update cumulative score.
      cumulative_score += get_score(i);
      // Candidate scores.
      float s1 = cumulative_score + score_structure[i+1].first;
      float s2 = score_structure[i+1].second;
      // Check if one is an optimal score (and therefore a valid starting point for an optimal binning).
      if (s1==optimal_score) {
        time = i;
        open = true;
        break;
      }
      if (s2==optimal_score) {
        time = i;
        open = false;
        break;
      }
    }
    // No partial opening time achieved an optimal score. The optimal score must start from the beginning (time=0).
    if (time==0 && score_structure[0].first>score_structure[0].second) open = true;
  }
  // Otherwise, there are restrictions on how whether the door was open or closed, and for how long it was in that state.
  else {
    open = is_open;
    time = macro_bin_size - micro_bins_in_current_state;
  }

  // Move forward in time, finding a schedule that would result in the optimal score.
  vector<pair<int, int> > binning;
  while (time<total_micro_bins-1) {
    if (open) {
      if (score_structure[time+1].first + get_score(time)==score_structure[time].first) ++time; // remain open
      // Close
      else {
        open = false;
        if (time>0) binning.push_back(std::make_pair(open_time, time-1));
        time += (macro_bin_size);
      }
    }
    // Door was closed.
    else {
      if (score_structure[time+1].second==score_structure[time].second) ++time; // Remain closed.
      else {
        open = true;
        open_time = time;
        time += (macro_bin_size);
      }
    }
  }
  // If we need to check what happens for the last bins.
  if (time==total_micro_bins-1) {
    if (open) {
      if (get_score(total_micro_bins-1)>=0) binning.push_back(std::make_pair(open_time, total_micro_bins-1));
      else binning.push_back(std::make_pair(open_time, total_micro_bins-2));
    }
    else if (get_score(total_micro_bins-1)>0) binning.push_back(std::make_pair(total_micro_bins-1, total_micro_bins-1));
  }
  // If bin ends open, set the "closing time" to be the total time.
  else if (open) binning.push_back(std::make_pair(open_time, total_micro_bins-1));

  // Return the optimal binning.
  return binning;
}

void Schedule::compute_micro_bin_scores(int micro_bins) {
  if (micro_bins<0) micro_bins = total_micro_bins;
  // Find a score for each micro bin.
  micro_bin_scores = vector<float>(micro_bins, 0);
  micro_bin_counts = vector<int>(micro_bins, 0);
  for (const auto pr : left_events) {
    float t = pr.first, sc = pr.second;
    int bin = static_cast<int>(floor((t-current_time)/micro_bin_time));
    if (0<=bin && bin<micro_bins) {
      micro_bin_scores[bin] += sc;
      ++micro_bin_counts[bin];
    }
  }
  for (const auto pr : right_events) {
    float t = pr.first, sc = pr.second;
    int bin = static_cast<int>(floor((t-current_time)/micro_bin_time));
    if (0<=bin && bin<micro_bins) {
      micro_bin_scores[bin] -= sc;
      --micro_bin_counts[bin];
    }
  }
}

inline void Schedule::compute_score_structure(int micro_bins) {
  if (micro_bins<0) micro_bins = total_micro_bins;
  else {
    score_structure.clear();
    micro_bin_scores.clear();
  }

  if (score_structure.empty()) score_structure = vector<point>(micro_bins, make_pair(0.f, 0.f));
  // Make sure micro bin scores have been computed.
  if (micro_bin_scores.empty()) compute_micro_bin_scores(micro_bins);
  
  // Compute window score.
  vector<float> window_score(micro_bins, 0.f);
  float cumulative_score = 0.f;
  for (int i=0; i<macro_bin_size; ++i) cumulative_score += get_score(i);
  window_score[0] = cumulative_score;
  for (int i=1; i<micro_bins - macro_bin_size; ++i) {
    // Loss from the score at that is now behind the bin.
    cumulative_score -= get_score(i-1);
    // Gain from the score at the end of the bin.
    if (i+macro_bin_size-1<micro_bins) cumulative_score += get_score(i+macro_bin_size-1);
    window_score[i] = cumulative_score;
  }

  // Initialize last [macro_bin_size] entries.
  int i = 1, index = micro_bins - 1;
  cumulative_score = get_score(index);
  // Last entry
  score_structure[index].first = std::max(cumulative_score, 0.f);
  score_structure[index].second = std::max(cumulative_score, 0.f);
  --index;
  for (; i<macro_bin_size && i<micro_bins; ++i, --index) {
    cumulative_score += get_score(index);
    score_structure[index].first = std::max(get_score(index) + score_structure[index+1].first, 0.f);
    score_structure[index].second = std::max(cumulative_score, score_structure[index + 1].second);
  }
  // Propogate solution backwards.
  index = micro_bins - 1 - macro_bin_size;
  for (; i<micro_bins; ++i) {
    score_structure[index].first = std::max(score_structure[index+1].first + get_score(index), score_structure[index+macro_bin_size].second);
    score_structure[index].second = std::max(score_structure[index+macro_bin_size].first + window_score[index], score_structure[index+1].second);
    --index;
  }
}

inline void Schedule::reset() {
  // Recalculate times and sizes.
  micro_bin_time = macro_bin_time / macro_bin_size;
  total_micro_bins = static_cast<int>(ceil(time_sight / micro_bin_time));
  // Clear optimal score and data structures.
  optimal_score = -1.f;
  micro_bin_scores.clear();
  micro_bin_counts.clear();
  score_structure.clear();
}
