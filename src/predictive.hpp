#ifndef __PREDICTIVE_HPP__
#define __PREDICTIVE_HPP__

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

#include <cmath>

template<typename T> inline T max(const T a, const T b) {
  return a>b ? a : b;
}

typedef pair<float,float> point;

class Schedule {
public:
  //! \brief Default constructor.
  Schedule() {};

  //! \brief Setting constructor.
  Schedule(float, int);

  //! \brief Find the max score for the events.
  float max_score();
  float test_max_score();

  //! \brief Find the (best) baseline score, using as small bins as possible.
  float base_line_score();

  //! \brief Get the score of a set of open intervals.
  float get_score(const vector<pair<int,int> >&);

  //! \brief Make sure a binning is valid, that the lengths of open and closed intervals is at least macro_bin_size.
  bool valid_binning(const vector<pair<int,int> >&);

  //! \brief Return a binning that yields the optimal score.
  //!
  //! Returns the periods where the door was open.
  vector<pair<int, int> > optimal_binning();

  // --- Mutators

  //! \brief Set the left and right events.
  void set_events(const vector<pair<float, float>>&, const vector<pair<float, float>>&);

  //! \brief Set the time sight.
  void set_time_sight(float);

  //! \brief Set the time length of a macro bin.
  void set_macro_bin_time(float);

  //! \brief Set the number of micro bins in a macro bin.
  void set_macro_bin_size(int);

  // --- Accessors

  //! \brief Get the vector of microbin scores.
  const vector<float>& get_micro_bin_scores() { return micro_bin_scores; }

  //! \brief Get the score structure.
  const auto& get_score_structure() { return score_structure; }

  //! \brief Get the macro bin time.
  float get_macro_bin_time() { return macro_bin_time; }

private:

  //! \brief Create and compute the micro bin score vector.
  inline void compute_micro_bin_scores();

  //! \brief Reset everything.
  inline void reset();

  //! \brief The left and right events: times when they occur, and their weights.
  vector<pair<float, float> > left_events, right_events;

  //! \brief The amount of time we can look ahead.
  float time_sight = 0;
  
  //! \brief The time length of a macro bin.
  float macro_bin_time = 0.01;

  //! \brief The minimun number of bins that can be used in a macro bin.
  int macro_bin_size = 5;

  //! \brief The micro bin length, in seconds.
  //!
  //! Must be greater than 0. Will be equal to macro_bin_time / macro_bin_size.
  float micro_bin_time = 0.002;

  //! \brief Total number of micro bins. This should be floor(time_sight / micro_bin_time).
  int total_micro_bins = 0;

  //! \brief The score of each microbin, micro_bin_score[t] = left_events[in bin t] - right_events[in bin t].
  vector<float> micro_bin_scores;

  //! \brief Record the last optimal score. An optimal score of -1 means that it has not been calculated yet.
  float optimal_score = -1.f;

  //! \brief The data structure for finding the best score.
  //!
  //! [bin]  [start time]    {open, closed}
  //!   T    0 - macro_bin   (pair of entries)
  vector<vector<pair<float, float> > > score_structure;

};

#endif // __PREDICTIVE_HPP__