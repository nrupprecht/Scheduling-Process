#ifndef __EXPERIMENTS_HPP__
#define __EXPERIMENTS_HPP__

#include "predictive.hpp"

//! \brief Typedef for point type.
typedef pair<float,float> point;

//! \brief This class encapsulates a number of experiments that can be done using the Scheduler.
class Experiments {
public:
  //! \brief Set up random number generator.
  Experiments();

  //! \brief Generate number-demon type events, where the score of each event is 1, and the events are a poisson point process.
  void generate_events(vector<point>&, vector<point>&, float, float, float);

  // --- Experiments

  //! \brief Run a single scheduling simulation.
  //!
  //! The events should be generated before this function is called, e.g. by generate_number_events.
  point score();

  //! \brief Run a simulation many times, fill the vector with the difference between max and baseline scores.
  void generate_scores(vector<float>&, int, float, float);

  //! \brief Grid search across system parameters (for the event densities of the right and left sides) returning the average max and basline
  //! scores for each (rhoL, rhoR) point.
  void create_system_grid(vector<vector<float> >&, float, float, float, float, int);

  //! \brief Fill a vector with data on how the max and baseline scores change as the total time sight varies.
  void score_vary_total_time(vector<vector<float> >&, float, float, float, float, int, int);

  //! \brief Get data that describes the ratio of average max / baseline scores as a function of number of microbins per macrobin.
  void score_vary_bin_granularity(vector<pair<int,float> >&, int, int, int, int, float, float);

  //! \brief Get data that describes the ratio of average max / baseline scores as a function of number of microbins per macrobin. 
  //! Also returns the average max score per unit time and average baseline scores as functions of macrobin granularity.
  void score_vary_bin_granularity(vector<pair<int,float> >&, int, int, int, int, float, float, vector<pair<int,float> >&, vector<pair<int,float> >&);

  void number_demon_vary_time_sight(vector<point>&, float, float, float, int, int, float, float);

  //! \brief Evaluate the score of a demon that can only see a fixed amount of time into the future, and must recalculate the best schedule whenever it can
  //! compared to a demon that can see all the events at once.
  //!
  //! The events should be generated before this function is called, e.g. by generate_number_events.
  float simulate_number_demon(float, float);

  // --- Accessors
  const auto& get_score_structure() { return schedule.get_score_structure(); }

  // --- Mutators

  //! \brief Generate number demon events, poisson point processes with the specified intensities, and each event has a score of 1.
  void generate_number_events(float, float);

  //! \brief Set the macro bin size.
  void setGranularity(int);

  //! \brief Set the time length of a macrobin.
  void setMacroBinTime(float);

  //! \brief Set the total time that we generate events for.
  void setTotalTime(float);

  //! \brief Set the time sight of the schedule.
  void setTimeSight(float);

  //! \brief Set time sight to be the total time, and sets the current time to be zero.
  void resetTimes();

  //! \brief Seed the random number generator
  void setSeed(unsigned);

private:
  // Create scheduler.
  Schedule schedule;

  //! \brief The time sight of the 
  float time_sight = 0.f;
  
  //! \brief The total time the schedule should include.
  float total_time = 0.f;

  //! \brief A seed for the random number generators.
  unsigned seed;

  //! \brief Mersenne twister generator.
  std::mt19937 generator;
};

#endif // __EXPERIMENTS_HPP__