#ifndef __EXPERIMENTS_HPP__
#define __EXPERIMENTS_HPP__

#include "predictive.hpp"

#include "noncentral_chi_squared_zero.hpp"

//! \brief Typedef for point type.
typedef pair<float,float> point;

//! \brief Enum classes for the operating modes of the demon.
enum class DemonMode { Number, Energy };

//! \brief This class encapsulates a number of experiments that can be done using the Scheduler.
class Experiments {
public:
  //! \brief Set up random number generator.
  Experiments();

  // --- Experiments

  //! \brief Run a single scheduling simulation.
  //!
  //! The events should be generated before this function is called, e.g. by generate_events.
  point score();

  //! \brief Run a simulation many times, fill the vector with the difference between max and baseline scores.
  void generate_scores(vector<float>&, int, float, float, bool=true);

  //! \brief Grid search across system parameters (for the event densities of the right and left sides) returning the average max and basline
  //! scores for each (rhoL, rhoR) point.
  void create_system_grid(vector<vector<float> >&, float, float, float, float, int);

  //! \brief Fill a vector with data on how the score changes as the response time is changed (sight is unlimited).
  void score_vary_tau(vector<point>&, float, float, float, float, int, int);

  //! \brief Fill a vector with data on how the max and baseline scores change as the total time sight varies.
  void score_vary_total_time(vector<vector<float> >&, float, float, float, float, int, int);

  //! \brief Fill a vector with data on how the max and baseline scores change as the right subsystem varies in density.
  void score_vary_system(vector<vector<float> >&, float, float, float, int, int);

  //! \brief Get data that describes the ratio of average max / baseline scores as a function of number of microbins per macrobin.
  void score_vary_bin_granularity(vector<pair<int,float> >&, int, int, int, int, float, float);

  //! \brief Get data that describes the ratio of average max / baseline scores as a function of number of microbins per macrobin. 
  //! Also returns the average max score per unit time and average baseline scores as functions of macrobin granularity.
  void score_vary_bin_granularity(vector<pair<int,float> >&, int, int, int, int, float, float, vector<pair<int,float> >&, vector<pair<int,float> >&);

  //! \brief Simulate a number demon with various time sights.
  void demon_vary_time_sight(vector<point>&, float, float, float, int, int, float, float);

  //! \brief Repeatedly simulate a demon, varying its response time.
  void demon_vary_tau(vector<point>&, float, float, float, float, int, int);

  //! \brief Evaluate the score of a demon that can only see a fixed amount of time into the future, and must recalculate the best schedule whenever it can
  //! compared to a demon that can see all the events at once.
  //!
  //! The events should be generated before this function is called, e.g. by generate_events.
  //!
  //! \param time The total amount of time in which we generate.
  //! \param t_sight How far the demon can see into the future.
  float simulate_demon(float, float);

  // --- Accessors
  const auto& get_score_structure() { return schedule.get_score_structure(); }

  // --- Mutators

  //! \brief Generate demon events, changing the densities to the specified values.
  void generate_events(const float, const float);

  //! \brief Generate demon type events and place them into the vectors.
  void generate_events(vector<point>&, vector<point>&, float);

  //! \brief Generate events based on the system parameters.
  void generate_events();

  //! \brief Set the macro bin size.
  void setGranularity(const int);

  //! \brief Set the time length of a macrobin.
  void setTau(const float);

  //! \brief Set the total time that we generate events for.
  void setTotalTime(const float);

  //! \brief Set the time sight of the schedule.
  void setTimeSight(const float);

  //! \brief Set the dimensionality of the simulation.
  void setDimensionality(const unsigned int);

  //! \brief Set the area of the demon's gate (if dimensionality>1).
  void setArea(const float);

  //! \brief Set time sight to be the total time, and sets the current time to be zero.
  void resetTimes();

  //! \brief Seed the random number generator
  void setSeed(unsigned);

  //! \brief Set the demon mode.
  void setDemonMode(DemonMode);

  //! \brief Set the system densities.
  void setDensities(float, float);

  //! \brief Set the system inverse temperatures.
  void setBetas(float, float);

private:
  //! \brief A helper function that fills a vector with events.
  inline void fill_event_vector(vector<point>&, const float, const float, const float);

  //! \brief A scheduler class.
  Schedule schedule;

  //! \brief The mode that the demon is operating in (for generating events).
  //!
  //! Modes: 0 (number demon), 1 (energy demon).
  DemonMode demon_mode = DemonMode::Number;

  //! \brief The dimensionality of the system. Used to simulate an energy demon.
  unsigned int dimensionality = 2;

  //! \brief The time sight of the the demon.
  float time_sight = 0.f;
  
  //! \brief The total time the schedule should include.
  float total_time = 0.f;

  //! \brief Number densities and (inverse) temperatures of the subsystems.
  float rho_left = 100.f, rho_right = 100.f, beta_left = 1.f, beta_right = 1.f;

  //! \brief The area of the gate. This will be 1 for a 1 dimensional simulation.
  float gate_area = 1.f;

  //! \brief The value we use for the particles' mass. For generating energy demons.
  float particle_mass = 1.f;

  //! \brief A seed for the random number generators.
  unsigned seed;

  //! \brief A standard normal distribution.
  std::normal_distribution<float> standard_normal;

  std::gamma_distribution<float> energy_distribution;

  //! \brief Mersenne twister generator.
  std::mt19937 generator;
};

#endif // __EXPERIMENTS_HPP__