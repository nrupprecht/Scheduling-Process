#include "../experiments.hpp"

int main(int argc, char **argv) {
  // Parameters.
  int trials = 100;
  int granularity = 10;
  float total_time = 10.f;
  float sight = 1.f;
  float macro_bin_time = 0.05;
  float rhoL = 50, rhoR = 50;
  string name = "data.csv";
  int seed = -1;

  // Parse command line args.
  ArgParse parser(argc, argv);
  parser.get("trials", trials);
  parser.get("bins", granularity);
  parser.get("time", total_time);
  parser.get("sight", sight);
  parser.get("tau", macro_bin_time);
  parser.get("rhoL", rhoL);
  parser.get("rhoR", rhoR);
  parser.get("name", name);
  parser.get("seed", seed);
  // Check for illegal options.
  try {
    parser.check();
  }
  catch (ArgParse::UncheckedToken illegal) {
    cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
    exit(1);
  }

  // Seed random number generator.
  srand48(time(0));
  // Create experiment handler.
  Experiments experiments;
  if (0<=seed) experiments.setSeed(seed);
  experiments.setGranularity(granularity);
  experiments.setTimeSight(total_time);
  experiments.setMacroBinTime(macro_bin_time);

  // ------ Experiments ------ //

  /*
  float max_value = 0.f, finite_value = 0.f;
  for (int i=0; i<trials; ++i) {
    experiments.generate_number_events(rhoL, rhoR);

    float fs = experiments.simulate_number_demon(total_time, sight); // << endl;
    
    experiments.resetTimes();
    auto [ms, bs] = experiments.score();
    max_value += ms; finite_value += fs;
  }
  max_value *= 1./trials; finite_value *= 1./trials;
  cout << finite_value / max_value << endl;
  */

  vector<point> data;
  experiments.number_demon_vary_time_sight(data, total_time, macro_bin_time, 15.f*macro_bin_time, 15.f, trials, rhoL, rhoR);

  // --> Look at the S_t/S^*_t stochastic process
  // experiments.generate_number_events(rhoL, rhoR);
  // experiments.score();
  // auto process = experiments.get_score_structure();
  // print_to_csv(process, "process-rhoL50-rhoR100.csv");

  // --> Generate many scores.
  // vector<float> score_record;
  // experiments.setTimeSight(total_time);
  // experiments.generate_scores(score_record, trials, rhoL, rhoR);
  // print_to_csv(score_record, name);
  
  // --> See how averege performance (per unit time) changes as total time increases.
  // vector<vector<float> > growth;
  // int divisions = 20;
  // float min_t = 0.5f, max_t = 20.f;
  // experiments.score_vary_total_time(growth, rhoL, rhoR, min_t, max_t, divisions, trials);
  // print_to_csv(growth, name);

  
  // --> See how max score varies as the macro bin granularity increases.
  // vector<pair<int,float> > data;
  // experiments.score_vary_bin_granularity(data, 1, granularity, 1, trials, rhoL, rhoR);
  // print_to_csv(data, name);
  

  // --> Look at the average structure record.
  /*
  schedule.set_macro_bin_size(granularity);
  // Run a first time to initialize.
  generate_events(left, right, rhoL, rhoR, total_time);
  schedule.set_events(left, right);
  float ave_score = schedule.max_score();
  auto structure_record = schedule.get_score_structure();

  for (int i=1; i<trials; ++i) {
    // Create new events.
    generate_events(left, right, rhoL, rhoR, total_time);
    schedule.set_events(left, right);
    ave_score += schedule.max_score();
    auto structure = schedule.get_score_structure();
    // Add to structure
    for (int t=0; t<structure_record.size(); ++t) 
      for (int j=0; j<structure_record[t].size(); ++j) {
        structure_record[t][j].first += structure[t][j].first;
        structure_record[t][j].second += structure[t][j].second;
      }
  }
  // Averaging
  for (int t=0; t<structure_record.size(); ++t) {
    for (int j=0; j<structure_record[t].size(); ++j) {
      structure_record[t][j].first /= static_cast<float>(trials);
      structure_record[t][j].second /= static_cast<float>(trials);

      cout << "(" << structure_record[t][j].first << ", " << structure_record[t][j].second << ") ";
    }
    cout << "\n";
  }
  cout << "\nAverage score: " << ave_score / trials << "\n";
  */


  // --> Compute a grid of max and baseline scores for different nl / nr.
  /*
  vector<vector<float> > data2;
  for (int nl=1; nl<40; ++nl)
    for (int nr=1; nr<40; ++nr) {
      vector<float> entry;
      // Choose some densities.
      rhoL = 10*nl;
      rhoR = 10*nr;
      // First two entries are expected number.
      entry.push_back(rhoL*macro_bin_time);
      entry.push_back(rhoR*macro_bin_time);
      // Vary number of bins
      for (int b=2; b<=30; b+=4) {
        float average_best = 0.f, average_base = 0.f;
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
          average_best += max_score;
          average_base += base_line_score;
        }
        average_best /= trials;
        average_base /= trials;
        // Push back bins and average ratio.
        entry.push_back(static_cast<float>(b));
        entry.push_back(average_best);
        entry.push_back(average_base);
      }
      // Add entry
      data2.push_back(entry);
    }
  print_to_csv(data2, "multidata.csv");
  */

  // vector<vector<float> > data2;
  // experiments.create_system_grid(data2, 10, 40, 10, 40, 40);
  // print_to_csv(data2, "multidata-B.csv");

  return 0;
}