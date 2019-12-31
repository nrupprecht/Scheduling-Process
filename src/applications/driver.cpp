#include "../experiments.hpp"

int main(int argc, char **argv) {
  // Parameters.
  int trials = 100;
  int bins = 25;
  int granularity = 50;
  float total_time = 10.f;
  float sight = 1.f;
  float macro_bin_time = 0.05;
  float rhoL = 50, rhoR = 50;
  string name = "data.csv";
  int seed = -1;
  string mode = "number";

  // Parse command line args.
  ArgParse parser(argc, argv);
  parser.get("trials", trials);
  parser.get("bins", bins);
  parser.get("granularity", granularity);
  parser.get("time", total_time);
  parser.get("sight", sight);
  parser.get("tau", macro_bin_time);
  parser.get("rhoL", rhoL);
  parser.get("rhoR", rhoR);
  parser.get("name", name);
  parser.get("seed", seed);
  parser.get("mode", mode);
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
  experiments.setTau(macro_bin_time);
  experiments.setDensities(rhoL, rhoR);
  // Set demon mode.
  if (mode=="number") experiments.setDemonMode(DemonMode::Number);
  if (mode=="energy") experiments.setDemonMode(DemonMode::Energy);

  // Vectors for data.
  vector<point> data;
  vector<vector<float> > multidata;

  // ------ Experiments ------ //
  float min_tau = 0.01, max_tau = 0.5;
  experiments.demon_vary_tau(data, total_time, sight, min_tau, max_tau, bins, trials);
  print_to_csv(data, name);
  return 0;
  
  // --> Vary system parameters (number density).
  // experiments.score_vary_system(multidata, rhoL, 0.2*rhoL, 5*rhoL, bins, trials);
  // print_to_csv(multidata, name);

  // --> Vary the time sight of a number demon.
  // experiments.number_demon_vary_time_sight(data, total_time, 2*macro_bin_time, 25*macro_bin_time, bins, trials, rhoL, rhoR);
  // print_to_csv(data, name);

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
  float min_t = 0.5f, max_t = 10.f;
  experiments.score_vary_total_time(multidata, rhoL, rhoR, min_t, max_t, bins, trials);
  print_to_csv(multidata, name);

  
  // --> See how max score varies as the macro bin granularity increases.
  // vector<pair<int,float> > p_data;
  // experiments.score_vary_bin_granularity(p_data, 1, granularity, 1, trials, rhoL, rhoR);
  // print_to_csv(p_data, name);
  

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


  // experiments.create_system_grid(multidata, 10, 40, 10, 40, 40);
  // print_to_csv(multidata, "multidata-B.csv");

  return 0;
}