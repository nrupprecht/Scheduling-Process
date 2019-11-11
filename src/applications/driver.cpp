#include "../experiments.hpp"
#include "../ArgParse.hpp"

#include <fstream>

void generate_events(vector<pair<float,float>>& left, vector<pair<float,float>>& right, float rhoL, float rhoR, float time) {
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

template<typename S, typename T> inline bool print_to_csv(const vector<pair<S, T> >& data, const string& fileName) {
  if (data.empty()) return true;
  std::ofstream fout(fileName);
  if (fout.fail()) {
    cout << "Could not open file [" << fileName << "].\n";
    return false;
  }
  for (const auto& pr : data) fout << pr.first << "," << pr.second << "\n";
  fout.close();
  return true;
}

template<typename T> inline bool print_to_csv(const vector<T>& data, const string& fileName) {
  if (data.empty()) return true;
  std::ofstream fout(fileName);
  if (fout.fail()) {
    cout << "Could not open file [" << fileName << "].\n";
    return false;
  }
  for (int i=0; i<data.size(); ++i) {
    fout << data[i];
    if (i!=data.size()-1) fout << ",";
  }
  fout.close();
  return true;
}

template<typename T> inline bool print_to_csv(const vector<vector<T> >& data, const string& fileName) {
  if (data.empty()) return true;
  std::ofstream fout(fileName);
  if (fout.fail()) {
    cout << "Could not open file [" << fileName << "].\n";
    return false;
  }
  for (const auto& entry : data) {
    for (int i=0; i<entry.size(); ++i) {
      fout << entry[i];
      if (i!=entry.size()-1) fout << ",";
    }
    fout << "\n";
  }
  fout.close();
  return true;
}

int main(int argc, char **argv) {
  // Parameters.
  int trials = 100;
  int macro_bin_size = 10;
  float total_time = 1.f;
  float macro_bin_time = 0.05;
  float rhoL = 50, rhoR = 50;
  string name = "data.csv";

  // Parse command line args.
  ArgParse parser(argc, argv);
  parser.get("trials", trials);
  parser.get("bins", macro_bin_size);
  parser.get("time", total_time);
  parser.get("tau", macro_bin_time);
  parser.get("rhoL", rhoL);
  parser.get("rhoR", rhoR);
  parser.get("name", name);
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
  experiments.setMacroBinSize(macro_bin_size);
  experiments.setTotalTime(total_time);
  experiments.setMacroBinTime(macro_bin_time);


  
  vector<float> score_record;
  experiments.setTotalTime(total_time);
  experiments.generate_scores(score_record, trials, rhoL, rhoR);
  print_to_csv(score_record, "histogram-35.csv");


  // vector<point> growth;
  // int divisions = 20;
  // experiments.score_vary_total_time(growth, rhoL, rhoR, 0.5, 20., divisions, trials);
  // print_to_csv(growth, name);

  
  // --> See how max score varies as the macro bin granularity increases.
  /*
  // Count events.
  vector<pair<int,float> > data;
  for (int b=2; b<=100; ++b) {
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
  // Print data to csv.
  print_to_csv(data, "data.csv");
  */
  

  // --> Look at the average structure record.
  /*
  schedule.set_macro_bin_size(macro_bin_size);
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