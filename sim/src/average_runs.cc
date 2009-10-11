/*! \file average_runs.cc
  \brief A program to average a given number of data files.
*/
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

void split(const std::string& str,const std::string delimiters,
           std::vector<std::string>& tokens)
{
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos) {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos     = str.find_first_of(delimiters, lastPos);
  }
}

void split_to_float(const std::string &str, const std::string delimiters,
                     std::vector<float>& values)
{
  std::vector<std::string> tokens;
  split(str, delimiters, tokens);
  std::vector<std::string>::iterator it;
  float f;

  for (it = tokens.begin(); it != tokens.end(); it++) {
    std::stringstream ss(*it);
    ss >> f;
    values.push_back(f);
  }
}


int main(int argc, char* argv[])
{

  float timeStep = 0.;
  float maxTime = 0.;
  
  unsigned int maxStep = 0;

  bool do_errors = false;
  bool alive = false;
  unsigned int verbose = 0;
  
  std::vector<std::string> inputFiles;
  std::string outputBase;
  
  po::options_description command_line_options
    ("Usage: average_runs -d arg -s arg -o arg [input-files]\n\nAllowed options");
  command_line_options.add_options()
    ("help", "produce help message")    
    ("time-step,d", po::value<float>(),
     "discretized time step")    
    ("output-file,o", po::value<std::string>(),
     "output file")    
    ("steps-file,s", 
     "write steps into a file")    
    ("init,i",
     "produce initial conditions file")    
    ("errors,e",
     "include errors in output file")    
    ("verbose,v",
     "produce verbose output")    
    ("very-verbose,V",
     "produce very verbose output")    
    ("max,m", po::value<float>(),
     "maximum time")    
    ("alive,a", 
     "condition on the process still going")    
    ;
    
  po::options_description hidden_options;
  
  hidden_options.add_options()
    ("input-files", po::value< std::vector<std::string> >(), "input files")
    ;

  po::options_description all_options;
  all_options.add(command_line_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);
        
  po::variables_map vm;
  store(po::command_line_parser(argc, argv).
        options(all_options).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("time-step")) {
    timeStep = vm["time-step"].as<float>();
  } else {
    std::cerr << "ERROR: no time-step specified" << std::endl;
    std::cerr << command_line_options << std::endl;
    return 1;
  }
  
  if (vm.count("output-file")) {
    outputBase = vm["output-file"].as<std::string>();
  } else {
    std::cerr << "ERROR: no output file specified" << std::endl;
    std::cerr << command_line_options << std::endl;
    return 1;
  }
  
  if (vm.count("input-files")) {
    inputFiles = vm["input-files"].as< std::vector <std::string> >();
  } else{
    std::cerr << "ERROR: no input files specified" << std::endl;
    std::cerr << command_line_options << std::endl;
    return 1;
  }

  if(vm.count("errors")) {
    do_errors = true;
  }

  if(vm.count("verbose")) {
    verbose = 1;
  }

  if(vm.count("very-verbose")) {
    verbose = 2;
  }

  if(vm.count("max")) {
    maxTime = vm["max"].as<float>();
  }

  if(vm.count("alive")) {
    alive = true;
  }

  if (maxTime > 0) {
    maxStep = static_cast<unsigned int>(maxTime / timeStep);
  } else {
    std::cerr << "ERROR: no --max specified" << std::endl;
    return 0;
  }
  
  std::vector< std::vector<float> > values(maxStep+1);
  std::vector< std::vector<float> > squares(maxStep+1);
  std::vector< unsigned int > no_files(maxStep+1, 0);

  unsigned int nFiles = inputFiles.size();
  unsigned int nColumns = 0;

  for (std::vector<std::string>::iterator it = inputFiles.begin();
       it != inputFiles.end(); it++) {
    std::string baseName = it->substr(0,it->rfind(".dat", it->size()));
    if (verbose) {
      std::cout << "Opening " << (*it) << " for reading" << std::endl;
    }
    std::ifstream ifs((*it).c_str());
    std::ofstream ofs_steps;
    if (vm.count("steps-file")) {
      ofs_steps.open((baseName + ".steps").c_str(), std::ios::out);
    }

    if (ifs.is_open()) {
      std::string line;
      unsigned int lineCount = 0;
      float currentTime = 0;
      unsigned int currentStep = 0;
      bool firstLine = true;
      
      std::vector<float> line_contents;
      std::vector<float> line_squares;
      std::vector<float> previous_line_contents;
      std::vector<float> previous_line_squares;

      while (!ifs.eof() && (currentStep <= maxStep)) {
        
        std::getline(ifs, line);
        lineCount++;
        line_contents.clear();
        split_to_float(line, " \t", line_contents);
        if (line_contents.size() > 0) {
          if (do_errors) {
            line_squares.clear();
            std::vector<float>::iterator it = line_contents.begin();
            line_squares.push_back(*it);
            for (it++;it != line_contents.end(); it++) {
              line_squares.push_back((*it)*(*it));
            }
          }

          if (firstLine) {
            previous_line_contents = std::vector<float>(line_contents);
            if (do_errors) {
              previous_line_squares = std::vector<float>(line_squares);
            }
          }

          if (nColumns == 0) {
            // set number of columns
            nColumns = line_contents.size();
            if (verbose) {
              std::cout << "Number of columns: " << nColumns << std::endl;
            }
            for (unsigned int i = 0; i < maxStep + 1; ++i) {
              values[i].resize(nColumns-1, 0.);
              squares[i].resize(nColumns-1, 0.);
            }
          } else if (nColumns != line_contents.size()) {
            std::cerr << "Warning: wrong number of columns in line "
                      << lineCount << ": " << line_contents.size()
                      << " instead of " << nColumns << std::endl;
          }

          // get time of next line
          currentTime = *(line_contents.begin());
          // need to write data?
          while (currentTime >= currentStep * timeStep &&
                 currentStep <= maxStep) {
            if (verbose >= 2) {
              std::cout << "processing data: current " << currentTime
                        << " time step " << currentStep << std::endl;
            }

            // interpolate values
            std::vector<float> interpol(nColumns-1);
            std::vector<float> interpol_sq(nColumns-1);
            if (currentTime == currentStep * timeStep) {
              interpol = std::vector<float>(line_contents.begin()+1,
                                            line_contents.end());
            } else {
              for (unsigned int i = 0; i < nColumns-1; ++i) {
                interpol[i] = previous_line_contents[i+1] +
                  (line_contents[i+1] - previous_line_contents[i+1]) *
                  (currentStep * timeStep - previous_line_contents[0]) /
                  (currentTime - previous_line_contents[0]);
              }
            }
            if (do_errors) {
              for (unsigned int i = 0; i < nColumns-1; ++i) {
                interpol_sq[i] = interpol[i] * interpol[i];
              }
            }
            
            for (unsigned int i=0; i < nColumns-1; i++) {
              values[currentStep][i] += interpol[i];
              if (do_errors) {
                squares[currentStep][i] += interpol_sq[i];
              }
            }
            if (verbose >= 2) {
              std::cout << currentStep*timeStep;
              for (unsigned int j = 1; j < nColumns ; j++) {
                std::cout << " " << values[currentStep][j-1];
              }
              std::cout << std::endl;
            }
            ++no_files[currentStep];
            ++currentStep;
	    firstLine = false;
          }
          if (ofs_steps.is_open()) {
            ofs_steps << previous_line_contents[0] << " " << currentStep - 1 << std::endl;
          }
          previous_line_contents = std::vector<float>(line_contents);
          if (do_errors) {
            previous_line_squares = std::vector<float>(line_squares);
          }
        }
      }
      if (!alive) {
        while (currentStep < maxStep) {
          ++currentStep;
          ++no_files[currentStep];
          for (unsigned int i=0; i < nColumns-1; i++) {
            values[currentStep][i] += line_contents[i+1];
            if (do_errors) {
              squares[currentStep][i] += line_squares[i+1];
            }
          }
        }
      }
      if (ofs_steps.is_open()) {
        ofs_steps.close();
      }
      ifs.close();
    } else {
      std::cerr << "Error reading " << *it << "." << std::endl;
      --nFiles;
    }
  }
  float time = 0.;
  if (verbose) {
    std::cout << "Opening output file " << outputBase << std::endl;
  }
  std::ofstream ofs(outputBase.c_str(), std::ios::out);
  
  if (verbose) {
    std::cout << "Number of lines: " << values.size() << std::endl;;
  }          
  if (ofs.is_open()) {
    for (unsigned int i = 0; i < values.size(); i++, time+=timeStep) {
      if (verbose >=2 ) {
        std::cout << "Writing timeStep " << time << ", "
                  << no_files[i] << " files" << std::endl;
      }
      ofs << time;
      for (unsigned int j = 0; j < values[i].size() ; j++) {
        ofs << " " << values[i][j]/no_files[i];
      }
      if (do_errors) {
        for (unsigned int j = 0; j < squares[i].size() ; j++) {
          double avg_value = values[i][j]/no_files[i];
          ofs << " " << sqrt(squares[i][j]/no_files[i] - avg_value*avg_value);
        }
      }
      ofs << std::endl;
    }
    ofs.close();
  } else {
    std::cerr << "Error writing to " << outputBase << std::endl;
    return 1;
  }

  if (vm.count("init")) {
    // write file with averaged initial conditions
    if (verbose) {
      std::cout << "Opening " << outputBase << ".init to write initial conditions"
                << std::endl;
    }
    ofs.open((outputBase+".init").c_str(), std::ios::out);
    if (ofs.is_open()) {
      std::vector<float> firstLine = (*values.begin());
      for (std::vector<float>::iterator it = firstLine.begin();
           it != firstLine.end(); it++) {
        ofs << (*it)/nFiles << std::endl;
      }
      ofs << std::endl;
      ofs.close();
    } else {
      std::cerr << "Error writing to " << outputBase << ".init" << std::endl;
      return 1;
    }
  }
  
  return 0;
}
