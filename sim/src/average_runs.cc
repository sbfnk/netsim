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
  float stopTime = 0.;

  bool do_errors = false;
  bool verbose = false;
  
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
    ("init,i",
     "produce initial conditions file")    
    ("errors,e",
     "include errors in output file")    
    ("verbose,v",
     "produce verbose output")    
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
    verbose = true;
  }

  std::vector< std::vector<float> > values;
  std::vector< std::vector<float> > squares;
  std::vector<unsigned int> no_files;

  unsigned int nFiles = inputFiles.size();
  unsigned int nColumns = 0;
  bool firstFile = true;

  for (std::vector<std::string>::iterator it = inputFiles.begin();
       it != inputFiles.end(); it++) {
    if (verbose) {
      std::cout << "Opening " << (*it) << " for reading" << std::endl;
    }
    std::ifstream ifs((*it).c_str());

    if (ifs.is_open()) {
      std::string line;
      std::getline(ifs, line);
      unsigned int lineCount = 0;
      float currentTime = 0;
      unsigned int currentStep = 0;
      bool firstLine = true;
      
      std::vector<float> line_contents;
      std::vector<float> line_squares;
      std::vector<float> previous_line_contents;
      std::vector<float> previous_line_squares;

      while (!ifs.eof()) {
        
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
            firstLine = false;
          }

          if (nColumns == 0) {
            // set number of columns
            nColumns = line_contents.size();
            if (verbose) {
              std::cout << "Number of columns: " << nColumns << std::endl;
            }
          } else if (nColumns != line_contents.size()) {
            std::cerr << "Warning: wrong number of columns in line "
                      << lineCount << ": " << line_contents.size()
                      << " instead of " << nColumns << std::endl;
          }

          // get time of next line
          currentTime = *(line_contents.begin());
          // need to write data?
          while ((currentTime) > (currentStep*timeStep)) {

            if ((firstFile) || (stopTime < currentStep*timeStep)) {
              values.push_back(std::vector<float>(previous_line_contents.begin()+1,
                                                  previous_line_contents.end()));
              if (do_errors) {
                squares.push_back(std::vector<float>(previous_line_squares.begin()+1,
                                                     previous_line_squares.end()));
              }
              no_files.push_back(1);
              stopTime = currentStep*timeStep;
            } else {
              for (unsigned int i=1; i < nColumns; i++) {
                values[currentStep][i-1] += previous_line_contents[i];
                if (do_errors) {
                  squares[currentStep][i-1] += previous_line_squares[i];
                }
              }
              no_files[currentStep]++;
            }
            ++currentStep;
          }
          previous_line_contents = std::vector<float>(line_contents);
          if (do_errors) {
            previous_line_squares = std::vector<float>(line_squares);
          }
          std::getline(ifs, line);
        }
      }
      firstFile = false;
      ifs.close();
    } else {
      std::cerr << "Error reading " << *it << "." << std::endl;
      --nFiles;
    }
  }
  if (verbose) {
    std::cout << "Stop time: " << stopTime << std::endl;;
  }          
  float time = 0.;
  if (verbose) {
    std::cout << "Opening output file " << outputBase << ".sim.dat" 
              << std::endl;
  }
  std::ofstream ofs((outputBase+".sim.dat").c_str(), std::ios::out);
  
  if (verbose) {
    std::cout << "Number of lines: " << values.size() << std::endl;;
  }          
  if (ofs.is_open()) {
    for (unsigned int i = 0; i < values.size(); i++, time+=timeStep) {
      if (no_files[i] > nFiles/2) {
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
    }
    ofs.close();
  } else {
    std::cerr << "Error writing to " << outputBase  << ".sim.dat" << std::endl;
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
