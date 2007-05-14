#include <iostream>
#include <fstream>
#include <string>

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
    ("errors,e",
     "include errors in output file")    
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

  std::vector< std::vector<float> > values;
  std::vector< std::vector<float> > squares;

  unsigned int nFiles = inputFiles.size();
  unsigned int nColumns = 0;
  bool firstFile = true;

  for (std::vector<std::string>::iterator it = inputFiles.begin();
       it != inputFiles.end(); it++) {
    std::ifstream ifs((*it).c_str());
    
    if (ifs.is_open()) {
      std::string line;
      std::getline(ifs, line);
      unsigned int lineCount = 1;
      float currentTime = 0;
      unsigned int currentStep = 0;
      
      std::vector<float> line_contents;
      std::vector<float> line_squares;
      std::vector<float> previous_line_contents;
      std::vector<float> previous_line_squares;

      while (!ifs.eof() && (firstFile || (currentTime <= stopTime))) {
        
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

          if (nColumns == 0) {
            // set number of columns
            nColumns = line_contents.size();
          } else if (nColumns != line_contents.size()) {
            std::cerr << "Warning: wrong number of columns in line "
                      << lineCount << ": " << line_contents.size()
                      << " instead of " << nColumns << std::endl;
          }
          
          currentTime = *(line_contents.begin());
          
          if (currentTime >= currentStep*timeStep) {
            while ((currentTime-timeStep) > (currentStep*timeStep)
                   && (firstFile || (stopTime > currentStep*timeStep))) {
              if (firstFile) {
                values.push_back(std::vector<float>(previous_line_contents.begin()+1,
                                                    previous_line_contents.end()));
                if (do_errors) {
                  squares.push_back(std::vector<float>(previous_line_squares.begin()+1,
                                                       previous_line_squares.end()));
                }
                stopTime = (currentStep+1)*timeStep;
              } else {
                for (unsigned int i=1; i < nColumns; i++) {
                  values[currentStep][i-1] += previous_line_contents[i];
                  if (do_errors) {
                    squares[currentStep][i-1] += previous_line_squares[i];
                  }
                }
              }
              ++currentStep;
            }
            
            // create line vector
            if (firstFile) {
              values.push_back(std::vector<float>(line_contents.begin()+1,
                                                  line_contents.end()));
              if (do_errors) {
                squares.push_back(std::vector<float>(line_squares.begin()+1,
                                                     line_squares.end()));
              }
              stopTime = (currentStep+1)*timeStep;
            } else if (stopTime > currentStep*timeStep) {
              for (unsigned int i=1; i < nColumns; i++) {
                values[currentStep][i-1] += line_contents[i];
                if (do_errors) {
                  squares[currentStep][i-1] += line_squares[i];
                }
              }
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
    }
  }
  float time = 0.;
  std::ofstream ofs((outputBase+".sim.dat").c_str(), std::ios::out);
  
  if (ofs.is_open()) {
    for (unsigned int i = 0; i < values.size(); i++, time+=timeStep) {
      ofs << time;
      for (unsigned int j = 0; j < values[i].size() ; j++) {
        ofs << " " << values[i][j]/nFiles;
      }
      if (do_errors) {
        for (unsigned int j = 0; j < squares[i].size() ; j++) {
          ofs << " " << (squares[i][j] - values[i][j])/nFiles;
        }
      }
      ofs << std::endl;
    }
    ofs.close();
  } else {
    std::cerr << "Error writing to " << outputBase  << ".sim.dat" << std::endl;
    return 1;
  }

  // write file with averaged initial conditions
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
    std::cout << "Error writing to " << outputBase << ".init" << std::endl;
    return 1;
  }

  return 0;
}
