#include <iostream>
#include <fstream>
#include <string>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

std::vector<std::string> split(const std::string& str,
                          const std::string delimiters)
{
  std::vector<std::string> tokens; 

  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos) {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos     = str.find_first_of(delimiters, lastPos);
  }

  return tokens;
}

std::vector<double> split_to_double(const std::string &str,
                                    const std::string delimiters ) {
  
  std::vector<std::string> tokens = split(str, delimiters);
  std::vector<double> values;
  std::vector<std::string>::iterator it;
  double f;

  for (it = tokens.begin(); it != tokens.end(); it++) {
    std::stringstream ss(*it);
    ss >> f;
    values.push_back(f);
  }
  return values;
}


int main(int argc, char* argv[])
{

  double timeStep;
  double stopTime;
  
  std::vector<std::string> inputFiles;
  std::string outputFile;
  
  po::options_description command_line_options
    ("Usage: average_runs -d arg -s arg -o arg [input-files]\n\nAllowed options");
  command_line_options.add_options()
    ("help", "produce help message")    
    ("time-step,d", po::value<double>(),
     "discretized time step")    
    ("stop,s", po::value<double>(),
     "stopping time")    
    ("output-file,o", po::value<std::string>(),
     "output file")    
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
    timeStep = vm["time-step"].as<double>();
  } else {
    std::cerr << "ERROR: no time-step specified" << std::endl;
    std::cerr << command_line_options << std::endl;
    return 1;
  }
  
  if (vm.count("stop")) {
    stopTime = vm["stop"].as<double>();
  } else {
    std::cerr << "ERROR: no stop time specified" << std::endl;
    std::cerr << command_line_options << std::endl;
    return 1;
  }
  
  if (vm.count("output-file")) {
    outputFile = vm["output-file"].as<std::string>();
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

  double nSteps = (stopTime / timeStep);
  
  std::vector< std::vector<double> > values;

  unsigned int nFiles = inputFiles.size();
  unsigned int nColumns = 0;

  for (std::vector<std::string>::iterator it = inputFiles.begin();
       it != inputFiles.end(); it++) {
    std::ifstream ifs((*it).c_str());
    if (ifs.is_open()) {
      std::string line;
      unsigned int lineCount = 0;
      unsigned int nextStep = 0;
      
      while (!ifs.eof() && !(nextStep*timeStep > stopTime)) {
        lineCount++;
        std::getline(ifs, line);
        
        std::vector<double> line_contents = split_to_double(line, " \t");
        
        if (line_contents.size() > 0) {
          if (nColumns == 0) {
            // set number of columns
            nColumns = line_contents.size()-1;
            // fill vector with zeros
            std::vector<double> v(nColumns,0.);
            values.resize(static_cast<unsigned int>(nSteps)+1, v);
          } else if (nColumns+1 != line_contents.size()) {
            std::cerr << "Warning: wrong number of columns in line "
                      << lineCount << ": " << line_contents.size()
                      << " instead of " << nColumns << std::endl;
          }
          std::vector<double>::iterator lineIt = line_contents.begin();
          for (double time = *lineIt;
               !(time < nextStep*timeStep) && !(nextStep*timeStep > stopTime);
               ++nextStep) {
            ++lineIt;
            for (unsigned int i=0; lineIt != line_contents.end();
                 lineIt++, i++) {
              values[nextStep][i] += *lineIt;
            }
          } 
        }
      }
      ifs.close();
    } else {
      std::cout << "Error reading " << *it << "." << std::endl;
    }
  }
  double time = 0.;
  std::ofstream ofs(outputFile.c_str(), std::ios::out);

  if (ofs.is_open()) {
    for (std::vector< std::vector<double> >::iterator it = values.begin();
         it != values.end(); it++, time+=timeStep) {
      ofs << time;
      for (std::vector<double>::iterator dit = (*it).begin();
           dit != (*it).end(); dit++) {
        ofs << " " << (*dit)/nFiles;
      }
      ofs << std::endl;
    }
  } else {
    std::cout << "Error writing to " << outputFile << "." << std::endl;
    return 1;
  }

  return 0;
}
