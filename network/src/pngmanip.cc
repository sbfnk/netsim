#define NO_FREETYPE // don't include freetype fonts
#include <pngwriter.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
enum pixelTypes {Susceptible, Infected, Recovered};

int main(int argc, char* argv[])
{

  unsigned int verbose = 0;

  std::string inputFile = "";
  std::string outputFile = "";

  po::options_description command_line_options
    ("\nUsage: pngmanip [options]... \n\nMain options");

  command_line_options.add_options()
    ("help,h",
     "produce help message")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ("png-file,f",po::value<std::string>(),
     "png file to manipulate ")
    ("output,o",po::value<std::string>(),
     "file to write result to")
    ("some-colour,s",
     "turn into black and white with some colour")
    ;

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(command_line_options).run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }
  
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << command_line_options << std::endl;
    return 0;
  }

  if (vm.count("verbose")) {
    verbose = 1;
  }
  if (vm.count("very-verbose")) {
    verbose = 2;
  }

  pngwriter* reader;
  pngwriter* writer;
  if (vm.count("png-file")) {
    inputFile = vm["png-file"].as<std::string>();
    reader = new pngwriter(1,1,0,inputFile.c_str());
    try {
      reader->readfromfile(inputFile.c_str());
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
      return 1;
    }
  } else {
    std::cerr << "ERROR: must specify input file using --png-file" << std::endl;
    std::cout << command_line_options << std::endl;
    return 0;
  }

  if (vm.count("output")) {
    outputFile = vm["output"].as<std::string>();
    writer = new pngwriter(reader->getwidth()*4, reader->getheight()*4, 0, outputFile.c_str());
//     writer = new pngwriter(reader->getwidth()*7, reader->getheight()*7, 0, outputFile.c_str());
  } else {
    std::cerr << "ERROR: must specify output file using --file" << std::endl;
    std::cout << command_line_options << std::endl;
    return 0;
  }

  for (int i = 0; i < reader->getwidth(); ++i) {
    for (int j = 0; j < reader->getheight(); ++j) {
      std::vector<int> pixelColour, bgColour, fgColour;
      pixelColour.push_back(reader->read(i,j,1));
      pixelColour.push_back(reader->read(i,j,2));
      pixelColour.push_back(reader->read(i,j,3));
      if (pixelColour[2] > 0) {
        //Susceptible
        bgColour.push_back(1.25*(65535-pixelColour[2]));
        bgColour.push_back(1.25*(65535-pixelColour[2]));
        bgColour.push_back(1.25*(65535-pixelColour[2]));
        if (vm.count("some-colour")) {
          fgColour.push_back(0);
          fgColour.push_back(0);
          fgColour.push_back(65535);
        } else {
//           unsigned int brightness = pixelColour[2] > 40000 ? 0 : 65535;
          unsigned int brightness = 0;
          fgColour.push_back(brightness);
          fgColour.push_back(brightness);
          fgColour.push_back(brightness);
        }
        writer->filledsquare(i*4,j*4,i*4+3,j*4+3,bgColour[0], bgColour[1], bgColour[2]);
//         writer->filledsquare(i*7,j*7,i*7+6,j*7+6,bgColour[0], bgColour[1], bgColour[2]);
//         writer->circle(i*7+3, j*7+3, 2, fgColour[0], fgColour[1], fgColour[2]);
      } else if (pixelColour[0] > 0) {
        //Infected
        bgColour.push_back(65535);
        bgColour.push_back(0);
        bgColour.push_back(0);
//         bgColour.push_back(1.25*(65535-pixelColour[0]));
//         bgColour.push_back(1.25*(65535-pixelColour[0]));
//         bgColour.push_back(1.25*(65535-pixelColour[0]));
        if (vm.count("some-colour")) {
          fgColour.push_back(65535);
          fgColour.push_back(0);
          fgColour.push_back(0);
        } else {
//           unsigned int brightness = pixelColour[0] > 40000 ? 0 : 65535;
          unsigned int brightness = 0;
          fgColour.push_back(brightness);
          fgColour.push_back(brightness);
          fgColour.push_back(brightness);
        }
        //paint symbol
        writer->filledsquare(i*4,j*4,i*4+3,j*4+3,bgColour[0], bgColour[1], bgColour[2]);
//         writer->filledsquare(i*7,j*7,i*7+6,j*7+6,bgColour[0], bgColour[1], bgColour[2]);
//         writer->cross(i*7+3, j*7+3, 7, 7, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7, j*7, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+1, j*7+1, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+2, j*7+2, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+3, j*7+3, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+4, j*7+4, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+5, j*7+5, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+6, j*7+6, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7, j*7+6, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+1, j*7+5, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+2, j*7+4, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+4, j*7+2, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+5, j*7+1, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+6, j*7, fgColour[0], fgColour[1], fgColour[2]);
//         writer->filledcircle(i*7+3, j*7+3, 2, fgColour[0], fgColour[1], fgColour[2]);
      } else if (pixelColour[1] > 0) {
        //Recovered
//         bgColour.push_back(1.25*(65535-pixelColour[1]));
//         bgColour.push_back(1.25*(65535-pixelColour[1]));
//         bgColour.push_back(1.25*(65535-pixelColour[1]));
        bgColour.push_back(40000);
        bgColour.push_back(0);
        bgColour.push_back(0);
        if (vm.count("some-colour")) {
          fgColour.push_back(65535);
          fgColour.push_back(0);
          fgColour.push_back(0);
        } else {
//           unsigned int brightness = pixelColour[1] > 40000 ? 0 : 65535;
          unsigned int brightness = 0;
          fgColour.push_back(brightness);
          fgColour.push_back(brightness);
          fgColour.push_back(brightness);
        }
        //paint cross
        writer->filledsquare(i*4,j*4,i*4+3,j*4+3,bgColour[0], bgColour[1], bgColour[2]);
//         writer->filledsquare(i*7,j*7,i*7+6,j*7+6,bgColour[0], bgColour[1], bgColour[2]);
//         writer->plot(i*7, j*7, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+1, j*7+1, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+2, j*7+2, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+3, j*7+3, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+4, j*7+4, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+5, j*7+5, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+6, j*7+6, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7, j*7+6, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+1, j*7+5, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+2, j*7+4, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+4, j*7+2, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+5, j*7+1, fgColour[0], fgColour[1], fgColour[2]);
//         writer->plot(i*7+6, j*7, fgColour[0], fgColour[1], fgColour[2]);
//         writer->cross(i*7+3, j*7+3, 7, 7, fgColour[0], fgColour[1], fgColour[2]);
//          writer->filledcircle(i*7+3, j*7+3, 2, fgColour[0], fgColour[1], fgColour[2]);
      } 
    }
  }

  writer->close();

  return 0;
}

