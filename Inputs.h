#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>

class Inputs {
  public:
    std::vector<std::string> sources;
    std::string outputname;
    bool doMC;
    bool doWeight;
    int  isoCut;

    int argc;
    char **argv;
    std::vector<std::string> indices;

    Inputs(int argc, char **argv);
    void ShowUsage(std::string argv);
    int ParseOptions();
    void ShowOptions();
};

Inputs::Inputs(int _argc, char **_argv) {
    argc = _argc;
    argv = _argv;
    doMC = false;
    doWeight = false;
    isoCut = 0;
    indices.push_back("-i");
    indices.push_back("-o");
    indices.push_back("-m");
    indices.push_back("-w");
    indices.push_back("-c");
    indices.push_back("-h");
}

void Inputs::ShowOptions() {
  std::cout << "Input options are given:\n"
            << " doMC\t\t\t" << doMC << std::endl
            << " doWeight\t\t" << doWeight << std::endl
            << " isoCut\t\t\t" << isoCut << std::endl
            << " Output file\t\t" << outputname << std::endl;
  std::cout << " Input files\n";
  for (unsigned int i=0; i<sources.size(); i++) {
    std::cout << "   " << sources[i] << std::endl;
  }
  std::cout << std::endl;
}

void Inputs::ShowUsage(std::string argv)
{
   std::cerr << "Usage: " << argv << " -i input1.root input2.root -o output.root\n"
             << "Options:\n"
             << "  -i\tName of input files\n"
             << "  -o\tName of output file\n"
             << "  -w\tApply weight (Default: 0)\n"
             << "  -m\tIs this MC file? (Default: 0)\n"
             << "  -c\tIsolation cut type (Default: 13)\n"
             << "  -h\tShow this help message\n"
             << std::endl;
}

int Inputs::ParseOptions() {
  if (argc<5) {
    ShowUsage(argv[0]);
    return 1;
  }

  bool multipleInput = false;
  for (unsigned int i=1; i<argc-1; i++) { // skip the program name(argv[0]) and the last option
    std::string argu = argv[i];
    std::string nextArgu = argv[i+1];

    if (argu=="-h") {
      ShowUsage(argv[0]);
      return 1;

    } else if (argu=="-i" || multipleInput) {
      // Are multiple input files given?
      if (std::find(indices.begin(), indices.end(), nextArgu) != indices.end())
        multipleInput = false;  // Next input argument is a different input option
      else multipleInput = true;
      if (i+1<argc && multipleInput) { // Make sure that this is not the end of argv!
        sources.push_back(nextArgu);
      }

    } else if (argu=="-o") {
      if (i+1<argc) { // Make sure that this is not the end of argv!
        outputname = nextArgu;
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "-o option requires 1 argument." << std::endl;
        return 1;
      }
    } else if (argu=="-w") {
      if (i+1<argc) { // Make sure that this is not the end of argv!
        doWeight = atoi(nextArgu.c_str());
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "-w option requires 1 argument." << std::endl;
        return 1;
      }
    } else if (argu=="-m") {
      if (i+1<argc) { // Make sure that this is not the end of argv!
        doMC = atoi(nextArgu.c_str());
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "-m option requires 1 argument." << std::endl;
        return 1;
      }
    } else if (argu=="-c") {
      if (i+1<argc) { // Make sure that this is not the end of argv!
        isoCut = atoi(nextArgu.c_str());
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "-c option requires 1 argument." << std::endl;
        return 1;
      }
    }

  }

  return 0;
}


