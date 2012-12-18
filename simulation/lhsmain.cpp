#include "lhs.h"
#include "params.h"
#include "xmlhandler.h"
#include <boost/program_options.hpp>
Rand g_Rand(1337);
namespace po = boost::program_options;

void printUsage(char* pArgv0, po::options_description& desc)
{
  std::cout << "Usage: " << pArgv0 << " [options]\n" << desc << std::endl;
}

void printVersion()
{
  //std::cout << "Version: " << GR_VERSION << std::endl;
}

int main(int argc, char** argv) {
  size_t seed = time(NULL);
  size_t nSamples = 0;
  std::string inputFileName;
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help,h", "Help message")
  ("input-file,i", po::value(&inputFileName), "Input file name")
  ("seed,s", po::value(&seed), "Seed")
  ("samples,n", po::value(&nSamples), "Number of samples")
  ("version,v", "Version number")
  ("log-scale", "Use log scale intervals for double-type parameters");
  po::variables_map vm;
  try
  {
    po::positional_options_description p;
    p.add("input-file", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);
  }
  catch(std::exception& e)
  {
    std::cerr << "Error processing arguments: " << e.what() << std::endl;
    printUsage(argv[0], desc);
    return 1;
  }
  if(vm.count("version"))
  {
    printVersion();
    return 0;
  }
  if(vm.count("help"))
  {
    printUsage(argv[0], desc);
    return 0;
  }
  if(nSamples == 0)
  {
    std::cerr<<"Required number of samples" << std::endl;
    printUsage(argv[0], desc);
    return 1;
  }
  if(inputFileName.empty())
  {
    std::cerr<<"Required input file" << std::endl;
    printUsage(argv[0], desc);
    return 1;
  }
  g_Rand.setSeed(seed);
  boost::property_tree::ptree pt;
  std::auto_ptr<ParamFileHandler> handler(new XMLHandler("GR"));
  std::ifstream _if(inputFileName.c_str());
  handler->read(_if, pt);
  assert(_if.good());
  _if.close();
  Params p;
  LHS<Params> lhs(nSamples, pt, p, handler.get(), vm.count("log-scale"));
  lhs.performLHS();
  return 0;
}
