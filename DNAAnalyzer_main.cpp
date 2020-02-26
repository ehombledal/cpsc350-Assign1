#include "DNAAnalyzer.h"

int main(int argc, char **argv)
{
  bool cont = 0;
  string newFile = argv[1];
  int ranBefore = 0;
  do
  {
    DNAAnalyzer *analyzer = new DNAAnalyzer();

    analyzer-> getInfo(newFile); //only works off command line params.
    analyzer-> generateStatistics(ranBefore);
    analyzer-> generateDNA();

    delete analyzer;
    cout << "Done with "<< newFile << "!" << endl; 

    cout << "Would you like to process another file? (true = 1/false = 0)" << endl;
    cin >> cont;
    ranBefore = 1;

    if (cont == true)
    {
      cout << "what is the name of the file you want to process?" << endl;
      cin >> newFile;
    }

  } while(cont == true);

    cout << "Thank you for using the DNA Analyzer!" << endl;
}
