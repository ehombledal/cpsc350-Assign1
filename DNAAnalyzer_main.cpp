#include "DNAAnalyzer.h"

int main(int argc, char **argv)
{
    DNAAnalyzer *analyzer = new DNAAnalyzer();

    
    analyzer-> getInfo(argv[1]); //only works off command line params.
    analyzer-> generateStatistics();
    analyzer-> generateDNA();

    delete analyzer;
    //DOESNT REPEAT BUT WE'LL FIX THAT
}
