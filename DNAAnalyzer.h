#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class DNAAnalyzer
{
  public:
          DNAAnalyzer(); //default constructor

          void initializeFileStream(string fileName);
          void getInfo(string fileName);
          void calculateStatistics();
          void calculateProbabilty();
          void printOutput();
          void generateDNA();

          //i know ints should be private, but i am NOT doing getters and setters for all of these.

          //holder for the line of DNA
          string m_DNAString;

          //Count for each nucleotide
          int m_ACount;
          int m_CCount;
          int m_TCount;
          int m_GCount;
          int m_totCount;

          //Bigrams starting with A
          int m_bigramAC;
          int m_bigramAT;
          int m_bigramAG;

          //Bigrams starting with C
          int m_bigramCA;
          int m_bigramCT;
          int m_bigramCG;

          //Bigrams starting with T
          int m_bigramTA;
          int m_bigramTC;
          int m_bigramTG;

          //Bigrams starting with G
          int m_bigramGA;
          int m_bigramGC;
          int m_bigramGT;

          //Holders for statistics
          int    m_sum;
          double m_mean;
          double m_variance;
          double m_stddev;

  private:
};
