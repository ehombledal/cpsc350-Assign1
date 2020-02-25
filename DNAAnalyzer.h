#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>

using namespace std;

class DNAAnalyzer
{
  public:
          DNAAnalyzer(); //default constructor

          int getInfo(string fileName);
          void generateStatistics();
          void generateDNA();

          //i know ints should be private, but i am NOT doing getters and setters for all of these.

          //holder for the line of DNA
          string m_DNAString;
          string m_StringLengthHolder;

          //count for how many strings were read (used for stats)
          int m_stringCount;

          //Count for each nucleotide
          double m_ACount;
          double m_CCount;
          double m_TCount;
          double m_GCount;
          double m_totCount;

          //Bigrams starting with A
          double m_bigramAA;
          double m_bigramAC;
          double m_bigramAT;
          double m_bigramAG;

          //Bigrams starting with C
          double m_bigramCC;
          double m_bigramCA;
          double m_bigramCT;
          double m_bigramCG;

          //Bigrams starting with T
          double m_bigramTT;
          double m_bigramTA;
          double m_bigramTC;
          double m_bigramTG;

          //Bigrams starting with G
          double m_bigramGG;
          double m_bigramGA;
          double m_bigramGC;
          double m_bigramGT;

          //Bigram count
          double m_bigramTot;

          //Holders for statistics
          double m_sum;
          double m_mean;
          double m_variance;
          double m_stddev;

  private:

        void valueCounter(char dna);
        void valueCounter(char dna1, char dna2);
        double calculateVariance();
};
