
#include "DNAAnalyzer.h"

DNAAnalyzer::DNAAnalyzer() //god is dead.
{
  DNAString = "null";

  m_ACount   = 0;
  m_CCount   = 0;
  m_TCount   = 0;
  m_GCount   = 0;
  m_totCount = 0;

  m_bigramAC = 0;
  m_bigramAG = 0;
  m_bigramAT = 0;

  m_bigramCA = 0;
  m_bigramCG = 0;
  m_bigramCT = 0;

  m_bigramTA = 0;
  m_bigramTC = 0;
  m_bigramTG = 0;

  m_bigramTA = 0;
  m_bigramTC = 0;
  m_bigramTG = 0;

  m_sum      = 0;
  m_mean     = 0;
  m_variance = 0;
  m_stddev   = 0;
}


DNAAnalyzer::getInfo(string fileName)
{
  ifstream inFS;     // Input file stream
  string DNAString;    // Data value from file. It is more convenient to have input stream as string, and convert it


  if (argc < 2)
  {
    cout << "invalid command line params" << endl;
    return 1; // 1 indicates error
  }

  // Try to open file
  cout << "Opening file " << fName << endl;
  inFS.open(fileName); //can also have this where numFile is a String, and you get the input from command line.

  if (!inFS.is_open())
  {
     cout << "Could not open file " <<fileName << endl;
     return 1; // 1 indicates error
  }


  //Print read file data to output
  cout <<"Reading and printing data" << endl;

  while (!inFS.eof()) //goes until end of file
   {
    inFS >> DNAString;
     if (!inFS.fail())
      {
        //PUT VALUES INTO VARIABLES HERE 
      }
   }
}

DNAAnalyzer::calculateStatistics()
{

}
