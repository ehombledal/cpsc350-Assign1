
#include "DNAAnalyzer.h"

DNAAnalyzer::DNAAnalyzer() //god is dead.
{
  m_DNAString = "null";
  m_StringLengthHolder;

  m_ACount   = 0;
  m_CCount   = 0;
  m_TCount   = 0;
  m_GCount   = 0;
  m_totCount = 0;

  m_bigramAA = 0;
  m_bigramAC = 0;
  m_bigramAG = 0;
  m_bigramAT = 0;

  m_bigramCA = 0;
  m_bigramCC = 0;
  m_bigramCG = 0;
  m_bigramCT = 0;

  m_bigramTA = 0;
  m_bigramTC = 0;
  m_bigramTT = 0;
  m_bigramTG = 0;

  m_bigramTA = 0;
  m_bigramTC = 0;
  m_bigramTG = 0;
  m_bigramGG = 0;

  m_bigramTot = 0;

  m_sum      = 0;
  m_mean     = 0;
  m_variance = 0;
  m_stddev   = 0;
}


void DNAAnalyzer::valueCounter(char dna)
{
  char toCheck = toupper(dna);
  if (toCheck == 'A')
  {
    m_ACount++;
    m_totCount++;
  } else if (toCheck == 'C')
  {
    m_CCount++;
    m_totCount++;
  } else if (toCheck == 'T')
  {
    m_TCount++;
    m_totCount++;
  } else if (toCheck == 'G')
  {
    m_GCount++;
    m_totCount++;
  } else {
    cout << toCheck << " is not a valid letter! Skipping..." << endl;
  }
}

void DNAAnalyzer::valueCounter(char dna1, char dna2)
{
  char first = toUpper(dna1);
  char second = toUpper(dna2);

  if (first == 'A') //ALL BIGRAMS THAT START WITH A
  {
    if (second == 'A')
    {
      m_bigramAA++;
      m_bigramTot++;
    }
    else if (second == 'C')
    {
      m_bigramAC++;
      m_bigramTot++;
    } else if (second == 'T')
    {
      m_bigramAT++;
      m_bigramTot++;
    } else if (second == 'G')
    {
      m_bigramAG++;
      m_bigramTot++;
    } else {
      cout << first << " " << second << " is not a valid bigram combo! Skipping..." << endl;
    }

  } else if (first == 'C') //ALL BIGRAMS THAT START WITH C
  {
    if (second == 'A')
    {
      m_bigramCA++;
      m_bigramTot++;
    }  else if (second == 'C')
    {
      m_bigramCC++;
      m_bigramTot++;
    } else if (second == 'T')
    {
      m_bigramCT++;
      m_bigramTot++;
    } else if (second == 'G')
    {
      m_bigramCG++;
      m_bigramTot++;
    } else {
      cout << first << " " << second << " is not a valid bigram combo! Skipping..." << endl;
    }

  } else if (first == 'T') //ALL BIGRAMS THAT START WITH T
  {
    if (second == 'A')
    {
      m_bigramTA++;
      m_bigramTot++;
    } else if (second == 'C')
    {
      m_bigramTC++;
      m_bigramTot++;
    } else if (second == 'T')
    {
      m_bigramTT++;
      m_bigramTot++;
    }else if (second == 'G')
    {
      m_bigramTG++;
      m_bigramTot++;
    } else {
      cout << first << " " << second << " is not a valid bigram combo! Skipping..." << endl;
    }
  } else if (first == 'G') //ALL BIGRAMS THAT START WITH G
  {
    if (second == 'A')
    {
      m_bigramGA++;
      m_bigramTot++;
    } else if (second == 'C')
    {
      m_bigramGC++;
      m_bigramTot++;
    } else if (second == 'T')
    {
      m_bigramGT++;
      m_bigramTot++;
    } else if (second == 'G')
    {
      m_bigramGG++;
      m_bigramTot++;
    }else {
      cout << first << " " << second << " is not a valid bigram combo! Skipping..." << endl;
    }
  } else {
    cout << first << " is not a valid letter! Skipping..." << endl;
  }
}

int DNAAnalyzer::getInfo(string fileName)
{
  ifstream inFS;     // Input file stream
  string DNAString;


  // Try to open file
  cout << "Opening file " << fileName << endl;
  inFS.open(fileName);

  if (!inFS.is_open())
  {
     cout << "Could not open file " <<fileName << endl;
     return 1; // 1 indicates error
  }

  //Print read file data to output
  cout <<"Reading data" << endl;

  while (!inFS.eof()) //goes until end of file
   {
     inFS >> DNAString;
     m_StringLengthHolder += to_string(DNAString.size()) + ",";

     if (!inFS.fail())
      {
        if (DNAString.size() % 2 == 0) //string is even, can count by 2.
        {
          for (int i = 0; i < DNAString.size(); i+=2)
          {
            valueCounter(DNAString[i]);
            valueCounter(DNAString[i+1]);
            valueCounter(DNAString[i], DNAString[i+1]);
          }
        } else { //STRING IS ODD, HAVE TO COUNT LAST INDEX SEPERATE
          for (int i = 0; i <DNAString.size() - 1; i+=2)
          {
            valueCounter(DNAString[i]);
            valueCounter(DNAString[i+1]);
            valueCounter(DNAString[i], DNAString[i+1]);
          }
           valueCounter(DNAString[DNAString.size() - 1]); //counts end value. Not part of a pair.
        }
        m_stringCount++;
      }
   }

   cout << "Closing file" << fileName << endl;
   inFS.close(); // Done with file, so close it
   return 0;
}


double DNAAnalyzer::calculateVariance() //cant use arrays, but strings are okay, so...
{
  double variance = 0;

   for (int n = 0; n < m_stringCount; ++n)
   {
     int found = m_StringLengthHolder.find(","); //stored values in string, split by comma. This just grabs them back out.
     string stringLength = m_StringLengthHolder.substr(0,found);

     m_StringLengthHolder = m_StringLengthHolder.substr(found + 1); //removes number just found by substring.
     int length = stoi(stringLength); //converts to int

     variance += (length - m_mean) * (length - m_mean); //actual math
   }

   variance = variance / m_stringCount;
   return variance;
}


void DNAAnalyzer::generateStatistics()
{
  ofstream outFS;    //output file stream

  //generate statistics
  outFS.open("Hombledal.txt");
  outFS << "Erik Hombledal, 2325523" << endl;
  outFS << endl;
  outFS << "STATISTICS" << endl;
  outFS<< endl;

  outFS << "The sum of the lengths of the DNA strings is: " << m_totCount << endl;

  m_mean = m_totCount / m_stringCount;
  outFS << "The mean of the DNA strings is: " << m_mean << endl;

  m_variance = calculateVariance();
  outFS << "The variance of the DNA strings is: " << m_variance << endl;

  m_stddev = sqrt(m_variance);
  outFS << "The standard deviation of the DNA strings is: " << m_stddev << endl;
  outFS << endl;


  //generate probabilities
  outFS << "PROBABILITIES" << endl;
  outFS << endl;

  double apercent = m_ACount / m_totCount;
  double tpercent = m_TCount / m_totCount;
  double cpercent = m_CCount / m_totCount;
  double gpercent = m_GCount / m_totCount;

  outFS << "Nucleotide Probabilities --  A: " << (apercent * 100)
                                  << "%  T: " << (tpercent * 100)
                                  << "%  C: " << (cpercent * 100)
                                  << "%  G: " << (gpercent * 100)
                                  << "%" << endl;

  outFS << "Bigram Probabilities - A" << endl;
  outFS << "AA: " << ((m_bigramAA / m_bigramTot)) << endl;
  outFS << "AT: " << ((m_bigramAT / m_bigramTot)) << endl;
  outFS << "AC: " << ((m_bigramAC / m_bigramTot)) << endl;
  outFS << "AG: " << ((m_bigramAG / m_bigramTot)) << endl;

  outFS << "Bigram Probabilities - T" << endl;
  outFS << "TA: " << ((m_bigramTA / m_bigramTot)) << endl;
  outFS << "TT: " << ((m_bigramTT / m_bigramTot)) << endl;
  outFS << "TC: " << ((m_bigramTC / m_bigramTot)) << endl;
  outFS << "TG: " << ((m_bigramTG / m_bigramTot)) << endl;

  outFS << "Bigram Probabilities - C" << endl;
  outFS << "CA: " << ((m_bigramCA / m_bigramTot)) << endl;
  outFS << "CT: " << ((m_bigramCT / m_bigramTot)) << endl;
  outFS << "CC: " << ((m_bigramCC / m_bigramTot)) << endl;
  outFS << "CG: " << ((m_bigramCG / m_bigramTot)) << endl;

  outFS << "Bigram Probabilities - G" << endl;
  outFS << "GA: " << ((m_bigramAT / m_bigramTot)) << endl;
  outFS << "GT: " << ((m_bigramAC / m_bigramTot)) << endl;
  outFS << "GC: " << ((m_bigramAG / m_bigramTot)) << endl;
  outFS << "GG: " << ((m_bigramGG / m_bigramTot)) << endl;


  outFS.close();
}

void DNAAnalyzer::generateDNA()
{
  //TODO: this.
}
