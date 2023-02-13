#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <iomanip>
#include <algorithm>

using namespace std;

void extendSpace (const vector<int> &equationVector, const int stEq, vector<mpq_class> myCoeffs);
void toLatex (const vector<int> &equationVector, const int stEq, const vector<mpq_class> &myCoeffs);

int main ()
{
  const vector<int> vanderPol
  {
    1, 0,  0, 0,  1, 0,
    1, 0,  1, 0,  2, 0,
    1, 0,  0, 0,  0, 1,
    0, 1,  0, 0,  1, 0
  };

  vector<mpq_class> vanderPolCoeffs{
   mpq_class(1), mpq_class(-1,3), mpq_class(-1), mpq_class(1)
  };


  cout << "\n The space extension for van der Pol ODE:";
  cout << endl;
  extendSpace (vanderPol, 2, vanderPolCoeffs);
  cout << "--------------------------------\n\n";


  const vector<int> quarticAnharmonicOscillator
  {
    1, 0,  0, 0,  0, 1,
    0, 1,  0, 0,  1, 0,
    0, 1,  1, 0,  2, 0
  };

  vector<mpq_class> quarticAnharmonicOscillatorCoeffs{
   mpq_class(1), mpq_class(-1), mpq_class(-1)
  };   // mu is 1, k1 is 1, k2 is 1


  cout << "\n The space extension for classical quartic anharmonic oscillator ODE:";
  cout << endl;
  extendSpace (quarticAnharmonicOscillator, 2, quarticAnharmonicOscillatorCoeffs);
  cout << "--------------------------------\n\n";


  const vector<int> highPowers
  {
    1, 0,  50, 50,  50, 50,
    0, 1,   0,  0,   1,  0,
    0, 1,   1,  0,   2,  0
  };

  vector<mpq_class> highPowersCoeffs{
   mpq_class(1), mpq_class(1), mpq_class(1)
  };


  cout << "\n The space extension for ODE with high powers:";
  cout << endl;
  extendSpace (highPowers, 2, highPowersCoeffs);
  cout << "--------------------------------\n\n";


  const vector<int> henonHeiles
  {
    1, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,
    0, 1, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0,
    0, 1, 0, 0,  0, 0, 1, 0,  1, 0, 0, 0,
    0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 1,
    0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 1, 0,
    0, 0, 0, 1,  1, 0, 0, 0,  1, 0, 0, 0,
    0, 0, 0, 1,  0, 0, 1, 0,  0, 0, 1, 0
  };

  vector<mpq_class> henonHeilesCoeffs{
   mpq_class(1), mpq_class(-1), mpq_class(-2), mpq_class(1), mpq_class(-1), mpq_class(-1), mpq_class(1)
  }; // lambda is 1


  cout << "\n The space extension for Henon-Heiles ODE:";
  cout << endl;
  extendSpace (henonHeiles, 4, henonHeilesCoeffs);
  cout << "--------------------------------\n\n";


  const vector<int> rabinovichFabrikant
  {
    1, 0, 0,  0, 0, 1,  0, 1, 0,
    1, 0, 0,  0, 0, 0,  0, 1, 0,
    1, 0, 0,  1, 0, 0,  1, 1, 0,
    1, 0, 0,  0, 0, 0,  1, 0, 0,
    0, 1, 0,  0, 0, 1,  1, 0, 0,
    0, 1, 0,  0, 0, 0,  1, 0, 0,
    0, 1, 0,  1, 0, 0,  2, 0, 0,
    0, 1, 0,  0, 0, 0,  0, 1, 0,
    0, 0, 1,  0, 0, 0,  0, 0, 1,
    0, 0, 1,  0, 1, 0,  1, 0, 1
  };

  vector<mpq_class> rabinovichFabrikantCoeffs{
   mpq_class(1), mpq_class(-1), mpq_class(1), mpq_class(1), mpq_class(3), mpq_class(1), mpq_class(-1), mpq_class(1), mpq_class(-2), mpq_class(-2)
  }; // gamma is 1, alpha is 1

  cout << "\n The space extension for Rabinovich-Fabrikant ODE:";
  cout << endl;
  extendSpace (rabinovichFabrikant, 3, rabinovichFabrikantCoeffs);
  cout << "--------------------------------\n\n";


  return 0;
}


void extendSpace (const vector<int> &equationVector, const int stEq, vector<mpq_class> myCoeffs)
{

  vector<int> runVec = equationVector;

  bool appears = 0;

  int i = 0;
  int numEqs = stEq;

  vector<int> rightHandSide (stEq, 0);

  while (i < runVec.size ())
  {
    for (int j = stEq; j <= 2 * stEq; j += stEq)
    {
      for (int k = 0; k < stEq; k++)
      {
        rightHandSide[k] = runVec[i + j + k];
      }
      for (int i = 0; i < runVec.size (); i += 3 * stEq)
      {
        vector<int> subvector(stEq, 0);
        subvector = {runVec.begin () + i, runVec.begin () + i + stEq};
        if (subvector == rightHandSide)
        {
          appears = 1;
          break;
        }
      }

      if (appears == 0)
      {
        numEqs++;
        vector<int> allZero(stEq,0);
        vector<int> zeroEntry(3*stEq, 0);

        if (rightHandSide == allZero)
        {
          runVec.insert (runVec.end (), zeroEntry.begin(), zeroEntry.end() );
          myCoeffs.push_back(mpq_class(0));
        }
        else
        {
          for (int i = 0; i < equationVector.size (); i += 3*stEq)
          {
            vector<int>::const_iterator it;
            it = find (equationVector.begin () + i, equationVector.begin () + i + 3*stEq, 1);
            int oneInd = (it - equationVector.begin ());

            bool x = 0;
            vector<int> myTemp(stEq, 0);

            for (int j = 0; j < stEq; j++)
            {
              if (oneInd == i+j && rightHandSide[j] > 0)
              {
                x = 1;
                for (int k = 0 ; k < stEq; k++)
                {
                  myTemp[k] =
                    equationVector[i + stEq + k] + equationVector[i + stEq + stEq + k] +
                    rightHandSide[k];
                  if (j == k)
                  {
                    myCoeffs.push_back(rightHandSide[j] * myCoeffs[i/(3*stEq)]);
                    myTemp[k] -= 1;
                  }
                }
              }
            }

            if (x == 0)
            {
              continue;
            }
            vector<int> lefts;
            int oneCount = 0;

            for (int counter = 0; counter < stEq; counter++)
              runVec.push_back (rightHandSide[counter]);

            for (int i = 0; i < myTemp.size (); i++)
            {
              if (myTemp[i] == 1)
              {
                oneCount++;
              }
              if (myTemp[i] != 1 || (myTemp[i]==1 && (oneCount%2==1)))
              {
                lefts.push_back(myTemp[i] / 2);
                runVec.push_back(myTemp[i] / 2);
              }
              else if (myTemp[i]==1 && oneCount%2==0)
              {
                lefts.push_back(1);
                runVec.push_back(1);
              }
            }
            for (int i = 0; i < lefts.size (); i++)
            {
              runVec.push_back(myTemp[i] - lefts[i]);
            }
          }
        }
      }
      appears = 0;
    }

    i += 3*stEq;
  }

  cout << "\nThe original ODE set is " << endl;
  toLatex (equationVector, stEq, myCoeffs);

  cout << "\nThe new ODE set is " << endl;
  toLatex (runVec, stEq, myCoeffs);

  cout << endl;

  return;
}

void toLatex (const vector<int> &equationVector, const int stEq, const vector<mpq_class> &myCoeffs)
{
  cout << endl;
  cout << "\\begin{eqnarray}" << endl;

  vector<int> prevLeftHandSide(stEq);
  vector<int> leftHandSide(stEq);
  for (auto i = 0; i < equationVector.size(); i+=3*stEq)
  {
    bool same = 0;
    leftHandSide.assign(equationVector.begin()+i, equationVector.begin()+i+stEq-1); 
    if ( i != 0 && leftHandSide == prevLeftHandSide )
    {
      same = 1;
    }


    if (!same)
    {
      if ( i != 0 )
      {
         cout << "\\\\" << endl;
      }
      cout << "  u^{(";
      for (auto j = 0; j < stEq; j++)
      {
        if (j != 0)
        {
          cout << ",";
        }
        cout << equationVector[i+j];
      }
      cout << ")} &=& ";
   }


    mpq_class myZero(0);
    bool isPos = 1;

    if ( myCoeffs[i/(3*stEq)] < myZero )
    {
       isPos = 0;
    }


    mpq_class coeffToPrint = myCoeffs[i/(3*stEq)];
    if ( same )
    {
      if ( isPos )
      {
        cout << " \\nonumber \\\\" << endl << " &+& ";
      }
      else
      {
        cout << " \\nonumber \\\\" << endl << " &-& ";
        coeffToPrint = -coeffToPrint;
      }
    }
    cout << coeffToPrint << "\\,";

    cout << " u^{(";
    for (auto k = i + stEq; k < i + stEq + stEq; k++)
    {
      if (k != i + stEq)
      {
        cout << ",";
      }
      cout << equationVector[k];
    }
    cout << ")}";

    cout << " u^{(";
    for (auto k = i + stEq + stEq; k < i + stEq + stEq + stEq; k++)
    {
      if (k != i + stEq + stEq)
      {
        cout << ",";
      }
      cout << equationVector[k];
    }
    cout << ")}";

    prevLeftHandSide = leftHandSide;
  }
  cout << endl << "\\end{eqnarray}" << endl;
}


// END OF PROGRAM
