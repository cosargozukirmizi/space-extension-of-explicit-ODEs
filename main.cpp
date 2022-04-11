#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

void extendSpace (const vector<int> equationVector, const int stEq);

int main ()
{
  const vector<int> vanderPol
  {
    1, 0,  0, 0,  1, 0,
    1, 0,  1, 0,  2, 0,
    1, 0,  0, 0,  0, 1,
    0, 1,  0, 0,  1, 0
  };

  cout << "\n The space extension for van der Pol ODE is";
  cout << endl;
  extendSpace (vanderPol, 2);
  cout << "--------------------------------\n\n";

  const vector<int> quarticAnharmonicOscillator
  {
    1, 0,  0, 0,  0, 1,
    0, 1,  0, 0,  1, 0,
    0, 1,  1, 0,  2, 0
  };

  cout << "\n The space extension for classical quartic anharmonic oscillator ODE is";
  cout << endl;
  extendSpace (quarticAnharmonicOscillator, 2);
  cout << "--------------------------------\n\n";

  const vector<int> highPowers
  {
    1, 0,  50, 50,  50, 50,
    0, 1,   0,  0,   1,  0,
    0, 1,   1,  0,   2,  0
  };

  cout << "\n The space extension for ODE with high powers is";
  cout << endl;
  extendSpace (highPowers, 2);
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

  cout << "\n The space extension for Henon-Heiles ODE is";
  cout << endl;
  extendSpace (henonHeiles, 4);
  cout << "--------------------------------\n\n";

  const vector<int> rabinovichFabrikant
  {
    1, 0, 0,  0, 0, 1,  0, 1, 0,
    1, 0, 0,  0, 0, 0,  0, 1, 0,
    1, 0, 0,  1, 0, 0,  1, 1, 0,
    1, 0, 0,  0, 0, 0,  1, 0, 0,
    0, 1, 0,  0, 1, 0,  1, 0, 1,
    0, 1, 0,  0, 0, 0,  1, 0, 0,
    0, 1, 0,  1, 0, 0,  2, 0, 0,
    0, 1, 0,  0, 0, 0,  0, 1, 0,
    0, 0, 1,  0, 0, 0,  0, 0, 1,
    0, 0, 1,  0, 0, 1,  1, 1, 0
  };

  cout << "\n The space extension for Rabinovich-Fabrikant ODE is";
  cout << endl;
  extendSpace (rabinovichFabrikant, 3);
  cout << "--------------------------------\n\n";

  return 0;
}


void extendSpace (const vector<int> equationVector, const int stEq)
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
                oneCount++;

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

  for (auto i = 0; i < runVec.size (); i++)
  {
    cout << setw(3);
    cout << runVec[i];

    if (i % stEq == stEq-1)
      cout << ' ';

    if (i % (3*stEq) == (3*stEq-1))
      cout << endl;
  }

  cout << " There are " << numEqs << " equations.\n";

  return;
}

// END OF PROGRAM
