#include <iostream>
#include <vector>
#include <map>
#include <gmpxx.h>
#include <iomanip>
#include <algorithm>
#include <cassert>

using namespace std;

void extendSpace (const vector<int> &equationVector, vector<mpq_class> myCoeffs, const vector<mpq_class> initVal);
void toLatex (const vector<int> &equationVector, const int stEq, const vector<mpq_class> &myCoeffs);
void printAsTripleF(const vector<int> rowInd, const vector<int> colInd, const vector<mpq_class> myValue);
void constructF(vector<int>& rowInd, vector<int>& colInd, vector<mpq_class>& myValue, const vector<int> runVec, map<vector<int>, int> myMap, const vector<mpq_class> myCoeffs, const int stEq);
void condensedKroneckerProduct(vector<mpq_class>& result, const vector<mpq_class> a, const vector<mpq_class> b);
void sparseMatTimesVec(vector<mpq_class>& result, const vector<int> rowInd, const vector<int> colInd, const vector<mpq_class> myCoeffs, const vector<mpq_class> x, const int resSize);
void constructAugmentedInitVal(vector<mpq_class>& runInitVal, const vector<mpq_class> initVal, map<vector<int>, int> myMap);

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

  vector<mpq_class> vanderPolInitVal{
   mpq_class(1,2), mpq_class(1,2)
  };

  cout << "\n The space extension for van der Pol ODE:";
  cout << endl;
  extendSpace (vanderPol, vanderPolCoeffs, vanderPolInitVal);
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

  vector<mpq_class> quarticAnharmonicOscillatorInitVal{
   mpq_class(1,2),mpq_class(1,2)
  };

  cout << "\n The space extension for classical quartic anharmonic oscillator ODE:";
  cout << endl;
  extendSpace (quarticAnharmonicOscillator, quarticAnharmonicOscillatorCoeffs, quarticAnharmonicOscillatorInitVal);
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

  vector<mpq_class> henonHeilesInitVal{
   mpq_class(1,2), mpq_class(1,2), mpq_class(1,2), mpq_class(1,2)
  };

  cout << "\n The space extension for Henon-Heiles ODE:";
  cout << endl;
  extendSpace (henonHeiles, henonHeilesCoeffs, henonHeilesInitVal);
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

  vector<mpq_class> rabinovichFabrikantInitVal{
   mpq_class(1,2), mpq_class(1,2), mpq_class(1,2)
  };

  cout << "\n The space extension for Rabinovich-Fabrikant ODE:";
  cout << endl;
  extendSpace (rabinovichFabrikant, rabinovichFabrikantCoeffs, rabinovichFabrikantInitVal);
  cout << "--------------------------------\n\n";

  return 0;
}


void extendSpace (const vector<int> &equationVector, vector<mpq_class> myCoeffs, const vector<mpq_class> initVal)
{

  const int stEq = initVal.size();
  vector<int> runVec = equationVector;
  map<vector<int>, int> myMap;

  bool appears = 0;

  int i = 0;
  int numEqs = stEq;

  vector<int> rightHandSide (stEq, 0);
  vector<int> tempLeft (stEq, 0);

  for (auto i = 0; i < runVec.size (); i+=3*stEq)
  {
    tempLeft.assign(equationVector.begin()+i, equationVector.begin()+i+stEq);
    myMap.insert({tempLeft, 0});
  }

  while (i < runVec.size ())
  {
    for (int j = stEq; j <= 2 * stEq; j += stEq)
    {
      for (int k = 0; k < stEq; k++)
      {
        rightHandSide[k] = runVec[i + j + k];
      }
      if (myMap.find(rightHandSide) == myMap.end())
      {
         myMap.insert({rightHandSide, 0});
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

  cout << "\nThe new ODE set is " << endl;
  toLatex (runVec, stEq, myCoeffs);

  cout << endl;

  map<vector<int>, int> ::iterator myMapIterator;
  vector<int> tempVec;


  int myDummyVariable = 0;

  for(myMapIterator=myMap.begin();myMapIterator !=myMap.end();++myMapIterator)
  {
    myMapIterator->second = myDummyVariable++;
  }

cout << "\nThe map for stacking the equations is " << endl;

  for(myMapIterator=myMap.begin();myMapIterator !=myMap.end();++myMapIterator)
  {
    tempVec = myMapIterator->first;
    cout << "[ ";

    for (auto dummy : tempVec)
      cout << dummy << " ";

    cout << "]:";
    cout << myMapIterator->second << ", " << endl;
  }

  vector<mpq_class> runInitVal(myMap.size());

  cout << "\n" << "Now we will try to construct initial vector." << endl;
  cout << "It is given as index value pairs below. " << endl;

  constructAugmentedInitVal(runInitVal, initVal, myMap);


  for(myMapIterator=myMap.begin();myMapIterator !=myMap.end();++myMapIterator)
  {
    cout << myMapIterator->second << " " << runInitVal[myMapIterator->second] << endl;
  }

  cout << "\n" << "Now we will try to construct F: " << endl;

  vector<int> rowInd;
  vector<int> colInd;
  vector<mpq_class> myValue;

  constructF(rowInd, colInd, myValue, runVec, myMap, myCoeffs, stEq);

  vector<mpq_class> resOfKronProd(runInitVal.size()*(runInitVal.size()+1)/2, 0);
  vector<mpq_class> resOfMatVecProd(runInitVal.size(), 0);
  vector<mpq_class> resOfKronProdAccumulator(runInitVal.size()*(runInitVal.size()+1)/2, 0);

  const int max_iter = 5;

  vector<vector<mpq_class>> rho;

  rho.push_back(runInitVal);

  vector<mpq_class> solutionIngredient(runInitVal.size(), 0);
  vector<mpq_class> solution(runInitVal.size(), 0);

  const mpq_class aHalf(1,2);

  int myconstant;

  for(int j=0; j < max_iter; j++)
  {
    if( j%2 == 0)
    {
       for(int k=0; k<(j/2-1); k++)
       {
         condensedKroneckerProduct(resOfKronProd, rho[k], rho[j-k]);
         std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());
       }
       condensedKroneckerProduct(resOfKronProd, rho[j/2], rho[j/2]);
       std::transform(resOfKronProd.begin(), resOfKronProd.end(), resOfKronProd.begin(), [&aHalf](auto& c){return c*aHalf;});
       std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());

    }
    else
    {
       for(int k=0; k<(j-1)/2; k++)
       {
         condensedKroneckerProduct(resOfKronProd, rho[k], rho[j-k]);
         std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());
       }
    }

    sparseMatTimesVec(resOfMatVecProd, rowInd, colInd, myValue, resOfKronProdAccumulator, runInitVal.size());
    myconstant = j + 1;
    std::transform(resOfMatVecProd.begin(), resOfMatVecProd.end(), resOfMatVecProd.begin(), [&myconstant](auto& c){return c*myconstant;});
    rho.push_back(resOfMatVecProd);
  }


  const mpq_class t(1, 10);    // time value under consideration
  mpq_class leftPart(1);


  for(int i = 0; i < max_iter; i++)
  {
    solutionIngredient = rho[i];
    std::transform(solutionIngredient.begin(), solutionIngredient.end(), solutionIngredient.begin(), [&leftPart](auto& c){return c*leftPart;});
    std::transform (solution.begin(), solution.end(), solutionIngredient.begin(), solution.begin(), std::plus<mpq_class>());
    leftPart = leftPart * t;
  }


  cout << "\nThe rho vectors are\n";
  for(int i = 0; i < max_iter; i++)
  {
    cout << "\nrho[" << i << "] :\n";
    for (auto element : rho[i])
    {
       cout << element << endl;
    }
  }

  cout << "\nThe solution of the initial value problem at t=" << t <<  " is\n";

  for(myMapIterator=myMap.begin();myMapIterator !=myMap.end();++myMapIterator)
  {

    tempVec = myMapIterator->first;
    cout << "[ ";

    for (auto dummy : tempVec)
      cout << dummy << " ";

    cout << "]: ";
    cout << solution[myMapIterator->second] << endl;
  }

  cout << endl; 

  return;
}

void toLatex (const vector<int> &equationVector, const int stEq, const vector<mpq_class> &myCoeffs)
{
  cout << endl;
  cout << "\\begin{eqnarray}" << endl;

  vector<int> prevLeftHandSide(stEq);
  vector<int> leftHandSide(stEq);
  vector<int> allZero(stEq,0);

  for (auto i = 0; i < equationVector.size(); i+=3*stEq)
  {
    bool same = 0;
    leftHandSide.assign(equationVector.begin()+i, equationVector.begin()+i+stEq); 
    if ( i != 0 && leftHandSide == prevLeftHandSide )
    {
      same = 1;
    }

    if( !same && (leftHandSide == allZero) )
    {
      if ( i != 0 )
      {
         cout << "\\\\" << endl;
      }
      cout << "  \\dot{u}^{(";
      for (auto j = 0; j < stEq; j++)
      {
        if (j != 0)
        {
          cout << ",";
        }
        cout << equationVector[i+j];
      }
      cout << ")} &=& 0";

      continue;
    }


    if (!same)
    {
      if ( i != 0 )
      {
         cout << "\\\\" << endl;
      }
      cout << "  \\dot{u}^{(";
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
  return;
}

int findPlace(const int a, const int b, const int n)
{
  return a*n+b-(a*(a+1))/2;
}


void constructAugmentedInitVal(vector<mpq_class>& runInitVal, const vector<mpq_class> initVal, map<vector<int>, int> myMap)
{
  map<vector<int>, int> ::iterator myMapIterator;
  vector<int> tempVec;
  mpq_class tempVal(1);
  mpq_class tempInside(1);



  for(myMapIterator=myMap.begin();myMapIterator !=myMap.end();++myMapIterator)
  {
    tempVec = myMapIterator->first;
    tempVal = 1;

    for(int i = 0; i<tempVec.size(); i++)
    {

      tempInside = 1;
      for(int j=1; j<=tempVec[i]; j++)
      {
        tempInside = tempInside * initVal[i];
      }
      tempVal = tempVal * tempInside;
    }

    runInitVal[myMapIterator->second] = tempVal;
  }
}

void constructF(vector<int>& rowInd, vector<int>& colInd, vector<mpq_class>& myValue, const vector<int> runVec, map<vector<int>, int> myMap, const vector<mpq_class> myCoeffs, const int stEq)
{

  vector<int> left(stEq);
  vector<int> middle(stEq);
  vector<int> right(stEq);

  for(auto i = 0; i < runVec.size(); i+=3*stEq)
  {
    left = { runVec.begin() + i, runVec.begin() + i + stEq };
    middle = { runVec.begin() + i + stEq, runVec.begin() + i + stEq +stEq };
    right = { runVec.begin() + i + stEq + stEq, runVec.begin() + i + stEq + stEq + stEq};

    const int leftAsInd = myMap[left];
    const int middleAsInd = myMap[middle];
    const int rightAsInd = myMap[right];

    rowInd.push_back(leftAsInd);

    colInd.push_back(findPlace(middleAsInd, rightAsInd, myMap.size()));
    if (middleAsInd == rightAsInd)
    {
      myValue.push_back(myCoeffs[i/(3*stEq)]);
    }
    else
    {
      myValue.push_back(myCoeffs[i/(3*stEq)]/2);
    }
  }

  printAsTripleF(rowInd, colInd, myValue);
  return;
}

void printAsTripleF(const vector<int> rowInd, const vector<int> colInd, const vector<mpq_class> myValue)
{
  cout << setw(5) << "row" << " " << setw(5) << "col" << " " << setw(5) << "val" << endl;

  for(auto i = 0; i < rowInd.size(); i++)
    cout << setw(5) << rowInd[i] << " " << setw(5) << colInd[i] << " " << setw(5) << myValue[i] << endl;

  return;
}

void condensedKroneckerProduct(vector<mpq_class>& result, const vector<mpq_class> a, const vector<mpq_class> b)
{
  assert(a.size() == b.size());
  const int n = a.size();
  result.clear();

  for (int i = 0; i < n; i++)
  {
    result.push_back(a[i]*b[i]);
    for (int j = i+1; j < n; j++)
    {
      result.push_back(a[i]*b[j]+a[j]*b[i]);
    }
  }
  return;
}

void sparseMatTimesVec(vector<mpq_class>& result, const vector<int> rowInd, const vector<int> colInd, const vector<mpq_class> myCoeffs, const vector<mpq_class> x, const int resSize)
{
  assert(rowInd.size()==colInd.size());
  assert(rowInd.size()==myCoeffs.size());

  result.clear();
  result.resize(resSize,0);

  for(int i=0; i < rowInd.size(); i++)
  {
    result[rowInd[i]] += myCoeffs[i] * x[colInd[i]];
  }
  return;
}

// END OF PROGRAM
