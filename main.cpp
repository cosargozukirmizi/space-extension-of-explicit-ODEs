#include <iostream>
#include <vector>
#include <map>
#include <gmpxx.h>
#include <iomanip>
#include <thread>
#include <semaphore>
#include <algorithm>
#include <cassert>
#include <chrono>

void extendSpace (const std::vector<int>& equationVector, std::vector<mpq_class>& myCoeffs, const std::vector<mpq_class>& initVal, const int max_iter = 10, const bool newHeuristic = 0);
void toLatex (const std::vector<int>& equationVector, const int stEq, const std::vector<mpq_class>& myCoeffs);
void printAsTriple(const std::vector<int>& rowInd, const std::vector<int>& colInd, const std::vector<mpq_class>& myValue);
void construct2F(std::vector<int>& rowInd, std::vector<int>& colInd, std::vector<mpq_class>& myValue, const std::vector<int>& runVec, const std::map<std::vector<int>, int>& myMap, const std::vector<mpq_class>& myCoeffs, const int stEq);
void condensedKroneckerProduct(std::vector<mpq_class>& result, const std::vector<mpq_class>& a, const std::vector<mpq_class>& b);
void sparseMatTimesVec(std::vector<mpq_class>& result, const std::vector<int>& rowInd, const std::vector<int>& colInd, const std::vector<mpq_class>& myCoeffs, const std::vector<mpq_class>& x);
void constructAugmentedInitVal(std::vector<mpq_class>& runInitVal, const std::vector<mpq_class>& initVal, const std::map<std::vector<int>, int>& myMap);


int main ()
{
  using namespace std;

  const int max_iter = 50;   // maximum truncation order is one less than max_iter
  const bool newHeuristic = 1; // if this is true, heuristic H4 is used; otherwise heuristic H1 is used
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

  const vector<mpq_class> vanderPolInitVal{
   mpq_class(1,2), mpq_class(1,2)
  };

  cout << "\n The space extension for van der Pol ODE:";
  cout << '\n';
  extendSpace (vanderPol, vanderPolCoeffs, vanderPolInitVal, max_iter,newHeuristic);
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

  const vector<mpq_class> quarticAnharmonicOscillatorInitVal{
   mpq_class(1,2),mpq_class(1,2)
  };

  cout << "\n The space extension for classical quartic anharmonic oscillator ODE:";
  cout << '\n';
  extendSpace (quarticAnharmonicOscillator, quarticAnharmonicOscillatorCoeffs, quarticAnharmonicOscillatorInitVal, max_iter, newHeuristic);
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

  const vector<mpq_class> henonHeilesInitVal{
   mpq_class(1,2), mpq_class(1,2), mpq_class(1,2), mpq_class(1,2)
  };

  cout << "\n The space extension for Henon-Heiles ODE:";
  cout << '\n';
  extendSpace (henonHeiles, henonHeilesCoeffs, henonHeilesInitVal, max_iter, newHeuristic);
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

  const vector<mpq_class> rabinovichFabrikantInitVal{
   mpq_class(1,2), mpq_class(1,2), mpq_class(1,2)
  };

  cout << "\n The space extension for Rabinovich-Fabrikant ODE:";
  cout << '\n';
  extendSpace (rabinovichFabrikant, rabinovichFabrikantCoeffs, rabinovichFabrikantInitVal, max_iter, newHeuristic);
  cout << "--------------------------------\n\n";

  return 0;
}


void extendSpace (const std::vector<int>& equationVector, std::vector<mpq_class>& myCoeffs, const std::vector<mpq_class>& initVal, const int max_iter, const bool newHeuristic)
{
  using namespace std;
  assert(max_iter > 1);
  const int stEq = initVal.size();
  vector<int> runVec = equationVector;
  map<vector<int>, int> myMap;

  bool appears = 0;

  vector<int>::size_type i = 0;
  auto numEqs = stEq;

  vector<int> rightHandSide (stEq, 0);
  vector<int> tempLeft (stEq, 0);

  for (vector<int>::size_type i = 0; i < runVec.size (); i+=3*stEq)
  {
    tempLeft.assign(equationVector.begin()+i, equationVector.begin()+i+stEq);
    myMap.insert({tempLeft, 0});
  }

  while (i < runVec.size ())
  {
    for (auto j = stEq; j <= 2 * stEq; j += stEq)
    {
      for (auto k = 0; k < stEq; k++)
      {
        rightHandSide[k] = runVec[i + j + k];
      }
      if (myMap.find(rightHandSide) == myMap.end())
      {
         myMap.insert({rightHandSide, 0});
      }

      for (vector<int>::size_type i = 0; i < runVec.size (); i += 3 * stEq)
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
        ++numEqs;
        vector<int> allZero(stEq,0);
        vector<int> zeroEntry(3*stEq, 0);

        if (rightHandSide == allZero)
        {
          runVec.insert (runVec.end (), zeroEntry.begin(), zeroEntry.end() );
          myCoeffs.push_back(mpq_class(0));
        }
        else
        {
          for (vector<int>::size_type i = 0; i < equationVector.size (); i += 3*stEq)
          {
            vector<int>::const_iterator it;
            it = find (equationVector.begin () + i, equationVector.begin () + i + 3*stEq, 1);
            vector<int>::size_type oneInd = (it - equationVector.begin ());

            bool x = 0;
            vector<int> myTemp(stEq, 0);

            for (auto j = 0; j < stEq; j++)
            {
              if (oneInd == i+j && rightHandSide[j] > 0)
              {
                x = 1;
                for (auto k = 0 ; k < stEq; k++)
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
            for (auto counter = 0; counter < stEq; counter++)
              runVec.push_back (rightHandSide[counter]);

            if (newHeuristic == 0)
            {
              int oneCount = 0;

              for (vector<int>::size_type i = 0; i < myTemp.size (); ++i)
              {
                if (myTemp[i] == 1)
                {
                  ++oneCount;
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
            }
            else
            {
              int oddCount = 0;

              for (vector<int>::size_type i = 0; i < myTemp.size (); ++i)
              {
                if (myTemp[i] % 2 == 1)
                {
                  ++oddCount;
                }

                if (myTemp[i] % 2 == 0 || (myTemp[i] % 2 == 1 && (oddCount%2 == 1)))
                {
                  lefts.push_back(myTemp[i] / 2);
                  runVec.push_back(myTemp[i] / 2);
                }
                else if (myTemp[i] % 2 == 1 && oddCount%2 == 0)
                {
                  lefts.push_back( (myTemp[i]+1) / 2 );
                  runVec.push_back( (myTemp[i]+1) / 2 );
                }
              }
            }

            for (vector<int>::size_type i = 0; i < lefts.size (); ++i)
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

  cout << "\nThe new ODE set is \n";
  toLatex (runVec, stEq, myCoeffs);

  cout << '\n';

  for(auto myDummyVariable = 0; auto& [key, value] : myMap)
  {
    value = myDummyVariable++;
  }

  cout << "\nThe map for stacking the equations is \n";

  for(const auto& [key, value] : myMap)
  {
    cout << "[ ";

    for (const auto& dummy : key)
      cout << dummy << " ";

    cout << "]:";
    cout << value<< ", " << '\n';
  }

  vector<mpq_class> runInitVal(myMap.size());

  cout << "\n" << "Now we will try to construct initial vector.\n";
  cout << "It is given as index value pairs below. \n";

  constructAugmentedInitVal(runInitVal, initVal, myMap);


  for(map<vector<int>, int>::const_iterator myMapIt=myMap.begin(); myMapIt != myMap.end(); ++myMapIt)
  {
    cout << myMapIt->second << " " << runInitVal[myMapIt->second] << '\n';
  }

  cout << "\n" << "Now we will try to construct 2*F: \n";

  vector<int> rowInd;
  vector<int> colInd;
  vector<mpq_class> myValue;

  construct2F(rowInd, colInd, myValue, runVec, myMap, myCoeffs, stEq);


  vector<vector<mpq_class>> rho(max_iter, vector<mpq_class>(runInitVal.size()));

  auto start = std::chrono::steady_clock::now();

  rho.at(0) = (runInitVal);

  vector<mpq_class> resOfKronProd(runInitVal.size()*(runInitVal.size()+1)/2, 0);

  condensedKroneckerProduct(resOfKronProd, rho[0], rho[0]);
  std::transform(resOfKronProd.begin(), resOfKronProd.end(), resOfKronProd.begin(), [](mpq_class& c){return c / 2;});

  sparseMatTimesVec(rho.at(1), rowInd, colInd, myValue, resOfKronProd);

  std::binary_semaphore smphSignalMainToThread{1}, smphSignalThreadToMain{0};  // oddSummer has green light, it is good to go

  auto evenSummer = [&]()
  {

    for(auto j=2; j < max_iter-1; j=j+2)
    {
       vector<mpq_class> resOfKronProd(runInitVal.size()*(runInitVal.size()+1)/2, 0);
       vector<mpq_class> resOfMatVecProd(runInitVal.size(), 0);
       vector<mpq_class> resOfKronProdAccumulator(runInitVal.size()*(runInitVal.size()+1)/2, 0);

       std::fill(resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), 0);
       condensedKroneckerProduct(resOfKronProd, rho[j/2], rho[j/2]);
       std::transform(resOfKronProd.begin(), resOfKronProd.end(), resOfKronProd.begin(), [](mpq_class& c){return c / 2;});
       std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());


       for(auto k=(j/2-1); k>0; --k)
       {
         condensedKroneckerProduct(resOfKronProd, rho[k], rho[j-k]);
         std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());
       }

       smphSignalThreadToMain.acquire();

       condensedKroneckerProduct(resOfKronProd, rho[0], rho[j]);
       std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());

       sparseMatTimesVec(resOfMatVecProd, rowInd, colInd, myValue, resOfKronProdAccumulator);

       mpq_class myConst(j+1);

       std::transform(resOfMatVecProd.begin(), resOfMatVecProd.end(), resOfMatVecProd.begin(), [&myConst](mpq_class& c){return c / myConst;});

       rho.at(j+1) = resOfMatVecProd;

       smphSignalMainToThread.release();
    }

  };


  auto oddSummer = [&]()
  {
    for(auto j=1; j < max_iter-1; j=j+2)
    {
       vector<mpq_class> resOfKronProd(runInitVal.size()*(runInitVal.size()+1)/2, 0);
       vector<mpq_class> resOfMatVecProd(runInitVal.size(), 0);
       vector<mpq_class> resOfKronProdAccumulator(runInitVal.size()*(runInitVal.size()+1)/2, 0);

       std::fill(resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), 0);

       for(auto k=(j-1)/2; k>0; --k)  // 0 not included for a reason
       {
         condensedKroneckerProduct(resOfKronProd, rho[k], rho[j-k]);
         std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());
       }

       smphSignalMainToThread.acquire();

       condensedKroneckerProduct(resOfKronProd, rho[0], rho[j]);
       std::transform (resOfKronProdAccumulator.begin(), resOfKronProdAccumulator.end(), resOfKronProd.begin(), resOfKronProdAccumulator.begin(), std::plus<mpq_class>());

       sparseMatTimesVec(resOfMatVecProd, rowInd, colInd, myValue, resOfKronProdAccumulator);

       mpq_class myConst(j+1);

       std::transform(resOfMatVecProd.begin(), resOfMatVecProd.end(), resOfMatVecProd.begin(), [&myConst](mpq_class& c){return c / myConst;});

       rho.at(j+1) = resOfMatVecProd;

       smphSignalThreadToMain.release();
    }
  };

  thread t(oddSummer);
  evenSummer();

  t.join();

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed wall clock time for computing rho vectors: " << elapsed_seconds.count() << "s\n";

  cout << "\nAfter "<< max_iter-1 << " iterations, the rho vectors are\n";

  for (auto i = 0; const auto& myRho : rho)
  {
    cout << "\nrho[" << i++ << "] :\n";
    for (const auto& element : myRho)
    {
       cout << element << '\n';
    }
  }

  cout << '\n';


  const mpq_class t_max(1, 10);    // time value under consideration
  const mpq_class t_step(t_max/10);


  for(mpq_class indQ=0; indQ<=t_max; indQ=indQ+t_step)
  {

    vector<mpq_class> solution(runInitVal.size(), 0);

    cout << '\n';
    mpq_class leftPart(1);

    for(auto i = 0; i < max_iter; ++i)
    {
      vector<mpq_class> solutionIngredient(rho[i].begin(), rho[i].end());

      std::transform (solutionIngredient.begin(), solutionIngredient.end(), solutionIngredient.begin(), [&leftPart](auto& c){return c*leftPart;});
      std::transform (solution.begin(), solution.end(), solutionIngredient.begin(), solution.begin(), std::plus<mpq_class>());

      cout << "\n\nAfter "<< i << " iteration" << (i==1 ? ", " : "s, ") << "the solution of the initial value problem at t=" << indQ <<  " is\n";

      // truncation order is the number of terms.
      // truncation order is also the number of iterations.

      for(const auto& [key, value] : myMap)
      {
        cout << "[ ";

        for (const auto& dummy : key)
          cout << dummy << " ";

        cout << "]: ";
        cout << solution[value] << '\n' << "= ";

        mpf_class f(solution[value], 500);
        gmp_printf ("%.*Ff\n", 50, f);
      }

      leftPart *= indQ;
    }

  }
  cout << '\n';

  return;
}

void toLatex (const std::vector<int>& equationVector, const int stEq, const std::vector<mpq_class>& myCoeffs)
{
  using namespace std;

  cout << '\n';
  cout << "\\begin{eqnarray}\n";

  vector<int> prevLeftHandSide(stEq);
  vector<int> leftHandSide(stEq);
  vector<int> allZero(stEq,0);

  for (vector<int>::size_type i = 0; i < equationVector.size(); i+=3*stEq)
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
         cout << "\\\\" << '\n';
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
         cout << "\\\\" << '\n';
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

    bool isPos = 1;

    if ( sgn(myCoeffs[i/(3*stEq)]) == -1 )
    {
       isPos = 0;
    }


    mpq_class coeffToPrint = myCoeffs[i/(3*stEq)];
    if ( same )
    {
      if ( isPos )
      {
        cout << " \\nonumber \\\\" << '\n' << " &+& ";
      }
      else
      {
        cout << " \\nonumber \\\\" << '\n' << " &-& ";
        coeffToPrint = -coeffToPrint;
      }
    }
    cout << coeffToPrint << "\\,";

    cout << " u^{(";
    for (auto k = i + stEq; k < i + stEq + stEq; ++k)
    {
      if (k != i + stEq)
      {
        cout << ",";
      }
      cout << equationVector[k];
    }
    cout << ")}";

    cout << " u^{(";
    for (auto k = i + stEq + stEq; k < i + stEq + stEq + stEq; ++k)
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
  cout << '\n' << "\\end{eqnarray}" << '\n';
  return;
}

int findPlace(const int a, const int b, const int n)
{
  return a*n+b-(a*(a+1))/2;
}


void constructAugmentedInitVal(std::vector<mpq_class>& runInitVal, const std::vector<mpq_class>& initVal, const std::map<std::vector<int>, int>& myMap)
{
  using namespace std;

  for(const auto& [key, value] : myMap)
  {
    mpq_class tempVal(1);

    for(vector<int>::size_type i = 0; i<key.size(); ++i)
    {
      mpq_class tempInside(1);

      for(auto j = 1; j<=key[i]; ++j)
      {
        tempInside = tempInside * initVal[i];
      }
      tempVal = tempVal * tempInside;
    }

    runInitVal[value] = tempVal;
  }
}


// The function below constructs 2*F.
void construct2F(std::vector<int>& rowInd, std::vector<int>& colInd, std::vector<mpq_class>& myValue, const std::vector<int>& runVec, const std::map<std::vector<int>, int>& myMap, const std::vector<mpq_class>& myCoeffs, const int stEq)
{
  using namespace std;

  vector<int> left(stEq);
  vector<int> middle(stEq);
  vector<int> right(stEq);

  for(vector<int>::size_type i = 0; i < runVec.size(); i+=3*stEq)
  {
    left = { runVec.begin() + i, runVec.begin() + i + stEq };
    middle = { runVec.begin() + i + stEq, runVec.begin() + i + stEq +stEq };
    right = { runVec.begin() + i + stEq + stEq, runVec.begin() + i + stEq + stEq + stEq};

    const int leftAsInd = myMap.at(left);
    const int middleAsInd = myMap.at(middle);
    const int rightAsInd = myMap.at(right);

    rowInd.push_back(leftAsInd);

    colInd.push_back(findPlace(middleAsInd, rightAsInd, myMap.size()));
    if (middleAsInd == rightAsInd)
    {
      myValue.push_back(myCoeffs[i/(3*stEq)]*2);
    }
    else
    {
      myValue.push_back(myCoeffs[i/(3*stEq)]);
    }
  }

  printAsTriple(rowInd, colInd, myValue);
  return;
}

void printAsTriple(const std::vector<int>& rowInd, const std::vector<int>& colInd, const std::vector<mpq_class>& myValue)
{
  using namespace std;

  cout << setw(5) << "row" << ' ' << setw(5) << "col" << ' ' << setw(5) << "val" << '\n';

  for(vector<int>::size_type i = 0; i < rowInd.size(); i++)
    cout << setw(5) << rowInd[i] << ' ' << setw(5) << colInd[i] << ' ' << setw(5) << myValue[i] << '\n';

  return;
}

void condensedKroneckerProduct(std::vector<mpq_class>& result, const std::vector<mpq_class>& a, const std::vector<mpq_class>& b)
{
  using namespace std;

  vector<mpq_class>::size_type myIndex = 0;

  for (vector<mpq_class>::size_type i = 0; i < a.size(); ++i)
  {
    result[myIndex++] = (a[i]*b[i]);
    for (vector<mpq_class>::size_type j = i+1; j < a.size(); ++j)
    {
      result[myIndex++] = (a[i]*b[j]+a[j]*b[i]);
    }
  }
  return;
}

void sparseMatTimesVec(std::vector<mpq_class>& result, const std::vector<int>& rowInd, const std::vector<int>& colInd, const std::vector<mpq_class>& myCoeffs, const std::vector<mpq_class>& x)
{
  using namespace std;

  std::fill(result.begin(), result.end(), mpq_class(0));

  for(vector<int>::size_type i = 0; i < rowInd.size(); ++i)
  {
    result[rowInd[i]] += myCoeffs[i] * x[colInd[i]];
  }

  return;
}

// END OF PROGRAM
