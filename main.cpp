#include <iostream>
#include <vector>

using namespace std;

int
main ()
{

  vector < int >equationVector
  {
  1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0};


  vector < int >runVec = equationVector;

  bool appears = 0;

  int i = 0;
  int temp = 0;
  int temp2 = 0;
  int temp3 = 0;
  int temp4 = 0;
  int temp5 = 0;
  int temp6 = 0;

  vector < int >rightHandSide
  {
  0, 0};

  while (i < runVec.size ())
  {
    for (int j = 2; j <= 4; j += 2)
    {
      rightHandSide[0] = runVec[i + j];
      rightHandSide[1] = runVec[i + j + 1];

      for (int i = 0; i < runVec.size (); i += 6)
      {
	if (rightHandSide[0] == runVec[i]
	    && rightHandSide[1] == runVec[i + 1])
	{
	  appears = 1;
	  break;
	}
      }

      if (appears == 0)
      {
	if (rightHandSide[0] == 0 && rightHandSide[1] == 0)
	{
	  runVec.insert (runVec.end (),
	                 {
	                 0, 0, 0, 0, 0, 0}
	  );
	}
	else
	{

	  for (int i = 0; i < equationVector.size (); i += 6)
	  {
	    if (equationVector[i] == 1 && equationVector[i + 1] == 0)
	    {
	      if (rightHandSide[0] > 0)
	      {
		temp =
		  equationVector[i + 2] + equationVector[i + 4] +
		  rightHandSide[0] - 1;
		temp2 =
		  equationVector[i + 3] + equationVector[i + 5] +
		  rightHandSide[1];

		if (temp == temp2)
		{
		  temp3 = temp / 2;
		  temp5 = temp - temp3;
		  runVec.push_back (rightHandSide[0]);
		  runVec.push_back (rightHandSide[1]);
		  runVec.push_back (temp3);
		  runVec.push_back (temp5);
		  runVec.push_back (temp5);
		  runVec.push_back (temp3);
		}
		else
		{
		  temp3 = temp / 2;
		  temp5 = temp - temp3;
		  temp4 = temp2 / 2;
		  temp6 = temp2 - temp4;
		  runVec.push_back (rightHandSide[0]);
		  runVec.push_back (rightHandSide[1]);
		  runVec.push_back (temp3);
		  runVec.push_back (temp4);
		  runVec.push_back (temp5);
		  runVec.push_back (temp6);
		}
	      }

	    }
	    else if (equationVector[i] == 0 && equationVector[i + 1] == 1)
	    {
	      if (rightHandSide[1] > 0)
	      {
		temp =
		  equationVector[i + 2] + equationVector[i + 4] +
		  rightHandSide[0];
		temp2 =
		  equationVector[i + 3] + equationVector[i + 5] +
		  rightHandSide[1] - 1;

		if (temp == temp2)
		{
		  temp3 = temp / 2;
		  temp5 = temp - temp3;
		  runVec.push_back (rightHandSide[0]);
		  runVec.push_back (rightHandSide[1]);
		  runVec.push_back (temp3);
		  runVec.push_back (temp5);
		  runVec.push_back (temp5);
		  runVec.push_back (temp3);
		}
		else
		{
		  temp3 = temp / 2;
		  temp5 = temp - temp3;
		  temp4 = temp2 / 2;
		  temp6 = temp2 - temp4;
		  runVec.push_back (rightHandSide[0]);
		  runVec.push_back (rightHandSide[1]);
		  runVec.push_back (temp3);
		  runVec.push_back (temp4);
		  runVec.push_back (temp5);
		  runVec.push_back (temp6);
		}
	      }
	    }
	    else
	    {
	      cout << "If the program gets here, there is sth wrong." << endl;
	    }

	  }
	}

      }
      appears = 0;

    }
    i += 6;
  }
  for (int i = 0; i < runVec.size (); i++)
  {
    cout << runVec[i] << ' ';

    if (i % 6 == 5)
      cout << endl;
  }

  return 0;
}

// END OF PROGRAM
