#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <assert.h>

int main(){
  // Number of cols in input data files
  const int numCols = 20;

  // Names of input data files, sorted according to Q2
  const char * genFile = "Amunu.txt";
  const char * fortFile = "fort.73";

  // Acceptable fractional difference between kinematics
  // (genie - fortran)/fortran < error
  double Q2Error = 0.01;
  double q0Error = 0.01;
  double dqError = 0.01;
  double rError = 0.01;

  std::string headers[21] = {"type","Q2","q0","dq","radius",
			     "A00","A03","A11","A12","A33",
			     "T0","T3","R00","R03","R11","R33",
			     "kFn","kFp","CN","CT","CL"};

  // Stream to write data for matching kinematics to file
  std::ofstream dataStream;
  dataStream.open("matching_kine.txt");
  for(int i=0; i<numCols+1; i++)
    dataStream << std::left << std::setw(12) << headers[i];
  dataStream << std::endl << std::endl;

  // Stream to test what the allowed errors should be
  std::ofstream testStream;
  testStream.open("error_test.txt");
  for(int i=0; i<5; i++)
    testStream << std::left << std::setw(12) << headers[i];
  testStream << std::endl << std::endl;

  // Streams to read from the data file
  // Assume that the files are sorted accoring to Q2,
  // so if the Q2 values don't match I can get the next line
  // from the one with a lower Q2 and the next line will
  // have a higher Q2
  std::ifstream genStream(genFile); // GENIE output
  std::ifstream fortStream(fortFile); // Fortran output

  bool readFile = false;
  std::string genLine;
  std::string fortLine;

  readFile = std::getline(genStream, genLine);
  readFile = std::getline(fortStream, fortLine);
  
  double genVals[numCols];
  double fortVals[numCols];

  while(readFile){
    // Move data to string stream to read
    std::istringstream genStrStream(genLine);
    std::istringstream fortStrStream(fortLine);

    // Check Q2
    assert(genStrStream >> genVals[0]);
    assert(fortStrStream >> fortVals[0]);
    if(std::abs(genVals[0]-fortVals[0])/std::abs(fortVals[0]) < Q2Error){
      // Check other kinematics
      assert(genStrStream >> genVals[1] >> genVals[2] >> genVals[3]);
      assert(fortStrStream >> fortVals[1] >> fortVals[2] >> fortVals[3]);

      testStream << std::left << std::setw(12) << "GENIE";
      for(int i=0; i<4; i++)
	testStream << std::left << std::setw(12) << genVals[i];
      testStream << std::endl;
      testStream << std::left << std::setw(12) << "Fortran";
      for(int i=0; i<4; i++)
	testStream << std::left << std::setw(12) << fortVals[i];
      testStream << std::endl << std::endl;

      if(std::abs(genVals[1]-fortVals[1])/std::abs(fortVals[1]) < q0Error &&
	 std::abs(genVals[2]-fortVals[2])/std::abs(fortVals[2]) < dqError &&
	 std::abs(genVals[3]-fortVals[3])/std::abs(fortVals[3]) < rError){
	// Get the values of the nucleon tensor and print all
	for(int i=4; i<numCols; i++){
	  assert(genStrStream >> genVals[i]);
	  assert(fortStrStream >> fortVals[i]);
	}

	// Print values from GENIE
	dataStream << std::left << std::setw(12) << "GENIE";
	for(int i=0; i<numCols; i++){
	  dataStream << std::left << std::setw(12) << genVals[i];
	}
	dataStream << std::endl;

	// Print values from Fortran code
	dataStream << std::left << std::setw(12) << "Fortran";
	for(int i=0; i<numCols; i++)
	  dataStream << std::left << std::setw(12) << fortVals[i];
	dataStream << std::endl;

	// Print difference
	dataStream << std::left << std::setw(12) << "Frac Diff";
	for(int i=0; i<numCols; i++)
	  dataStream << std::left << std::setw(12) 
		     << std::abs(fortVals[i]-genVals[i])/std::abs(fortVals[i]);
	dataStream << std::endl << std::endl;

      } // Check q0, dq, and r match
    } // Check if q2 matches

    // Try to read the next line from the file with lower Q2
    if(genVals[0] <= fortVals[0]){
      readFile = std::getline(genStream, genLine);
    }else{
      readFile = std::getline(fortStream, fortLine);
    }
  }// While(readFile)

  dataStream.close();
  testStream.close();
  std::cout << "read_sorted_data finished\n";
  return 0;
} 
