#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main(int argc, char * argv[]){

  fstream myfile("test.skf", ios_base::in);

  double a;
  while (myfile >> a)
  {
      printf("%f ", a);
  }

  getchar();

  return 0;
}
