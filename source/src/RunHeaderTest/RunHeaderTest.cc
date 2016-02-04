#include <iostream>
#include <iomanip>

int main(int argc, char** argv)
{
  for (int i=0; i<argc; ++i)
    std::cout << std::setw(4) << i << " : "  << argv[i] << std::endl;
  return 0;
}
