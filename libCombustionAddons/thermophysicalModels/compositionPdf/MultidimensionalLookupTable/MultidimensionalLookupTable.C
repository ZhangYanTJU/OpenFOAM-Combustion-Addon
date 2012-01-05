#include "MultidimensionalLookupTable.H"

#include "/usr/include/math.h"

// 4. Änderung: Zusammenfassen
double Round(double Zahl, int Stellen)
{
  return floor(Zahl * pow( 10, Stellen) + 0.5) * pow(10, -Stellen);
}
     
