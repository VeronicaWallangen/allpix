/**
 *  Author:
 *    John Idarraga <idarraga@cern.ch>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#ifndef AllPixLETCalculatorDigit_h
#define AllPixLETCalculatorDigit_h 1

#include "AllPixDigitInterface.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/**
 *  Digit AllPixLETCalculator class.
 */
class AllPixLETCalculatorDigit : public AllPixDigitInterface {

public:
  AllPixLETCalculatorDigit();
  ~AllPixLETCalculatorDigit();

  AllPixLETCalculatorDigit(const AllPixLETCalculatorDigit&);
  const AllPixLETCalculatorDigit& operator=(const AllPixLETCalculatorDigit&);
  int operator==(const AllPixLETCalculatorDigit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4int pixelIDX;    
  G4int pixelIDY;    
  G4int pixelCounts;
  G4double pixelEnergyDep; // MC // Corrected MC charge (detector effects included, at Digi step)
  G4ThreeVector primaryVertex;
  
public:
  
  inline void SetPixelIDX(G4int pidX)   {pixelIDX = pidX;};
  inline void SetPixelIDY(G4int pidY)   {pixelIDY = pidY;};
  inline void SetPixelCounts(G4int pc)  {pixelCounts = pc;};
  inline void SetPixelEnergyDep(G4double ed)  {pixelEnergyDep = ed;}; // MC // Corrected MC charge (detector effects included, at Digi step)
  inline void SetPrimaryVertex(G4ThreeVector pv)  {primaryVertex = pv;}; // MC vertex //
  inline void IncreasePixelCounts()     {pixelCounts++;};

  inline G4int GetPixelIDX()   {return pixelIDX;};
  inline G4int GetPixelIDY()   {return pixelIDY;};
  inline G4int GetPixelCounts()  {return pixelCounts;};
  inline G4double GetPixelEnergyDep()  {return pixelEnergyDep;}; // MC //
  inline G4ThreeVector GetPrimaryVertex()  {return primaryVertex;}; // MC //

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<AllPixLETCalculatorDigit> AllPixLETCalculatorDigitsCollection;

extern G4Allocator<AllPixLETCalculatorDigit> AllPixLETCalculatorDigitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* AllPixLETCalculatorDigit::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) AllPixLETCalculatorDigitAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void AllPixLETCalculatorDigit::operator delete(void* aDigi)
{
  AllPixLETCalculatorDigitAllocator.FreeSingle((AllPixLETCalculatorDigit*) aDigi);
}

#endif

