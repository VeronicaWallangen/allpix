/*
 *  Authors:
 *    Veronica Wallangen veronica.wallangen@cern.ch
 *    Benjamin Nachman bnachman@cern.ch
 *    Gilberto Giugiarelli gilberto.giugiarelli@cern.ch
 *
 *  allpix Authors:
 *    John Idarraga <idarraga@cern.ch>
 *    Mathieu Benoit <benoit@lal.in2p3.fr>
 *
 *  Brief Description:
 *    This is a digitization model for radiation damage effects in 3D pixel sensors.
 *
 */


#include "TMath.h"
#include "CLHEP/Random/RandGauss.h"

//temporary to write to textfiles
#include <iostream>
#include <fstream>
#include <vector>
//#include <algorithm>

void CLHEPtest() {
    
    for (int i=0; i>100; i++) {
        cout << CLHEP::RandFlat::shoot(0.,1.) << "\n";
    }
    
}
