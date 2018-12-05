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

#include "AllPixFEI4RadDamage3DDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"
#include "TMath.h"

//Included for RadDamage
#include "TFile.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"
#include "TROOT.h"

//temporary to write to textfiles
#include <iostream>
#include <fstream>
#include <vector>
//#include <algorithm>

AllPixFEI4RadDamage3DDigitizer::AllPixFEI4RadDamage3DDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	//////////////////////////////////////////////////////////////////
	/////Setup all the inputs/////////////////////////////////////////
	//////////////////////////////////////////////////////////////////

	// Physics switches
	ITkGeometry = false; //use ITk dimensions (50umx50um), if false use IBL dimensions (250umx50um)
	doDiff = true; //include diffusion
    doRamo = true; //include induced charge (set to false to neglect induced charge and only record charge from electrons that make it all the way to the collecting electrode)
	doRamoNeighbors = false; //include induced charge contribution from nearest neighboring pixels
	doChunkCorrection = false; //include unsmearing

	// Fetching geometry info from pixeldetector.xml
	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	detectorThickness = gD->GetSensorZ();
	pitchX = gD->GetPixelX(); 	// total length of pixel in x
	pitchY = gD->GetPixelY();	// total length of pixel in y
	nPixX= gD->GetNPixelsX(); 	// total number of pixels in x
	nPixY= gD->GetNPixelsY();	// total number of pixels in y

	// Phenomenological parameters
	precision = 50; //this is the number of charges to divide the G4 hit into

	// Constants
	betaElectrons = 3.0E-16; //cm2/ns; The charge-trapping probability is t = -tau*ln(u), for u ~ Uniform(0,1) and tau^{-1}=beta*fluence.  The value of beta might be slightly higher for holes than electrons, but it is hard to say (we ignore this effect here).  See e.g. https://cds.cern.ch/record/685542/files/indet-2003-014.pdf.
	//Conditions
	elec = 3.64*eV;//average energy required to produce a e/h pair in Si.  This has been known for a long time - see for instance http://journals.aps.org/prb/pdf/10.1103/PhysRevB.1.2945.  It depends a little on temperature, but we are ignoring that effect here.
	fluence = 1e14; // neq/cm^2
	temperature = 263.5;// ITk: 258 K, IBL: 263.5 K
	threshold = 600*elec; //This is the threshold for charge collection. Min. stable threshold RD53A: 600e
	//tuning = 5./(20000*elec); //for X ToT @ Y e, this is X/Y so that a deposited energy of Y gives a ToT of X.  Typical values are 5 ToT @ 20ke.
	tuning = 5./(11000*elec);

	//////////////////////////////////////////////////////////////////
	//Maps: E-field, charge collection time, Ramo potential
	//////////////////////////////////////////////////////////////////

	std::cout << "Load the input maps " << std::endl;

	// Get electric field data
	eFieldMap=0;
	TFile* inputfile;
	if (!ITkGeometry) {
		inputfile=new TFile("/allpix/allpix/TCADmaps/phi_1e14_20V.root");
	} else {
		inputfile=new TFile("/allpix/allpix/TCADmaps/TCADmaps/phi_1e14_20V.root");
	}
	eFieldMap=(TH2F*)inputfile->Get("efield");
	if (eFieldMap == 0) G4cout << "Unsuccessful picking up histogram: eFieldMap" << G4endl;
	
	// Get charge collection time data (expected time to reach electrode)
	timeMap_e=0;
	timeMap_h=0;
	timeMap_e=(TH2F*)inputfile->Get("etimes");
	timeMap_h=(TH2F*)inputfile->Get("htimes");
	if (timeMap_e == 0 || timeMap_h == 0) G4cout << "Unsuccessful picking up histogram: timeMap" << G4endl;
	
	// Get charge position data
	xPositionMap_e=0;
	xPositionMap_h=0;
	yPositionMap_e=0;
	yPositionMap_h=0;
	xPositionMap_e=(TH3F*)inputfile->Get("xPosition_e");
	xPositionMap_h=(TH3F*)inputfile->Get("xPosition_h");
	yPositionMap_e=(TH3F*)inputfile->Get("yPosition_e");
	yPositionMap_h=(TH3F*)inputfile->Get("yPosition_h");
	if (xPositionMap_e == 0 || xPositionMap_h == 0 || yPositionMap_e ==0 || yPositionMap_h == 0) G4cout << "Unsuccessful picking up histogram: positionMap" << G4endl;
	
	// Get Ramo data
	ramoPotentialMap=0;
	ramoPotentialMap=(TH2F*)inputfile->Get("ramo");
	if (ramoPotentialMap == 0){
	  G4cout << "Unsuccessful picking up histogram: ramoPotentialMap" << G4endl;
        }

	// Get average charge data (for charge chunk effect correction)
	avgChargeMap_e=0;
	avgChargeMap_h=0;
	avgChargeMap_e=(TH2F*)inputfile->Get("avgCharge_e");
	avgChargeMap_h=(TH2F*)inputfile->Get("avgCharge_h");
	if (avgChargeMap_e == 0 || avgChargeMap_h == 0) {
	  G4cout << "Unsuccessful picking up histogram: avgChargeMap" << G4endl;
        }

	//Compute the trapping time.
	if(fluence!=0.0)
	  {
	    trappingTimeElectrons = 1.0/(betaElectrons*fluence); //[ns]
	  }
	else { //fluence = 0 so do not trap!
	  trappingTimeElectrons = 1000;//*s;
	}
	trappingTimeHoles = trappingTimeElectrons;

}


AllPixFEI4RadDamage3DDigitizer::~AllPixFEI4RadDamage3DDigitizer(){

}


///////////////////////////////////
// Radiation Damage Calculations //
///////////////////////////////////

G4double AllPixFEI4RadDamage3DDigitizer::GetElectricField(G4double x, G4double y){ 

	int n_binx = eFieldMap->GetXaxis()->FindBin(x*1000); //position coordinates in um to return electric field in V/cm
	int n_biny = eFieldMap->GetYaxis()->FindBin(y*1000);
	G4double electricField = eFieldMap->GetBinContent(n_binx,n_biny);
	return electricField*1.0E-7; //return efield in MV/mm (for mobility calculation)
}

G4double AllPixFEI4RadDamage3DDigitizer::GetMobility(G4double electricField, G4bool isHoleBit){ 
	
	G4double vsat = 0;
	G4double ecrit = 0;
	G4double beta = 0;	

	//These parameterizations come from C. Jacoboni et al., Solid‐State Electronics 20 (1977) 77‐89. (see also https://cds.cern.ch/record/684187/files/indet-2001-004.pdf).
	if(!isHoleBit){
		vsat = 15.3*pow(temperature,-0.87);	// mm/ns
		ecrit = 1.01E-7*pow(temperature,1.55);	// MV/mm
		beta = 2.57E-2*pow(temperature,0.66);
	}
	if(isHoleBit){
		vsat = 1.62*pow(temperature,-0.52);	// mm/ns
		ecrit = 1.24E-7*(temperature,1.68);	// MV/mm
		beta = 0.46*pow(temperature,0.17);
	}

	G4double mobility = (vsat/ecrit)/pow(1+pow((electricField/ecrit),beta),(1/beta));	
	return mobility;  					// mm^2/(MV*ns)
}

G4double AllPixFEI4RadDamage3DDigitizer::GetDriftTime(G4bool isHoleBit){
	G4double u = CLHEP::RandFlat::shoot(0.,1.); // 
	G4double driftTime = 0;

	if(!isHoleBit) driftTime = (-1.)*trappingTimeElectrons*TMath::Log(u); // ns
	if(isHoleBit) driftTime = (-1.)*trappingTimeHoles*TMath::Log(u); // ns
	return driftTime;
}

G4double AllPixFEI4RadDamage3DDigitizer::GetTimeToElectrode(G4double x, G4double y, G4bool isHoleBit){
	// Uses (x,y) position in um to return time to electrode in ns
	G4double timeToElectrode = 0;
	if(!isHoleBit)
	  {
	    int n_binx = timeMap_e->GetXaxis()->FindBin(x*1000); //convert from mm to um
	    int n_biny = timeMap_e->GetYaxis()->FindBin(y*1000); 
	    timeToElectrode = timeMap_e->GetBinContent(n_binx,n_biny);
	  }
	if(isHoleBit)
	  {
	    int n_binx = timeMap_h->GetXaxis()->FindBin(x*1000);
	    int n_biny = timeMap_h->GetYaxis()->FindBin(y*1000);  
	    timeToElectrode = timeMap_h->GetBinContent(n_binx,n_biny);
	  }
	return timeToElectrode; //[ns]
}

G4double AllPixFEI4RadDamage3DDigitizer::GetTrappingPositionX(G4double initX, G4double initY, G4double driftTime, G4bool isHoleBit){
	G4double finalX=initX;
	if (!isHoleBit)	{
		int n_binx = xPositionMap_e->GetXaxis()->FindBin(initX*1000);
		int n_biny = xPositionMap_e->GetYaxis()->FindBin(initY*1000);
		int n_bint = xPositionMap_e->GetZaxis()->FindBin(driftTime);
		finalX = xPositionMap_e->GetBinContent(n_binx,n_biny,n_bint);
	} else {
		int n_binx = xPositionMap_h->GetXaxis()->FindBin(initX*1000);
		int n_biny = xPositionMap_h->GetYaxis()->FindBin(initY*1000);
		int n_bint = xPositionMap_h->GetZaxis()->FindBin(driftTime);
		finalX = xPositionMap_h->GetBinContent(n_binx,n_biny,n_bint);
	}
	return finalX*1e-3; //[mm]
}

G4double AllPixFEI4RadDamage3DDigitizer::GetTrappingPositionY(G4double initX, G4double initY, G4double driftTime, G4bool isHoleBit){
	G4double finalY=initY;
	if (!isHoleBit)	{
		int n_binx = yPositionMap_e->GetXaxis()->FindBin(initX*1000);
		int n_biny = yPositionMap_e->GetYaxis()->FindBin(initY*1000);
		int n_bint = yPositionMap_e->GetZaxis()->FindBin(driftTime);
		finalY = yPositionMap_e->GetBinContent(n_binx,n_biny,n_bint);
	} else {
		int n_binx = yPositionMap_h->GetXaxis()->FindBin(initX*1000);
		int n_biny = yPositionMap_h->GetYaxis()->FindBin(initY*1000);
		int n_bint = yPositionMap_h->GetZaxis()->FindBin(driftTime);
		finalY = yPositionMap_h->GetBinContent(n_binx,n_biny,n_bint);
	}
	return finalY*1e-3; //[mm]
}

G4double AllPixFEI4RadDamage3DDigitizer::GetRamoPotential(G4double x, G4double y){

	int n_binx = ramoPotentialMap->GetXaxis()->FindBin(x*1000);
	int n_biny = ramoPotentialMap->GetYaxis()->FindBin(y*1000);
	G4double ramoPotential = ramoPotentialMap->GetBinContent(n_binx,n_biny);
	return ramoPotential;
}

// ***** Process for each event *****
void AllPixFEI4RadDamage3DDigitizer::Digitize(){

	// create the digits collection
	m_digitsCollection = new AllPixFEI4RadDamage3DDigitsCollection("AllPixFEI4RadDamage3DDigitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > pixelsContent;		// stored energy per (countx,county) pixel	
	pair<G4int, G4int> tempPixel;					// (countx,county) which pixel

	// Loop over the whole Hits Collection - this spreads out the effect of the hit over its path through the pixel 
	G4int nEntries = hitsCollection->entries();
    //G4int nEntries = 13; //23;
    
    ofstream outfileE;
    outfileE.open("/allpix/allpix/outputdata/entries.txt", std::ofstream::app);
    outfileE << nEntries << "\n";
    outfileE.close();

	std::vector<pair<int, int>> hitPixelPair;

	double inducedCharge = 0;
	double depositedCharge = 0;
	//***temporary for cluster size calculations
	//double clusterThreshold = threshold;

	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		//Get pixel coordinates (note that here the origin has been adjusted to the bottom left corner NOT the middle of the pixel)
		G4double xpos = 110*1e-3;
		G4double ypos = 45*1e-3;
		//G4double zpos = 10*itr+10;//detectorThickness/2. - (*hitsCollection)[itr]->GetPosWithRespectToPixel().z();

		//G4double xpos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().x() + pitchX/2;
		//G4double ypos = (*hitsCollection)[itr]->GetPosWithRespectToPixel().y() + pitchY/2;

		G4double efield = GetElectricField(xpos,ypos);

		//G4double eHitTotal = (*hitsCollection)[itr]->GetEdep()/(1.0*eV);
		//G4double eHitTotal = (*hitsCollection)[itr]->GetEdep(); // Energy deposition for the hit, (MeV)
		G4double eHitTotal = 2.91*1e-3;//800*3.64;

	    	/*ofstream outfile;	   	
	    	outfile.open("/tmp/vewallan/protontrack.txt", std::ofstream::app);
	    	outfile << xpos*1e3 << "	" << ypos*1e3 << "	" << zpos*1e3 << "	" << eHitTotal << "\n";
	    	outfile.close();*/
	
		depositedCharge += eHitTotal; 

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();

		//if (tempPixel.first > 30 && tempPixel.first < 51 && tempPixel.second > 100 && tempPixel.second < 121) {
	    		/*ofstream outfile;	   	
	    		outfile.open("/tmp/vewallan/sampletrack_basic.txt", std::ofstream::app);
	    		outfile << ((tempPixel.first-1)*pitchX+xpos)*1e3 << "	" << ((tempPixel.second-1)*pitchY+ypos)*1e3 << "	" << zpos*1e3 << "	" << eHitTotal*1e-3 << "	" << tempPixel.first << "	" << tempPixel.second << "\n";
	    		outfile.close();*/
		//}

		// In case the charge moves into a neighboring pixel
		pair<G4int, G4int> extraPixel;
		extraPixel = tempPixel;

		//only process hits which are not on the electrodes (E-field zero)
		if (efield !=0) {

			// Split the charge into subcharges (# = precision) that are diffused separately to the electrode
			for(G4int nQ  = 0 ; nQ < precision ; nQ++) {
				
		  	    G4double eHit = G4double(eHitTotal)/precision; // energy of one electron-hole pair "chunk" [eV];

		  	    //Need to determine how many elementary charges this charge chunk represents.
			    double chunk_size = eHit/elec; //number of electrons/holes
			    //set minimum limit to prevent dividing into smaller subcharges than one fundamental charge
			    if (chunk_size < 1) chunk_size = 1;
			    double kappa = 1./sqrt(chunk_size);
		  
		  	    // Loop over everything twice: once for electrons and once for holes
		  	    for(G4int eholes=0 ; eholes<2 ; eholes++) {

				isHole = false; // Set a condition to keep track of electron/hole-specific functions
		    	if (eholes == 1) isHole = true;

		    	// Reset extraPixel coordinates each time through loop
		    	extraPixel = tempPixel;

		    	G4double timeToElectrode = GetTimeToElectrode(xpos,ypos,isHole); //ns
				if (timeToElectrode == 0) continue;
				G4double driftTime = GetDriftTime(isHole);
				//determine if subcharge trapped
				bool isTrapped = false;
				if (driftTime < timeToElectrode) isTrapped = true;

				//initialize diffusion coordinates
		    		G4double xposDiff=xpos;
		    		G4double yposDiff=ypos;

		    		if (doDiff) { 

                        // if subcharge trapped use drift time to determine diffusion distribution, else use time-to-electrode
                        G4double Dt = GetMobility(efield,isHole)*(0.024)*min(driftTime,timeToElectrode)*temperature/273.;
                        G4double rdif=sqrt(Dt)/1000; //in mm
                        xposDiff=xpos+rdif*CLHEP::RandGauss::shoot(0,1);
                        yposDiff=ypos+rdif*CLHEP::RandGauss::shoot(0,1);

		    			// Account for drifting into another pixel 
		    			while (xposDiff > pitchX){
		      				extraPixel.first = extraPixel.first + 1;               // increments or decrements pixel count in x
		      				xposDiff = xposDiff - pitchX;                              // moves xpos coordinate 1 pixel over in x
		    			}
		    			while (xposDiff < 0){
		      				extraPixel.first = extraPixel.first - 1;
		      				xposDiff = xposDiff + pitchX;
		    			}
		    			while (yposDiff > pitchY){
		      				extraPixel.second = extraPixel.second + 1;               // increments or decrements pixel count in y
		      				yposDiff = yposDiff - pitchY;                              // moves xpos coordinate 1 pixel over in y
		    			}
		    			while (yposDiff < 0){
		      				extraPixel.second = extraPixel.second - 1;
		      				yposDiff = yposDiff + pitchY;
		    			}


		    		} //doDiff

				G4double drift_time_constant = trappingTimeElectrons;
				if (isHole) drift_time_constant = trappingTimeHoles;
				G4double charge_correction_exp = exp(-timeToElectrode/drift_time_constant);
				G4double average_charge = avgChargeMap_e->GetBinContent(avgChargeMap_e->GetXaxis()->FindBin(xpos*1000),avgChargeMap_e->GetYaxis()->FindBin(ypos*1000));
				if (isHole) average_charge = avgChargeMap_h->GetBinContent(avgChargeMap_h->GetXaxis()->FindBin(xpos*1000),avgChargeMap_h->GetYaxis()->FindBin(ypos*1000));

	    			//for decoupling
		    		int locX = tempPixel.first;
                		int locY = tempPixel.second;

				G4double xposFinal = GetTrappingPositionX(xposDiff, yposDiff, min(driftTime,timeToElectrode), isHole);
				G4double yposFinal = GetTrappingPositionY(xposDiff, yposDiff, min(driftTime,timeToElectrode), isHole);
				//G4double xposFinal = GetTrappingPositionX(xpos, ypos, min(driftTime,timeToElectrode), isHole)+(xposDiff-xpos);
				//G4double yposFinal = GetTrappingPositionY(xpos, ypos, min(driftTime,timeToElectrode), isHole)+(yposDiff-ypos);

				if (doRamo){ //account for Ramo in hit pixel plus nearest neighbors

					// -- loop in the x-coordinate
					for (int i=-1; i<=1; i++){

			  			// -- loop in the y-coordinate
			  			for (int j=-1; j<=1; j++){
//int i =0;
//int j = 0;

							if (!doRamoNeighbors) {

								if (i!=0 || j!=0) continue; //only consider primary pixel, but exclude nearest neighbors		

							}
							
						extraPixel.first = locX + i;
			    			extraPixel.second = locY + j;

						//Ramo init different if charge diffused into neighboring pixel -> change primary pixel!!
						//IBL dimensions (Ramo map over 500umx350um pixel area)
                        G4double ramoInit = GetRamoPotential(xpos+pitchX*(0.5+i),ypos+pitchY*(3+j));
						G4double ramoFinal = GetRamoPotential(xposFinal+pitchX*(0.5+i),yposFinal+pitchY*(3+j));
						if (ITkGeometry) {//ITk dimensions (Ramo map over 200umx200um pixel area)
                            ramoInit = GetRamoPotential(xpos+pitchX*(1.5+i),ypos+pitchY*(1.5+j));
							ramoFinal = GetRamoPotential(xposFinal+pitchX*(1.5+i),yposFinal+pitchY*(1.5+j));
						}

/*ofstream outfile;	   	
	    outfile.open("/tmp/vewallan/positions.txt", std::ofstream::app);
	    outfile << (xpos+pitchX*(0.5+i))*1e3 << "	" << (ypos+pitchY*(3+j))*1e3 << "	" << (xposFinal+pitchX*(0.5+i))*1e3 << "	" << (yposFinal+pitchY*(3+j))*1e3 << "\n";
	    outfile.close();*/

			    			// Record deposit
						G4double eHitRamo = (1-2*isHole)*eHit*(ramoFinal - ramoInit);
						inducedCharge+=eHitRamo;

						if (!isHole) {
							ofstream outfile;
                            outfile.open("/allpix/allpix/outputdata/track_electrons_phi1e14_Diff.txt", std::ofstream::app);
                            outfile << xpos*1e3 << "	" << ypos*1e3 << "	" << xposFinal*1e3 << "	" << yposFinal*1e3 << "	" << eHitRamo/eHit << "\n";
                             outfile.close();
						} /*else {
							ofstream outfile1;	   	
	    						outfile1.open("/tmp/vewallan/track_holes.txt", std::ofstream::app);
	    						outfile1 << xpos*1e3 << "	" << xposFinal*1e3 << "	" << ypos*1e3 << "	" << yposFinal*1e3 << "	" << eHitRamo/eHit << "\n";
	    						outfile1.close();
						}*/


							if (doChunkCorrection) {

								eHitRamo = eHit*average_charge + kappa*(eHitRamo - eHit*average_charge);

							}

							/*if (!isHole) {
								ofstream outfile1;	   	
	    							outfile1.open("/tmp/vewallan/corrections_electrons.txt", std::ofstream::app);
	    							outfile1 << GetTimeToElectrode(xpos,ypos,false) << "	" << ((1-2*isHole)*eHit*(ramoFinal - ramoInit))/eHit << "	" << eHitRamo/eHit << "	" << average_charge << "\n";
	    						outfile1.close();
							} else {
								ofstream outfile2;	   	
	    							outfile2.open("/tmp/vewallan/corrections_holeTtE.txt", std::ofstream::app);
	    							outfile2 << GetTimeToElectrode(xpos,ypos,true) << "	" << ((1-2*isHole)*eHit*(ramoFinal - ramoInit))/eHit << "	" << eHitRamo/eHit << "	" << average_charge << "\n";
							}*/
		
							//*******temporary for debugging*******
							/*if (!isHole) {
								inducedCharge_e += eHit*(ramoFinal - ramoInit);
								inducedCharge_e_corr += eHitRamo;
							} else {
								inducedCharge_h += -eHit*(ramoFinal - ramoInit);
								inducedCharge_h_corr += eHitRamo;
							}*/
							//*************************************

							pixelsContent[extraPixel] += eHitRamo;


			  			} //loop over y

					} //loop over x	
				

		      		} else { //no Ramo


					/*if (!isHole) {

						G4double eHitCollected = 0;
					
						if (!isTrapped) eHitCollected = eHit;

						if (doChunkCorrection) {

							eHitCollected = eHit*(charge_correction_exp+kappa*(1.-charge_correction_exp));
							if (isTrapped) eHitCollected = eHit*(charge_correction_exp+kappa*(0.-charge_correction_exp));

						} 

						//*******temporary for debugging*******
						if (!isTrapped) collectedCharge += eHit;
						collectedCharge_corr += eHitCollected;
						//*************************************

						pixelsContent[extraPixel] += eHitCollected;

					}*/

				} //doRamo
		    
			  } // loop over holes/electrons

		  	} // end loop over nQ charges

	   	} //efield zero

	    		/*ofstream outfile;	   	
	    		outfile.open("/tmp/vewallan/sampletrack_basic.txt", std::ofstream::app);
	    		outfile << ((tempPixel.first-1)*pitchX+xpos)*1e3 << "	" << ((tempPixel.second-1)*pitchY+ypos)*1e3 << "	" << zpos*1e3 << "	" << eHitTotal*1e-3 << "	" << tempPixel.first << "	" << tempPixel.second << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second-1}]*1e-3 << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second}]*1e-3 << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second+1}]*1e-3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second-1}]*1e-3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second}]*1e-3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second+1}]*1e-3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second-1}]*1e-3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second}]*1e-3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second+1}]*1e-3 << "\n";
	    		outfile.close();*/

		/*if (tempPixel.first > 30 && tempPixel.first < 51 && tempPixel.second > 100 && tempPixel.second < 121) {
	    		ofstream outfile;	   	
	    		outfile.open("/tmp/vewallan/sampletrack.txt", std::ofstream::app);
	    		outfile << pixelsContent[{tempPixel.first-1,tempPixel.second-1}]*1e-3 << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second}]*1e-3 << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second+1}]*1e-3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second-1}]*1e-3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second}]*1e-3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second+1}]*1e-3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second-1}]*1e-3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second}]*1e-3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second+1}]*1e-3 << "\n";
	    		outfile.close();
		}*/

	    /*ofstream outfile;	   	
	    outfile.open("/tmp/vewallan/tracks50_theta0eta0.txt", std::ofstream::app);
	    outfile << xpos*1e3 << "	" << ypos*1e3 << "	" << zpos << "	" << eHitTotal*1e3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second+1}]*1e3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second}]*1e3 << "	" << pixelsContent[{tempPixel.first+1,tempPixel.second-1}]*1e3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second+1}]*1e3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second}]*1e3 << "	" << pixelsContent[{tempPixel.first,tempPixel.second-1}]*1e3 << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second+1}]*1e3 << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second}]*1e3 << "	" << pixelsContent[{tempPixel.first-1,tempPixel.second-1}]*1e3 << "\n";
	    outfile.close();*/

	} // end loop over nEntries

	    /*int ToTValue = TMath::FloorNint(inducedCharge*tuning);
	    if (ToTValue >=15) ToTValue = 15;
	    if (inducedCharge < threshold) ToTValue = 0;*/

	    /*ofstream outfile;	   	
	    outfile.open("/tmp/vewallan/3DITk_ToT_phi6e15_Vbias140.txt", std::ofstream::app);
	    outfile << inducedCharge << "	" << depositedCharge << "	" << ToTValue << "\n";
	    outfile.close();*/

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

	// Now that pixelContent is filled, create one digit per pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

	//***temporary for cluster properties
	G4int clusterSize = 0; 
	G4double clusterCharge = 0; 
	G4double clusterToT = 0;
	//***
	
	for( ; pCItr != pixelsContent.end() ; pCItr++) {
	    G4double deposited_energy = (*pCItr).second;
	    int ToT = TMath::FloorNint(deposited_energy*tuning);
	    if (ToT >=15) ToT = 15; //FEI4 is 4-bit.
	    AllPixFEI4RadDamage3DDigit * digit = new AllPixFEI4RadDamage3DDigit;
	    digit->SetPixelIDX((*pCItr).first.first);
	    digit->SetPixelIDY((*pCItr).first.second);
	    digit->SetPixelCounts(ToT);
	    if (deposited_energy < threshold){
	      digit->SetPixelCounts(0);
	    }
	    if (deposited_energy >= threshold && ToT > 0){
	      m_digitsCollection->insert(digit);
	      //***temporary for cluster properties
	      clusterCharge += deposited_energy;
	      clusterToT += ToT;
	      clusterSize++;
	      //***
	      //***to visualize clusters
	      /*ofstream outfile1;	 
	      outfile1.open("/tmp/vewallan/pixelIDcharge_phi1e16_theta40_phi20_1000evts.txt", std::ofstream::app); //theta = 90 degrees - traditional theta
	      outfile1 << (*pCItr).first.first << " " << (*pCItr).first.second << "	" << deposited_energy << "	" << ToT << "\n";
	      outfile1.close();*/
	    }
	}
	//***

	//***temporary for cluster properties
	/*if (clusterSize > 0) {
		ofstream outfile;	 
		outfile.open("/tmp/vewallan/clusterInfo_phi1e16_theta80_10000evts.txt", std::ofstream::app); //theta = 90 degrees - traditional theta
		outfile << clusterSize << " " << clusterCharge << "	" << clusterToT << "\n";
		outfile.close();
	}*/
	//***
	
	G4int dc_entries = m_digitsCollection->entries();
	if(dc_entries > 0)
	  {
	    G4cout << "--------> Digits Collection : " << collectionName[0]
		   << "(" << m_hitsColName[0] << ")"
		   << " contains " << dc_entries
		   << " digits" << G4endl;
	  }
	
	StoreDigiCollection(m_digitsCollection);
	
} // end Digitize function
