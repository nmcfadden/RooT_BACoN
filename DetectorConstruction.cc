//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
// use of stepping action to set the accounting volume
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"

#include "G4Torus.hh"
#include "G4UnionSolid.hh"

#include "G4String.hh"
#include "math.h"

#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

using namespace std;
//using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const G4double DetectorConstruction::LambdaE = 2.0*TMath::Pi()*1.973269602e-16 * m * GeV;
DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
  checkOverlaps = true;
  //***//
  //TPB//
  //***//    
  G4String pathFile = "External_data/tpbGhemann.root";
  TFile *tpbFile = new TFile(pathFile.data());
  if (!tpbFile ) 
    G4cout<<" DetectorConstruction ERROR:: file " << pathFile << " not found " << G4endl;
  else
    G4cout<<" DetectorConstruction INFO:: file " << pathFile << " opened " << G4endl;
  fTPBspec=NULL;
  tpbFile->GetObject("tpbGhemann",fTPBspec);
  if (!fTPBspec ) 
    G4cout<<" DetectorConstruction ERROR:: not graph tpbBhemann in file " << pathFile <<G4endl;
  else 
    G4cout<<" DetectorConstruction info tpbBhemann graph found " <<G4endl;
    
  // create directory 
  fDir = LegendAnalysis::Instance()->topDir()->mkdir("detec");
  fDir->cd();
  G4cout<<" DetectorAction working root directory  is  " << G4endl;  
  gDirectory->pwd();
  G4cout << " ... " << G4endl;
  G4double LowE =1.4*eV;//885.6013 2.4796*eV;//500 nm
  G4double HighE = 12.3984*eV;//100 nm
  G4double ArHighE = LambdaE /(115*nanometer);
  G4double ArLowE = LambdaE /(650*nanometer);
  hWLSPhotonE = new TH1F("WLSPhotonE"," photon energy in WLS",1000,LowE,HighE);
  hWLSPhotonWavelength = new TH1F("WLSPhotonWavelength"," photon Wavelength in WLS",1000,LambdaE/HighE,LambdaE/LowE);
  hArPhotonE = new TH1F("ArPhotonE"," photon energy in LAr",100,ArLowE,ArHighE);
  hArPhotonWavelength = new TH1F("ArPhotonWavelength"," photon Wavelength in LAr",1000,LambdaE/ArHighE,LambdaE/ArLowE);  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  nist = G4NistManager::Instance();
  innerVessel_FillMaterial = "ArgonLiquid";//"NitrogenGas";
  // used in this include 
	#include "LegendDetectorMaterials.icc"
  ArgonOpticalProperties();
  //mat_fill->GetMaterialPropertiesTable()->DumpTable();
  
  WLSOpticalProperties();
  
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4Material* GeMaterial = nist->FindOrBuildMaterial("G4_Ge");

  ////////////////////////////////////////////////////////////////////////////////////////
	//
  // World
  //
  ////////////////////////////////////////////////////////////////////////////////////////
  solidWorld = new G4Box("sol_World",50*m,50*m,50*m);
  logicalWorld = new G4LogicalVolume(solidWorld,mat_air,"log_World");
	logicalWorld->SetVisAttributes (G4VisAttributes::Invisible);
  physicalWorld = new G4PVPlacement(0,G4ThreeVector(),logicalWorld,"phy_World",0,false,0,checkOverlaps);

  /*
	G4Box* solid_Rock = new G4Box("sol_Rock",50*m,50*m,30*m);
	G4Box* solid_Lab = new G4Box("sol_Lab",35*m,10*m,4*m);
	G4SubtractionSolid *solid_Rock2 = new G4SubtractionSolid("sol_Rock2", solid_Rock, solid_Lab ,0 , G4ThreeVector(-25*m,0,10.5*m));
	G4Tubs* solid_CutOut = new G4Tubs("sol_CutOut",0, 6.50001*m ,6.50001*m, 0, 2*M_PI);
	G4SubtractionSolid *solid_Rock3 = new G4SubtractionSolid("sol_Rock2", solid_Rock2, solid_CutOut ,0 , G4ThreeVector(0,0,0));


  G4LogicalVolume* logical_Rock = new G4LogicalVolume(solid_Rock3,mat_Rock,"log_Rock");
	logical_Rock->SetVisAttributes ( new G4VisAttributes(G4Colour(0.1,0.1,0.7,0.5) ) );//(0.7, 0.7, 0.7, 0.5) )); //grey 50% transparent
  new G4PVPlacement(0,G4ThreeVector(),logical_Rock,"phy_Rock",logicalWorld,false,0,checkOverlaps);
  */

  vector<string> Detectors;
  string input;
  string s;

  G4int Detector_number;
  G4int Detector_slices;
  vector<G4double> Detector_r;
  vector<G4double> Detector_z;

  G4double sum_x1 = 0;
  G4double sum_y1 = 0;
  G4double sum_z1 = 0;
  G4double sum_x2 = 0;
  G4double sum_y2 = 0;
  G4double sum_z2 = 0;
  G4double group1zmax= 0;
  G4double group1rmax= 0;
  G4double group2zmax= 0;
  G4double group2rmax= 0;


  int j =0;

  std::vector<G4ThreeVector> detRout;
  std::vector<G4ThreeVector> detZhalf;
  std::vector<G4ThreeVector> detPositions;
  std::vector<G4int> detNumbers;
  std::vector<G4String> detNames;
  std::vector<G4String> detNamesPhys;
  std::vector<G4String> detNamesLog;

  ifstream infile("Detectorposition.txt",ios_base::in);
  if (! infile.is_open()) G4cout<< "DetectorConstruction unable to open Detectorposition.txt" << G4endl;
  else  G4cout<< "DetectorConstruction opened Detectorposition.txt" << G4endl;

  while(getline(infile,input)){
    Detectors.push_back(input);
  }
  infile.close();
  G4cout << "DetectorConstruction Detector size = " << Detectors.size() <<G4endl;
  char mess[200]; 
  G4int nhalfDet = Detectors.size()/2;

  // here I do the detector placement and some approximate calculations 
  // for the surrounding cylinders. (MG)

  for (int i =0;i< Detectors.size() ;i++){
    G4String Detector_name;
    G4String Detector_name_sol;
    G4String Detector_name_log;
    
    //G4cout << " i= " << i << " Detectors " << Detectors[i] << endl;
    j = 0;
    Detector_number = 0;
    Detector_r.clear();
    Detector_z.clear();
    G4ThreeVector Detector_position;
    stringstream ss(Detectors[i]);
    while (ss >> s) {
     j++;
     if (j==1) Detector_number += atoi(s.c_str())*100;
     else if (j==2) Detector_number += atoi(s.c_str())*10;
     else if (j==3) Detector_number += atoi(s.c_str())*1;
     else if (j==4) Detector_position.setX(atof(s.c_str())*mm);
     else if (j==5) Detector_position.setY(atof(s.c_str())*mm);
     else if (j==6) Detector_position.setZ(atof(s.c_str())*mm);
     else if (j==7) Detector_name = s;
     else if (j==8) Detector_slices = atoi(s.c_str());
     else if (j==9) Detector_r.push_back((G4double) atof(s.c_str())*mm);
     else if (j==10) Detector_z.push_back((G4double) atof(s.c_str())*mm);
     else if (j==11) Detector_r.push_back((G4double) atof(s.c_str())*mm);
     else if (j==12) Detector_z.push_back((G4double) atof(s.c_str())*mm);
     else if (j==13) Detector_r.push_back((G4double) atof(s.c_str())*mm);
     else if (j==14) Detector_z.push_back((G4double) atof(s.c_str())*mm);

    }
    //G4cout << "\t " << i << " j=" << j << " number " << Detector_number << " position " << Detector_position << " slices " << Detector_slices << G4endl;         
    Detector_name_sol = Detector_name + "_sol";
    Detector_name_log = Detector_name + "_log";
    if(Detector_slices==2){
      Detector_r.push_back(Detector_r[1]+0.00001*mm);
      Detector_z.push_back(Detector_z[1]+0.00001*mm);
    }
    detNumbers.push_back(Detector_number);
    detNames.push_back(Detector_name);
    detNamesPhys.push_back(Detector_name_sol);
    detNamesLog.push_back(Detector_name_log);

    G4ThreeVector r_i(0,0,0);
    G4ThreeVector r(Detector_r[0],Detector_r[1],Detector_r[2]);
    G4ThreeVector z(Detector_z[0],Detector_z[1],Detector_z[2]);

    /*
    sprintf(mess," placement of detector %i number %i at (%f,%f,%f) \n",i,Detector_number,Detector_position.x(),Detector_position.y(),Detector_position.z());
    G4cout<< mess;
    sprintf(mess," zPlane %f %f %f \n",z[0],z[1],z[2]);
    G4cout<< mess;
    sprintf(mess," r in   %f %f %f \n",r_i[0],r_i[1],r_i[2]);
    G4cout<< mess;
    sprintf(mess," r out  %f %f %f \n",r[0],r[1],r[2]);
    G4cout<< mess;
    */

    // mean positions and maximum sizes 
    G4double maxz = max(-z[0],z[2]);
    if(i<nhalfDet) {
      sum_x1 += Detector_position.x() /G4double(nhalfDet);
      sum_y1 += Detector_position.y() /G4double(nhalfDet);
      sum_z1 += Detector_position.z() /G4double(nhalfDet);
      if(maxz>group1zmax) group1zmax=maxz;
      for(int ir=0 ; ir<3 ; ++ir) if(r[ir]>group1rmax) group1rmax=r[ir];
    } else {
     sum_x2 += Detector_position.x() /G4double(nhalfDet);
     sum_y2 += Detector_position.y() /G4double(nhalfDet);
     sum_z2 += Detector_position.z() /G4double(nhalfDet);
     if(maxz>group2zmax) group2zmax=maxz;
     for(int ir=0 ; ir<3 ; ++ir) if(r[ir]>group2rmax) group2rmax=r[ir];
    }
    detPositions.push_back(Detector_position);
    detRout.push_back(r);
    detZhalf.push_back(z);
  }

 
  G4ThreeVector group1pos(sum_x1,sum_y1,sum_z1);
  G4ThreeVector group2pos(sum_x2,sum_y2,sum_z2);
  
  // calculate cylinder radii.
  G4double maxRho1 =0; 
  G4double maxZee1 =0; 
  G4double maxRho2 =0; 
  G4double maxZee2 =0; 
  
  std::vector<G4ThreeVector> detRelPositions;
  for(unsigned idet=0; idet < detPositions.size(); ++ idet) {
    G4ThreeVector rel;
    if(idet<detPositions.size()/2) rel = detPositions[idet] - group1pos;
    else rel = detPositions[idet] - group2pos;
    detRelPositions.push_back(rel);
    G4double rho = rel.getRho();
    if(idet<detPositions.size()/2) {
      if(rho>maxRho1) maxRho1=rho;
      if(rel.z()>maxZee1) maxZee1 = rel.z();
    } else {
      if(rho>maxRho2) maxRho2=rho;
      if(rel.z()>maxZee2) maxZee2 = rel.z();
    }
  }
  group1rmax  += maxRho1;
  group1zmax  += maxZee1;
  // inflate by 10%
  group1zmax *= 1.1;
  group1rmax *= 1.1;

  group2rmax  += maxRho2;
  group2zmax  += maxZee2;
  // inflate by 2%
  group2zmax *= 1.1;
  group2rmax *= 1.1;

  // make both tubes the same
  grouprmax = max(group1rmax,group2rmax);
  groupzmax = max(group1zmax,group2zmax);



  G4cout<< "\t *******************************************" <<G4endl;
  G4cout<< "\t DetectorConstruction -- placement info " <<G4endl;
  sprintf(mess,"\t position group 1  (%f,%f,%f) rmax %f  zmax %f\n",sum_x1,sum_y1,sum_z1,group1rmax,group1zmax); G4cout<< mess;
  sprintf(mess,"\t position group 2  (%f,%f,%f) rmax %f  zmax %f\n",sum_x2,sum_y2,sum_z2,group2rmax,group2zmax); G4cout<< mess;
  sprintf(mess,"\t group radius  %f  half z %f \n",grouprmax,groupzmax); G4cout<< mess;
  G4cout<< "\t *******************************************" <<G4endl;

  //First Place World volume of LAr...different than LAr cylinders
  G4int sourceRadius = 4*groupzmax;
  G4Sphere* larSolid = new G4Sphere("source",0,sourceRadius,0,2*M_PI,0,2*M_PI);
  
  //mat_fill defined in DetectorConstruction.icc from innerVessel_FillMaterial or command file
  larSourceLogical = new G4LogicalVolume(larSolid,mat_fill,"source_log"); 
  larSourceLogical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1, 0.1, 0.9,0.5)));
  
  G4ThreeVector source_center(0,0,0);
  larPhysical = new G4PVPlacement (0,source_center,larSourceLogical,"larPhysical",logicalWorld,false,0,checkOverlaps);

  /*
  // construct cylindrical tubes surrounding detectors
  // G4Tubs(Name,RMin,RMax,Half-z,StartPhi,DeltaPhi)
  */
  //G4double tubeWall = 1*mm;
  //G4Tubs* groupTube = new G4Tubs("groupTube",0,grouprmax,grouprmax+tubeWall,groupzmax,twopi);
  
  //G4LogicalVolume* group1Logical = new G4LogicalVolume(groupTube,fTPB,"group1Logical" );
  //G4LogicalVolume* group2Logical = new G4LogicalVolume(groupTube,fTPB,"group2Logical" );
  
  //Added by Neil, previous declaration set starting angle at groupzmax
  G4Tubs* groupTube = new G4Tubs("groupTube",0,grouprmax,groupzmax,0,twopi);
  //Added by Neil, Note the change in material...Previous the detectors were submerged in TPB not LAr!!! 
  G4LogicalVolume* group1Logical = new G4LogicalVolume(groupTube,mat_fill,"group1Logical" );
  G4LogicalVolume* group2Logical = new G4LogicalVolume(groupTube,mat_fill,"group2Logical" );

  group1Logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1,0.5,0.7)));
  group2Logical->SetVisAttributes(new G4VisAttributes(G4Colour(0.1,0.5,0.7)));

  
  // place detectors in tubes
  for(unsigned idet=0; idet < detPositions.size(); ++ idet) {
    G4ThreeVector r_i(0,0,0);
    G4ThreeVector r=detRout[idet];
    G4ThreeVector z=detZhalf[idet];
    G4Polycone* Det_solid = new G4Polycone(detNamesPhys[idet],0,2*M_PI,3,&z[0],&r_i[0],&r[0]);
    G4LogicalVolume* Det_logical = new G4LogicalVolume(Det_solid,GeMaterial,detNamesLog[idet]);
    if(idet<detPositions.size()/2) 
      new G4PVPlacement (0,detRelPositions[idet],Det_logical,detNames[idet],group1Logical,false,detNumbers[idet],checkOverlaps);
     else 
      new G4PVPlacement (0,detRelPositions[idet],Det_logical,detNames[idet],group2Logical,false,detNumbers[idet],checkOverlaps);     
  }
    
  G4VPhysicalVolume* group1Physical = new G4PVPlacement(0,group1pos,group1Logical,"group1Physical",larSourceLogical,false,0,checkOverlaps);
  G4VPhysicalVolume* group2Physical = new G4PVPlacement(0,group2pos,group2Logical,"group2Physical",larSourceLogical,false,1,checkOverlaps);

  /////////////WLS Cylinder around groupTubs, does not cover top or bottom, PMTs will do that below/////////////
  WLSHalfThickness = 0.05*mm;  // half thickness
  G4Tubs* WLSgroupTube = new G4Tubs("PMTgroupTube",grouprmax,grouprmax+WLSHalfThickness,groupzmax,0,twopi);
  //I think we need only one logical volume and just place it twice. There might be a nuance to this when using Sensitive Detector
  G4LogicalVolume* WLSgroupTubeLogical = new G4LogicalVolume(WLSgroupTube,fTPB,"WLSgroup1TubeLogical"); 
  
  WLSgroupTubeLogical->SetVisAttributes ( new G4VisAttributes(G4Colour(0.6,0.1,0.7) ) );
  
  //Create a physical placement to create rough surfaces for tpb/LAr interface
  G4double roughness = 0.5;
  G4OpticalSurface* WLSoptSurf = new G4OpticalSurface("WLS_rough_surf",glisur,ground,dielectric_dielectric,roughness);

  G4PVPlacement* phys_PMTWLS =  new G4PVPlacement (0,group1pos,WLSgroupTubeLogical,"WLSgroup1TubeLogical", larSourceLogical,false,0,checkOverlaps);
  new G4LogicalBorderSurface("Phys_PMT_WLS_cylinder_1",group1Physical,phys_PMTWLS,WLSoptSurf);
  new G4LogicalBorderSurface("Phys_PMT_WLS_cylinder_1",phys_PMTWLS,group1Physical,WLSoptSurf);

  phys_PMTWLS = new G4PVPlacement (0,group2pos,WLSgroupTubeLogical,"WLSgroup2TubeLogical", larSourceLogical,false,1,checkOverlaps);  
  new G4LogicalBorderSurface("Phys_PMT_WLS_cylinder_2",group2Physical,phys_PMTWLS,WLSoptSurf);
  new G4LogicalBorderSurface("Phys_PMT_WLS_cylinder_2",phys_PMTWLS,group2Physical,WLSoptSurf);
  

  /////////////PMT coated in WLS/////////////
  //WLSHalfThickness defined above
  G4double glassHalfThickness = 10*mm;  // half thickness
  G4double housingHalfThickness = 10*mm;  // half thickness
  G4Tubs* PMTDiskTubs = new G4Tubs("PMTDiskTubs",0.,grouprmax,housingHalfThickness,0,twopi);
  
  //Metal housing, Kovar is a Ni Co alloy
  G4Material* materialPMTHousing = G4Material::GetMaterial("Kovar");
  logicalPmtHousing = new G4LogicalVolume(PMTDiskTubs,materialPMTHousing,"logicalPmtHousing");  
  //logicalPmtHousing->SetVisAttributes ( new G4VisAttributes(G4Colour(0.9,0.1,0.1) ) );
  
  G4Tubs* PMTGlassTubs = new G4Tubs("PMTGlassTubs",0,grouprmax,glassHalfThickness,0,twopi);
  G4Material* materialPMTGlass = G4Material::GetMaterial("Quartz"); 
  logicalPmtGlass = new G4LogicalVolume(PMTGlassTubs,materialPMTGlass,"logicalPmtGlass");            
  //logicalPmtGlass->SetVisAttributes ( new G4VisAttributes(G4Colour(0.1,0.9,0.1) ) );
  
  G4Tubs* PMTWlsTubs = new G4Tubs("PMTWlsTubs",0,grouprmax,WLSHalfThickness,0,twopi);
  
  logicalPMTWLS = new G4LogicalVolume(PMTWlsTubs,fTPB,"logicalPmtGlassWLS");   
  logicalPMTWLS->SetVisAttributes ( new G4VisAttributes(G4Colour(0.6,0.1,0.7) ) );
  //fPMTGlassOptSurface defined in LegendDetectorMaterials.icc
  new G4LogicalSkinSurface("PMTGlass_surf",logicalPmtGlass,fPMTGlassOptSurface);
  
  G4cout<< "\t *******************************************" <<G4endl;
  // construcnt and put into larSourceLogical needs larPhysical for LogicalBorderSuface
  G4double pmtZOffset = 2*(glassHalfThickness+WLSHalfThickness)+housingHalfThickness;//+tubeWall;
  G4ThreeVector rPMT1top = group1pos + G4ThreeVector(0,0,groupzmax+pmtZOffset);
  G4ThreeVector rPMT1bot = group1pos + G4ThreeVector(0,0,-groupzmax-pmtZOffset);
  G4ThreeVector rPMT2top = group2pos + G4ThreeVector(0,0,groupzmax+pmtZOffset);
  G4ThreeVector rPMT2bot = group2pos + G4ThreeVector(0,0,-groupzmax-pmtZOffset);
  PlacePMT(rPMT1top,-1,0);
  PlacePMT(rPMT1bot,1,1);
  PlacePMT(rPMT2top,-1,2);
  PlacePMT(rPMT2bot,1,3);

  G4cout << "\t done constructing detector  " << G4endl;
  return physicalWorld;
}


void DetectorConstruction::PlacePMT(G4ThreeVector rhousing,double top_or_bot,int num)
{
  //PMT Housing 
  new G4PVPlacement(0,rhousing,logicalPmtHousing,"phys_PMTHousing_"+std::to_string(num),larSourceLogical,false,num,checkOverlaps);
  //PMT GLASS
  G4double zhalf = dynamic_cast<G4Tubs*>(logicalPmtHousing->GetSolid())->GetZHalfLength() +  
    dynamic_cast<G4Tubs*>(logicalPmtGlass->GetSolid())->GetZHalfLength();
  zhalf *= top_or_bot;
  G4ThreeVector rglass = rhousing + G4ThreeVector(0,0,zhalf);
  new G4PVPlacement(0,rglass,logicalPmtGlass,"phys_PMTGlass_"+std::to_string(num),larSourceLogical,false,num,checkOverlaps);
  //WLS PMT 
  zhalf += top_or_bot*( dynamic_cast<G4Tubs*>(logicalPmtGlass->GetSolid())->GetZHalfLength()
      +dynamic_cast<G4Tubs*>(logicalPMTWLS->GetSolid())->GetZHalfLength());
  G4ThreeVector rwls = rhousing + G4ThreeVector(0,0,zhalf);
  G4PVPlacement* phys_PMTWLS = new G4PVPlacement(0,rwls,logicalPMTWLS,"phys_WLSGlassPmt_"+std::to_string(num),larSourceLogical,false,num,checkOverlaps);
  /* border between WLS and larPhysicsl crossing in both directions WLS is rough on both sides */
  G4double roughness = 0.5;
  G4OpticalSurface* WLSoptSurf = new G4OpticalSurface("WLS_rough_surf",glisur,ground,dielectric_dielectric,roughness);
  new G4LogicalBorderSurface("Phys_PMT_WLS_"+ std::to_string(num),larPhysical,phys_PMTWLS,WLSoptSurf);
  new G4LogicalBorderSurface("Phys_PMT_WLS_"+ std::to_string(num),phys_PMTWLS,larPhysical,WLSoptSurf);
}

//
/// optical properties of lar in several places
void DetectorConstruction::ArgonOpticalProperties()
{
	//Taken from home/gold/MaGe-master/munichteststand/src/GEMPIKLArGe.cc
  G4int num = 500;
  static const G4double temp = 88.5*kelvin;
 // static const G4double LambdaE = twopi *1.973269602e-16 * m * GeV;
  G4double scint_yield = 28120./MeV;//40000./MeV; //sono 40000

  G4int ji;
  G4double e;
  G4double ee;

  G4double PPCKOVHighE = LambdaE / (115*nanometer);
  G4double PPCKOVLowE = LambdaE / (650*nanometer); 
  G4double de = ((PPCKOVHighE - PPCKOVLowE) / ((G4double)(num-1)));
  // liquid argon (LAr)  
  G4double LAr_PPCK[(num)];
  G4double LAr_RIND[(num)];
  G4double LAr_RAYL[(num)];
  G4double LAr_ABSL[(num)];
  
  G4double lar_absl_xuv = 60*cm;
	G4double lar_absl_vis = 1000*m;
 
  for (ji = 0; ji < num; ji++) {
    e = PPCKOVLowE + ((G4double)ji) * de;
    LAr_PPCK[ji] = e;
    LAr_RIND[ji] = LArRefIndex((LambdaE / e));
    LAr_RAYL[ji] = LArRayLength((LambdaE / e), temp);

    if (((LambdaE / e)/nm) < 200.0) {
      LAr_ABSL[ji] = lar_absl_xuv;
    } else {
      LAr_ABSL[ji] = lar_absl_vis;
    }
  }

  G4double PPSCHighE = LambdaE /(115*nanometer);
  G4double PPSCLowE = LambdaE /(650*nanometer);
  G4double dee = ((PPSCHighE - PPSCLowE) / ((G4double)(num-1)));
  G4double LAr_SCIN[num];
  G4double LAr_SCPP[num];
  for (ji = 0; ji < num; ji++) {
    ee=PPSCLowE+ ((G4double)ji) * dee;
    LAr_SCPP[ji]=ee;
    LAr_SCIN[ji]=ArScintillationSpectrum(ee);
    // plot spectrum
    hArPhotonE->SetBinContent(  hArPhotonE->FindBin(ee), LAr_SCIN[ji]);
    hArPhotonWavelength->SetBinContent(hArPhotonWavelength->FindBin(LambdaE/ee) ,LAr_SCIN[ji]);
  }


  G4MaterialPropertiesTable* LAr_mt = new G4MaterialPropertiesTable();

  LAr_mt->AddProperty("RINDEX",        LAr_PPCK, LAr_RIND, num);
  LAr_mt->AddProperty("RAYLEIGH",      LAr_PPCK, LAr_RAYL, num);
  LAr_mt->AddProperty("ABSLENGTH",     LAr_PPCK, LAr_ABSL, num);
  if ( (LAr_SCPP[0] >= PPCKOVLowE) &&
       (LAr_SCPP[(sizeof(LAr_SCPP)/sizeof(G4double) - 1)] <= PPCKOVHighE) )
    {
      LAr_mt->AddProperty("FASTCOMPONENT",LAr_SCPP,LAr_SCIN,num);
      LAr_mt->AddProperty("SLOWCOMPONENT",LAr_SCPP,LAr_SCIN,num);
    }
  LAr_mt->AddConstProperty("SCINTILLATIONYIELD",scint_yield);
  LAr_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  LAr_mt->AddConstProperty("FASTTIMECONSTANT", 5.95*ns);//6.*ns);
  LAr_mt->AddConstProperty("SLOWTIMECONSTANT",922*ns);//1590.*ns);
  LAr_mt->AddConstProperty("YIELDRATIO",0.23);
  // G4cout<<LAr_SCIN<<G4endl;
  G4double fano = 0.11;
  LAr_mt->AddConstProperty("RESOLUTIONSCALE",fano); 
  mat_ArLiq->SetMaterialPropertiesTable(LAr_mt); // G4Material defined in Detector_Materials.icc
  mat_ArLiq->GetIonisation()->SetBirksConstant(5.1748e-4*cm/MeV);//0.0144*mm/MeV);
}

G4double DetectorConstruction::LArEpsilon(const G4double lambda)
{
  G4double epsilon;
  if (lambda < 110*nanometer) return 1.0e4; // lambda MUST be > 110.0 nm
  epsilon = lambda / micrometer; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  G4double LArRho = 1.396*g/cm3;
  G4double GArRho = 1.66e-03*g/cm3;
  epsilon *= (LArRho / GArRho); // density correction (Ar gas -> LAr liquid)
  if (epsilon < 0.0 || epsilon > 0.999999) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti
  return epsilon;
}

G4double DetectorConstruction::LArRefIndex(const G4double lambda)
{
 //G4cout<< ( sqrt(LArEpsilon(lambda)))<<G4endl;
 return ( sqrt(LArEpsilon(lambda)) ); // square root of dielectric constant
}
G4double DetectorConstruction::LArRayLength(const G4double lambda,const
				   G4double temp)
{
  G4double dyne = 1.0e-5*newton;
  static const G4double LArKT = 2.18e-10 * cm2/dyne; // LAr isothermal compressibility
  static const G4double k = 1.380658e-23 * joule/kelvin; // the Boltzmann constant
  G4double h;
  h = LArEpsilon(lambda);
  if (h < 1.00000001) h = 1.00000001; // just a precaution
  h = (h - 1.0) * (h + 2.0); // the "dielectric constant" dependance
  h *= h; // take the square
  h *= LArKT * temp * k; // compressibility * temp * Boltzmann constant
  h /= lambda * lambda * lambda * lambda; // (lambda)^4
  h *= 9.18704494231105429; // (2 * Pi / 3)^3
  if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a precaution
  if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 * nanometer); // just a precaution
  return ( 1.0 / h );
}

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager* SDman   =  G4SDManager::GetSDMpointer();  
  PMTSD* sd = new PMTSD("PhotoCathode",1,"PhCathodeHC" );    
  SDman->AddNewDetector(sd); 
  logicalPmtGlass->SetSensitiveDetector(sd);
}

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetOverlapsCheck(G4bool f_check)
{
	//checkOverlaps = f_check;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetFillMaterial(G4String smaterial)
{
	innerVessel_FillMaterial = smaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetShieldStyle(G4String f_type)
{
	//detector_type = f_type;
}


void DetectorConstruction::WLSOpticalProperties()
{
  // add WLS 
  // -- WLS: TPB (Tetraphenyl butadiene)
  // --M.Gold from Gehmann et al plot

  fTPB = G4Material::GetMaterial("TPB", false);
  if (fTPB == 0) {
    G4Element* elementH = nist->FindOrBuildElement("H");
    G4Element* elementC = nist->FindOrBuildElement("C");
    fTPB= new G4Material("TPB", 1*g/cm3, 2, kStateSolid);
    fTPB->AddElement(elementH, 22);
    fTPB->AddElement(elementC, 28);
  }
   // Now attach the optical properties to it.
   // Build table with photon energies
   
   const G4int numTPB =100;// 63;;
   G4double HighETPB = LambdaE /(115*nanometer);
   G4double LowETPB = LambdaE /(650*nanometer);//(650*nanometer); //598
   G4double deeTPB = ((HighETPB - LowETPB) / ((G4double)(numTPB-1)));
   G4double LAr_SCPPTPB[numTPB];
   for (G4int ji = 0; ji < numTPB; ji++)  LAr_SCPPTPB[ji]=LowETPB+ ((G4double)ji) * deeTPB;
   
   G4double WLS_absorption[numTPB];
   G4double WLS_emission[numTPB];
   G4double Refraction[numTPB];
   
   // make new table 
   tpbTable = new G4MaterialPropertiesTable();
   for (G4int ji=0;ji < numTPB; ji++) {
     Refraction[ji] = 1.6; //this is just a guess
     // Should the TPB shift the Cherenkov light?
     // This makes a tail starting at 128 until the visible.
     if (LAr_SCPPTPB[ji] > 3.31*eV){// < 374.57 nm 
       // For the moment set it to always absorb photons
       WLS_absorption[ji] = 0.001*nm; //absorbs UV (always)
     } else {
       // < 350 nm
       WLS_absorption[ji] = 10.5*m; //otherwise transparent
     }
     WLS_emission[ji] = TPBEmissionSpectrum(LAr_SCPPTPB[ji]);
     hWLSPhotonE->SetBinContent(ji,WLS_emission[ji]);
     hWLSPhotonWavelength->SetBinContent(numTPB-1-ji,WLS_emission[ji]);
     //G4cout<<" WLS Emmsion "<<WLS_emission[ji]<<" LAr Energy "<<LAr_SCPPTPB[ji]<<G4endl;
     //G4cout<<" WLS Absorption Length "<<WLS_absorption[ji]<<" LAr Energy "<<LAr_SCPPTPB[ji]<<G4endl;
   }

   G4cout << " to here " << G4endl;
   tpbTable->AddProperty("RINDEX",LAr_SCPPTPB,Refraction,numTPB);
   tpbTable->AddProperty("WLSABSLENGTH",LAr_SCPPTPB,WLS_absorption,numTPB);
   tpbTable->AddProperty("WLSCOMPONENT",LAr_SCPPTPB,WLS_emission,numTPB);
   // From WArP
   tpbTable->AddConstProperty("WLSTIMECONSTANT", 0.01*ns);
   G4double WLSyield = 1.2;// should be integral of TGraph!MG!!
   tpbTable->AddConstProperty("WLSMEANNUMBERPHOTONS",WLSyield);
   fTPB->SetMaterialPropertiesTable(tpbTable);
}

// arg is energy , returns probability 
//  FWHM at 80K from J Chem Phys vol 91 (1989) 1469 E Morikawa et al
G4double DetectorConstruction::ArScintillationSpectrum(const G4double ee)
{
  G4double meanWave = 128.0*nm; //nm
  G4double meanE = LambdaE/meanWave;
  G4double sigmaE =  (0.56/2.355)*eV;  // sigma(ruler) from FWHM and convert ruler to electron volts
  G4double x = (ee-meanE)/sigmaE;
  G4double emit = exp(-0.5*x*x);
  return emit;
}

