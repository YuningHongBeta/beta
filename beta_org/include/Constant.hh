// Constant.hh
#ifndef Constants_h
#define Constants_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

constexpr G4double M_proton = 938.2720813; // [MeV/c2]  
constexpr G4double M_neutron = 939.5654133; // [MeV/c2]  
constexpr G4double M_e = 0.51099895000; // [MeV/c2]  
constexpr G4double M_Kaon = 493.677; // [MeV/c2]  
constexpr G4double M_pion0 = 134.9766; // [MeV/c2]  
constexpr G4double M_pion = 139.57061; // [MeV/c2]  
constexpr G4double M_Lambda = 1115.683; // [MeV/c2]  
constexpr G4double M_Omega = 1672.45; // [MeV/c2]  
constexpr G4double M_neue = 0.; // [MeV/c2]  

constexpr G4int shapeFlag = 1; //0 -> Box, 1 -> Shere 
// constexpr G4int nLayer = 20;//30
// constexpr G4int nSector = 8;//60
constexpr G4int nLayer = 15;
constexpr G4int nSector = 15;
//constexpr G4int nb_cryst = 50; //For Box: Comment out for shapeFlag=1
constexpr G4int nb_cryst = nLayer*nSector; //For Sphere: Comment out for shapeFlag=0
constexpr G4int phi = 360./double(nSector);
constexpr G4double absoThickness = 20.*cm;
constexpr G4double ScintiThickness = 1.*cm;
constexpr G4double calorSizeXY   = 40.*cm; //only used in Box shape
// constexpr G4double cellAngle = 160./double(nLayer);

// 目標の穴（半頂角）[deg]
constexpr G4double thetaHoleDown = 9.698;  // +z側（downstream） 89.8 msr
constexpr G4double thetaHoleUp   = 5.666;  // -z側（upstream）   30.7 msr（または 6.45deg）
// constexpr G4double thetaHoleUp   = 9.59;  // -z側（upstream）   87.9 msr for K1.8

// 結晶で覆う θ 範囲 [deg]
constexpr G4double thetaMin = thetaHoleUp;
constexpr G4double thetaMax = 180.0 - thetaHoleDown;

// 各層の幅（今までの cellAngle の定義を置き換え）
constexpr G4double cellAngle = (thetaMax - thetaMin) / double(nLayer);

constexpr G4double Rmin = 30*cm;
constexpr G4double Rmax = absoThickness + Rmin;
constexpr G4double TargetSizeX  = 4.*cm; 
constexpr G4double TargetSizeY  = 4.*cm; 
constexpr G4double TargetThickness  = 8.*cm; 
constexpr G4double TargetRadius = 1.5*cm;
constexpr G4double TargetLength = 30.*cm;

constexpr G4double TLength = 60.*cm;
constexpr G4double ThRmin = 1.5*cm;
constexpr G4double ThRmax = 2.1*cm;
constexpr G4double TlcRmin = 2.2*cm;
constexpr G4double TlcRmax = 2.8*cm;
constexpr G4int nSegTH = 30;
constexpr G4int nSegTLC = 30;
constexpr G4double dPhiTH = 360.*deg/double(nSegTH);
constexpr G4double dPhiTLC = 360.*deg/double(nSegTLC);

// -----------------------------
// physics switch
// 1: FTFP_BERT
// 2: INCLXX & ABLA++ (coupling in RunAction)
// 3: INCLXX & (default de-excitation = G4ExcitationHandler/GEM系)
// 4: (3) + StackingAction(mode=3)
// 5: (4) + BiasingOperator(mode=3)  ※BGO内のみ attach 想定
// -----------------------------
constexpr G4int PhysicsFlag = 1;

// bias settings (必要なら調整)
constexpr G4double NeutronScale  = 2.0; // copy size of enutron tracks
constexpr G4double InelasticBias = 3.0;  // Bias factor for pi- inelastic

#endif
