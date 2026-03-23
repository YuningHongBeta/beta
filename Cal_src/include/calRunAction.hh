// calRunAction.hh
#ifndef calRunAction_h
#define calRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class TH1D;
class TTree;
class TFile;

/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit 
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following 
/// physics quantities:
/// - Edep in absorber
/// - Edep in gap
/// - Track length in absorber
/// - Track length in gap
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in calAnalysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed 
/// dispersion is printed.
///

class calRunAction : public G4UserRunAction
{
private:
  TFile* hfile;
  TH1D* Cal_dE;
  TTree* tree;
  int RunSS_Time,RunE_Time;

  public:
    calRunAction();
    virtual ~calRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

