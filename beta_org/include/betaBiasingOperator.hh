#ifndef betaBiasingOperator_hh
#define betaBiasingOperator_hh 1

#include "G4VBiasingOperator.hh"
#include "G4BOptnChangeCrossSection.hh"

class betaBiasingOperator : public G4VBiasingOperator {
public:
    betaBiasingOperator(G4String name, G4double factor);
    virtual ~betaBiasingOperator();

    // 必須の実装メソッド
    virtual G4VBiasingOperation* ProposeOccurenceBiasingOperation(
        const G4Track* track, const G4BiasingProcessInterface* callingProcess) override;

    // 今回は使わないが実装必須なメソッド
    virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation(
        const G4Track*, const G4BiasingProcessInterface*) override { return 0; }
    virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation(
        const G4Track*, const G4BiasingProcessInterface*) override { return 0; }

private:
    G4BOptnChangeCrossSection* fChangeXS; // 断面積を変更する実働部隊
    G4double fXSFactor;                   // 倍率
    G4int fProcessID;                     // プロセス識別のためのID
};
#endif
