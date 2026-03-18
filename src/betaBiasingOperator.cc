#include "betaBiasingOperator.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4PionMinus.hh"

betaBiasingOperator::betaBiasingOperator(G4String name, G4double factor)
    : G4VBiasingOperator(name), fXSFactor(factor)
{
    // "pi-Inelastic" という名前の変更操作オブジェクトを作っておく
    fChangeXS = new G4BOptnChangeCrossSection("ChangeXS");
}

betaBiasingOperator::~betaBiasingOperator()
{
    delete fChangeXS;
}

G4VBiasingOperation *betaBiasingOperator::ProposeOccurenceBiasingOperation(
    const G4Track *track, const G4BiasingProcessInterface *callingProcess)
{
    // 1. 対象が pi- かどうかチェック
    if (track->GetDefinition() != G4PionMinus::PionMinus())
        return 0;

    // 2. 呼んできたプロセスが「非弾性散乱」かチェック
    //    (Geant4ではバイアス用にラップされたプロセス名を見る必要があります)
    const G4VProcess *physicsProcess = callingProcess->GetWrappedProcess();
    if (!physicsProcess)
        return 0;
    if (physicsProcess->GetProcessName().find("Inelastic") == std::string::npos)
        return 0;

    // --- ここで断面積を変更する ---

    // 本来の断面積を取得
    G4double interactionLength = callingProcess->GetCurrentInteractionLength();

    // 断面積がある（反応する可能性がある）場合のみ処理
    if (interactionLength > 0.0)
    {
        // 断面積を N倍 に設定
        fChangeXS->SetBiasedCrossSection((1.0 / interactionLength) * fXSFactor);
        fChangeXS->Sample();

        // ★重要: Geant4は自動で「確率を上げた分、重み(Weight)を下げる」処理をしますが、
        // もしあなたが「物理的に断面積がでかい世界」として扱いたいなら、
        // 重みは解析時に無視する（あるいは1.0とみなす）必要があります。

        return fChangeXS;
    }

    return 0;
}
