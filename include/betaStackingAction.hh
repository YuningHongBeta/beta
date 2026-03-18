#ifndef betaStackingAction_hh
#define betaStackingAction_hh

#include <unordered_map>
#include <unordered_set>
#include <string>

#include "G4UserStackingAction.hh"
#include "G4VUserTrackInformation.hh"
#include "G4ClassificationOfNewTrack.hh"

class G4Track;

// クローン判定用の最小 TrackInfo
// クローン判定用の最小 TrackInfo
class betaTrackInfo : public G4VUserTrackInformation {
public:
  enum OriginType {
    kUnknown = 0,
    kEvaporation = 1,
    kFission = 2
  };

  betaTrackInfo(bool fromPiMinus = false, bool isClone = false, OriginType origin = kUnknown)
  : fromPiMinus_(fromPiMinus), isClone_(isClone), origin_(origin) {}

  bool FromPiMinus() const { return fromPiMinus_; }
  bool IsClone()     const { return isClone_; }
  OriginType GetOrigin() const { return origin_; }

private:
  bool fromPiMinus_;
  bool isClone_;
  OriginType origin_;
};

class betaStackingAction : public G4UserStackingAction {
public:
  explicit betaStackingAction(int nClones = 1)
  : nClones_(nClones) {}

  ~betaStackingAction() override;

  void SetNClones(int n) { nClones_ = n; }

  // イベント開始時に呼ばれるので、ここで per-event 情報をクリア
  void PrepareNewEvent() override;

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override;

private:
  bool IsPiMinusSeedProcess(const G4Track* t) const;

  int nClones_ = 5;
  double neutronCloneEkinMin_ = 20.0;

  // trackID -> “π-起源フラグ”
  std::unordered_map<int, bool> fromPiMinus_;

  // 同一trackを複数回スプリットしないため（保険）
  std::unordered_set<int> alreadySplit_;
};

#endif