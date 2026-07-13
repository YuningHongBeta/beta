# BGO + TH/TLC classifier v2

## 目的と統計

`physicsFlag=4`のelectron/π− sampleを各100,000 event使い、BGO threshold
3, 4, 5 MeVとsegmentationを比較する。固定detector operating pointは100k全event、
学習scoreだけをeven eventID=train / odd eventID=testで評価する。

## 入力制約とreadout mode

- 必須tree: `runmeta`, `calarr`, `th`, `tlc`
- 常時使用: `calarr.{eventID,dE_MeV}`, `th.{eventID,dE_MeV,time_ns}`,
  `tlc.{eventID,cherenkovExpectedPhotons,cherenkovTime_ns}`
- `strict-timing`だけ: `th.zReco_mm`
- 禁止: `target`, `pdg`, `pid`, `th.chargedPath_truth_mm`,
  `tlc.dE_truth_MeV`, `tlc.chargedPath_truth_mm`, `tlc.cherenkovPath_mm`
- treeはentry番号でなく`eventID`で結合する。duplicate、missing、extraはdefaultで
  errorにする。

modeはsidecarの`readoutMode`に保存する。

- `legacy-readout`: timing追加前のclean geometry scan用。明示的に
  `--allow-legacy-hodo`を付ける。`th.zReco_mm`は読まず、`timingResponse=null`、
  Δz/dE/dx timing結果はoptimizerで`null`になる。
- `strict-timing`: default。TH左右time/difference/zRecoとrunmetaのTH半径・z範囲・
  有効光速・timing modelを全て要求する。fallbackは禁止する。

## Buildと抽出

```bash
g++ -O3 -std=c++17 analysis/bgo_extract_features_v2.cc \
  $(root-config --cflags --libs) -o analysis/bgo_extract_features_v2
analysis/bgo_extract_features_v2 --self-test-geometry

# legacy geometry scan
analysis/bgo_extract_features_v2 INPUT_E.root GEOM_e_T3.bgo2 3 \
  --allow-legacy-hodo
analysis/bgo_extract_features_v2 INPUT_PIM.root GEOM_pim_T3.bgo2 3 \
  --allow-legacy-hodo
python3 analysis/bgo_features_v2.py GEOM_e_T3.bgo2 --expect-rows 100000
```

strict timingでは`--allow-legacy-hodo`を付けない。T4/T5も同様に抽出する。
現在のschemaは130 float features。equal-solidではcos(theta) bin edgeの中点を
cell方向に用いる。

結果JSONの`source_root_files`はbeta repository root相対の再生成元で、こちらを
authoritative provenanceとする。`electron`/`pion`/`pi0`および
`feature_files_ephemeral`の`tmp/*.bgo2`は一時抽出物で、repositoryには保存しない。

QA-only truth検証はclassifierと分離する。

```bash
g++ -O2 -std=c++17 analysis/validate_hodo_root_schema.cc \
  $(root-config --cflags --libs) -o analysis/validate_hodo_root_schema
analysis/validate_hodo_root_schema STRICT.root
```

このvalidatorだけはtruth pathを読み、vector長、path関係、timing sentinel、
`zReco=0.5*v_eff*(tLeft-tRight)`を確認する。classifierには渡さない。

## BGOとTLC

BGOは4/8-neighbor cluster、seed=1.5×threshold、cluster sum/size/RMS、isolated
fraction、local 3×3、multiplicity、theta/phi spanを出す。any-hit vetoと
`nCluster==1`もbaselineとして保存する。

TLCはanalytic expected photonsにcollection×PDE=1,2,5,10%を掛け、eventIDと
efficiencyだけをseedにPoisson Npeを生成する。主operating pointは

`hit = any of 30 TLC segments has Npe >= fixed threshold (1,2,3)`

で、optimizerにNpe thresholdを再選択させない。100k全eventのdirect効率を
`fixed_tlc_binary_response`に保存する。

- `matched_pm1_ablation`: BGO方向に対応するTH segment±1内のTLCだけ
- `event_total_npe_comparison`: 全segment Npe和の比較だけ
- continuous expected photons、total/max Npe、自由cutは
  `diagnostic_response_scan`だけで最終推奨には使わない。

## strict TH dE/dxとΔz

BGO leading cell中心をtarget中心からTHへ外挿し、predicted phi segment±1のdE最大barを
matchする。truth pathは使わない。

- path: `(thRMax-thRMin)/sin(theta_BGO)`
- `zPred=0.5*(thRMin+thRMax)*cot(theta_BGO)`
- Δzの主比較はvalid event条件付きで`abs(zReco-zPred)<90 mm`のpass/rejectを出す。
  `valid && absDeltaZ<90 mm`を全eventへ要求したfull selectionは別欄にし、
  missing readoutによるrejectとΔz shapeによるrejectを混ぜない。

source timingはno-smearingを要求する。per-end独立Gaussian
σt=0,0.2,0.5,1.0 nsをevent/segment hashで加え、固定90 mm cutをscanする。
`v_eff`はrunmeta値（現在150 mm/ns）。実測calibration未入力なのでσtは最適化しない。
TH dE/dxも0.2,0.3,0.4,0.5,0.75,1.0 MeV/mmの事前固定gridを100k全eventで出し、
best thresholdは自動選択しない。

`fixed_combined_response`は主TLC point（5%, Npe≥1 any-segment）単独、
TH dE/dx<0.3/0.4とのAND、full Δz selectionとのANDを100k全eventで比較する。
主operating pointはTLC単独で、追加cutはelectron保持とのtrade-offとして扱う。

## Optimizer

```bash
python3 analysis/bgo_optimize_v2.py \
  --electron GEOM_e_T3.bgo2 --pion GEOM_pim_T3.bgo2 \
  --output GEOM_T3.json --groups bgo --methods fisher
```

- `legacy_response_scan`: BGO + TH dE + 固定binary TLC
- `final_response_scan`: strictだけ。BGO + geometry TH dE/dx + Δz + 固定binary TLC
- `diagnostic_response_scan`: continuous/free Npe systematic診断
- `--pi0 PI0.bgo2`: strict固定Δz/dE/dx表へπ0効率を追加
- learned scoreだけπ− reject 93/94/95%をtrainで合わせ、testにWilson 95% CIを付ける。

geometry、segmentation、physicsFlag、threshold、readout mode不一致はerrorにする。

## 検証済み事項

- equal-solid 10×20とuniform-theta 15×15の方向self-test
- TH/TLC entry順を逆転してもBGO2 byte一致
- clean exit0 a10x20 e/π−各100k、missing/extra=0、T3/T4/T5
- legacy sidecarは`th.zReco_mm`を記載せずtiming/Δzを生成しない
- fixed TLC binaryはefficiency×Npe thresholdを迂回せず直接評価する