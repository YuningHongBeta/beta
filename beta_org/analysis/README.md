# BGO + TH/TLC classifier v2

## BGOegg photon-counter baseline and exact-frustum design scan

The focused ROOT macro below reports the historical four-neighbour
`nCluster==1` fraction and the additional scalar photon-counter veto. It reads
the segmentation from `runmeta` and uses defaults of 1 MeV per hit, 1 MeV per
cluster, and 0.5 MeV summed photon-counter plastic energy:

```bash
root -l -b -q 'analysis/summarize_bgo_pc.C("INPUT.root")'
```

This is only a reproducible baseline. The BGOegg feasibility program should
replace the rectangular-neighbour cluster definition with the exact crystal
adjacency and test hit-pattern/dE classifiers.

The exact-frustum Photon Counter study uses runmanager schema
`beta-bgo-th-tlc-pc-design-v5`. It fixes BGOegg to 31x60, z=-10 cm and the BGO
per-segment threshold to 3 MeV. The endcap is a conical annulus or disk with
configurable inner/outer polar angles held constant through every Pb/plastic layer. The
extractor joins `evt` by `eventID` and adds `pcSumE_MeV`, `pcDownE_MeV`, and
`pcUpE_MeV`. New ROOT files also provide unique gamma-entry count and energy
branches for the first Pb layer, allowing geometry and intrinsic response to be
separated. BGO-only features and the frozen score definition are unchanged.
The PC analysis threshold is scanned separately and is not the BGO threshold.

Photon Counter (PC) は BGO score と結合した分類器にはせず、独立した
hard-hit veto として評価する。`bgoegg_pc_detector_evaluate.py` は修論と同じ
固定 cut で BGO、TH、TLC、delta-z、PC を順に評価する。
`bgoegg_pc_aperture_evaluate.py` は PC について、gamma が Pb 第1層へ入る
幾何学的入射率、入射事象に対する plastic hit 率、両者を含む実効 PC veto
率を分離する。全生成、BGO `Ncluster=1` 後、全 pre-PC cut 後の分母を混ぜない。
BGO threshold は常に 3 MeV で、PC plastic threshold は独立に走査する。
独立seedの200万π0確認は`bgoegg_pc_confirm_evaluate.py`でgeometryとPC cutを
固定したまま集計する。

## 目的と統計

`physicsFlag=4`のelectron/π− sampleを各100,000 event使い、BGO threshold
3, 4, 5 MeVとsegmentationを比較する。固定detector operating pointは100k全event、
学習scoreだけをeven eventID=train / odd eventID=testで評価する。

## 3 MeV固定のBGOegg比較

`bgo_threshold3_analysis_v1.py`はper-segment thresholdを3 MeVに固定し、
current ball、BGOegg 22x60、端部ring追加案を共通の背景除去条件で比較する。
3 MeV未満のdepositはextractorでzero化され、classifier入力には使わない。

追加したhurdle likelihoodは`nHit=0,1,2,3,4,5+`を離散成分として扱い、
各bin内でtotal energyとshower topologyの対角Gaussian likelihoodを作る。
空または低統計binはspecies全体のmomentへpseudo-count 50でshrinkする。
feature、method、cutは`eventID mod 4 == 0/2`だけで決め、odd eventは
candidate evaluationである。複数候補の比較にもodd eventを使っているため、
study全体のblind holdoutではない。

```bash
python3 analysis/bgo_threshold3_analysis_v1.py

# 新seed ROOTを3 MeVで抽出後、fit/cutを変えず独立確認する。
python3 analysis/bgo_threshold3_confirm_v1.py \
  --confirm-cache tmp/bgo_threshold3_confirm_s7302026
```

比較の最低条件は1 MeV current-ball基準のπ−除去92.212%、
π0除去98.894%である。`guarded` operating pointは25k validation背景sampleの
一側95%程度の統計ゆらぎを見込み、validation目標だけを92.50%、99.01%へ
厳しくしたものを表す。

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
g++ -O2 -std=c++17 analysis/check_bgoegg_geometry.cc \
  -o analysis/check_bgoegg_geometry
analysis/check_bgoegg_geometry

# legacy geometry scan
analysis/bgo_extract_features_v2 INPUT_E.root GEOM_e_T3.bgo2 3 \
  --allow-legacy-hodo
analysis/bgo_extract_features_v2 INPUT_PIM.root GEOM_pim_T3.bgo2 3 \
  --allow-legacy-hodo
python3 analysis/bgo_features_v2.py GEOM_e_T3.bgo2 --expect-rows 100000
```

strict timingでは`--allow-legacy-hodo`を付けない。T4/T5も同様に抽出する。
現在のschemaは130 float features。equal-solidではcos(theta) bin edgeの中点を
cell方向に用いる。`bgoegg_published`ではAnnual Reportの非一様ring tableと
各frustumの8頂点平均を使い、`bgoZOffset_cm`を加えたworld座標からtarget方向を
計算する。uniform/equal-solid geometryでもz offsetをcell方向へ反映する。

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
