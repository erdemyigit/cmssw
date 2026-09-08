# Branch `econ-v4-production` — notes

## Purpose

Configuration and concentrator fixes needed to run the **v4-trained** ECON autoencoder models
(AE / CAE / FiLM, one graph per eLink class 2/3/4/5) in the HGCAL L1 trigger-primitive
simulation. The branch is based on `CMSSW_12_5_2_patch1` at
`b46a31701452e7ded5a11db754cca48f735e2f76`, which is a verbatim snapshot of the LPC working
area that produced the July-2026 ntuples. **That branch is the historical record and is not
modified.** All changes live on this branch only.

Evidence lives in the private `ECON_CAE` repository; only paths are cited here, never content.

## Changes

### 1. `L1Trigger/L1THGCal/src/concentrator/HGCalConcentratorAutoEncoderImpl.cc` — conditioning inputs

Replaced with the May-2026 corrected file (`cmssw/corrected/HGCalConcentratorAutoEncoderImpl_fixed.cc`).
Verified before applying: the fork's file was byte-identical (md5 `f41e58c7…`) to the production
copy `cmssw/production/HGCalConcentratorAutoEncoderImpl.cc`, and `diff production corrected`
contains exactly three hunks:

| hunk | what | why |
|---|---|---|
| `originalCALQsum` loop | skip `(u,v) == (0,0)` | training (`process_data.py`) sums CALQ over indices 1..63; production included CALQ0 |
| decoder condition `[16]` (eta) | module-centre eta via `triggerTools_.getTriggerGeometry()->getModulePosition(getModuleFromTriggerCell(detId))` | training uses `wafer.eta` (module centre); production used the first trigger cell's eta |
| comment near the `log(sumCALQ)` step | text only | — |

The header (`interface/concentrator/HGCalConcentratorAutoEncoderImpl.h`) already declares
`HGCalTriggerTools triggerTools_` and `int verbose_`; `AEinputUtil.h` is untouched. No
interface change was required.

### 2. Same file — debug prints gated on `verbose_`

Five groups of `printf` ran unconditionally (`nInputs_` in the constructor; the
`originalCALQsum/ADCsum/INPUTsum/modSum` block; `"Conditional Info"`; the decoder-input
dump + `"END OF Conditionals"`; `"Past Decoder"`). They are now inside `if (verbose_) { … }`,
using the flag the class already reads from `conf.getParameter<int>("verbose")`. Nothing is
removed. Note that `verbose` is supplied only by `L1THGCalUtilities/python/concentrator.py::CreateAutoencoder(verbose=…)`;
the bare `autoEncoder_conc_proc` PSet in the cfi has no `verbose` key (pre-existing, unchanged).
The production test config sets `verbose=True`, so its stdout is unchanged by this hunk.

### 3. `L1Trigger/L1THGCalUtilities/test/test_CAE.py` — v4 model configuration

The July-2026 production block is retained as a comment; a new v4 block follows it.

* **`linkToGraphMap` override removed.** Upstream default (fork cfi, unchanged):
  `linkToGraphMapping = [0,0,0,1,2,3,3,3,3,3,3,3,3,3,3]`, i.e. nLinks 2→graph 0, 3→1, 4→2,
  ≥5→3. The v4 preprocessing groups wafers with CMSSW's own link table and the same
  `class_from_nlinks` map (`reports/PREPROC_V4.md`), so training grouping == router only under
  this default. The production override `[0,0,0,1,1,2,2,3,…]` sent nLinks 4 to the eLink-3 graph
  and nLinks 5,6 to the eLink-4 graph (`reports/ALLOC_COND_AUDIT.md` §4.3; `docs/DEPLOYMENT.md` §3).
* **`bitsPerLink` set explicitly:** `[0, 1, 3, 5, 7, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]` (indexed by
  nLinks). Since `bitsPerOutput = bitsPerLink.at(nLinks)` and `graphIndex = linkToGraphMap.at(nLinks)`
  are looked up independently, this gives graph 0/1/2/3 exactly 3/5/7/9 bits at every nLinks
  routed to it under the upstream map. It equals the upstream default — the production defect
  was purely the map override desynchronising the two vectors (32.7 % of wafers,
  `docs/DEPLOYMENT.md` §4).
* **nLinks == 1 — DECISION REQUIRED, not taken here.** 2.51 % of wafers (485,731; planes 7-47)
  have nLinks = 1 under the CMSSW table (`docs/DATA.md` §10; `reports/PREPROC_V4.md`). They are
  routed to the eLink-2 graph (3-bit trained) but truncated at `bitsPerLink[1] = 1`. The config
  keeps the upstream value 1 and carries a `# DECISION REQUIRED` block with the options
  (leave 1 / set 3 / train a 1-bit graph).
* **`decoderShape`:** `[1,16]` for the AE arm, `[1,24]` for CAE and FiLM (16 latent ++ 8
  conditions; the FiLM graph is exported as a single-input graph that slices internally, so it
  uses the CAE code path). `encoderShape` stays `[1,8,8,1]`.
* **Model paths are PLACEHOLDERS** (`ECON_V4_MODEL_ROOT = '/PATH/TO/econ_v4_models'`). The files
  are produced by `training/final_run/export_cmssw.py` (ECON_CAE repo) as
  `<ARM>/eLink<N>/{encoder,decoder}_<ARM>_model.pb`; copy them under
  `L1Trigger/L1THGCal/data/models/<tag>/` and set the root. `cms.FileInPath` fails at
  configuration time until this is done — intentionally loud.
* `preserveModuleSum=True` (explicit; also the upstream default). A `FILM` concentrator is
  registered but `standard_concentrators = ['CAE','AE','Threshold0']` is unchanged from
  production; append `'FILM'` to run it. Everything else in the file is as production.

## Known, unfixed: B20 (V16 detector + V11 link map)

The production `Phase2C17I13M9` config pairs a V16 detector (47 planes) with the V11 trigger
link map `hgcal_trigger_link_mapping_120links_v1.json` and a V11-shaped `DisconnectedLayers`
list (`docs/BUG_HISTORY.md` B20). This branch does **not** touch geometry/link-map files —
out of scope. Note the v4 training grouping was derived from this same V11 table on purpose so
that training and router agree; fixing B20 later means regrouping the training data too.

## Not verified

* **Not compiled.** No CMSSW release is available on the machine that made these edits. The
  C++ hunks are small and use only members already declared in the header, but `scram b` has
  not been run.
* The Python config was checked with `python3 -m py_compile` / `ast.parse` only — not loaded
  under `cmsRun` (it cannot load until the model-path placeholders are real files).
* No event was processed; no comparison of ntuple output vs the July production exists yet.

## How to test on LPC

```bash
cmsrel CMSSW_12_5_2_patch1 && cd CMSSW_12_5_2_patch1/src && cmsenv
git cms-init   # or: git init; then
git remote add erdemyigit https://github.com/erdemyigit/cmssw && git fetch erdemyigit econ-v4-production
git checkout erdemyigit/econ-v4-production -- L1Trigger/L1THGCal L1Trigger/L1THGCalUtilities
# drop the exported v4 .pb files under L1Trigger/L1THGCal/data/models/<tag>/{AE,CAE,FILM}/eLink{2,3,4,5}/
# edit ECON_V4_MODEL_ROOT in L1Trigger/L1THGCalUtilities/test/test_CAE.py accordingly
scram b -j 8
cd L1Trigger/L1THGCalUtilities/test
# point process.source.fileNames at the staged 160-event DoubleElectron PU200 file, then
cmsRun test_CAE.py > run.log 2>&1
```

Expected in `run.log` (with `verbose=True`): one `nInputs 64` line per concentrator instance
at construction; per encoded module the `----` separator, `tc (u, v) has ADC …` lines,
`originalCALQsum` / `originalADCsum` / `originalINPUTsum` / `modSum` values, `INPUT` and
`CALQ INPUT` 8×8 dumps, `bitsPerOutput` = 3/5/7/9 according to the module's nLinks (1 for
nLinks = 1 until the decision above is taken), for CAE/FiLM a `Conditional Info … END OF
Conditionals` block with 24 decoder inputs (index 16 = module-centre eta / 3.1), then
`Past Decoder` and the `OUTPUTS` 8×8 dump. No `BadInitialization` exception. The ntuple
`ntuple.root` should contain the `CAE`, `AE` and `Threshold0` trigger-cell collections.
