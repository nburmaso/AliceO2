// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Digitizer.cxx
/// \brief Implementation of the ITS/MFT digitizer

#include "DataFormatsITSMFT/Digit.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "FCTSimulation/DPLDigitizerParam.h"
#include "FCTSimulation/Digitizer.h"
#include "MathUtils/Cartesian.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DetectorsRaw/HBFUtils.h"

#include <TRandom.h>
#include <climits>
#include <vector>
#include <numeric>
#include <fairlogger/Logger.h> // for LOG

using o2::itsmft::Digit;
using o2::itsmft::Hit;
using Segmentation = o2::itsmft::SegmentationAlpide;

using namespace o2::fct;

//_______________________________________________________________________
void Digitizer::init()
{
  mNumberOfChips = mGeometry->getNumberOfChips();
  mChips.resize(mNumberOfChips);
  mChipIDOffsets.resize(mGeometry->mNumberOfLayers); // starting chipID for a layer

  LOGP(info, "mNumberOfChips={}", mNumberOfChips);

  int32_t start = 0;
  for (int32_t i = 0; i < mGeometry->mNumberOfLayers; ++i) {
    mChipIDOffsets[i] = start;
    LOGP(info, "mChipIDOffsets[{}]={}", i, start);
    start += mGeometry->mRowsX[i] * mGeometry->mRowsY[i];
  }
  LOGP(info, "end={}", start);

  // todo: no noise, no dead channels
  // for (int i = mNumberOfChips; i--;) {
  //   mChips[i].setChipIndex(i);
  //   if (mNoiseMap) {
  //     mChips[i].setNoiseMap(mNoiseMap);
  //   }
  //   if (mDeadChanMap) {
  //     mChips[i].disable(mDeadChanMap->isFullChipMasked(i));
  //     mChips[i].setDeadChanMap(mDeadChanMap);
  //   }
  // }

  mParams.print();
  mIRFirstSampledTF = o2::raw::HBFUtils::Instance().getFirstSampledTFIR();
}

auto Digitizer::getChipResponse(int chipID)
{
  return 0; // mAlpSimRespMFT;
}

//_______________________________________________________________________
void Digitizer::process(const std::vector<Hit>* hits, int evID, int srcID)
{
  // digitize single event, the time must have been set beforehand

  //  LOG(info) << "Digitizing " << mGeometry->getName() << " hits of entry " << evID << " from source "
  //            << srcID << " at time " << mEventTime << " ROFrame= " << mNewROFrame << ")"
  //            << " cont.mode: " << isContinuous()
  //            << " Min/Max ROFrames " << mROFrameMin << "/" << mROFrameMax;

  int nHits = hits->size();
  std::vector<int> hitIdx(nHits);
  std::iota(std::begin(hitIdx), std::end(hitIdx), 0);
  // sort hits to improve memory access
  std::sort(hitIdx.begin(), hitIdx.end(),
            [hits](auto lhs, auto rhs) {
              return (*hits)[lhs].GetPosStart().z() < (*hits)[rhs].GetPosStart().z();
            });

  for (int i : hitIdx) {
    processHit((*hits)[i], mROFrameMax, evID, srcID);
  }
}

//_______________________________________________________________________
void Digitizer::setEventTime(const o2::InteractionTimeRecord& irt)
{
  // assign event time in ns
  mEventTime = irt;
  if (!mParams.isContinuous()) {
    mROFrameMin = 0; // in triggered mode reset the frame counters
    mROFrameMax = 0;
  }
  // RO frame corresponding to provided time
  mCollisionTimeWrtROF = mEventTime.timeInBCNS; // in triggered mode the ROF starts at BC (is there a delay?)
  if (mParams.isContinuous()) {
    auto nbc = mEventTime.differenceInBC(mIRFirstSampledTF);
    if (mCollisionTimeWrtROF < 0 && nbc > 0) {
      nbc--;
    }
    mNewROFrame = nbc / mParams.getROFrameLengthInBC();
    // in continuous mode depends on starts of periodic readout frame
    mCollisionTimeWrtROF += (nbc % mParams.getROFrameLengthInBC()) * o2::constants::lhc::LHCBunchSpacingNS;
  } else {
    mNewROFrame = 0;
  }

  if (mNewROFrame < mROFrameMin) {
    LOG(error) << "New ROFrame " << mNewROFrame << " (" << irt << ") precedes currently cashed " << mROFrameMin;
    throw std::runtime_error("deduced ROFrame precedes already processed one");
  }

  if (mParams.isContinuous() && mROFrameMax < mNewROFrame) {
    mROFrameMax = mNewROFrame - 1; // all frames up to this are finished
  }
}

//_______________________________________________________________________
void Digitizer::fillOutputContainer(uint32_t frameLast)
{
  // fill output with digits from min.cached up to requested frame, generating the noise beforehand
  if (frameLast > mROFrameMax) {
    frameLast = mROFrameMax;
  }
  // make sure all buffers for extra digits are created up to the maxFrame
  getExtraDigBuffer(mROFrameMax);

  LOG(info) << "Filling " << mGeometry->getName() << " digits output for RO frames " << mROFrameMin << ":"
            << frameLast;

  o2::itsmft::ROFRecord rcROF;

  // we have to write chips in RO increasing order, therefore have to loop over the frames here
  for (; mROFrameMin <= frameLast; mROFrameMin++) {
    rcROF.setROFrame(mROFrameMin);
    rcROF.setFirstEntry(mDigits->size()); // start of current ROF in digits

    auto& extra = *(mExtraBuff.front().get());
    for (auto& chip : mChips) {
      if (chip.isDisabled()) {
        continue;
      }
      // todo: add noise
      // chip.addNoise(mROFrameMin, mROFrameMin, &mParams);
      auto& buffer = chip.getPreDigits();
      if (buffer.empty()) {
        continue;
      }
      auto itBeg = buffer.begin();
      auto iter = itBeg;
      ULong64_t maxKey = chip.getOrderingKey(mROFrameMin + 1, 0, 0) - 1; // fetch digits with key below that
      for (; iter != buffer.end(); ++iter) {
        if (iter->first > maxKey) {
          break; // is the digit ROFrame from the key > the max requested frame
        }
        auto& preDig = iter->second; // preDigit
        if (preDig.charge >= mParams.getChargeThreshold()) {
          int digID = mDigits->size();
          mDigits->emplace_back(chip.getChipIndex(), preDig.row, preDig.col, preDig.charge);
          mMCLabels->addElement(digID, preDig.labelRef.label);
          auto& nextRef = preDig.labelRef; // extra contributors are in extra array
          while (nextRef.next >= 0) {
            nextRef = extra[nextRef.next];
            mMCLabels->addElement(digID, nextRef.label);
          }
        }
      }
      buffer.erase(itBeg, iter);
    }
    // finalize ROF record
    rcROF.setNEntries(mDigits->size() - rcROF.getFirstEntry()); // number of digits
    if (isContinuous()) {
      rcROF.getBCData().setFromLong(mIRFirstSampledTF.toLong() + mROFrameMin * mParams.getROFrameLengthInBC());
    } else {
      rcROF.getBCData() = mEventTime; // RSTODO do we need to add trigger delay?
    }
    if (mROFRecords) {
      mROFRecords->push_back(rcROF);
    }
    extra.clear(); // clear container for extra digits of the mROFrameMin ROFrame
    // and move it as a new slot in the end
    mExtraBuff.emplace_back(mExtraBuff.front().release());
    mExtraBuff.pop_front();
  }
}

//_______________________________________________________________________
void Digitizer::processHit(const o2::itsmft::Hit& hit, uint32_t& maxFr, int evID, int srcID)
{
  // convert single hit to digits
  // todo: implement digitization based on ALPIDE response simulation (see MFT simulation)

  float timeInROF = hit.GetTime() * sec2ns;
  if (timeInROF > 20e3) {
    const int maxWarn = 10;
    static int warnNo = 0;
    if (warnNo < maxWarn) {
      LOG(warning) << "Ignoring hit with time_in_event = " << timeInROF << " ns"
                   << ((++warnNo < maxWarn) ? "" : " (suppressing further warnings)");
    }
    return;
  }

  if (isContinuous()) {
    timeInROF += mCollisionTimeWrtROF;
  }
  // calculate RO Frame for this hit
  if (timeInROF < 0) {
    timeInROF = 0.;
  }
  float tTot = mParams.getSignalShape().getMaxDuration();
  // frame of the hit signal start wrt event ROFrame
  int roFrameRel = int(timeInROF * mParams.getROFrameLengthInv());
  // frame of the hit signal end  wrt event ROFrame: in the triggered mode we read just 1 frame
  uint32_t roFrameRelMax = mParams.isContinuous() ? (timeInROF + tTot) * mParams.getROFrameLengthInv() : roFrameRel;
  int nFrames = roFrameRelMax + 1 - roFrameRel;
  uint32_t roFrameMax = mNewROFrame + roFrameRelMax;
  if (roFrameMax > maxFr) {
    maxFr = roFrameMax; // if signal extends beyond current maxFrame, increase the latter
  }

  // global position
  math_utils::Vector3D<float> xyzGlobS(hit.GetPosStart()); // start position

  // fixme: hit.GetDetectorID() is incorrect for FCT -> using this workaround to get layer ID
  int16_t layer = -1;
  for (int32_t i = 0; i < mGeometry->mBotLeft.size(); ++i) {
    auto& layerBL = mGeometry->mBotLeft[i];
    if (std::abs(xyzGlobS.z() - layerBL.z()) < 1.) {
      layer = i;
      break;
    }
  }

  // local sensor position
  math_utils::Vector3D<float> xyzBotLeft(mGeometry->mBotLeft[layer]);
  math_utils::Vector3D<float> xyzDetS(xyzGlobS - xyzBotLeft);
  auto rowX = static_cast<int32_t>(std::floor(xyzDetS.x() / mGeometry->mPadSizeX));
  auto rowY = static_cast<int32_t>(std::floor(xyzDetS.y() / mGeometry->mPadSizeY));
  // from bottom left corner (looking from IP)
  int32_t chipID = mChipIDOffsets[layer] + rowX + rowY * mGeometry->mRowsY[layer];
  LOGP(debug, "x={}, y={}, z={}, detZ={}, rowX={}, rowY={}, chipID={}, layer={}",
       xyzDetS.x(), xyzDetS.y(), xyzDetS.z(), xyzBotLeft.z(), rowX, rowY, chipID, layer);
  auto& chip = mChips[chipID];
  chip.setChipIndex(chipID);

  // here we start stepping in the depth of the sensor to generate charge diffusion
  float nStepsInv = mParams.getNSimStepsInv();

  float nElectrons = hit.GetEnergyLoss() * mParams.getEnergyToNElectrons(); // total number of deposited electrons
  nElectrons *= nStepsInv;                                                  // N electrons injected per step

  o2::MCCompLabel lbl(hit.GetTrackID(), evID, srcID, false);
  auto roFrameAbs = mNewROFrame + roFrameRel;

  // 1 digit per sensor
  int32_t nEle = gRandom->Poisson(nElectrons);

  registerDigits(chip, roFrameAbs, timeInROF, nFrames, rowX, rowY, nEle, lbl);
}

//________________________________________________________________________________
void Digitizer::registerDigits(itsmft::ChipDigitsContainer& chip, uint32_t roFrame, float tInROF, int nROF,
                               int32_t row, int32_t col, int nEle, o2::MCCompLabel& lbl)
{
  // Register digits for given pixel, accounting for the possible signal contribution to
  // multiple ROFrame. The signal starts at time tInROF wrt the start of provided roFrame
  // In every ROFrame we check the collected signal during strobe

  float tStrobe = mParams.getStrobeDelay() - tInROF; // strobe start wrt signal start
  for (int i = 0; i < nROF; i++) {
    uint32_t roFr = roFrame + i;
    int nEleROF = mParams.getSignalShape().getCollectedCharge(nEle, tStrobe, tStrobe + mParams.getStrobeLength());
    tStrobe += mParams.getROFrameLength(); // for the next ROF

    // discard too small contributions, they have no chance to produce a digit
    if (nEleROF < mParams.getMinChargeToAccount()) {
      continue;
    }
    if (roFr > mEventROFrameMax) {
      mEventROFrameMax = roFr;
    }
    if (roFr < mEventROFrameMin) {
      mEventROFrameMin = roFr;
    }
    auto key = chip.getOrderingKey(roFr, row, col);
    itsmft::PreDigit* pd = chip.findDigit(key);
    if (!pd) {
      chip.addDigit(key, roFr, row, col, nEleROF, lbl);
    } else { // there is already a digit at this slot, account as PreDigitExtra contribution
      pd->charge += nEleROF;
      if (pd->labelRef.label == lbl) { // don't store the same label twice
        continue;
      }
      ExtraDig* extra = getExtraDigBuffer(roFr);
      int& nxt = pd->labelRef.next;
      bool skip = false;
      while (nxt >= 0) {
        if ((*extra)[nxt].label == lbl) { // don't store the same label twice
          skip = true;
          break;
        }
        nxt = (*extra)[nxt].next;
      }
      if (skip) {
        continue;
      }
      // new predigit will be added in the end of the chain
      nxt = extra->size();
      extra->emplace_back(lbl);
    }
  }
}
