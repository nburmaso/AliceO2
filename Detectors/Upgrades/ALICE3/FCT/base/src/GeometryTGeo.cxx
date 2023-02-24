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

/// \file GeometryTGeo.cxx
/// \brief Implementation of the GeometryTGeo class
/// \author cvetan.cheshkov@cern.ch - 15/02/2007
/// \author ruben.shahoyan@cern.ch - adapted to ITSupg 18/07/2012
/// \author rafael.pezzi@cern.ch - adapted to ALICE 3 EndCaps 14/02/2021

// ATTENTION: In opposite to old AliITSgeomTGeo, all indices start from 0, not from 1!!!

#include "FCTBase/GeometryTGeo.h"
#include "DetectorsBase/GeometryManager.h"
#include "MathUtils/Cartesian.h"

#include <fairlogger/Logger.h> // for LOG

#include <TGeoBBox.h> // for TGeoBBox
#include <TGeoCompositeShape.h>
#include <TGeoManager.h>      // for gGeoManager, TGeoManager
#include <TGeoPhysicalNode.h> // for TGeoPNEntry, TGeoPhysicalNode
#include <TGeoTube.h>
#include <TGeoShape.h>  // for TGeoShape
#include <TMath.h>      // for Nint, ATan2, RadToDeg
#include <TString.h>    // for TString, Form
#include "TClass.h"     // for TClass
#include "TGeoMatrix.h" // for TGeoHMatrix
#include "TGeoNode.h"   // for TGeoNode, TGeoNodeMatrix
#include "TGeoVolume.h" // for TGeoVolume
#include "TMathBase.h"  // for Max

#include "TObjArray.h" // for TObjArray
#include "TObject.h"   // for TObject
#include <cctype>      // for isdigit
#include <cstdio>      // for snprintf, NULL, printf
#include <cstring>     // for strstr, strlen

using namespace TMath;
using namespace o2::fct;
using namespace o2::detectors;

ClassImp(o2::fct::GeometryTGeo);

std::unique_ptr<o2::fct::GeometryTGeo> GeometryTGeo::sInstance;

std::string GeometryTGeo::sVolumeName = "FCTV";          ///< Mother volume name
std::string GeometryTGeo::sInnerVolumeName = "FCTInner"; ///< Mother inner volume name
std::string GeometryTGeo::sLayerName = "FCTLayer";       ///< Layer name
std::string GeometryTGeo::sChipName = "FCTChip";         ///< Sensor name
std::string GeometryTGeo::sSensorName = "FCTSensor";     ///< Sensor name

//__________________________________________________________________________
GeometryTGeo::GeometryTGeo(bool build, int loadTrans) : o2::itsmft::GeometryTGeo(DetID::FCT)
{
  // default c-tor, if build is true, the structures will be filled and the transform matrices
  // will be cached
  if (sInstance) {
    LOG(fatal) << "Invalid use of public constructor: o2::fct::GeometryTGeo instance exists";
    // throw std::runtime_error("Invalid use of public constructor: o2::fct::GeometryTGeo instance exists");
  }

  if (build) {
    Build(loadTrans);
  }
}

//__________________________________________________________________________
void GeometryTGeo::Build(int loadTrans)
{
  if (isBuilt()) {
    LOG(warning) << "Already built";
    return; // already initialized
  }

  if (!gGeoManager) {
    // RSTODO: in future there will be a method to load matrices from the CDB
    LOG(fatal) << "Geometry is not loaded";
  }

  for (auto vol : *gGeoManager->GetListOfVolumes()) {
    LOGP(info, "vol name: {}", vol->GetName());
  }

  mNumberOfLayers = extractNumberOfLayers();

  LOGP(info, "mNumberOfLayers={}", mNumberOfLayers);

  mRowsX.resize(mNumberOfLayers);
  mRowsY.resize(mNumberOfLayers);
  mBotLeft.resize(mNumberOfLayers);

  mSize = 0;

  // fixme: assuming every layer to be sensitive
  for (int32_t i = 0; i < mNumberOfLayers; ++i) {
    std::string layerName = GeometryTGeo::getFCTLayerPattern() + std::string("_") + std::to_string(i);
    TGeoVolume* volLay = gGeoManager->GetVolume(layerName.c_str());
    if (!volLay) {
      LOG(fatal) << "can't find " << layerName << " volume";
      return;
    }
    std::string path = std::string("/cave_1/barrel_1/FCTV_2/") + layerName + std::string("_1");
    gGeoManager->cd(path.c_str());
    auto* trans = gGeoManager->GetCurrentMatrix()->GetTranslation();
    float origX = trans[0];
    float origY = trans[1];
    float origZ = trans[2];
    float blX, blY, blZ;
    int32_t rowsX, rowsY;
    auto* shape = volLay->GetShape();
    if (shape->IsA() == TGeoTube::Class()) {
      auto* tube = (TGeoTube*)shape;
      float maxR = tube->GetDX();
      rowsX = static_cast<int32_t>(std::ceil(2.f * maxR / mPadSizeX));
      rowsY = static_cast<int32_t>(std::ceil(2.f * maxR / mPadSizeY));
      blX = origX - maxR;
      blY = origY - maxR;
      blZ = origZ;
    }
    if (shape->IsA() == TGeoCompositeShape::Class()) {
      auto* square = (TGeoCompositeShape*)shape;
      float dx = square->GetDX();
      float dy = square->GetDY();
      rowsX = static_cast<int32_t>(std::ceil(2.f * dx / mPadSizeX));
      rowsY = static_cast<int32_t>(std::ceil(2.f * dy / mPadSizeY));
      float origX = square->GetOrigin()[0];
      float origY = square->GetOrigin()[1];
      float origZ = square->GetOrigin()[2];
      blX = origX - dx;
      blY = origY - dy;
      blZ = origZ;
    }
    mRowsX[i] = rowsX;
    mRowsY[i] = rowsY;
    mSize += mRowsX[i] * mRowsY[i];
    mBotLeft[i] = math_utils::Vector3D<float>(blX, blY, blZ);
    LOGP(info, "bottom left={}, {}, {}", blX, blY, blZ);
    LOGP(info, "mRowsX[{}]={}, mRowsY[{}]={}, mSize={}", i, mRowsX[i], i, mRowsY[i], mSize);
    gGeoManager->cd();
  }

  fillMatrixCache(loadTrans);
}

//__________________________________________________________________________
int GeometryTGeo::extractVolumeCopy(const char* name, const char* prefix) const
{
  TString nms = name;
  if (!nms.BeginsWith(prefix)) {
    return -1;
  }
  nms.Remove(0, strlen(prefix) + 1);
  if (!isdigit(nms.Data()[0])) {
    return -1;
  }

  return nms.Atoi();
}

//__________________________________________________________________________
Int_t GeometryTGeo::extractNumberOfLayers()
{
  Int_t numberOfLayers = 0;

  TGeoVolume* volFCT = gGeoManager->GetVolume(getFCTVolPattern());
  if (!volFCT) {
    LOG(fatal) << "FCT volume " << getFCTVolPattern() << " is not in the geometry";
  }

  TObjArray* nodes = volFCT->GetNodes();
  int nNodes = nodes->GetEntriesFast();

  for (int j = 0; j < nNodes; j++) {
    Int_t layID = -1;
    auto* nd = (TGeoNode*)nodes->At(j);
    const Char_t* name = nd->GetName();

    if (strstr(name, getFCTLayerPattern())) {
      numberOfLayers++;
      if ((layID = extractVolumeCopy(name, getFCTLayerPattern())) < 0) {
        LOG(fatal) << "Failed to extract layer ID from the " << name;
        exit(1);
      }
    }
  }

  return numberOfLayers;
}

//__________________________________________________________________________
const char* GeometryTGeo::composeSymNameLayer(Int_t d, Int_t lr)
{
  return Form("%s/%s%d", composeSymNameFCT(d), getFCTLayerPattern(), lr);
}

//__________________________________________________________________________
const char* GeometryTGeo::composeSymNameChip(Int_t d, Int_t lr)
{
  return Form("%s/%s%d", composeSymNameLayer(d, lr), getFCTChipPattern(), lr);
}

//__________________________________________________________________________
const char* GeometryTGeo::composeSymNameSensor(Int_t d, Int_t lr)
{
  return Form("%s/%s%d", composeSymNameChip(d, lr), getFCTSensorPattern(), lr);
}

//__________________________________________________________________________
void GeometryTGeo::fillMatrixCache(int mask)
{
  // populate matrix cache for requested transformations
  //
}
