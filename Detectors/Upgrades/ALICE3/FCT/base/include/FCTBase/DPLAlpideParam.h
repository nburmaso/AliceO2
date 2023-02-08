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

#ifndef ALICEO2_FCTDPLBASEPARAM_H_
#define ALICEO2_FCTDPLBASEPARAM_H_

#include "DetectorsCommonDataFormats/DetID.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "CommonConstants/LHCConstants.h"

#include <string_view>

namespace o2
{
namespace fct
{

constexpr float DEFStrobeDelay = o2::constants::lhc::LHCBunchSpacingNS * 4; // ~100 ns delay

struct DPLAlpideParam : public o2::conf::ConfigurableParamHelper<DPLAlpideParam> {

  static constexpr std::string_view paramName{"FCTDigitizerParam"};
  static constexpr int roFrameLengthInBC{o2::constants::lhc::LHCMaxBunches / 18}; ///< ROF length in BC for continuos mode
  static constexpr float roFrameLengthTrig{6000.};                                ///< length of RO frame in ns for triggered mode
  float strobeDelay{DEFStrobeDelay};                                              ///< strobe start (in ns) wrt ROF start
  float strobeLengthCont{-1.};                                                    ///< if < 0, full ROF length - delay
  float strobeLengthTrig{100.};                                                   ///< length of the strobe in ns (sig. over threshold checked in this window only)
  int roFrameBiasInBC{0};                                                         ///< bias of the start of ROF wrt orbit start: t_irof = (irof*roFrameLengthInBC + roFrameBiasInBC)*BClengthMUS

  static_assert(o2::constants::lhc::LHCMaxBunches % roFrameLengthInBC == 0); // make sure ROF length is divisor of the orbit

  // boilerplate stuff + make principal key
  O2ParamDef(DPLAlpideParam, paramName.data());
};

DPLAlpideParam DPLAlpideParam::sInstance;

} // namespace fct

namespace framework
{
template <typename T>
struct is_messageable;
template <>
struct is_messageable<o2::fct::DPLAlpideParam> : std::true_type {
};

} // namespace framework

} // namespace o2

#endif
