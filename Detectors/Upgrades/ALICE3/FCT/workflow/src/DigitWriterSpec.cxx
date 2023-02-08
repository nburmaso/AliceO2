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

/// @brief  Processor spec for a ROOT file writer for FCT digits

#include <algorithm>
#include <string>
#include <vector>

#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "DataFormatsITSMFT/Digit.h"
#include "DataFormatsITSMFT/GBTCalibData.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "FCTWorkflow/DigitWriterSpec.h"
#include "Headers/DataHeader.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/IOMCTruthContainerView.h"
#include "SimulationDataFormat/MCCompLabel.h"

using namespace o2::framework;
using SubSpecificationType = o2::framework::DataAllocator::SubSpecificationType;

namespace o2::fct
{

template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;
using MCCont = o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>;

/// create the processor spec
/// describing a processor receiving digits for ITS/MFT and writing them to file
DataProcessorSpec getDigitWriterSpec(bool mctruth, bool dec, o2::header::DataOrigin detOrig)
{
  std::string detStr = o2::detectors::DetID::getName(o2::detectors::DetID::FCT);
  std::string detStrL = dec ? "o2_" : ""; // for decoded digits prepend by o2
  detStrL += detStr;
  std::transform(detStrL.begin(), detStrL.end(), detStrL.begin(), ::tolower);
  auto logger = [](std::vector<o2::itsmft::Digit> const& inDigits) {
    LOG(info) << "RECEIVED DIGITS SIZE " << inDigits.size();
  };

  // the callback to be set as hook for custom action when the writer is closed
  auto finishWriting = [](TFile* outputfile, TTree* outputtree) {
    const auto* brArr = outputtree->GetListOfBranches();
    int64_t nent = 0;
    for (const auto* brc : *brArr) {
      int64_t n = ((const TBranch*)brc)->GetEntries();
      if (nent && (nent != n)) {
        LOG(error) << "Branches have different number of entries";
      }
      nent = n;
    }
    outputtree->SetEntries(nent);
    // do not use TTree::Write .. as this writes to default directory (not the associated file)
    // instead of outputtree->Write("", TObject::kOverwrite)
    // --> better use TFile::Write or TFile::WriteObject
    outputfile->Write("", TObject::kOverwrite);
    outputfile->Close();
  };

  // handler for labels
  // This is necessary since we can't store the original label buffer in a ROOT entry -- as is -- if it exceeds a certain size.
  // We therefore convert it to a special split class.
  auto fillLabels = [](TBranch& branch, std::vector<char> const& labelbuffer, DataRef const& /*ref*/) {
    o2::dataformats::ConstMCTruthContainerView<o2::MCCompLabel> labels(labelbuffer);
    LOG(info) << "WRITING " << labels.getNElements() << " LABELS ";

    o2::dataformats::IOMCTruthContainerView outputcontainer;
    auto ptr = &outputcontainer;
    auto br = framework::RootTreeWriter::remapBranch(branch, &ptr);
    outputcontainer.adopt(labelbuffer);
    br->Fill();
    br->ResetAddress();
  };

  return MakeRootTreeWriterSpec((detStr + "DigitWriter" + (dec ? "_dec" : "")).c_str(),
                                (detStrL + "digits.root").c_str(),
                                MakeRootTreeWriterSpec::TreeAttributes{"o2sim", "Digits tree"},
                                MakeRootTreeWriterSpec::CustomClose(finishWriting),
                                // in case of labels we first read them as std::vector<char> and process them correctly in the fillLabels hook
                                BranchDefinition<std::vector<char>>{InputSpec{detStr + "_digitsMCTR", detOrig, "DIGITSMCTR", 0},
                                                                    detStr + "DigitMCTruth",
                                                                    (mctruth ? 1 : 0), fillLabels},
                                BranchDefinition<std::vector<itsmft::MC2ROFRecord>>{InputSpec{detStr + "_digitsMC2ROF", detOrig, "DIGITSMC2ROF", 0},
                                                                                    detStr + "DigitMC2ROF",
                                                                                    (mctruth ? 1 : 0)},
                                BranchDefinition<std::vector<itsmft::Digit>>{InputSpec{detStr + "digits", detOrig, "DIGITS", 0},
                                                                             detStr + "Digit",
                                                                             logger},
                                BranchDefinition<std::vector<itsmft::ROFRecord>>{InputSpec{detStr + "digitsROF", detOrig, "DIGITSROF", 0},
                                                                                 detStr + "DigitROF"})();
}

DataProcessorSpec getFCTDigitWriterSpec(bool mctruth, bool dec, bool calib)
{
  return getDigitWriterSpec(mctruth, dec, o2::header::gDataOriginFCT);
}

} // namespace o2::fct
