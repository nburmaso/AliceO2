# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

add_subdirectory(FemtoDream)

o2_add_library(PWGCFCore
               SOURCES  AnalysisConfigurableCuts.cxx
               PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel)

o2_target_root_dictionary(PWGCFCore
              HEADERS AnalysisConfigurableCuts.h
              LINKDEF PWGCFCoreLinkDef.h)

o2_add_dpl_workflow(dptdptcorrelations
                    SOURCES dptdptcorrelations.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisCore O2::AnalysisDataModel O2::PWGCFCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(correlations
                    SOURCES correlations.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(filter-cf
                    SOURCES filterCF.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(correlations-filtered
                    SOURCES correlationsFiltered.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisCore
                    COMPONENT_NAME Analysis)

