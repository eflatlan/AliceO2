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

/// \file   MIDWorkflow/RawDumpSpec.h
/// \brief  Device to dump decoded raw data
/// \author Diego Stocco <Diego.Stocco at cern.ch>
/// \date   17 February 2022

#ifndef O2_MID_RAWDUMPSPEC_H
#define O2_MID_RAWDUMPSPEC_H

#include "Framework/DataProcessorSpec.h"

namespace o2
{
namespace mid
{
framework::DataProcessorSpec getRawDumpSpec();
} // namespace mid
} // namespace o2

#endif // O2_MID_RAWDUMPSPEC_H
