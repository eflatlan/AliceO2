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

#include "O2Version.h"

namespace o2
{

std::string fullVersion()
{
  return "1.2.0";
}

std::string gitRevision()
{
  return "ae9ad5d22";
}

/// get information about build environment (for example OS and alibuild/alidist release when used)
std::string getBuildInfo()
{
  return "Built by --unavailable-- on OS:Linux-5.15.0-92-generic";
}

} // namespace o2
