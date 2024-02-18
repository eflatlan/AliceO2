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

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// ef ::added
#pragma link C++ class o2::hmpid::ClusterCandidate + ;
#pragma link C++ class std::vector < o2::hmpid::ClusterCandidate> + ;

// ef ::added

#pragma link C++ function o2::hmpid::Cluster::setClusterTopology;

#pragma link C++ class o2::hmpid::Digit + ;
#pragma link C++ class o2::hmpid::Topology + ;
#pragma link C++ class vector < o2::hmpid::Topology> + ;

#pragma link C++ class vector < int> + ;
#pragma link C++ class vector < o2::hmpid::Digit> + ;
#pragma link C++ class o2::hmpid::HitType + ;
#pragma link C++ class vector < o2::hmpid::HitType> + ;
#pragma link C++ class o2::hmpid::Cluster + ;
#pragma link C++ class vector < o2::hmpid::Cluster> + ;
#pragma link C++ class o2::hmpid::Trigger + ;
#pragma link C++ class vector < o2::hmpid::Trigger> + ;

#pragma link C++ struct o2::hmpid::CTFHeader + ;
#pragma link C++ struct o2::hmpid::CTF + ;
#pragma link C++ class o2::ctf::EncodedBlocks < o2::hmpid::CTFHeader, 8, uint32_t> + ;

#endif
