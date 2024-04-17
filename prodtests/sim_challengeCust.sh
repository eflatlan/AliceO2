
#!/bin/bash

# A simple chain of algorithms from MC to reco (and analysis)

# ------------ LOAD UTILITY FUNCTIONS ----------------------------
. ${O2_ROOT}/share/scripts/jobutils.sh
# ----------- START WITH ACTUAL SCRIPT ---------------------------


if [ -z "$SHMSIZE" ]; then export SHMSIZE=10000000000; fi

# default run number
# (for now set to a pilot beam run until we have all CCDB objects for default unanchored MC)
runNumDef=300000

# default time stamp --> will be determined from run number during the sim stage
# startTimeDef=$(($(date +%s%N)/1000000))

# default number of events
nevPP=10
nevPbPb=10

# default interaction rates in kHz
intRatePP=400
intRatePbPb=50

# default collision system
collSyst="pp"

generPP="pythia8pp"
generPbPb="pythia8hi"

# default sim engine
engine="TGeant3"

# options to pass to every workflow
gloOpt=" -b --run --shm-segment-size $SHMSIZE"

# ITS reco options depends on pp or pbpb
ITSRecOpt=""

# option to set the number of sim workers
simWorker=""

# option to set the number of tpc-lanes
tpcLanes=""


detectorsToskip=""

#
Usage()
{
  echo "Usage: ${0##*/} [-s system /pp[Def] or pbpb/] [-r IR(kHz) /Def = $intRatePP(pp)/$intRatePbPb(pbpb)] [-n Number of events /Def = $nevPP(pp) or $nevPbPb(pbpb)/] [-e TGeant3|TGeant4] [-t startTime/Def = $startTimeDef] [-run runNumber/Def = $runNumDef] [-f fromstage sim|digi|reco /Def = sim]"
  exit
}

fromstage="sim"
while [ $# -gt 0 ] ; do
  case $1 in
    -n) nev=$2;  shift 2 ;;
    -s) collSyst=$2; shift 2 ;;
    -r) intRate=$2; shift 2 ;;
    -e) engine=$2; shift 2 ;;
    -f) fromstage=$2; shift 2 ;;
    -j) simWorker="-j $2"; shift 2 ;;
    -l) tpcLanes="--tpc-lanes $2"; shift 2 ;;
    -t) startTime=$2; shift 2 ;;
    -p) pdg=$2; shift 2 ;;
    -run) runNumber=$2; shift 2 ;;
    -h) Usage ;;
    *) echo "Wrong input"; Usage;
  esac
done

# convert to lower case (the bash construct ${collSyst,,} is less portable)
collSyst=`echo "$collSyst" | awk '{print tolower($0)}'`
if [ "$collSyst" == "pp" ]; then
    gener="$generPP"
    ITSRecOpt=" --configKeyValues \"ITSVertexerParam.phiCut=0.5;ITSVertexerParam.clusterContributorsCut=3;ITSVertexerParam.tanLambdaCut=0.2\""
    [[ "nev" -lt "1"  ]] && nev="$nevPP"
    [[ "intRate" -lt "1"  ]] && intRate="$intRatePP"
elif [ "$collSyst" == "pbpb" ]; then
    gener="$generPbPb"
    [[ "nev" -lt "1"  ]] && nev="$nevPbPb"
    [[ "intRate" -lt "1"  ]] && intRate="$intRatePbPb"
else
    echo "Wrong collision system $collSyst provided, should be pp or pbpb"

fi

[[ -z $startTime ]] && startTime=$startTimeDef
[[ -z $runNumber ]] && runNumber=$runNumDef
dosim="0"
dosimp="0"
dosimk="0"
dosimpr="0"
dodigi="0"
dotrdtrap="0"
doreco="0"
# convert to lowercase
fromstage=`echo "$fromstage" | awk '{print tolower($0)}'`
if [ "$fromstage" == "simp" ]; then
  dosimp="1"
  dodigi="1"
  dotrdtrap="1"
  doreco="1"
elif [ "$fromstage" == "simk" ]; then
  dosimk="1"
  dodigi="1"
  dotrdtrap="1"
  doreco="1"
elif [ "$fromstage" == "sim" ]; then
  dosim="1"
  dodigi="1"
  dotrdtrap="1"
  doreco="1"  
elif [ "$fromstage" == "simpr" ]; then
  dosimpr="1"
  dodigi="1"
  dotrdtrap="1"
  doreco="1"
elif [ "$fromstage" == "digi" ]; then
  dodigi="1"
  dotrdtrap="1"
  doreco="1"
elif [ "$fromstage" == "reco" ]; then
  doreco="1"
else
  echo "Wrong stage string $fromstage provided, should be sim or digi or reco"
fi



if [ "$dosim" == "1" ]; then
  #---- GRP creation ------
  echo "Creating GRPs ... and publishing in local CCDB overwrite"
  taskwrapper grp.log o2-grp-simgrp-tool createGRPs --run ${runNumber} --publishto GRP -o mcGRP

  #---------------------------------------------------
  echo "Running simulation for $nev $collSyst events with $gener generator and engine $engine and run number $runNumber"
  taskwrapper sim.log o2-sim -n"$nev"  -m PIPE ITS TPC FT0 HMP TRD TOF CTP --configKeyValues "Diamond.width[2]=6." -g "$gener" -e "$engine" $simWorker --run ${runNumber}

  ##------ extract number of hits
  taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C
fi



if [ "$dosimp" == "1" ]; then
  #---- GRP creation ------
  echo "Creating GRPs ... and publishing in local CCDB overwrite"
  taskwrapper grp.log o2-grp-simgrp-tool createGRPs --run ${runNumber} --publishto GRP -o mcGRP

  #---------------------------------------------------
  echo "Running Pion -211 gun simulation for $nev Events "
  #taskwrapper sim.log o2-sim -n"$nev" --configKeyValues "Diamond.width[2]=6." -g "$gener" -e "$engine" $simWorker --run ${runNumber}

  o2-sim -n"$nev" -e TGeant3 -g boxgen --configKeyValues "BoxGun.pdg=-211; BoxGun.phirange[0]=-5; BoxGun.phirange[1]=60; BoxGun.number=2; BoxGun.eta[0]=-0.5 ; BoxGun.eta[1]=0.5; BoxGun.prange[0]=5.8; BoxGun.prange[1]=5.83;" $simWorker --run ${runNumber}

  ##------ extract number of hits
  taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C
fi

if [ "$dosimk" == "1" ]; then
  #---- GRP creation ------
  echo "Creating GRPs ... and publishing in local CCDB overwrite"
  taskwrapper grp.log o2-grp-simgrp-tool createGRPs --run ${runNumber} --publishto GRP -o mcGRP 

  #---------------------------------------------------
  echo "Running Ka gun simulation for $nev  "
  #taskwrapper sim.log o2-sim -n"$nev" --configKeyValues "Diamond.width[2]=6." -g "$gener" -e "$engine" $simWorker --run ${runNumber}

  o2-sim -n"$nev" -e TGeant3 -g boxgen -m PIPE ITS TPC FT0 HMP TRD TOF CTP --configKeyValues "BoxGun.pdg=321; BoxGun.phirange[0]=-5; BoxGun.phirange[1]=60; BoxGun.number=2; BoxGun.eta[0]=-0.5 ; BoxGun.eta[1]=0.5; BoxGun.prange[0]=2.8; BoxGun.prange[1]=2.83;" $simWorker --run ${runNumber}

  ##------ extract number of hits
  taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C
fi



if [ "$dosimpr" == "1" ]; then
  #---- GRP creation ------
  echo "Creating GRPs ... and publishing in local CCDB overwrite"
  taskwrapper grp.log o2-grp-simgrp-tool createGRPs --run ${runNumber} --publishto GRP -o mcGRP

  #---------------------------------------------------
  echo "Running Proton gun simulation for $nev  "
  #taskwrapper sim.log o2-sim -n"$nev" --configKeyValues "Diamond.width[2]=6." -g "$gener" -e "$engine" $simWorker --run ${runNumber}

  o2-sim -n"$nev" -e TGeant3 -g boxgen --configKeyValues "BoxGun.pdg=2212; BoxGun.phirange[0]=-5; BoxGun.phirange[1]=60; BoxGun.number=2; BoxGun.eta[0]=-0.5 ; BoxGun.eta[1]=0.5; BoxGun.prange[0]=2.8; BoxGun.prange[1]=2.83;" $simWorker --run ${runNumber}



  ##------ extract number of hits
  taskwrapper hitstats.log root -q -b -l ${O2_ROOT}/share/macro/analyzeHits.C
fi

if [ "$dodigi" == "1" ]; then
  echo "Running digitization for $intRate kHz interaction rate"
  intRate=$((1000*(intRate)));
  if [[ "$dotrdtrap" == "1" ]]; then trddigioption="--disable-trd-trapsim"; fi # no need to run the TRAP simulation twice
  taskwrapper digi.log o2-sim-digitizer-workflow --skipDet EMC,MFT,FDD,FV0,MID,MCH,CPV,PHS,ZDC  $gloOpt $trddigioption --interactionRate $intRate $tpcLanes --configKeyValues "HBFUtils.runNumber=${runNumber}" --early-forward-policy always --combine-devices
  echo "Return status of digitization: $?"
  # existing checks
  #root -b -q O2/Detectors/ITSMFT/ITS/macros/test/CheckDigits.C+
fi


#add CTP, CPV, PHS, ZDC

# CTP trengs av tpc

# --skipDet istedet --onlyDet HMP,TOF,ITS,TPC,TRD,FT0
#--skipDet EMC,MFT,FDD,FV0,MID,MCH 
if [ "$dotrdtrap" == "1" ]; then
  echo "Running TRD trap simulator"
  taskwrapper trdtrap.log o2-trd-trap-sim $gloOpt
  echo "Return status of trd trap sim: $?"
fi


if [ "$doreco" == "1" ]; then

  echo "Running TPC reco flow"
  #needs TPC digitized data
  taskwrapper tpcreco.log o2-tpc-reco-workflow $gloOpt --input-type digits --output-type clusters,tracks,send-clusters-per-sector  --configKeyValues "GPU_rec.maxTrackQPtB5=20"
  echo "Return status of tpcreco: $?"

  echo "Running ITS reco flow"
  taskwrapper itsreco.log  o2-its-reco-workflow --trackerCA --tracking-mode async $gloOpt $ITSRecOpt
  echo "Return status of itsreco: $?"

  # existing checks
  # root -b -q O2/Detectors/ITSMFT/ITS/macros/test/CheckClusters.C+
  # root -b -q O2/Detectors/ITSMFT/ITS/macros/test/CheckTracks.C+

  #echo "Running MFT reco flow"
  #needs MFT digitized data
  #MFTRecOpt=" --configKeyValues \"MFTTracking.forceZeroField=false;MFTTracking.LTFclsRCut=0.0100;\""
  #taskwrapper mftreco.log  o2-mft-reco-workflow  --disable-mc $gloOpt $MFTRecOpt # provde denne .. OK
  #echo "Return status of mftreco: $?"

  #echo "Running MCH reco flow"
  #taskwrapper mchreco.log o2-mch-reco-workflow $gloOpt
  #echo "Return status of mchreco: $?"

  echo "Running FT0 reco flow"
  #needs FT0 digitized data
  taskwrapper ft0reco.log o2-ft0-reco-workflow --disable-mc $gloOpt
  echo "Return status of ft0reco: $?"

  #echo "Running FDD reco flow"
  #needs FDD digitized data  
  #taskwrapper fddreco.log o2-fdd-reco-workflow --disable-mc $gloOpt
  #echo "Return status of fddreco: $?"

  #echo "Running FV0 reco flow"
  #needs FV0 digitized data
  #taskwrapper fv0reco.log o2-fv0-reco-workflow  --disable-mc $gloOpt
  #echo "Return status of fv0reco: $?"

  #echo "Running MID reco flow"
  #needs MID digitized data
  #taskwrapper midreco.log "o2-mid-digits-reader-workflow | o2-mid-reco-workflow $gloOpt"
  #echo "Return status of midreco: $?"

  echo "Running HMPID reco flow to produce clusters"
  #needs HMPID digitized data
  taskwrapper hmpreco.log "o2-hmpid-digits-to-clusters-workflow $gloOpt" # can not disable --disable-mc
  echo "Return status of hmpid cluster reco: $?"

  #echo "Running MCH-MID matching flow"
  #taskwrapper mchmidMatch.log "o2-muon-tracks-matcher-workflow $gloOpt"
  #echo "Return status of mchmidmatch: $?"

  echo "Running ITS-TPC matching flow"
  #needs results of o2-tpc-reco-workflow, o2-its-reco-workflow and o2-fit-reco-workflow
  taskwrapper itstpcMatch.log o2-tpcits-match-workflow --use-ft0  $gloOpt # can not disable --disable-mc
  echo "Return status of itstpcMatch: $?"

  echo "Running TRD matching to ITS-TPC and TPC"
  #needs results of o2-tpc-reco-workflow, o2-tpcits-match-workflow and o2-trd-tracklet-transformer
  taskwrapper trdTrkltTransf.log o2-trd-tracklet-transformer --disable-mc $gloOpt
  echo "Return status of trdTrkltTransf: $?"
  taskwrapper trdMatch.log o2-trd-global-tracking $gloOpt
  echo "Return status of trdTracker: $?"


  echo "Running TOF reco flow to produce clusters"
  #needs results of TOF digitized data and results of o2-tpcits-match-workflow
  taskwrapper tofReco.log o2-tof-reco-workflow  $gloOpt # can not disable --disable-mc
  echo "Return status of tof cluster reco : $?"

  echo "Running Track-TOF macthing flow"
  #needs results of TOF clusters data from o2-tof-reco-workflow and results of o2-tpc-reco-workflow and ITS-TPC matching
  taskwrapper tofMatchTracks.log o2-tof-matcher-workflow $gloOpt  # can not disable --disable-mc
  echo "Return status of o2-tof-matcher-workflow: $?"

  echo "Running Track-HMPID macthing flow"
  #needs results of HMPID clusters data from o2-hmpid-digits-to-clusters-workflow
  taskwrapper hmpidMatchTracks.log o2-hmpid-matcher-workflow $gloOpt
  echo "Return status of o2-hmpid-matcher-workflow: $?"
fi
