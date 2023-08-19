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
    Usage
fi

[[ -z $startTime ]] && startTime=$startTimeDef
[[ -z $runNumber ]] && runNumber=$runNumDef

dosim="0"
dodigi="0"
dotrdtrap="0"
doreco="0"
# convert to lowercase
fromstage=`echo "$fromstage" | awk '{print tolower($0)}'`
if [ "$fromstage" == "sim" ]; then
  dosim="1"
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
  Usage
fi


if [ "$doreco" == "1"  ]; then


  echo "Running Track-HMPID macthing flow"
  #needs results of HMPID clusters data from o2-hmpid-digits-to-clusters-workflow
  taskwrapper hmpidMatchTracks.log o2-hmpid-matcher-workflow $gloOpt
  echo "Return status of o2-hmpid-matcher-workflow: $?"


  # the strangeness trackin is now called from the secondary-vertexing. To enable it as a standalone workflow
  # one should run the previous o2-secondary-vertexing-workflow with options
  # --configKeyValues "svertexer.createFullV0s=true;svertexer.createFullCascades=true;svertexer.createFull3Bodies=true" --disable-strangeness-tracker
  #  echo "Running strangeness tracking flow"
  #  #needs results of S.Vertexer + ITS reco
  #  taskwrapper sttracking.log o2-strangeness-tracking-workflow $gloOpt
  #  echo "Return status of strangeness tracking: $?"

  echo "Producing AOD"
  taskwrapper aod.log o2-aod-producer-workflow $gloOpt --aod-writer-keep dangling --aod-writer-resfile "AO2D" --aod-writer-resmode UPDATE --aod-timeframe-id 1 --run-number 300000
  echo "Return status of AOD production: $?"

  # let's do some very basic analysis tests (mainly to enlarge coverage in full CI) and enabled when SIM_CHALLENGE_ANATESTING=ON
  if [[ ${O2DPG_ROOT} && ${SIM_CHALLENGE_ANATESTING} ]]; then
    # to be added again: Efficiency
    for t in ${ANATESTLIST:-MCHistograms Validation PIDTOF PIDTPC EventTrackQA WeakDecayTutorial}; do
      ${O2DPG_ROOT}/MC/analysis_testing/analysis_test.sh ${t}
      echo "Return status of ${t}: ${?}"
    done
  fi
fi
