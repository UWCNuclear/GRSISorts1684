#!/bin/bash 
#./runPb_Sr82.sh

DATADIR=/home/nico/Nikita/S1684
FRAGDIR=FragmentTrees
ANALDIR=AnalysisTrees
CALFILE=S1684new.cal
#CALFILE=S1684safe.cal

DLEN=${#DATADIR}

# Add $DATADIR/runXXXXX_*.mid to list below, where XXXXX is the run number
 for f in $DATADIR/run44471_*.mid

do
 g=${f:DLEN+4}
 h=${g:0:${#g}-4} 
 i=${h:6:${#h}-6}
 FFILE=$FRAGDIR/fragment$h.root
 AFILE=$ANALDIR/analysis$h.root

 echo "Processing run$g"
    
#if [ "$AFILE" -ot "$f" ];
#		then
#			echo "File $AFILE exists but is older than $f"
#			echo "grsisort -laq $f $CALFILE --sort-depth 1000000000"
#		 	grsisort -laq $f $CALFILE --sort-depth 1000000000 --suppress_errors
#			mv -f analysis$h.root $ANALDIR
#			mv -f fragment$h.root $FRAGDIR
#fi

if [ ! -f $AFILE ];  
		then
			echo "File $AFILE does not exist."
			echo "grsisort -laq $f $CALFILE --sort-depth 1000000000"
		 	grsisort -laq $f $CALFILE --sort-depth 1000000000 --suppress_errors
			mv -f analysis$h.root $ANALDIR
			mv -f fragment$h.root $FRAGDIR

	fi
	
	if [ -f $AFILE ];
		then
			SortCodes/SortData $ANALDIR/analysis$h.root $CALFILE Histograms.root 82Sr 208Pb
  	  mv Histograms.root HistFiles/Histograms_NoLasers_82Sr_Pb_$h.root
	fi 

done

gadd -f Sr82_Pb_NoLasers_Summed.root HistFiles/Histograms_NoLasers_82Sr_Pb_*.root



