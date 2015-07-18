#!/bin/csh -f
#Keeps running and automatically deletes useless, extraneous result information from the C3 server.
#This is to prevent server storage from overflowing. (As of the time of creation, the C3 server has only 255 GB of storage)
#
#Note: this csh script must be run from the /home/user081/ESP_permutation_analysis folder
#
#Author: Jonathan Ma
#Created: 12/25/2014 (Merry Christmas!)

while(1)
 foreach cc (`cat allcancers.txt`)
  echo "now cleaning $cc"
  set parentdir="Permuted_Results_SAM/$cc"
  @ j = 1
  while( $j <= 1000)
   #caution: the following directories are relative to /home/user081/ESP_permutation_analysis
   set maledir = "$parentdir/male${j}/"
   set femaledir = "$parentdir/female${j}/"
   set malefiledir = "${maledir}summary/1/result_NA_summly_n7503.txt"
   set femalefiledir = "${femaledir}summary/1/result_NA_summly_n7503.txt"
   set newmaledir="$parentdir/male${j}_result/"
   set newfemaledir="$parentdir/female${j}_result/"

   set garbagemalefiledir="${maledir}/summary/1/oDataSet.json"
   set garbagefemalefiledir="${femaledir}/summary/1/oDataSet.json"

   if( -e $malefiledir ) then
    #get rid of other useless garbage

    rm -rf "${maledir}/config.yaml"
    rm -rf "${maledir}/index.html"
    rm -rf "${maledir}/query_result"
    rm -rf "${maledir}/summary/rpt_full.mat"
    rm -rf "${maledir}/summary/summly_rpt.mat"
    rm -rf "${maledir}/summary/query_result_n1x69761.gctx"
    rm -rf "${maledir}/summary/sig_summly_tool_params.txt"
    rm -rf "${maledir}/summary/config.yaml"
    rm -rf "${maledir}/summary/index.html"
    rm -rf "${maledir}/summary/1/index.html"
    rm -rf "${maledir}/summary/1/oDataSet.json"
    rm -rf "${maledir}/summary/1/oFigures.json"
    rm -rf "${maledir}/summary/1/query_info.txt"
    echo "${cc} male${j} cleaned up"
   endif

   if( -e $femalefiledir ) then
    #get rid of other useless garbage

    rm -rf "${femaledir}/config.yaml"
    rm -rf "${femaledir}/index.html"
    rm -rf "${femaledir}/query_result"
    rm -rf "${femaledir}/summary/rpt_full.mat"
    rm -rf "${femaledir}/summary/summly_rpt.mat"
    rm -rf "${femaledir}/summary/query_result_n1x69761.gctx"
    rm -rf "${femaledir}/summary/sig_summly_tool_params.txt"
    rm -rf "${femaledir}/summary/config.yaml"
    rm -rf "${femaledir}/summary/index.html"
    rm -rf "${femaledir}/summary/1/index.html"
    rm -rf "${femaledir}/summary/1/oDataSet.json"
    rm -rf "${femaledir}/summary/1/oFigures.json"
    rm -rf "${femaledir}/summary/1/query_info.txt"
    echo "${cc} female${j} cleaned up"

   endif

   @ j = $j + 1
  end
 end
end
