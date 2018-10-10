#!/bin/bash

inname=$1
outname=${2:-eventList.txt}
treename=${3:-eventTree}
directory=${4:-ttHtaus}

root -b $inname <<EOF
$directory->cd();
$treename->SetScanField(0);
.> $outname
$treename->Scan("run:ls:nEvent","","colsize=32");
.>
EOF
