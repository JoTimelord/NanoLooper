#!/bin/bash


rm -f .jobs.txt
rm -r ~/NanoLooper/logfiles
rm -r ~/NanoLooper/outputs

mkdir ~/NanoLooper/logfiles
mkdir ~/NanoLooper/outputs

cp ~/ranscripts/.jobs.txt ~/NanoLooper/.jobs.txt

xargs.sh .jobs.txt
