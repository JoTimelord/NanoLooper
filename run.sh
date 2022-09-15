#!/bin/bash


rm -f .jobs.txt
cp ~/ranscripts/.jobs.txt ~/NanoLooper/.jobs.txt

xargs.sh .jobs.txt
