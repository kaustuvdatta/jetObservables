#!/bin/sh
while true
do
    nohup ./resub_crab.sh > nohup.out &
    sleep 900
done

