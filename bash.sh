#!/bin/bash

for ((i=1; i<=2000; i++))
do
    ps aux | grep "kuangguopeng" | awk '{print $2}' | xargs kill -9
done
