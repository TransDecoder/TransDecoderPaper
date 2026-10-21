#!/bin/bash
# usage: download.sh <outdir> <url> <md5>
set -u
out=$1; url=$2; md5=$3; f=$out/$(basename $url)
for i in 1 2 3 4 5; do
  curl -sSL --retry 5 -C - -o $f "$url" && break
  sleep 30
done
got=$(md5sum $f | cut -d' ' -f1)
if [ "$got" == "$md5" ]; then echo "MD5 OK $f $got"; else echo "MD5 FAIL $f got=$got expected=$md5"; fi
