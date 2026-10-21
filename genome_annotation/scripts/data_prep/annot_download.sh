#!/bin/bash
# Download genome + annotation for each organism listed in sources.tsv 
cd "$(dirname "$0")"
tail -n +2 sources.tsv | while IFS=$'\t' read -r org rel g a; do
  mkdir -p "$org"
  for u in "$g" "$a"; do
    f=$org/$(basename "${u//%20/_}")
    [ -s "$f" ] && continue
    if curl -sfL --retry 3 -o "$f.part" "$u" && mv "$f.part" "$f"; then echo "ok $f"; else echo "FAIL $u"; fi
  done
done
date > download.done
