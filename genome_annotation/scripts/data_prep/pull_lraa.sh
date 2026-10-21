#!/bin/bash
set -uo pipefail
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER/data/gm12878/longread_lraa
export PATH=/broad/hptmp/bhaas/tools/apptainer-1.5.3/bin:$PATH
export APPTAINER_CACHEDIR=$P/apptainer_cache APPTAINER_TMPDIR=$P/tmp TMPDIR=$P/tmp
export APPTAINER_DOCKER_USERNAME=oauth2accesstoken
export APPTAINER_DOCKER_PASSWORD="$(gcloud auth print-access-token 2>/dev/null)"
URI=us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:0.44.3
apptainer pull --force $P/image/lraa-core_0.44.3.sif docker://$URI && echo PULL_DONE
gcloud artifacts docker images describe "$URI" --format='value(image_summary.digest)'
sha256sum $P/image/lraa-core_0.44.3.sif
