#!/bin/bash

cat <<EOF \
| grep -r . \
    --include '*.h' \
    --include '*.cpp' \
    --exclude-dir '.git' \
    --exclude-dir 'Tests' \
    --exclude 'ConfigTreeNew*.*' \
    -f - -r -w -o
getConfParam.*)
getConfSubtree.*)
peekConfParam.*)
checkConfParam.*)
ignoreConfParam.*)
EOF

