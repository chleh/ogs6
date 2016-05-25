#!/bin/bash

#color="--color=always"
color=""

cat <<"EOF" \
| grep -r . \
    --include '*.h' \
    --include '*.cpp' \
    --exclude-dir '.git' \
    --exclude-dir 'Tests' \
    --exclude 'ConfigTree*.*' \
    -f - -r -n -o $color \
| sed -e 's_::_@@_g' -e's_:\s\+_:_' | column -t -s: | sed -e 's_@@_::_g'
^\s*//! \\ogs_project_file_parameter{[a-z_0-9]\+}\( \\todo project_file_docu\)\?$
checkConfParam.*)
getConfAttribute.*)
getConfParam.*)
getConfSubtree.*)
ignoreConfParam.*)
peekConfParam.*)
EOF
