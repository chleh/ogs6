#!/bin/bash

if [ $# -ne 1 ]; then
    echo "USAGE: ${0##*/} SRCDIR" >&2
    exit 1
fi

srcdir="`readlink -f "$1"`"

#color="--color=always"
color=""

cat <<"EOF" \
| grep -r $srcdir \
    --include '*.h' \
    --include '*.cpp' \
    --exclude-dir '.git' \
    --exclude-dir 'Tests' \
    --exclude 'ConfigTree*.*' \
    -f - -r -n -o $color \
| cut -c $((`expr length "$srcdir"` + 2))-
^\s*//! \\ogs_file_\(param\|attr\){[A-Za-z_0-9]\+}\( \\todo .*\)\?$
^\s*//! \\ogs_file_special$
^\s*//! \\ogs_file_\(param\|attr\)_special{[A-Za-z_0-9]\+}\( \\todo .*\)\?$
checkConfParam.*)
getConfAttribute.*)
getConfParam.*)
getConfSubtree.*)
ignoreConfAttribute.*)
ignoreConfParam.*)
peekConfParam.*)
EOF

# format as table:
# | sed -e 's_::_@@_g' -e's_:\s\+_:_' | column -t -s: | sed -e 's_@@_::_g'
