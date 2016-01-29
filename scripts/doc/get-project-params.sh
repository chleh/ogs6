#!/bin/bash

if [ $# -ne 2 ]; then
    echo "USAGE: ${0##*/} SRCDIR DOCAUXDIR" >&2
    exit 1
fi

further_script="`dirname "$0"`/check-project-params.py"
srcdir="`readlink -f "$1"`"
docauxdir="`readlink -f "$2"`"

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
| cut -c $((`expr length "$srcdir"` + 2))- \
| "$further_script" "$docauxdir"
^\s*//! \\ogs_project_file_parameter{[A-Za-z_0-9]\+}\( \\todo project_file_docu\)\?$
checkConfParam.*)
getConfAttribute.*)
getConfParam.*)
getConfSubtree.*)
ignoreConfAttribute.*)
ignoreConfParam.*)
peekConfParam.*)
EOF

# exits with the status of "$further_script"
exit $?

# format as table:
# | sed -e 's_::_@@_g' -e's_:\s\+_:_' | column -t -s: | sed -e 's_@@_::_g'
