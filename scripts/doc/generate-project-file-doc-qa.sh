#!/bin/bash

echo "======== $@"

if [ $# -ne 3 ]; then
    echo "USAGE: $0 SRCDIR BUILDDIR DATADIR" >&2
    exit 1
fi

srcdir="$1"
builddir="$2"
datadir="$3"

doxdir="$builddir/DocAux/dox"
toolsdir="$srcdir/scripts/doc"

mkdir -p "$doxdir"

qafile="$doxdir/project-file-doc-qa.dox"

cat <<"EOF" >"$qafile"
/*! \page project_file_doc_qa Project File Parameters&mdash;Quality Assurance

This is the QA page

EOF

"$toolsdir/get-project-params.sh" "$srcdir" "$doxdir/ProjectFile" >>"$qafile"

cat <<EOF >>"$qafile"

*/
EOF

"$toolsdir/append-xml-tags.py" prj "$datadir" "$doxdir/ProjectFile"
