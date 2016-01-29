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

qafile="$doxdir/project-file-doc-qa.dox"

cat <<"EOF" | tee /dev/stderr >"$qafile"
/*! \page project_file_doc_qa ProjectFile Documentation Quality Assurance

This is the QA page

EOF

"$toolsdir/get-project-params.sh" "$srcdir" "$doxdir" >>"$qafile"

cat <<EOF >>"$qafile"

*/
EOF

