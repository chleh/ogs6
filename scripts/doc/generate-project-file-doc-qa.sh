#!/bin/bash

echo "======== $@"

if [ $# -ne 3 ]; then
    echo "USAGE: $0 SRCDIR BUILDDIR DATADIR" >&2
    exit 1
fi

srcdir="$1"
builddir="$2"
datadir="$3"

docauxdir="$builddir/DocAux"
doxdir="$docauxdir/dox"
toolsdir="$srcdir/scripts/doc"

param_cache="$docauxdir/documented-parameters-cache.txt"

qafile="$doxdir/project-file-doc-qa.dox"
check_quality_script="`dirname "$0"`/check-project-params.py"

mkdir -p "$doxdir"

# gather information about documented parameters
"$toolsdir/get-project-params.sh" "$srcdir" >"$param_cache"

# write QA information
cat <<"EOF" >"$qafile"
/*! \page project_file_doc_qa Project File Parameters&mdash;Quality Assurance

This is the QA page

EOF

cat "$param_cache" | "$check_quality_script" "$doxdir/ProjectFile" >>"$qafile"

cat <<EOF >>"$qafile"

*/
EOF

# finish parameter documentation dox files
"$toolsdir/append-xml-tags.py" prj "$datadir" "$docauxdir"
