#!/bin/sh

# expect input from get-project-params.sh

base="Documentation/ProjectFile"

while read -r fn lno content; do
    [ "$content" = "${content#//!}" ] && continue
    tag_name="$(echo "$content" | sed -n -e 'sX^//! \\ogs_project_file_parameter{\([A-Za-z_0-9]\+\)}$X\1Xp')"
    [ -z "$tag_name" ] && continue
    tag_name="${tag_name//__/\/}"
    echo "$base/$tag_name"
done \
| sort -r \
| while read path; do
    dn="`dirname "$path"`"
    if [ ! -d "$dn" ]; then
        mkdir -p "$dn"

        bn="`basename "$path"`"
        echo "creating $path/_$bn.md"
        echo '\todo document' >"$path/_$bn.md"
    fi

    if [ -d "$path" ]; then
        bn="`basename "$path"`"
        if ! [ -f "$path/_$bn.md" ]; then
            echo "creating $path/_$bn.md"
            echo '\todo document' >"$path/_$bn.md"
        fi
    elif ! [ -f "$path.md" ]; then
        echo "creating $path.md"
        echo '\todo document' >"$path.md"
    fi

    if [ -d "$path" ] && [ -f "$path.md" ]; then
        echo "ERROR: both $path and $path.md exist!" >&2
    fi
done
