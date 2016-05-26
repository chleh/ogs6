#!/usr/bin/python

import sys
import re
import os.path

def debug(msg):
    sys.stderr.write(msg+"\n")

def write_out(*args):
    print("@@@".join([str(a) for a in args]))

# capture #1 is the parameter path
comment = re.compile(r"^//! \\ogs_project_file_parameter\{([A-Za-z_0-9]+)\}( \\todo .*)?$")

# capture #5 is the parameter name
getter = re.compile(r'^(get|check|ignore|peek)Conf(Param|Attribute|Subtree)(List|Optional|All)?'
                   +r'(<.*>)?'
                   +r'\("([a-zA-Z_0-9:]+)"[,)]')

state = "getter"
path = ""
lineno = 0
line = ""
tag_path_comment = ""

for inline in sys.stdin:
    oldpath = path; oldlineno = lineno; oldline = line

    path, lineno, line = inline.split(":", 2)

    if path != oldpath: debug(path)

    line = line.strip()
    lineno = int(lineno)

    m = comment.fullmatch(line)
    if m:
        if state != "getter":
            write_out("UNNEEDED", oldpath, oldlineno, oldline)
        state = "comment"

        tag_path_comment = m.group(1).replace("__", ".")
        debug(" {:>5}  //! {}".format(lineno, tag_path_comment))
        tag_name_comment = tag_path_comment.split(".")[-1]

        continue

    m = getter.match(line)
    if m:
        param = m.group(5)
        paramtype = m.group(4)[1:-1] if m.group(4) else ""
        method = m.group(1) + "Conf" + m.group(2) + (m.group(3) or "")

        if state != "comment" or oldpath != path:
            write_out("NODOC", path, lineno, "NONE", param, paramtype, method)
        else:
            debug(" {:>5}  {} {} ".format(lineno, param, paramtype))

            if param != tag_name_comment:
                debug("error: parameter name from comment and code do not match: "
                        + tag_name_comment + " vs. " + param)
                write_out("NODOC", path, lineno, tag_path_comment, param, paramtype, method)
            elif lineno != oldlineno+1:
                debug("error: the associated comment is not on the line preceding this one."
                        + " line numbers {} vs. {}".format(oldlineno, lineno))
                write_out("NODOC", path, lineno, tag_path_comment, param, paramtype, method)
            else:
                write_out("OK", path, lineno, tag_path_comment, param, paramtype, method)

        state = "getter"
        continue

    write_out("WRONGIN", path, lineno, line.strip())
    state = "getter" # reset state in order to avoid warnings
