#!/usr/bin/python

import sys
import re
import os.path

github_src_url = "https://github.com/ufz/ogs/tree/master"

def debug(msg):
    sys.stderr.write(msg+"\n")

if len(sys.argv) != 2:
    print("USAGE: {} DOCAUXDIR".format(sys.argv[0]))
    sys.exit(1)

docauxdir = sys.argv[1]
if not os.path.isdir(docauxdir):
    print("error: `{}' is not a directory".format(docauxdir))
    sys.exit(1)

excluded = set(("MathLib/LinAlg/LinearSolverOptions.cpp", ))

# capture #1 is the parameter path
comment = re.compile(r"^//! \\ogs_project_file_parameter\{([A-Za-z_0-9]+)\}$")

# capture #4 is the parameter name
getter = re.compile(r'^(get|check|ignore|peek)Conf(Param|Attribute|Subtree)(List|Optional)?'
                   +r'(<.*>)?'
                   +r'\("([a-zA-Z_0-9:]+)"[,)]')


undocumented = []
unneeded_comments = []
wrong_input = []

no_doc_page = []


state = "getter"
path = ""
lineno = 0
line = ""
tag_path_comment = ""

for inline in sys.stdin:
    oldpath = path; oldlineno = lineno; oldline = line

    path, lineno, line = inline.split(":", 2)
    if path in excluded: continue

    if path != oldpath: debug(path)

    line = line.strip()
    lineno = int(lineno)

    m = comment.fullmatch(line)
    if m:
        if state != "getter":
            unneeded_comments.append((oldpath, oldlineno, oldline))
        state = "comment"

        tag_path_comment = m.group(1).replace("__", ".")
        debug(" {:>5}  //! {}".format(lineno, tag_path_comment))
        tag_name_comment = tag_path_comment.split(".")[-1]

        dirs = tag_path_comment.split(".")[:-1]
        p = os.path.join(docauxdir, *dirs, )
        if     (not os.path.isfile(os.path.join(p,                   "t_" + tag_name_comment + ".dox"))) \
           and (not os.path.isfile(os.path.join(p, tag_name_comment, "i_" + tag_name_comment + ".dox"))) \
           and (not os.path.isfile(os.path.join(p, tag_name_comment, "c_" + tag_name_comment + ".dox"))) :
            no_doc_page.append((tag_path_comment, oldpath, oldlineno))

        continue
    
    m = getter.match(line)
    if m:
        param = m.group(5)
        paramtype = m.group(4)[1:-1] if m.group(4) else ""
        method = m.group(1) + "Conf" + m.group(2) + (m.group(3) or "")

        if state != "comment" or oldpath != path:
            undocumented.append((path, lineno, param, paramtype, method))
        else:
            debug(" {:>5}  {} {} ".format(lineno, param, paramtype))

            if param != tag_name_comment:
                debug("error: parameter name from comment and code do not match:",
                        tag_name_comment, "vs.", param)
                undocumented.append((path, lineno, param, paramtype))
            elif lineno != oldlineno+1:
                debug("error: the associated comment is not on the line preceding this one."
                        + " line numbers {} vs. {}".format(oldlineno, lineno))
                undocumented.append((path, lineno, param, paramtype))

        state = "getter"
        continue

    wrong_input.append(inline.strip())


if (undocumented):
    print()
    print("# Undocumented parameters")
    print("| File | Line | Parameter | Type | Method | Link |")
    print("| ---- | ---- | --------- | ---- | ------ | ---- |")
    for u in sorted(undocumented):
        print(("| {0} | {1} | {2} | <tt>{3}</tt> | <tt>{4}</tt> "
                + "| [&rarr; ufz/ogs/master]({5}/{0}#L{1})").format(*u, github_src_url))

if (unneeded_comments):
    print()
    print("# Comments not documenting anything")
    print("| File | Line | Comment |")
    print("| ---- | ---- | ------- |")
    for u in sorted(unneeded_comments):
        print("| {} | {} | {} |".format(*u))

if (wrong_input):
    print()
    print("# Lines of input to that script that have not been recognized")
    for u in sorted(wrong_input):
        print(" - ", u)

if (no_doc_page):
    print()
    print("# No documentation page")
    print("| Parameter | File | Line |")
    print("| --------- | ---- | ---- |")
    for n in sorted(no_doc_page):
        print("| {} | {} | {} |".format(*n))


# exit with error status if something was not documented.
if (not not undocumented) or (not not unneeded_comments) \
        or (not not wrong_input) or (not not no_doc_page):
            sys.exit(1)

sys.exit(0)

