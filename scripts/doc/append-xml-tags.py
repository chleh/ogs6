#!/usr/bin/python

# prevent broken pipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

import os

import xml.etree.cElementTree as ET

import argparse

parser = argparse.ArgumentParser(description="Print XML tags")

parser.add_argument("ext",     help="Extension of files to consider")
parser.add_argument("datadir", help="data directory")
parser.add_argument("docdir",  help="doc output directory")

args = parser.parse_args()
extension = '.' + args.ext
datadir = os.path.abspath(args.datadir)
docdir  = os.path.abspath(args.docdir)

# maps tags to the set of xml files they appear in
dict_tag_files = dict()

def dict_of_set_append(dict_, key, value):
    if key in dict_:
        dict_[key].add(value)
    else:
        dict_[key] = set((value,))


def print_tags(node, path, level, filepath):
    global dict_tag_files

    tag = node.tag
    if level>1: # skip root node
        tagpath = path + "." + tag
    else:
        tagpath = tag

    if level>0: # skip root node
        dict_of_set_append(dict_tag_files, (True, tagpath), filepath)
        for k in node.attrib:
            dict_of_set_append(dict_tag_files, (False, tagpath + "." + k), filepath)

    for child in node:
        print_tags(child, tagpath, level + 1, filepath)


print("data dir", datadir)
for (dirpath, _, filenames) in os.walk(datadir):
    for f in filenames:
        if not f.endswith(extension): continue

        filepath = os.path.join(dirpath, f)
        xmlroot = ET.parse(filepath).getroot()
        print_tags(xmlroot, "", 0, filepath[len(datadir)+1:])

if False:
    first = True
    for (tag, files) in sorted(dict_tag_files.items()):
        if first:
            first = False
        else:
            print()

        print("T |" if tag[0] else "A |", tag[1])
        for f in sorted(files):
            print("   ", f)


for (dirpath, _, filenames) in os.walk(docdir):
    reldirpath = dirpath[len(docdir)+1:]
    istag = True

    for f in filenames:
        if not f.endswith(".dox"): continue

        if f.startswith("__"):
            tagpath = reldirpath
        elif f.startswith("t_"):
            tagpath = os.path.join(reldirpath, f[2:-len(".dox")])
            istag = True
        elif f.startswith("a_"):
            tagpath = os.path.join(reldirpath, f[2:-len(".dox")])
            istag = False

        # TODO make work for IC etc, too
        tagpath = tagpath.replace(os.sep, ".")
        tagpath = ".".join(tagpath.split(".")[1:])

        path = os.path.join(dirpath, f)
        with open(path, "a") as fh:
            if tagpath:
                fh.write("\n# Used in the following test data files\n\n")
                try:
                    datafiles = dict_tag_files[(istag, tagpath)]

                    for df in sorted(datafiles):
                        fh.write("- {}\n".format(df))
                except KeyError:
                    fh.write("Used in no data files.\n")
            else:
                # no additional output for the main doc page
                pass

            fh.write("\n*/\n")
