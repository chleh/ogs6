/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TestTools.h"

boost::property_tree::ptree
readXml(const char xml[])
{
    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree,
             boost::property_tree::xml_parser::no_comments |
             boost::property_tree::xml_parser::trim_whitespace);
    return ptree;
}
