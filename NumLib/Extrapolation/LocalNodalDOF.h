#pragma once


namespace NumLib {


class LocalNodalDOF
{
    virtual std::vector<double> getElementNodalValues() = 0;
    virtual std::vector<double> getElementNodalValues(unsigned component) = 0;
};

}
