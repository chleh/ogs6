#pragma once


namespace NumLib {


class LocalNodalDOF
{
public:
    virtual std::vector<double> getElementNodalValues() = 0;
    virtual std::vector<double> getElementNodalValues(unsigned component) = 0;
};

}
