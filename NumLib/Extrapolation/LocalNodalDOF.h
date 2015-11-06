#pragma once


namespace NumLib {


class LocalNodalDOF
{
public:
    virtual std::vector<double> getElementNodalValues() const = 0;
    virtual std::vector<double> getElementNodalValues(unsigned component) const = 0;
};

}
