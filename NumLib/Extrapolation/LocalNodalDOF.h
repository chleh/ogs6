#pragma once


namespace NumLib {


class LocalNodalDOF
{
public:
    virtual std::vector<double> const& getElementNodalValues() = 0;
    virtual std::vector<double> const& getElementNodalValues(unsigned component) = 0;
};

}
