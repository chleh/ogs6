#pragma once

#include <cassert>
#include <cmath>

namespace MathLib
{

namespace Nonlinear
{

namespace detail
{

inline
bool same_sign(const double a, const double b)
{
    return (a < 0.0) == (b < 0.0);
}

}

class Bisect
{
public:
    typedef double (*const Function)(double);

    Bisect(Function f, double a, double b)
        : _f(f), _a(a), _b(b), _fa(f(a)), _fb(f(b))
    {
        if (_fa == 0.0)
            _b = _a;
        else if (_fb == 0.0)
            _a = _b;
        else
            assert(!detail::same_sign(_fa, _fb));
    }

    void step(const unsigned num_steps)
    {
        for (unsigned i=0; i<num_steps; ++i)
        {
            const double c = 0.5*(_a+_b);
            const double fc = _f(c);

            if (fc == 0.0) {
                _a = _b = c;
                return;
            } else if (detail::same_sign(_fa, fc)) {
                _a = c;
                _fa = fc;
            } else {
                _b = c;
                _fb = fc;
            }
        }
    }

    double get_result() const { return 0.5*(_a+_b); }

    double get_range() const { return std::fabs(_a-_b); }


private:
    const Function _f;
    double _a, _b, _fa, _fb;
};


template<typename T>
class RegulaFalsi : private T
{
public:
    typedef double (*const Function)(double);

    RegulaFalsi(Function f, double a, double b)
        : _f(f), _a(a), _b(b), _fa(f(a)), _fb(f(b))
    {
        if (_fa == 0.0) {
            _b = _a;
        } else if (_fb == 0.0) {
            _a = _b;
        } else {
            assert(!detail::same_sign(_fa, _fb));
        }
    }

    void step(const unsigned num_steps)
    {
        for (unsigned i=0; i<num_steps; ++i)
        {
            if (_a == _b) return;

            const double s = (_fb - _fa)/(_b - _a);
            const double c = _a - _fa/s;
            const double fc = _f(c);

            if (fc == 0.0) {
                _a = _b = c;
                return;
            } else if (!detail::same_sign(fc, _fb)) {
                _a = _b;
                _fa = _fb;
                _b = c;
                _fb = fc;
            } else {
                const double m = T::get_m(_fa, _fb, fc);
                _fa *= m;
                _b = c;
                _fb = fc;
            }
        }
    }

    double get_result() const
    {
        if (_a == _b) return _a;

        const double s = (_fb - _fa)/(_b - _a);
        const double c = _a - _fa/s;

        return c;
    }

    double get_range() const { return std::fabs(_a - _b); }

private:
    const Function _f;
    double _a, _b, _fa, _fb;
};

class Unmodified {
protected:
    double get_m(const double /*fa*/, const double /*fb*/, const double /*fc*/) const { return 1.0; }
};

class Illinois {
protected:
    double get_m(const double /*fa*/, const double /*fb*/, const double /*fc*/) const { return 0.5; }
};

class Pegasus {
protected:
    double get_m(const double /*fa*/, const double fb, const double fc) const { return fb / (fb+fc); }
};

class AndersonBjorck {
protected:
    double get_m(const double /*fa*/, const double fb, const double fc) const {
        const double f = 1.0 - fc / fb;
        return (f >= 0.0) ? f : 0.5;
    }
};

}

} // namespace MathLib

