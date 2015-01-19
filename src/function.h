#ifndef __yLADSKJH_FUNCTION_H
#define __yLADSKJH_FUNCTION_H
#include <vector>
#include <cassert>
#include <iostream>

typedef std::vector<std::vector<double>> Function;

inline std::ostream& operator<<(std::ostream& os, const Function& a) {
    os << "(";
    for (size_t d = 0; d < a.size(); d++) {
        os << "[";
        for (size_t i = 0; i < a[d].size(); i++) {
            os << a[d][i] << ", ";
        }
        os << "], ";
    }
    os << ")";
    return os;
}

//TODO: Maybe get rid of this. Interesting macro experiment. What about python template?
#define __FunctionLoopBegin \
    assert(a.size() == b.size());\
    for (size_t d = 0; d < a.size(); d++) {\
        assert(a[d].size() == b[d].size());\
        for (size_t i = 0; i < a[d].size(); i++) {

#define __FunctionLoopEnd \
        } \
    }

inline Function& operator+=(Function& a, const Function& b) {
    __FunctionLoopBegin;
    a[d][i] += b[d][i];
    __FunctionLoopEnd;
    return a;
}


inline Function& operator-=(Function& a, const Function& b) {
    __FunctionLoopBegin;
    a[d][i] -= b[d][i];
    __FunctionLoopEnd;
    return a;
}

inline Function& operator*=(Function& a, const double b) {
    for (size_t d = 0; d < a.size(); d++) {\
        for (size_t i = 0; i < a[d].size(); i++) {
            a[d][i] *= b;
        }
    }
    return a;
}

inline Function operator*(const Function& a, double b) {
    Function out = a;
    out *= b;
    return out;
}

inline Function constant_function(size_t components, size_t dofs, double value) {
    return Function(components, std::vector<double>(dofs, value));
}

inline Function operator-(const Function& a) {
    Function out(a.size());
    for (size_t d = 0; d < a.size(); d++) {
        out[d].resize(a[d].size());
        for (size_t i = 0; i < a[d].size(); i++) {
            out[d][i] = -a[d][i];  
        }
    }
    return out;
}

inline bool operator==(Function& a, const Function& b) {
    __FunctionLoopBegin;
    if (a[d][i] != b[d][i]) {
        return false;
    }
    __FunctionLoopEnd;
    return true;
}

#endif
