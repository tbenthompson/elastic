#include "UnitTest++.h"
#include "function.h"

struct Data {
    Function a{ {1,2}, {3,4} };
    Function b{ {-1,-2}, {-3,-4} };
    Function c{ {0,0}, {0,0} };
};

TEST_FIXTURE(Data, FunctionAdd) {
    a += b;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, FunctionSub) {
    a -= -b;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, FunctionMul) {
    a *= 2;
    a += b * 2;
    CHECK_EQUAL(a, c);
}
