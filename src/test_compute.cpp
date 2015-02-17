#include "UnitTest++.h"
#include "compute.h"

using namespace tbem;

TEST(DiagBlock) {
    auto empty = BlockOperator{0,0,{}};
    ComputedEquation lhs{
        {empty, "A", "B", "C"},
        {empty, "A", "A", "E"},
        {empty, "A", "C", "D"}
    };
    BlockFunction rhs{};
    LinearSystem ls{lhs, rhs};
    CHECK_EQUAL(ls.get_diag_block().function, "E");
}

TEST(DiagBlockException) {
    ComputedEquation lhs{};
    BlockFunction rhs{};
    LinearSystem ls{lhs, rhs};
    CHECK_THROW(ls.get_diag_block(), NoDiagonalException);
}
