#include "UnitTest++.h"
#include <iostream>
#include "vec.h"
#include "elastic.h"

TEST(SimpleLoadGood) {
    parse_json(load_file("test_data/good.in"));
}

TEST(SimpleLoadBad) {
    CHECK_THROW(parse_json(load_file("test_data/bad.in")), std::invalid_argument);
    CHECK_THROW(parse_json(load_file("doesnotexist")), std::invalid_argument);
}

int main() {
    return UnitTest::RunAllTests();
}
