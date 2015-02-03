#include "filenames.h"
#include "UnitTest++.h"

TEST(RemoveExtension) {
    CHECK_EQUAL(remove_extension("abc.in"), "abc");
}
