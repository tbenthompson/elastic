#include "filenames.h"
#include "UnitTest++.h"

TEST(RemoveExtension) {
    CHECK_EQUAL(remove_extension("abc.in"), "abc");
}

TEST(DispOutFilename) {
    CHECK_EQUAL(interior_disp_out_filename("def.in"), "def.disp_out_interior");
}
