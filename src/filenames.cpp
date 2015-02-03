#include "filenames.h"

std::string remove_extension(const std::string& filename) 
{
    size_t last_dot = filename.find_last_of(".");
    if (last_dot == std::string::npos) {
        return filename;
    }
    return filename.substr(0, last_dot); 
}

std::string disp_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".disp_out";
}

std::string interior_disp_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".disp_out_interior";
}

std::string trac_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".trac_out";
}

std::string interior_trac_out_filename(const std::string& filename) {
    auto in_filename_root = remove_extension(filename);
    return in_filename_root + ".trac_out_interior";
}

