#ifndef __TKJLQWENFjas12_FILENAMES_H
#define __TKJLQWENFjas12_FILENAMES_H

#include <string>

std::string remove_extension(const std::string& filename); 
std::string disp_out_filename(const std::string& filename);
std::string interior_disp_out_filename(const std::string& filename); 
std::string trac_out_filename(const std::string& filename); 
std::string interior_stress_out_filename(const std::string& filename);
std::string slip_out_filename(const std::string& filename);

#endif
