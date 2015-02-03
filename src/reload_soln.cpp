#include "reload_soln.h"
#include <hdf5.h>
#include <iostream>

std::string get_dataset_name(hid_t gid, size_t index) {
    const size_t max_len = 1024;
    char dataset_name[max_len];
    H5Gget_objname_by_idx(gid, index, dataset_name, max_len);
    return std::string(dataset_name);
}

size_t n_elements(hid_t dataset_id) 
{
    hsize_t shape[2];
    auto filespace = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace, shape, NULL);
    return shape[0];
}

std::vector<double> load_function(hid_t file_id, const std::string& dataset_name_str) 
{
    auto dataset_name_with_group = "/" + dataset_name_str;
    auto dataset_id = H5Dopen(file_id, dataset_name_with_group.c_str(), H5P_DEFAULT);

    std::vector<double> out(n_elements(dataset_id));
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    return out;
}

bool startswith(const std::string& str, const std::string& prefix) {
    return str.substr(0, prefix.size()) != prefix;
}

size_t dim_from_dataset_name(const std::string& name, const std::string& prefix) {
    std::string dim_str = name.substr(prefix.size());
    return std::stoi(dim_str);
}

//TODO: Move this to the library. Maybe in the IO file
//TODO: Write a test that saves and reloads a file and compares the values. Could this
//be done without actually writing to disk?
std::vector<std::vector<double>> load_surface(const std::string& filename) 
{
    const std::string dataset_prefix = "values";

    std::vector<std::vector<double>> out;

    auto file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t gid = H5Gopen(file_id, "/", H5P_DEFAULT);

    hsize_t n_objs;
    H5Gget_num_objs(gid, &n_objs);

    for (size_t i = 0; i < n_objs; i++) {
        auto dataset_name = get_dataset_name(gid, i);
        if (startswith(dataset_name, dataset_prefix)) {
            continue;
        }
        auto dim = dim_from_dataset_name(dataset_name, dataset_prefix);
        out.resize(std::max(out.size(), dim + 1));
        out[dim] = load_function(file_id, dataset_name);
    }

    return out;
};
