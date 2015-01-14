#ifndef __YYYCYYCXYZYZYYZI_SPEC_H
#define __YYYCYYCXYZYZYYZI_SPEC_H

#include <vector>
#include <string>

struct OperatorSpec 
{
    const std::string src_mesh;
    const std::string obs_mesh;
    const std::string kernel;
};

struct BlockOperatorSpec
{
    const std::vector<OperatorSpec> specs;
    const int n_rows;
    const int n_cols;
};

struct RHSSpec
{
    OperatorSpec operator_spec;
    const std::string bc_function_idx;  
};

struct BlockRHSSpec 
{
    const std::vector<RHSSpec> specs;
};

#endif
