#pragma once

#include "hnswlib.h"

namespace hnswlib
{
    class CloudySpace : public SpaceInterface<double>
    {
        DISTFUNC<double> fstdistfunc_;
        size_t data_size_;
        size_t dim_;

    public:
        CloudySpace(size_t dim, DISTFUNC<double> dist_func)
        {
            fstdistfunc_ = dist_func;

            dim_ = dim;
            data_size_ = dim * sizeof(double);
        }

        size_t get_data_size() override { return data_size_; }

        DISTFUNC<double> get_dist_func() override { return fstdistfunc_; }

        void* get_dist_func_param() override { return &dim_; }

        ~CloudySpace() {}
    };
};
