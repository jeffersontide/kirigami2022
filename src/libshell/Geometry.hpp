//
//  Geometry.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Geometry_h
#define Geometry_h

#include "common.hpp"

#include <functional>

#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/upsample.h>

#include "ReadVTK.hpp"
#include "GeometryBase.hpp"


class Geometry_Dummy : public PlateGeometry
{
    Eigen::MatrixXd my_vertices;
    Eigen::MatrixXi my_face2vertices;

    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        vertices = my_vertices;
        face2vertices = my_face2vertices;
        initVertexBoundaries(vertices.rows(), vertices_bc);
    }

public:
    Geometry_Dummy(const Eigen::Ref<const Eigen::MatrixXd> my_vertices_in, const Eigen::Ref<const Eigen::MatrixXi> my_face2vertices_in)
    {
        my_vertices = my_vertices_in;
        my_face2vertices = my_face2vertices_in;
    }
};


#endif /* Geometry_h */
