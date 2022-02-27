//
//  Geometry_Kirigami.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Geometry_Kirigami_h
#define Geometry_Kirigami_h

#include "common.hpp"

#include <functional>

#include "ReadVTK.hpp"
#include "GeometryBase.hpp"


class RandomKirigamiPlate : public PlateGeometry
{
protected:
    const Real radius;
    std::string filename;
    Eigen::MatrixXd cuts;
    const Real maxRes;
    const Real minRes;
    const Real decayRes;
    const std::array<Real, 2> clampParameters; // Assumes two clamped regions
    const bool fixedBoundary;
    const bool roundedEdges;

    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const int nPointsAlongPerimeter = (int)(2.0 * M_PI * radius / maxRes) + 1;
        const Real dTheta = 2.0 * M_PI / nPointsAlongPerimeter;

        // Matrix to store cut information: cutX, cutY, cutLength, cutWidth, cutAngle (in cm/radians)
        Eigen::MatrixXd cuts = getDataFromFile(filename, 5);
        int nCuts = cuts.rows();

        // Semi-circular caps will be added at each end of the cut, so the widest
        // part of the cut will actually be longer than the specified length.
        // Adjust number of vertices that lie on each semicircular cap:
        // const int nPointsPerCap = roundedEdges ? (int)(cuts(0, 3) / minRes) + 1 : 0;
        const int nPointsPerCap = roundedEdges ? 10 : 0;
        const int nPointsPerCut = 4 + 2 * nPointsPerCap;
        const Real dThetaCap = M_PI / (nPointsPerCap + 1);

        // Outer circle, cut boundaries, and clamped region
        Eigen::MatrixXd polygon_vertices(nPointsAlongPerimeter + nPointsPerCut * nCuts + 2 * 4, 2);
        Eigen::MatrixXi polygon_edges(nPointsAlongPerimeter + nPointsPerCut * nCuts + 2 * 3, 2);

        // Outer circle
        polygon_edges(0, 0) = 0;
        polygon_edges(nPointsAlongPerimeter - 1, 1) = 0;
        for (int i = 0; i < nPointsAlongPerimeter; i++) {
            const Real theta = dTheta * i;
            polygon_vertices(i, 0) = std::cos(theta) * radius;
            polygon_vertices(i, 1) = std::sin(theta) * radius;
            if (i < nPointsAlongPerimeter - 1) {
                polygon_edges(i, 1) = i + 1;
                polygon_edges(i + 1, 0) = i + 1;
            }
        }

        // Cut boundaries
        for (int i = 0; i < nCuts; i++) {
            int pos = nPointsAlongPerimeter + nPointsPerCut * i;
            Real cutX = cuts(i, 0);
            Real cutY = cuts(i, 1);
            Real cutLength = cuts(i, 2);
            Real cutWidth = cuts(i, 3);
            Real cutAngle = cuts(i, 4);

            // NE --> NW
            polygon_vertices(pos, 0) = cutX + cutLength * std::cos(cutAngle) / 2.0 + cutWidth * std::cos(cutAngle + M_PI / 2.0) / 2.0;
            polygon_vertices(pos, 1) = cutY + cutLength * std::sin(cutAngle) / 2.0 + cutWidth * std::sin(cutAngle + M_PI / 2.0) / 2.0;
            polygon_vertices(pos + 1, 0) = cutX - cutLength * std::cos(cutAngle) / 2.0 + cutWidth * std::cos(cutAngle + M_PI / 2.0) / 2.0;
            polygon_vertices(pos + 1, 1) = cutY - cutLength * std::sin(cutAngle) / 2.0 + cutWidth * std::sin(cutAngle + M_PI / 2.0) / 2.0;

            // W cap
            for (int j = 0; j < nPointsPerCap; j++) {
                polygon_vertices(pos + 2 + j, 0) = cutX - cutLength * std::cos(cutAngle) / 2.0 + cutWidth * std::cos(cutAngle + M_PI / 2.0 + dThetaCap * (j + 1)) / 2.0;
                polygon_vertices(pos + 2 + j, 1) = cutY - cutLength * std::sin(cutAngle) / 2.0 + cutWidth * std::sin(cutAngle + M_PI / 2.0 + dThetaCap * (j + 1)) / 2.0;
            }

            // SW --> SE
            polygon_vertices(pos + 2 + nPointsPerCap, 0) = cutX - cutLength * std::cos(cutAngle) / 2.0 - cutWidth * std::cos(cutAngle + M_PI / 2.0) / 2.0;
            polygon_vertices(pos + 2 + nPointsPerCap, 1) = cutY - cutLength * std::sin(cutAngle) / 2.0 - cutWidth * std::sin(cutAngle + M_PI / 2.0) / 2.0;
            polygon_vertices(pos + 3 + nPointsPerCap, 0) = cutX + cutLength * std::cos(cutAngle) / 2.0 - cutWidth * std::cos(cutAngle + M_PI / 2.0) / 2.0;
            polygon_vertices(pos + 3 + nPointsPerCap, 1) = cutY + cutLength * std::sin(cutAngle) / 2.0 - cutWidth * std::sin(cutAngle + M_PI / 2.0) / 2.0;

            // E cap
            for (int j = 0; j < nPointsPerCap; j++) {
                polygon_vertices(pos + 4 + nPointsPerCap + j, 0) = cutX + cutLength * std::cos(cutAngle) / 2.0 + cutWidth * std::cos(cutAngle - M_PI / 2.0 + dThetaCap * (j + 1)) / 2.0;
                polygon_vertices(pos + 4 + nPointsPerCap + j, 1) = cutY + cutLength * std::sin(cutAngle) / 2.0 + cutWidth * std::sin(cutAngle - M_PI / 2.0 + dThetaCap * (j + 1)) / 2.0;
            }

            // Connect all vertices sequentially
            for (int j = 0; j < nPointsPerCut; j++) {
                polygon_edges(pos + j, 0) = pos + j;
                polygon_edges(pos + j, 1) = (j == nPointsPerCut - 1) ? pos : pos + j + 1;
            }
        }

        // Clamped regions' vertices and edges (shrink them very slightly by
        // clampEps so that they will be caught by the bounds-checking in the
        // main script)
        {
            Real clampEps = 0.001;
            Real clampedWidth = clampParameters[0];
            Real clampedDepth = clampParameters[1];

            int pos = nPointsAlongPerimeter + nPointsPerCut * nCuts;

            Real tA = std::sqrt(radius * radius - clampedWidth * clampedWidth / 4.0); // x_max
            Real tB = clampedWidth / 2.0 - clampEps; // y
            Real tC = radius - clampedDepth + clampEps; // x_min

            // Right (+x) clamped region
            polygon_vertices(pos + 0, 0) = tA;
            polygon_vertices(pos + 0, 1) = tB;
            polygon_vertices(pos + 1, 0) = tC;
            polygon_vertices(pos + 1, 1) = tB;
            polygon_vertices(pos + 2, 0) = tC;
            polygon_vertices(pos + 2, 1) = -tB;
            polygon_vertices(pos + 3, 0) = tA;
            polygon_vertices(pos + 3, 1) = -tB;
            polygon_edges(pos + 0, 0) = pos + 0;
            polygon_edges(pos + 0, 1) = pos + 1;
            polygon_edges(pos + 1, 0) = pos + 1;
            polygon_edges(pos + 1, 1) = pos + 2;
            polygon_edges(pos + 2, 0) = pos + 2;
            polygon_edges(pos + 2, 1) = pos + 3;

            // Left (-x) clamped region
            polygon_vertices(pos + 4, 0) = -tA;
            polygon_vertices(pos + 4, 1) = tB;
            polygon_vertices(pos + 5, 0) = -tC;
            polygon_vertices(pos + 5, 1) = tB;
            polygon_vertices(pos + 6, 0) = -tC;
            polygon_vertices(pos + 6, 1) = -tB;
            polygon_vertices(pos + 7, 0) = -tA;
            polygon_vertices(pos + 7, 1) = -tB;
            polygon_edges(pos + 3, 0) = pos + 4;
            polygon_edges(pos + 3, 1) = pos + 5;
            polygon_edges(pos + 4, 0) = pos + 5;
            polygon_edges(pos + 4, 1) = pos + 6;
            polygon_edges(pos + 5, 0) = pos + 6;
            polygon_edges(pos + 5, 1) = pos + 7;
        }

        // Triangulation flags
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*maxRes*maxRes)+(verbose ? "" : "Q");

        // Set triangulation resolution function
        Real args[nCuts][4 + nPointsPerCut][3 + 4];
        for (int i = 0; i < nCuts; i++) {
            for (int j = 0; j < nPointsPerCut; j++) {
                int pos = nPointsAlongPerimeter + nPointsPerCut * i + j;
                int next_pos = (j < nPointsPerCut - 1) ? pos + 1 : pos - (nPointsPerCut - 1);
                args[i][j][0] = maxRes;
                args[i][j][1] = minRes;
                args[i][j][2] = decayRes;
                args[i][j][3] = polygon_vertices(pos, 0);       // x1
                args[i][j][4] = polygon_vertices(pos, 1);       // y1
                args[i][j][5] = polygon_vertices(next_pos, 0);  // x2
                args[i][j][6] = polygon_vertices(next_pos, 1);  // y2

                ::setTriangleResolutionFunction([](Real x, Real y, Real* arg) {
                    // https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
                    Real a = x - arg[3];
                    Real b = y - arg[4];
                    Real c = arg[5] - arg[3];
                    Real d = arg[6] - arg[4];

                    Real dot = a * c + b * d;
                    Real len_sq = c * c + d * d;
                    Real temp = -1.0;
                    if (len_sq != 0.0) { temp = dot / len_sq; }

                    // Find the location on the segment that the point is closest to.
                    Real xx;
                    Real yy;

                    if (temp < 0.0) {
                        xx = arg[3];
                        yy = arg[4];
                    }
                    else if (temp > 1.0) {
                        xx = arg[5];
                        yy = arg[6];
                    }
                    else {
                        xx = arg[3] + temp * c;
                        yy = arg[4] + temp * d;
                    }

                    Real dx = x - xx;
                    Real dy = y - yy;
                    Real dist = std::sqrt(dx * dx + dy * dy);

                    return arg[0] - (arg[0] - arg[1]) * exp(-dist / arg[2]);
                }, args[i][j]);
            }
        }

        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);

        if (fixedBoundary) {
            const int nV = vertices.rows();
            for(int i=0;i<nV;++i)
            {
                const Real radSq = vertices(i,0)*vertices(i,0) + vertices(i,1)*vertices(i,1);
                const Real diff = std::abs(radSq - radius*radius);
                if(diff < 1e-6*radius)
                {
                    vertices_bc(i,0) = true;
                    vertices_bc(i,1) = true;
                    vertices_bc(i,2) = true;
                }
            }
        }

        // Remove cuts from mesh
        for (int i = 0; i < nCuts; i++) {
            removePolygonFromMesh(vertices, face2vertices, vertices_bc,
                                  polygon_vertices.block(nPointsAlongPerimeter + nPointsPerCut * i, 0,
                                                         nPointsPerCut, 2));
        }

        // Remove any disjoint parts of the mesh
        removeUnconnectedComponents(vertices, face2vertices, vertices_bc);
    }

public:
    RandomKirigamiPlate(const Real radius, std::string filename,
                        const Real maxRes, const Real minRes, const Real decayRes,
                        const std::array<Real, 2> clampParameters,
                        const bool fixedBoundary, const bool roundedEdges = true):
    PlateGeometry(),
    radius(radius),
    filename(filename),
    maxRes(maxRes),
    minRes(minRes),
    decayRes(decayRes),
    clampParameters(clampParameters),
    fixedBoundary(fixedBoundary),
    roundedEdges(roundedEdges)
    {};

    template <typename tMesh>
    void fixClampedRegions(tMesh & mesh) {
        // Set boundary conditions for clamped areas of the mesh for pulling, based on rest positions
        auto rvertices = mesh.getRestConfiguration().getVertices();
        auto boundaryConditions = mesh.getBoundaryConditions().getVertexBoundaryConditions();
        for (int i = 0; i < rvertices.rows(); i++) {
            // Check if the vertex lies outside the range of clamping
            if (-radius + clampParameters[1] < rvertices(i, 0) && rvertices(i, 0) < radius - clampParameters[1]) { continue; }
            if (rvertices(i, 1) < -clampParameters[0] / 2.0 || clampParameters[0] / 2.0 < rvertices(i, 1)) { continue; }

            // Otherwise: specify that the vertex should remain fixed during energy minimization
            boundaryConditions(i, 0) = true;
            boundaryConditions(i, 1) = true;
            boundaryConditions(i, 2) = true;
        }
    }
};


#endif /* Geometry_Kirigami_h */
