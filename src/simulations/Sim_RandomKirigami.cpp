#include "Sim_RandomKirigami.hpp"
#include "Geometry_Kirigami.hpp"
#include "MaterialProperties.hpp"
#include "CombinedOperator_Parametric.hpp"
#include "EnergyOperatorList.hpp"


void Sim_RandomKirigami::init()
{}


void Sim_RandomKirigami::run()
{
    runRandomKirigami();
}


void Sim_RandomKirigami::runRandomKirigami()
{
    tag = "kirigami";

    // All lengths in units of cm, angles in radians
    const Real cmPerInch = 2.54;

    // Geometric parameters for mesh
    const Real radius = parser.parse<Real>("-radius", 3.5 * cmPerInch);

    // Geometric parameters for simulation
    // Extension angle is always 0 (sheet gets stretched in the x direction)
    const Real clampedWidth = parser.parse<Real>("-clampedWidth", 2.5);
    const Real clampedDepth = parser.parse<Real>("-clampedDepth", 0.5 * cmPerInch);

    // Resolution parameters for triangulation
    const Real maxRes = parser.parse<Real>("-maxRes", 0.5);
    const Real minRes = parser.parse<Real>("-minRes", 0.05);
    const Real decayRes = parser.parse<Real>("-decayRes", 1.0);

    // Read cut information from file
    // File must be in the following format:
    // cut1_X cut1_Y cut1_length cut1_width cut1_angle
    // cut2_X cut2_Y cut2_length cut2_width cut2_angle
    // ...
    const std::string filename = parser.parse<std::string>("-filename", "cuts.txt");

    // Initialize geometry/mesh
    RandomKirigamiPlate geometry(radius, filename, maxRes, minRes, decayRes, {clampedWidth, clampedDepth}, false);
    mesh.init(geometry);

    // Set boundary conditions for clamped areas of the mesh for pulling, based on rest positions
    geometry.fixClampedRegions(mesh);

    // Define material parameters
    const Real YoungModulus = parser.parse<Real>("-E", 1.0);
    const Real PoissonRatio = parser.parse<Real>("-nu", 0.4);
    const Real thickness = parser.parse<Real>("-h", 0.01);
    MaterialProperties_Iso_Constant matprop(YoungModulus, PoissonRatio, thickness);
    CombinedOperator_Parametric<tMesh> engOp(matprop);

    // Dump initial mesh
    dump(tag + "_init");

    // Define incremental strains (e.g. 0.1 would be a 10% strain)
    Eigen::VectorXd strains = getDataFromFile("strains.txt", 1);

    for (int step = 1; step < strains.size(); step++) {
        // Set strains
        Real prev_factor = strains(step - 1);
        Real factor = strains(step);

        auto cvertices = mesh.getCurrentConfiguration().getVertices();
        for (int i = 0; i < cvertices.rows(); i++) {
            cvertices(i, 0) *= (1.0 + factor) / (1.0 + prev_factor);
        }

        // Add noise to break symmetries (in z direction only)
        addNoiseToVertices_c<2>(0.01 * thickness);

        // Minimize mesh energy
        Real eps = 1e-2;
        minimizeEnergy(engOp, eps);

        // Dump results
        dump(tag + "_strained" + std::to_string(strains(step)));
    }
}
