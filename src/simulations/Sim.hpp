
#ifndef Sim_hpp
#define Sim_hpp

#include <iostream>

#include "common.hpp"
#include "ArgumentParser.hpp"

#ifdef USEVTK
#include "ReadVTK.hpp"
#include "WriteVTK.hpp"
#endif

#include "WriteSTL.hpp"
#include "ComputeCurvatures.hpp"

#include "EnergyOperator.hpp"

#ifdef USELIBLBFGS
#include "LBFGS_Wrapper.hpp"
#endif

#ifdef USEHLBFGS
#include "HLBFGS_Wrapper.hpp"
#endif

#include <igl/writeOBJ.h>


/*! \class Sim
 * \brief Base class for simulations.
 *
 * This class is called from main, and performs the main simulation. Every class that derives from here can implement a simulation case.
 */
class BaseSim
{
protected:
    ArgumentParser & parser;

public:
    BaseSim(ArgumentParser & parser_in):
    parser(parser_in)
    {
        parser.save_defaults(); // make sure all default values are printed as well from now on
        parser.save_options();
    }

    virtual void init() = 0;
    virtual void run() = 0;
    virtual int optimize()
    {
        std::cout << "Optimize not (yet) implemented for this class " << std::endl;
        return -1;
    };

    virtual ~BaseSim()
    {}
};


template<typename tMesh>
class Sim : public BaseSim
{
public:
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;

protected:
    std::string tag;
    tMesh mesh;
    std::vector< std::pair<std::string, Real> > currentEnergies;
    const bool verbose = parser.parse<bool>("-verbose", true);

    WriteVTK writer;

    virtual void writeSTL(const std::string filename)
    {
        WriteSTL::write(mesh.getTopology(), mesh.getCurrentConfiguration(), filename);
    }

    virtual std::string getArguments()
    {
        return parser.options_string();
    }

    virtual void addMesh()
    {
        writer.updateMesh(mesh.getCurrentConfiguration().getVertices(), mesh.getTopology().getFace2Vertices());
    }

    virtual void dump(const std::string filename)
    {
        if (!writer.initialized()) {
            addMesh();
        }

        // add arguments
        writer.addStringAttribute(getArguments(), "arguments");

        // add energies
        for (const auto & eng : currentEnergies) {
            writer.addScalarAttribute(eng.second, eng.first);
        }

        // add displacements
        const Eigen::MatrixXd displacements = mesh.getCurrentConfiguration().getVertices() - mesh.getRestConfiguration().getVertices();
        writer.addVectorFieldToVertices(displacements, "displacements");

        // add basic curvatures
        const int nFaces = mesh.getNumberOfFaces();
        Eigen::VectorXd gauss(nFaces);
        Eigen::VectorXd mean(nFaces);
        ComputeCurvatures<tMesh> computeCurvatures;
        computeCurvatures.compute(mesh, gauss, mean);
        writer.addScalarFieldToFaces(mean, "mean");
        writer.addScalarFieldToFaces(gauss, "gauss");

        // write
        writer.write(filename);

        writer.reset();
    }

    template<int component>
    void addNoiseToVertices_c(const Real ampl)
    {
        std::mt19937 gen;
        gen.seed(42);
        std::uniform_real_distribution<Real> distV(-ampl, ampl);
        auto perturb_v = [&](Eigen::Vector3d in)
        {
            Eigen::Vector3d retval;
            for(int d=0;d<3;++d)
                retval(d) = in(d) + (d == component ? distV(gen) : 0.0);
            return retval;
        };
        mesh.changeVertices(perturb_v);
    }

    template<typename tMeshOperator, bool verbose = true>
    int minimizeEnergy(const tMeshOperator & op, Real & eps, const Real epsMin=std::numeric_limits<Real>::epsilon(), const bool stepWise = false)
    {
#ifdef USELIBLBFGS
        // use the LBFGS_Energy class to directly minimize the energy on the mesh with these operators
        // LBFGS does not use the hessian
        LBFGS::LBFGS_Energy<Real, tMesh, tMeshOperator, true> lbfgs_energy(mesh, op);
        int retval = 0;
        while(retval == 0 && eps > epsMin)
        {
            eps *= 0.1;
            retval = lbfgs_energy.minimize(eps);
        }
#else
#ifdef USEHLBFGS

        HLBFGS_Methods::HLBFGS_Energy<tMesh, tMeshOperator, true> hlbfgs_wrapper(mesh, op);
        int retval = 0;
        if(stepWise)
        {
            while(retval == 0 && eps > epsMin)
            {
                eps *= 0.1;
                retval = hlbfgs_wrapper.minimize(tag+"_diagnostics.dat", eps);
            }
        }
        else
        {
            retval = hlbfgs_wrapper.minimize(tag+"_diagnostics.dat", epsMin);
            eps = hlbfgs_wrapper.get_lastnorm();
        }
#else
        std::cout << "should use liblbfgs or hlbfgs\n";
#endif
#endif

        // store energies
        {
            std::vector<std::pair<std::string, Real>> energies;
            op.addEnergy(energies);
            currentEnergies = energies;
            FILE * f = fopen((tag+"_energies.dat").c_str(), "a");
            for(const auto & eng : energies)
            {
                fprintf(f, "%s \t\t %10.10e\n", eng.first.c_str(), eng.second);
            }
            fclose(f);
        }
        return retval;
    }

public:
    Sim(ArgumentParser & parser_in):
    BaseSim(parser_in),
    tag("Sim")
    {
        writer = WriteVTK();
    }

    virtual ~Sim()
    {}
};

#endif /* Sim_hpp */
