#ifndef Sim_RandomKirigami_hpp
#define Sim_RandomKirigami_hpp

#include "Sim.hpp"
#include "Mesh.hpp"

class Sim_RandomKirigami : public Sim<Mesh>
{
    typedef Mesh tMesh;

protected:
    void runRandomKirigami();
    void runRectangularRandomKirigami();
    void runPentagonRandomKirigami();
    void runCurvedRandomKirigami();

public:
    Sim_RandomKirigami(ArgumentParser & parser):
    Sim<tMesh>(parser)
    {}

    void init() override;
    void run() override;
};

#endif /* Sim_RandomKirigami_hpp */
