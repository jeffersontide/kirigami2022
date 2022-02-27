#include "common.hpp"
#include "ArgumentParser.hpp"

#include "Sim.hpp"
#include "Sim_RandomKirigami.hpp"

#ifdef USETBB
#include "tbb/task_scheduler_init.h"
#endif

int main(int argc,  const char ** argv)
{
#ifdef USETBB
    const int num_threads = tbb::task_scheduler_init::default_num_threads();
    tbb::task_scheduler_init init(num_threads);
    std::cout << "Starting TBB with " << num_threads << " threads " << std::endl;
#endif

    BaseSim * sim = nullptr;

    ArgumentParser parser(argc,argv);

    parser.initialize();

    sim = new Sim_RandomKirigami(parser);

    sim->init();
    sim->run();

    delete sim;

    parser.finalize();

    return 0;
}
