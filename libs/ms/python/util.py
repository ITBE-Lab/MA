from MS import *

def promise_me(module, *args):
    arg_vector = VectorPledge()
    for arg in args:
        arg_vector.append(arg)
    if isinstance(module, Module):
        return libMS._util.ModulePledge(module, arg_vector)
    elif isinstance(module, VolatileModule):
        return libMS._util.VolatileModulePledge(module, arg_vector)
    else:
        raise Exception(
            "module must be an instance of Module or VolatileModule")

##
# @brief call this to set up a parallel section of a computational graph
# @details
# set up the parallel section in fSetUpGraph.
#
def parallel_graph(num_threads):
    for idx in range(num_threads):
        BasePledge.current_graph_thread = idx + 1
        yield idx + 1
    BasePledge.current_graph_thread = BasePledge.default_graph_thread