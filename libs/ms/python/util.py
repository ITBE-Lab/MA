from MS import *

def promise_me(module, *args):
    arg_vector = VectorPledge()
    for arg in args:
        arg_vector.append(arg)
    if isinstance(module, Module):
        return ModulePledge(module, arg_vector)
    elif isinstance(module, VolatileModule):
        return VolatileModulePledge(module, arg_vector)
    else:
        raise Exception(
            "module must be an instance of Module or VolatileModule")