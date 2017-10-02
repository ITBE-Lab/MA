from LAuS import *

class Module(__Module):
    def get_input_type(self):
        return [ContainerType.nothing]

    #override
    def get_output_type(self):
        return [ContainerType.nothing]

    #override
    def execute(self, input):
        return None

    def save_execute(self):
        print("python save execute disabled!")
        return self.execute()

    def promise_me(self, input):
        print("python save execute disabled!")
        return Pledge(self, self.get_output_type(), input)


class SweepAllReturnBest(Module):
    """executes the given module on all elements.
    """

    def __init__(self):
        self.liesweep = LineSweep()

    #override
    def get_input_type(self):
        return [ContainerType.seedsVector, ContainerType.query, ContainerType.ref_seq]

    #override
    def get_output_type(self):
        return [ContainerType.seeds]

    #override
    def execute(self, input):
        strips, query, ref_seq = input
        best_strip = []
        for strip in strips:
            best_strip.append(liesweep.execute((query, ref_seq, strip)))

        best = 0
        for index, strip in enumerate(best_strip):
            if strip.get_score() > best_strip[best].get_score():
                best = index

        return best_strip[best]

