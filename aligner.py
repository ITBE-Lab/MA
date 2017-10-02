from LAuS import *

class Module(Module__):
    def get_input_type(self):
        return [ContainerType.nothing]

    #override
    def get_output_type(self):
        return ContainerType.nothing

    #override
    def execute(self, input):
        return Module__.execute(input)

    def promise_me(self, input):
        return Pledge(self, self.get_output_type(), input)


class SweepAllReturnBest(Module):
    """executes the given module on all elements.
    """

    def __init__(self):
        self.linesweep = LineSweep()

    #override
    def get_input_type(self):
        return [ContainerType.seedsVector, ContainerType.query, ContainerType.ref_seq]

    #override
    def get_output_type(self):
        return ContainerType.seeds

    #override
    def execute(self, input):
        strips, query, ref_seq = input
        best_strip = []
        for strip in strips:
            app = self.linesweep.execute((query, ref_seq, strip))
            best_strip.append(app)

        best = 0
        for index, strip in enumerate(best_strip):
            if strip.get_score() > best_strip[best].get_score():
                best = index

        return best_strip[best]

