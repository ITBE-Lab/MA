from .printColumns import print_columns
import math

class AnalyzeRuntimes:
    def reset(self):
        self.times = {}
        self.counter = 0

    def __init__(self):
        self.times = {}
        self.counter = 0

    def register(self, name, pledge, average=False, func=lambda x: x.exec_time):
        if not name in self.times:
            self.times[name] = (self.counter, average, [])
            self.counter += 1
        self.times[name][1].append( (pledge, func) )

    def analyze(self, out_file=None):
        print("runtime analysis:")
        if not out_file is None:
            out_file.write("runtime analysis:\n")
        data = []
        max_before_dot = 0
        total_runtime = 0 
        for name, (counter, average, pledges) in self.times.items():
            total_runtime += sum(func(pledge) for pledge, func in pledges) / (len(pledges) if average else 1)
        for name, (counter, average, pledges) in self.times.items():
            seconds = round(sum(func(pledge) for pledge, func in pledges) / (len(pledges) if average else 1), 3)
            percentage = (100*seconds)//total_runtime
            percentage_str = str(percentage) + "%"
            if percentage < 100:
                percentage_str = " "+percentage_str
            if percentage < 10:
                percentage_str = " "+percentage_str
            max_before_dot = max(max_before_dot, int(math.log10(max(1,seconds))) )
            counter_str = str(counter)
            if counter < 10:
                counter_str = " " + counter_str
            if counter < 100:
                counter_str = " " + counter_str
            data.append(["[" + counter_str + "] " + name, seconds, percentage_str, average])
        data = [(x, str(" "*int(max_before_dot-int(math.log10(max(1,y))))) + str(y), z ) for x,y,z in data]
        data.sort()
        data.insert(0, ["Module name", "runtime [s]", "ratio [%]", "averaged?"])
        data.append([">>>>> Total", str(round(total_runtime, 3)), ""])
        print_columns(data, out_file)
        if not out_file is None:
            out_file.write("\n")