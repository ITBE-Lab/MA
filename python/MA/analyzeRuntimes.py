from .printColumns import print_columns
import math

class AnalyzeRuntimes:
    def reset(self):
        self.times = {}
        self.counter = 0

    def __init__(self):
        self.times = {}
        self.counter = 0

    def register(self, name, pledge, average=False, func=lambda x: x.exec_time,
                 w_func=lambda x: x.wait_on_lock_time() if "wait_on_lock_time" in dir(x) else 0):
        if not name in self.times:
            self.times[name] = (self.counter, average, [])
            self.counter += 1
        self.times[name][2].append( (pledge, func, w_func) )
    
    def prepend_space(self, max_len, val):
        return str(" "*int(max_len-int(math.log10(max(1,val))))) + str(val)

    def analyze(self, out_file=None):
        print("runtime analysis:")
        if not out_file is None:
            out_file.write("runtime analysis:\n")
        data = []
        total_runtime = 0 
        for name, (counter, average, pledges) in self.times.items():
            total_runtime += sum(func(pledge) for pledge, func, _ in pledges) / (len(pledges) if average else 1)
        max_before_dot = int(math.log10(max(1,total_runtime)))
        for name, (counter, average, pledges) in self.times.items():
            wait_sec = round(sum(w_func(pledge) for pledge, _, w_func in pledges) / (len(pledges) if average else 1), 3)
            seconds = round(sum(func(pledge) for pledge, func, _ in pledges) / (len(pledges) if average else 1), 3)
            percentage = (100*seconds)//total_runtime
            percentage_str = str(percentage) + "%"
            if percentage < 100:
                percentage_str = " "+percentage_str
            if percentage < 10:
                percentage_str = " "+percentage_str
            max_before_dot = max(max_before_dot, int(math.log10(max(1,seconds))), int(math.log10(max(1,wait_sec))) )
            counter_str = str(counter)
            if counter < 10:
                counter_str = " " + counter_str
            if counter < 100:
                counter_str = " " + counter_str
            data.append(["[" + counter_str + "] " + name, seconds, percentage_str, average, wait_sec])
        data_2 = [(x, self.prepend_space(max_before_dot, y), z, self.prepend_space(max_before_dot, b),
                   str(a) ) for x,y,z,a,b in data]
        data_2.sort()
        data_2.insert(0, ("Module name", "runtime [s]", "ratio [%]", "waittime [s]", "averaged?"))
        data_2.append([">>>>> Total", str(round(total_runtime, 3)), "", "", ""])
        print_columns(data_2, out_file)
        if not out_file is None:
            out_file.write("\n")