import sys

import datetime

class Unbuffered:
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

def progress_meter(starting_time, chrm, pos, bp_processed, total_bp_in_dataset):
    
    #sys.stdout=Unbuffered(sys.stdout) # make sure writing to std out isn't buffered

    ct = datetime.datetime.now()
    elapsed = (datetime.datetime.now() - starting_time)

    proportion_processed = bp_processed/float(total_bp_in_dataset)

    if elapsed.seconds % 30 == 0 and elapsed.seconds > 10:

        "INFO  19:29:31,683 ProgressMeter - GL343193.1:820101\t1.79e+09\t2.3 h\t4.6 s\t99.2%\t2.3 h\t62.5 s "
        status = "INFO  {} ProgressMeter - {}:{:,} {:.3} {:.2%}\n".format(ct.time(), chrm, pos, elapsed, proportion_processed)
        sys.stdout.write(status)
        return ct

    else:
        return starting_time