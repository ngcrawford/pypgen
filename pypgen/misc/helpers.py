import sys
import gzip
import datetime
import numpy


def float_2_string(value, places):
    """Convert value to string truncting the number
        of decimals at the number of places given by
        the places arguement. """

    float_types = [float, numpy.float, numpy.float128, numpy.float16, numpy.float32, numpy.float64]
    if type(value) in float_types:
        value = str(round(value, places))

    else:
        value = str(value)

    return value


class Unbuffered:
    """Process STOUT so that it is printed as it is produced.

        Note: This slows the program donw, a bit, but allows
        its progress to be easily tracked. I feel the tradeoff
        is worth it.
    """

    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)

        try:
            self.stream.flush()
        except IOError:  # Ensure exit is clean if piping into head (or similar)
            sys.exit()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


def open_vcf(args):
    if args.input.endswith('.gz') == True:  # To Do: This is hacky.
        fin = gzip.open(args.input, 'rb')
    else:
        fin = open(args.input, 'rU')

    return fin


def progress_meter(starting_time, chrm, pos, bp_processed, total_bp_in_dataset):

    #sys.stdout=Unbuffered(sys.stdout) # make sure writing to std out isn't buffered

    ct = datetime.datetime.now()
    elapsed = (datetime.datetime.now() - starting_time)

    proportion_processed = bp_processed / float(total_bp_in_dataset)

    if elapsed.seconds % 30 == 0 and elapsed.seconds > 10:

        "INFO  19:29:31,683 ProgressMeter - GL343193.1:820101\t1.79e+09\t2.3 h\t4.6 s\t99.2%\t2.3 h\t62.5 s "
        status = "INFO  {} ProgressMeter - {}:{:,} {:.3} {:.2%}\n".format(ct.time(), chrm, pos, elapsed, proportion_processed)
        sys.stdout.write(status)
        return ct

    else:
        return starting_time
