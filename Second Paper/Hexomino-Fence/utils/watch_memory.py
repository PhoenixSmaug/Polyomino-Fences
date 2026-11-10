import os
import sys
import time

LOGGING = sys.platform.startswith("linux")
PAGE = os.sysconf("SC_PAGE_SIZE")
MiB = 2**20

def rss_bytes():
    """Current resident set size."""
    with open("/proc/self/statm") as f:
        return int(f.read().split()[1]) * PAGE


def peak_bytes():
    """Kernel-tracked peak RSS (VmHWM). Catches spikes between samples."""
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024

def memory_checkpoint(only_time=False):
    tim = "%.3fh" % (time.process_time() / 3600) # time in hours
    if LOGGING and not only_time:
        return "%s rss=%.1f MiB peak=%.1f MiB" % (tim, rss_bytes() / MiB, peak_bytes() / MiB)
    return tim
try:
    memory_checkpoint()
except RuntimeError:
    LOGGING = False




