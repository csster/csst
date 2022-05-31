import random as rd
import subprocess
import sys
from multiprocessing import Process, Manager, Pool


def execute_thread(nobj, nthread, fun_thread, args):
    threads = []
    if nthread > nobj:
        nthread = nobj

    step = nobj // nthread
    nrest = nobj % nthread
    index_threads = [0]
    for i in range(nthread):
        if i < nrest:
            inext = index_threads[i] + step + 1
        else:
            inext = index_threads[i] + step
        index_threads.append(inext)
    manager = Manager()
    dict_manage = manager.dict()
    args = list(args)
    for i in range(nthread):
        index_slc = slice(index_threads[i], index_threads[i + 1])
        args.append(dict_manage)
        args.append(index_slc)
        args.append(i)
        t = Process(target=fun_thread, args=args)
        args.pop()
        args.pop()
        args.pop()
        threads.append(t)
    for i in range(nthread):
        threads[i].start()
    for i in range(nthread):
        threads[i].join()
    return dict_manage


def pool_nthread(fun_thread, nthread, args, rand_index=False, rand_seed=None, **kwd):
    """
    create multiprocessing pool and run
    """
    if nthread <= 0: nthread = multiprocessing.cpu_count()
    nargs = len(args)
    if len(args) <= 0: raise ValueError("args should have at least one elements")

    nrow = len(args[0])
    rangeii = range(nrow)
    if rand_index:
        rangeii = rd.sample(rangeii, nrow)

    pool = Pool(nthread)
    res = {}
    for ii in rangeii:
        argv = [iargv[ii] for iargv in args]
        res[ii] = pool.apply_async(fun_thread, argv, kwd)
    return res


def cmd_exists(cmd):
    """
    check wether a shell command exists
    """
    return subprocess.call("type " + cmd, shell=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0


def progressbar(count, total, bar_length=10, barchar='#', prefix='[', suffix=']', numshow=False, onlypercent=False,
                onlynum=False):
    percent = float(count) / total
    if len(barchar) > 1:
        barchar = '#'

    bar = barchar * int(round(percent * bar_length))
    spaces = ' ' * (bar_length - len(bar))
    percent = percent * 100
    perc = "%3.1f%%" % percent
    progbar = prefix + bar + spaces + suffix + perc
    if numshow:
        numperc = "%d/%d" % (count, total)
        progbar = prefix + bar + suffix + numperc
    if onlynum:
        numperc = "%d/%d" % (count, total)
        progbar = prefix + numperc + suffix
    if onlypercent:
        progbar = prefix + perc + suffix

    progbar = "\r" + progbar
    if count == total:
        progbar = progbar + '\n'
    sys.stdout.write(progbar)
    sys.stdout.flush()
