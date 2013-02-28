import datetime, os, time, argparse, multiprocessing, subprocess

def work(cmd):
    return subprocess.call(cmd, shell=False)

IR_PATH  = "IR"
VIS_PATH = "VIS"


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
parser.add_argument("square", type=int,
                    help="display a square of a given number")
    if not os.path.isdir(IR_PATH):
        os.mkdir(IR_PATH)
    if not os.path.isdir(VIS_PATH):
        os.mkdir(VIS_PATH)
    capture_delay = 5
    capture_num   = 5

    count = multiprocessing.cpu_count()  #get num of multiprocessors

    for i in range(capture_num):
        dt = datetime.datetime.now()
        dtstr = dt.strftime("%Y-%m-%d_%H_%M_%S")
        fn = "%simg%03d.png" % (dtstr,i)
        fpath0 = os.sep.join((IR_PATH,fn))
        fpath1 = os.sep.join((VIS_PATH,fn))
        cmd0 = "fswebcam -r 1280x1024 --png 9 -D 1 -S 3 --save %s -d /dev/video0" % fpath0 
        cmd1 = "fswebcam -r 1280x1024 --png 9 -D 1 -S 3 --save %s -d /dev/video1" % fpath1 
        cmd0 = cmd0.split()
        cmd1 = cmd1.split()
        pool = multiprocessing.Pool(processes=count)
        print pool.map(work, [cmd0,cmd1])
        time.sleep(capture_delay)
    

