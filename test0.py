import datetime
import multiprocessing
import subprocess

def work(cmd):
    return subprocess.call(cmd, shell=False)

WEBCAM_CAPTURE_FMT = "fswebcam -r 1280x1024 --png 9 -D 1 -S 3 --save %s.png -d /dev/video%d"

if __name__ == "__main__":
    capture_delay = 5
    capture_num   = 5
    for i in range(capt
    dt = datetime.datetime.now()
    dt.strftime("img_%Y-%m-%d_%H_%M_%S") 
count = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=count)
print pool.map(work, [w1,w2])

