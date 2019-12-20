from time import time, time_ns, sleep

for i in range(10):
    t_secs1 = time()
    sleep(10**(- i))
    t_secs2 = time()
    print(t_secs2 - t_secs1)