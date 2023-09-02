from multiprocessing import Pool
import multiprocessing as multiprocessing
from tqdm import tqdm

multiprocessing.set_start_method('fork')

def parallel_run(self, loop, iterable, tqdm_desc = None):
    pool = Pool(processes=self.n_procs)
    for _ in tqdm(pool.imap(loop, iterable), total = len(iterable), desc = tqdm_desc):
        pass
