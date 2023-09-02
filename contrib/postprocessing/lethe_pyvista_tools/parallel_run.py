from multiprocessing import Pool
import multiprocessing as multiprocessing
from tqdm import tqdm

multiprocessing.set_start_method('fork')


def parallel_run(self, loop, iterable, tqdm_desc = None):
    """
     Loop through data in parallel using multiprocessing Pool().imap().

     Number of processors in self.n_procs.

    :param loop -> name of inner function to call with iterator.
    :param iterable -> iterable object. Typically range(len(...)).
    :param tqdm_desc: String to appear by progress bar. If None, shows no text by progress bar.
    """
    pool = Pool(processes=self.n_procs)
    for _ in tqdm(pool.imap(loop, iterable), total = len(iterable), desc = tqdm_desc):
        pass
