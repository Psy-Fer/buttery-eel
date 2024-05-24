import os
import pyslow5
from itertools import chain
import time

import cProfile, pstats, io

def _get_slow5_batch(args, slow5_obj, reads, size=4096, slow5_filename=None, header_array=None):
    """
    re-batchify slow5 output
    """
    batch = []
    for read in reads:
        # if args.seq_sum:
        # get header once for each read group
        read_group = read["read_group"]
        # if read_group not in header_array:
        #     header_array[read_group] = slow5_obj.get_all_headers(read_group=read_group)
        # get aux data for ead read
        

        aux_data = {"channel_number": read["channel_number"], 
                            "start_mux": read["start_mux"],
                            "start_time": read["start_time"],
                            "read_number": read["read_number"],
                            "end_reason": read["end_reason"],
                            "median_before": read["median_before"],
                            "end_reason_labels": slow5_obj.get_aux_enum_labels('end_reason')
                            }
        read["aux_data"] = aux_data
        read["header_array"] = header_array[read_group]
        read["slow5_filename"] = slow5_filename
        
        batch.append(read)
        if len(batch) >= size:
            yield batch
            batch = []
    if len(batch) > 0:
        yield batch


def read_worker(args, iq):
    '''
    single threaded worker to read slow5 (with multithreading)
    '''
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    header_array = {}
    # is dir, so reading recursivley
    if os.path.isdir(args.input):
        # this adds a limit to how many reads it will load into memory so we
        # don't blow the ram up
        max_limit = int(args.max_read_queue_size / args.slow5_batchsize)
        for dirpath, _, files in os.walk(args.input):
            for sfile in files:
                if sfile.endswith(('.blow5', '.slow5')):
                    s5 = pyslow5.Open(os.path.join(dirpath, sfile), 'r')
                    reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')
                    # if args.seq_sum:
                    num_read_groups = s5.get_num_read_groups()
                    for read_group in range(num_read_groups):
                        header_array[read_group] = s5.get_all_headers(read_group=read_group)
                    batches = _get_slow5_batch(args, s5, reads, size=args.slow5_batchsize, slow5_filename=sfile, header_array=header_array)
                    # put batches of reads onto the queue
                    for batch in chain(batches):
                        # print(iq.qsize())
                        if iq.qsize() < max_limit:
                            iq.put(batch)
                        else:
                            while iq.qsize() >= max_limit:
                                time.sleep(0.01)
                            iq.put(batch)
        
    else:
        s5 = pyslow5.Open(args.input, 'r')
        filename_slow5 = args.input.split("/")[-1]
        reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')
        # if args.seq_sum:
        num_read_groups = s5.get_num_read_groups()
        for read_group in range(num_read_groups):
            header_array[read_group] = s5.get_all_headers(read_group=read_group)
        batches = _get_slow5_batch(args, s5, reads, size=args.slow5_batchsize, slow5_filename=filename_slow5, header_array=header_array)
        # this adds a limit to how many reads it will load into memory so we
        # don't blow the ram up
        max_limit = int(args.max_read_queue_size / args.slow5_batchsize)
        # put batches of reads onto the queue
        for batch in chain(batches):
            # print(iq.qsize())
            if iq.qsize() < max_limit:
                iq.put(batch)
            else:
                while iq.qsize() >= max_limit:
                    time.sleep(0.01)
                iq.put(batch)
    for _ in range(args.procs):
        iq.put(None)
    
    # if profiling, dump info into log files in current dir
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("read_worker.log", 'w') as f:
            print(s.getvalue(), file=f)

