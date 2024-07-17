import os, sys
import pyslow5
from itertools import chain
import time

import cProfile, pstats, io

def get_data_by_channel(args, dq):
    """
    Go through the the slow5 and group readIDs by channel and read_number
    Then duplex can be done by channel grouping where 1 client handles 1 channel
    duplex queue (dq) contains 1 full channel of data per entry
    """
    # {channel: sorted[[readID, read_number]], ...}
    duplex = {}
    if os.path.isdir(args.input):
        print()
        print("Please merge your blow5 files into a single blow5 with slow5tools merge")
        print("Duplex calling does not currently support directory based reading")
        print()
        sys.exit(1)
        # for dirpath, _, files in os.walk(args.input):
        #     for sfile in files:
        #         if sfile.endswith(('.blow5', '.slow5')):
        #             s5 = pyslow5.Open(os.path.join(dirpath, sfile), 'r')
        #             reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')

        #             # get all the channel data and stick it into duplex dic
        #             for read in reads:
        #                 readID = read['read_id']
        #                 channel = int(read['channel_number'])
        #                 # mux = read['start_mux']
        #                 read_num = int(read['read_number'])
        #                 if channel in duplex.keys():
        #                     duplex[channel].append([readID, read_num])
        #                 else:
        #                     duplex[channel] = [[readID, read_num]]

        # now sort based on read_num so it's in time order per channel
        # sorted keys so queue is first in first out
        # for ch in sorted(duplex.keys()):
        #     duplex[ch].sort(key = lambda x: x[1])
        #     print(duplex[ch])
        #     # send channel to the queue
        #     dq.put([ch, duplex[ch]])
        # # to break the reader worker
        # dq.put(None)
    else:
        s5 = pyslow5.Open(args.input, 'r')
        reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')

        # get all the channel data and stick it into duplex dic
        for read in reads:
            readID = read['read_id']
            channel = int(read['channel_number'])
            # mux = read['start_mux']
            read_num = int(read['read_number'])
            if channel in duplex.keys():
                duplex[channel].append([readID, read_num])
            else:
                duplex[channel] = [[readID, read_num]]
        print("Number of channels:", len(list(duplex.keys())))
        # now sort based on read_num so it's in time order per channel
        for ch in sorted(duplex.keys()):
            duplex[ch].sort(key = lambda x: x[1])
            # print(duplex[ch])
            # send channel to the queue
            dq.put([ch, duplex[ch]])
        # to break the reader worker
        dq.put(None)


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

def duplex_read_worker(args, dq, pre_dq):
    '''
    single threaded worker to read slow5 (with multithreading)
    organises data by channel for duplex calling
    
    2 methods, both use the same workers.
    1. memory method 
    - Read file once.
    - split and sort by channel
    - group channels into args.proc lists of lists
    - feed each one to it's own queue
    - worker makes it's way through it's queue till finished

    2. file method
    - read file and write each channel to it's own file
    - keep list of files and split them into args.procs lists
    - spawn args.procs readers from this read_worker to read each file in its list
    - each reader feeds its own queue qith an attached worker
    - worker keeps consuming from queue as reader makes its way through files/channels
    
    Method 2 is missing a read_numnber sort step so will take a while
    Perhaps read once and split and sort by chanel, then random access the file to write the individual
    files in the proper order.
    '''
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    dq_names = dq.keys()
    free_names = [i for i in dq_names]
    taken_names = []

    # read the file once, and get all the channels into the pre_dq queue
    get_data_by_channel(args, pre_dq)

    s5 = pyslow5.Open(args.input, 'r')
    filename_slow5 = args.input.split("/")[-1]
    header_array = {}
    num_read_groups = s5.get_num_read_groups()
    for read_group in range(num_read_groups):
        header_array[read_group] = s5.get_all_headers(read_group=read_group)
    # reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')

    readers = {}
    # break call
    ending = False
    # pull from the pre_dq
    while True:
        if len(free_names) == 0 or ending:
            if ending and len(taken_names) == 0:
                break
            # push a batch for each generator till queues are full
            for qname in readers.keys():
                reads = readers[qname]
                # TODO: make this related to max queue size
                if dq[qname].qsize() < 5:
                    batch = next(reads, None)
                    if batch is None:
                        free_names.append(qname)
                        taken_names.remove(qname)
                    else:
                        dq[qname].put(batch)
                else:
                    continue
            # if no more channels and all queues done, break
            # remove free items from readers so we don't double up in ending state
            for qn in free_names:
                if len(readers.keys()) > 0:
                    if qn in readers.keys():
                        readers.pop(qn)
        elif not ending:
            # populate the read generators with matched queues
            ch = pre_dq.get()
            if ch is None:
                ending = True
                continue
            channel = ch[0]
            print("processing channel: {}".format(channel))
            data = ch[1]
            read_list = [i for i, _ in data]
            q = free_names[0]
            free_names.pop(0)
            taken_names.append(q)
            reads = s5.get_read_list_multi(read_list, threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')
            batches = _get_slow5_batch(args, s5, reads, size=args.slow5_batchsize, slow5_filename=filename_slow5, header_array=header_array)
            readers[q] = batches
        else:
            print("Some end state not known!")
    
    for qname in dq_names:
        dq[qname].put(None)

    # if profiling, dump info into log files in current dir
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("read_worker.log", 'w') as f:
            print(s.getvalue(), file=f)


def duplex_read_worker_single(args, dq, pre_dq):
    '''
    Single proc method
    '''
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    # read the file once, and get all the channels into the pre_dq queue
    get_data_by_channel(args, pre_dq)

    s5 = pyslow5.Open(args.input, 'r')
    filename_slow5 = args.input.split("/")[-1]
    header_array = {}
    num_read_groups = s5.get_num_read_groups()
    for read_group in range(num_read_groups):
        header_array[read_group] = s5.get_all_headers(read_group=read_group)
    # reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')

    # readers = {}
    # break call
    # ending = False
    # pull from the pre_dq
    while True:
        ch = pre_dq.get()
        if ch is None:
            break
        channel = ch[0]
        print("[READER] - processing channel: {}".format(channel))
        data = ch[1]
        read_list = [i for i, _ in data]
        reads = s5.get_read_list_multi(read_list, threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')
        batches = _get_slow5_batch(args, s5, reads, size=args.slow5_batchsize, slow5_filename=filename_slow5, header_array=header_array)
        for batch in chain(batches):
            if dq.qsize() < 5:
                dq.put(batch)
            else:
                time.sleep(0.1)

    dq.put(None)

    # if profiling, dump info into log files in current dir
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("read_worker.log", 'w') as f:
            print(s.getvalue(), file=f)