# Thread Model

Here I will outline how Butter-eel is designed in terms of using cores (procs) and threads.

On the python side, I use a multiprocessing method, and on the C side in pyslow5 we use a threading method.

## Overview

- 1 proc is used for reading data, and puts that data into an input queue.
- 1 proc is used for writing data, and pulls that data from the output queue
    - Multiple threads are used when reading the input data
- multiple procs are used to pull batches of reads from the input queue, basecall them, and put the results into the output queue
- some procs (~2) are used by guppy/dorado to handle some interproccess communication, but usually minimal load.

## More detail

### Reading data

The 1 proc used for reading data, will exclusively use the `--slow5_threads` argument to set the number of theads used within th pyslow5 library to read and decompress the blow5 data.
This is generally quite fast, so the default of 4 threads is usually enough for single consumer level GPU workloads.
For HPC or multi-GPU setups, the basecalling will consume the data faster than it can be read. So increasing `--slow5_threads` will help scale with compute.
There is a limit on how much data will be stored in the input queue, so ram doesn't get out of hand. This is controlled by `--max_read_queue_size`. So in some cases this might need to be changed, either to limit ram usage, or to increase the input queue so more batches can be stored, for more worker procs to access.

The data is stored in the input queue as batches of reads, set by `--slow5_batchsize`, which can also be tweaked to make the reading and processing more efficient depending on the systems being used. As of dorado-server v7.4.12, the value of procs x slow5_batchsize > dorado-server gpu batch size (found in the basecalling logs). When this rule isn't met, there will be a pause of 30s for every batch to be processed, resulting in a large increase in basecalling time. A batch size of 4000 (default) is usually fine, though may need to be increased for GPUs with large VRAM (>40gb)

### Writing data

A single proc is used to read batches of reads from the output queue, and write them to the appropriate file/s.

The 1 proc used to write data should be sufficient to keep up with compute. If this is ever not the case, please create an issue and let me know,and i'll add an argument to increase this.

This 1 proc will spawn multiple threads, controlled by `--slow5_threads`, to decompress the batch of reads fetched.

### Processing data

By default, 4 procs are used for processing the data. This creates 4 separate client connections to the server which then each take a batch of N reads (set by `--slow5_batchsize` in the reading data step), packages each read into a data structure to then be sent to the basecalling server.

Once all reads are submitted, each process waits for the reads to be returned as basecalled, and handles things like sequencing summary and barcode summary output, then pushes the basecalled data into the output queue for writing.

When using multiple GPU systems, increasing the number of procs used with `--procs` should help scale out the data processing. However, it has to be scaled with `--slow5_threads` so enough batches of data are being populated into the input queue to fully utilise each proc and the full compute power of the GPUs. This will be different on different systems, as it depends on the speed of the storage system, type and memory size of GPU, how many GPUs, CPU speed, and system RAM.

An example of getting full utilisation out of a HPC node, with 40 cores and 4xTesla V100 (16GB) GPUs, reading from spinning disk storage (not nvme) would be:
```
–slow5_threads 10 –slow5_batchsize 4000 –procs 20
```
This will use 10 threads to get 4000 reads, where each thread decompresses reads for each batch, then submits that to the input queue.
The input queue has a default max size of 20,000 reads. 20 procs will consume 80,000 reads at a time (20 batches of 4000). This means the queue will be saturated, and allow for efficient feeding of the GPUs as the reading of the data should be faster than HAC or SUP basecalling. For systems with more GPUs or with larger VRAM, the `--max_read_queue_size` may need to be increased. Note that this value is used for both the GPU queue, and the buttery-eel queue limit. This is helpful to know when optimising RAM usage, as you essentially have double the amount of reads in this argument in flight in the system, plus whatever is being processed by the GPU and in the write queue, if the system can saturate the queues.

So we have:
- 1 proc for read
- 10 threads used in the reading
- 20 procs for working
- 1 proc for writing

That means we are using 32 of the 40 procs available on the node, with 8 left over for guppy/dorado to use for various basecalling/methylation calling/interprocesses communication type work.

