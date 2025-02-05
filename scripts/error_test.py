import pyslow5
import sys
import multiprocessing as mp
import time

from io import StringIO
from contextlib import contextmanager, redirect_stdout

try:
    from pybasecall_client_lib.pyclient import PyBasecallClient as pclient
    from pybasecall_client_lib import helper_functions

except ImportError:
    # maybe i can do a version check or something? hard to do from this side of the software
    print("Could not load pybasecall, trying for version earlier versions <=7.2.15 pyguppy lib")
    try:
        from pyguppy_client_lib.pyclient import PyGuppyClient as pclient
        from pyguppy_client_lib import helper_functions
    except ImportError:
        print("Can't import pybasecall_client_lib or pyguppy_client_lib, please check environment and try again.")
        sys.exit(1)



def reader_proc(fname):
    print("before")
    s5 = pyslow5.Open(fname, 'r', DEBUG=1)
    print("after")
    print(s5.get_num_read_groups())


def worker_proc(address, config, fname):
    print("worker_proc")
    print("worker_proc before")
    sys.exit(1)


    connected = False
    for t in range(3):
        print("t:", t)
        if connected:
            continue
        client_sub = pclient(address=address, config=config)
        connected = True
        with client_sub as client:
            for i in range(15):
                print("i:", i)
                print("status: {}".format(client.get_status()))
                print(type(client.get_status()))
                if client.get_status() == client.status.connected:
                    print("client is connected!")
                else:
                    print("client is not connected!")
                    connected = False
                # print("error_message: {}".format(client.get_error_message()))
                # print("internal state: {}".format(client.get_internal_state()))
                # print("server information: {}".format(client.get_server_information(address, 2000)))
                # if i == 3 and t < 2:
                #     result = client.disconnect()
                #     print("disconnect result:", result)
                #     connected = False
                if not connected:
                    break
                time.sleep(i)
                # s5 = pyslow5.Open(fname, 'r', DEBUG=1)
                # print(s5.get_num_read_groups())
            print("worker_proc after forloop")
    if not connected:
        raise ConnectionAbortedError
    print("worker_proc with")




# region start basecaller
@contextmanager
def start_guppy_server_and_client(basecaller_bin):
    """
    Thanks to Alex Payne for basis of this code to appropriately handle the server and client connections
    https://gist.github.com/alexomics/043bb120c74161e5b93e1b68fb00206c

    Starts server and connects client
    TODO: allow a connection to existing guppy server
    """
    server_args = []

    server_args.extend(["--log_path", "./logs",
                        "--config", "dna_r10.4.1_e8.2_400bps_5khz_fast.cfg",
                        "--port", "auto",
                        "--use_tcp",
                        "--device", "cuda:all",
                        # "--max_queued_reads", args.max_queued_reads,
                        # "--chunk_size", args.chunk_size,
                        ])
    # the high priority queue uses a different batch size which alters the basecalls when called with dorado
    # leaving this on default should set it to medium and give 'correct' results
    # funny enough, it will call R9.4.1 data at higher accuracy, and the opposite impact on R10.4.1
    # params = {"priority": PyGuppyClient.high_priority}
    params = {}

    # This function has it's own prints that may want to be suppressed
    with redirect_stdout(StringIO()) as fh:
        server, port = helper_functions.run_server(server_args, bin_path=basecaller_bin)

    if port == "ERROR":
        raise RuntimeError("Server couldn't be started")

    if port.startswith("ipc"):
        address = "{}".format(port)
    else:
        address = "localhost:{}".format(port)
    client = pclient(address=address, config="dna_r10.4.1_e8.2_400bps_5khz_fast.cfg")


    print("Setting params...")
    client.set_params(params)
    print("Connecting...")
    try:
        with client:
            yield [client, address, "dna_r10.4.1_e8.2_400bps_5khz_fast.cfg", params]
    finally:
        server.terminate()


def main():

    mp.set_start_method('spawn')

    # filename = sys.argv[1]

    # reader = mp.Process(target=reader_proc, args=(filename,), name='read_worker')
    # reader.start()

    # reader.join()
    # print("reader exitcode:", reader.exitcode)
    # if reader.exitcode < 0:
    #     sys.exit(1)
    # print("done")

    # filename = sys.argv[1]
    basecaller_bin = sys.argv[1]
    
    with start_guppy_server_and_client(basecaller_bin) as client_one:
        client, address, config, params = client_one
        print(client)
        print("guppy/dorado started...")
        print("basecaller version:", "{}".format(".".join([str(i) for i in client.get_software_version()])))
        print()

        worker = mp.Process(target=worker_proc, args=(address, config, "t.blow5",), name='worker_proc')
        worker.start()

        worker.join()
        print("worker exitcode:", worker.exitcode)
        if worker.exitcode != 0:
            print("killing program")
            sys.exit(1)
        print("last print before closing server")
    print("Server closed")


if __name__ == '__main__':
    main()
