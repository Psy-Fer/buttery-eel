#!/usr/bin/env python3

import argparse
import logging
import sys
import os
import time
from pathlib import Path
import numpy as np

import pyslow5

from pyguppy_client_lib.pyclient import PyGuppyClient
from pyguppy_client_lib import helper_functions

# import json

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)



def start_guppy_server(log, config, port, guppy_bin, max_queued_reads, chunk_size):
    """
    Start guppy server in separate process
    """
    server_args = ["--log_path", log,
                   "--config", config,
                   "--port", port,
                   "--use_tcp",
                   "-x", "auto",
                   "--max_queued_reads", max_queued_reads,
                   "--chunk_size", chunk_size]
    ret = helper_functions.run_server(server_args, bin_path=guppy_bin)
    return ret


def calibration(digitisation, range):
    """
    input:
        digitisation: float
        range: float
    output:
        scale: float
    """
    return range / digitisation


def write_fastq(fq, header, seq, qscore):
    fq.write("{}\n".format(header))
    fq.write("{}\n".format(seq))
    fq.write("+\n")
    fq.write("{}\n".format(qscore))


# def get_model_info(config, guppy_bin):
#     config = os.path.join(guppy_bin,"../data/", config)
#     model = ""
#     with open(config, 'r') as f:
#         for line in f:
#             line = line.strip('\n')
#             line = line.replace(" ", "")
#             line = line.split("=")
#             if len(line) > 1:
#                 if line[0] == "model_file":
#                     model = line[1]
#                     break
#
#     if len(model) == 0:
#         logger.warning("could not deduce model for fastq writing, falling back to default, model_version_id=conf")
#         model_version_id = config
#     else:
#         model_json = os.path.join(guppy_bin,"../data/", model)
#         with open(model_json, 'r') as f:
#             jdata = json.load(f)
#             model_version_id = jdata["version"]["id"]
#
#     return model_version_id


def main():

    VERSION = "v0.0.1"

    parser = MyParser(description="buttery-eel - wrapping guppy for file agnostic basecalling",
    epilog="Citation:...")

    """
    Example:
    ./buttery-eel.py --guppy_bin /home/jamfer/Downloads/ont-guppy-6.1.3/bin --chunk_size 200 --max_queued_reads 1000 --port 5558 -i ~/Data/bench/1_slow5/PAF25452_pass_bfdfd1d8_11.blow5 -o ~/Data/bench/buttery_test/test.fastq
    """

    # Args for the wrapper, and then probably best to just have free form args for guppy
    parser.add_argument("-i", "--input",
                        help="input blow5 file for basecalling")
    parser.add_argument("-o", "--output",
                        help="output .fastq file to write")
    parser.add_argument("--guppy_bin", type=Path,
                        help="path to ont_guppy/bin folder")
    parser.add_argument("--config", default="dna_r9.4.1_450bps_fast.cfg",
                        help="basecalling model config")
    parser.add_argument("--port", default="5558",
                        help="port to use between server/client")
    parser.add_argument("--log", default="buttery_guppy_logs",
                        help="guppy log folder path")
    parser.add_argument("--max_queued_reads", default="2000",
                        help="Number of reads to send to guppy server queue")
    parser.add_argument("--chunk_size", default="2000",
                        help="signal chunk size, lower this for lower VRAM GPUs")
    parser.add_argument("-v", "--version", action='version', version="buttery-eel - wraping guppy for file agnostic basecalling version: {}".format(VERSION),
                        help="Prints version")
    parser.add_argument("--debug", action="store_true",
                        help="Set logging to debug mode")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # turn on debug mode
    logger = logging.getLogger(__name__)
    if args.debug == 1:
        lev = logging.DEBUG
    else:
        lev = logging.WARNING

    logging.basicConfig(format='%(asctime)s - [%(levelname)s]: %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S', level=lev)

    logger.info("args: {}".format(args))

    logger.debug("starting: guppy_basecall_server")
    serv = start_guppy_server(args.log, args.config, args.port, args.guppy_bin, args.max_queued_reads, args.chunk_size)
    print(serv)
    server = serv[0]
    p = serv[1]
    logger.debug("guppy_basecall_server started...")

    logger.debug("Connecting to server with client...")
    client = PyGuppyClient(
    "127.0.0.1:{}".format(args.port),
    "{}".format(args.config))

    logger.debug("Setting params...")
    client.set_params({"priority": PyGuppyClient.high_priority})

    logger.debug("Connecting...")
    client.connect()
    logger.debug("Connected, testing connection...")
    logger.debug("Connection status:")
    print(client.get_status())

    # print(client.get_barcode_kits("127.0.0.1:{}".format(args.port), 10))
    # print(client.get_protocol_version())
    # print(client.get_server_information("127.0.0.1:{}".format(args.port), 10))
    # print(client.get_software_version())

    fq = open(args.output, 'w')
    s5 = pyslow5.Open(args.input, 'r')
    reads = s5.seq_reads()

    get_raw = True
    total_reads = 0
    read_counter = 0
    done = 0
    skipped = []
    for read in reads:
        # sys.exit()
        read_id = read['read_id']
        # calculate scale
        scale = calibration(read['digitisation'], read['range'])
        result = False
        tries = 0
        while not result:
            result = client.pass_read(
                    helper_functions.package_read(
                        read_id=read_id,
                        raw_data=np.frombuffer(read['signal'], np.int16),
                        daq_offset=read['offset'],
                        daq_scaling=scale,
                    )
                )
            if tries > 1:
                time.sleep(1)
            tries += 1
            if tries >= 5:
                if not result:
                    logging.warning("Skipped a read: {}".format(read_id))
                    skipped.append(read_id)
                    break

        if not result:
            continue

        else:
            read_counter += 1
            total_reads += 1

        if read_counter >= 1000:
            while done < read_counter:
                bcalled = client.get_completed_reads()
                if not bcalled:
                    time.sleep(0.5)
                    continue
                else:
                    for call in bcalled:
                        done += 1
                        if len(call) != 1:
                            # possible split reads?
                            logging.warning("Call is longer than 1: {}".format(len(call)))


                        # @read_id runid=bf... sampleid=NA12878_SRE read=476 ch=38 start_time=2020-10-26T19:58:23Z model_version_id=2021-05-17_dna_r9.4.1_minion_96_29d8704b
                        # model_version_id = get_model_info(args.config, args.guppy_bin)
                        header = "@{} model_version_id={}".format(call[0]['metadata']['read_id'], call[0]['metadata']['model_version_id'])
                        sequence = call[0]['datasets']['sequence']
                        qscore = call[0]['datasets']['qstring']
                        write_fastq(fq, header, sequence, qscore)
            done = 0
            read_counter = 0

    if read_counter > 0:
        while done < read_counter:
            bcalled = client.get_completed_reads()
            if not bcalled:
                time.sleep(0.5)
                continue
            else:
                for call in bcalled:
                    done += 1
                    if len(call) != 1:
                        # possible split reads?
                        logging.warning("Call is longer than 1: {}".format(len(call)))

                    # model_version_id = get_model_info(args.config, args.guppy_bin)
                    header = "@{} model_version_id={}".format(call[0]['metadata']['read_id'], call[0]['metadata']['model_version_id'])
                    sequence = call[0]['datasets']['sequence']
                    qscore = call[0]['datasets']['qstring']
                    write_fastq(fq, header, sequence, qscore)
        done = 0
        read_counter = 0


    print("Processed {} reads".format(total_reads))
    print("skipped {} reads".format(len(skipped)))
    # close file
    fq.close()
    s5.close()

    # terminate guppy_basecall_server
    client.disconnect()
    server.terminate()


if __name__ == '__main__':
    main()
