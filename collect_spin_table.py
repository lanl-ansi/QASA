import os, argparse, math, csv
from dwave.cloud import Client
from collections import namedtuple

CALL_NUM_READS = 10000
CALLS_PER_ROUND = 10
SPIN_TABLE_FILENAME='spin_table.csv'
Range = namedtuple('Range', ['lb', 'ub'])

def main(args):
    h_i = -1*args.h_range
    h_list = []
    while round(h_i, 4) <= args.h_range:
        h_list.append(round(h_i, 4))
        h_i+=args.h_step
    #get spin id's
    with Client.from_config(config_file=os.getenv("HOME")+"/dwave.conf", profile=args.profile, connection_close=True) as client:
        solver = client.get_solver()
        spin_ids = sorted(solver.nodes)
        if args.spin_set != None:
            spin_set = set(args.spin_set)
            filtered_spin_ids = []
            for s in spin_ids:
                if s in spin_set:
                    filtered_spin_ids.append(s)
            spin_ids = filtered_spin_ids
        site_range = Range(*solver.properties['h_range'])
        scaling_factor = 1.0

        for h_i in h_list:
            if h_i != 0:
                if h_i < site_range.lb:
                    scaling_factor = min(scaling_factor, site_range.lb/float(h_i))
                if h_i > site_range.ub:
                    scaling_factor = min(scaling_factor, site_range.ub/float(h_i))
        if scaling_factor < 1.0:
            print('info: rescaling field to {} with scaling factor {}'.format(site_range, scaling_factor))
            h_list = [v*scaling_factor for v in h_list]

    if os.path.isdir(args.directory) is False:
        os.mkdir(args.directory)

    spin_table_file = os.path.join(args.directory, SPIN_TABLE_FILENAME)
    h_visited = []
    if os.path.isfile(spin_table_file):
        with open(spin_table_file) as f:
            reader = csv.DictReader(f) # read rows into a dictionary format
            h_visited = [float(row['h']) for row in reader]
        print('previously collected h: ', h_visited)
    else:
        with open(spin_table_file, 'w') as f:
            header = ['h', 'samples'] + ['spin_{}'.format(s) for s in spin_ids]
            f.write('{}\n'.format(','.join(header)))

    h_list = [el for el in h_list if el not in h_visited]
    print('remaining h values to collect: ', h_list)
    with open(spin_table_file, 'a') as spin_table:
        for h_i in h_list:
            with Client.from_config(config_file=os.getenv("HOME")+"/dwave.conf", profile=args.profile, connection_close=True) as client:
                solver = client.get_solver()
                h = {}
                for i in spin_ids:
                    h[i] = h_i
                params = {
                'auto_scale': False,
                'num_reads': CALL_NUM_READS,
                'flux_drift_compensation': False,
                'annealing_time': args.annealing_time
                }
                print('total num reads: {}'.format(args.num_reads))
                print('')
                print('starting collection:')
                num_reads_remaining = args.num_reads
                num_reads = min(CALL_NUM_READS, num_reads_remaining)

                rounds = int(math.ceil(num_reads_remaining/(CALL_NUM_READS*CALLS_PER_ROUND)))
                solutions_all = {}
                spin_down = [0 for i in range(len(spin_ids))]
                iteration = 1
                retries = 0
                total_collected = 0
                while num_reads_remaining > 0:
                    try:
                        print('')
                        print('  collection round {} of {} (sample_ising calls per round {})'.format(iteration, rounds, CALLS_PER_ROUND))

                        submitted_problems = []
                        num_submitted_reads = 0 # number of reads submitted in this round
                        for i in range(CALLS_PER_ROUND):
                            num_reads = min(CALL_NUM_READS, num_reads_remaining - num_submitted_reads)
                            params['num_reads'] = num_reads
                            if args.spin_reversal_transform_rate != None:
                                params['num_spin_reversal_transforms'] = int(num_reads/args.spin_reversal_transform_rate)


                            print('    submit {} of {} remaining'.format(num_reads, num_reads_remaining-num_submitted_reads))

                            submitted_problems.append(solver.sample_ising(h, {}, **params))
                            num_submitted_reads += num_reads
                            if num_reads_remaining - num_submitted_reads <= 0:
                                break
                        print('    waiting...')
                        solutions_list = []
                        for i, problem in enumerate(submitted_problems):
                            if problem.wait(timeout = args.timeout) is False:
                                raise TimeoutError('    timed out after {} seconds while waiting for response from submitted problem'.format(args.timeout))


                            print('    collect {} of {} calls'.format(i+1, len(submitted_problems)))
                            answers = problem.result()
                            solutions_list.append(answers)
                    except Exception as error:
                        retries += 1
                        print(error)
                        print('    resubmitting round (retries: {})'.format(retries))
                    else:
                        retries = 0
                        num_reads_remaining -= num_submitted_reads
                        for s in solutions_list:
                            for i, solution in enumerate(s['solutions']):
                                spin_down = [sum(x) for x in zip(spin_down, [s['num_occurrences'][i] if solution[j] == -1 else 0 for j in spin_ids])]
                                total_collected += s['num_occurrences'][i]
                        print('    round complete')
                        #print('    num_reads_remaining = {}'.format(num_reads_remaining))
                        iteration += 1
            spin_down_string = ['{:d}'.format(int(v)) for v in spin_down]
            spin_down_string = ','.join(spin_down_string)
            spin_table.write('{:.5f},{:d},{}\n'.format(h_i, total_collected, spin_down_string))
            spin_table.flush()

            print('')
            print('total collected: {}'.format(total_collected))
            assert(total_collected == args.num_reads)
def build_cli_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', '--profile', help='connection details to load from dwave.conf', required=True)

    parser.add_argument('-hr', '--h-range', help='the maximum magnitude of h values to sweep', type=float, default=2.0)
    parser.add_argument('-hs', '--h-step', help='step size between consecutive h values', type=float, default=0.025)
    parser.add_argument('-d', '--directory', help='working directory', required=True)
    parser.add_argument('-nr', '--num-reads', help='number of samples to take for each h value', type=int, default=100000)
    parser.add_argument('-ss', '--spin-set', help='a set of spins that is used to filter the hardware graph', nargs='+', type=int)
    parser.add_argument('-at', '--annealing-time', help='annealing time in microseconds', type=int, default=1)
    parser.add_argument('-srtr', '--spin-reversal-transform-rate', help='the number of reads to take before each spin reversal transform', type=int)
    parser.add_argument('-to', '--timeout', help='the number of seconds before request to D-wave is resubmitted', type=int, default=3000)
    return parser


if __name__ == '__main__':
    parser = build_cli_parser()
    main(parser.parse_args())
