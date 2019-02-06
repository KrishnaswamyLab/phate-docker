import sys
import argparse
import tasklogger

from utils import check_filetype
from run_phate import run_phate_from_file


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run PHATE for visualization of '
        'high-dimensional data.',
        epilog='For help, visit phate.readthedocs.io or '
        'krishnaswamylab.org/get-help',
        add_help=True, allow_abbrev=True)

    io_group = parser.add_argument_group('Data IO')
    filename = io_group.add_mutually_exclusive_group(required=True)
    filename.add_argument('--filename', type=str, default=None,
                          help='Input data. Allowed types: csv, tsv, mtx, '
                          'hdf5/h5 (10X format), directory/zip (10X format)')
    filename.add_argument('--validate', action='store_true', default=False,
                          help='Run PHATE on a test dataset to ensure '
                          'output is correct.')
    sparse = io_group.add_mutually_exclusive_group()
    sparse.add_argument('--sparse', action='store_true',
                        help='Use sparse data format', dest='sparse',
                        default=None)
    sparse.add_argument('--dense', action='store_false',
                        help='Use dense data format', dest='sparse',
                        default=None)
    gene_names = io_group.add_mutually_exclusive_group()
    gene_names.add_argument('--gene-names', action='store_true',
                            help='Use gene name headers in data file'
                            ' (csv, tsv, fcs)',
                            dest='gene_names', default=True)
    gene_names.add_argument('--no-gene-names', action='store_false',
                            help='Do not use gene names'
                            ' (csv, tsv, fcs, mtx)', dest='gene_names',
                            default=True)
    gene_names.add_argument('--gene-name-file', type=str,
                            help='Use gene name headers in FILE'
                            ' (csv, tsv, fcs, mtx)', metavar='FILE',
                            dest='gene_names', default=True)
    cell_names = io_group.add_mutually_exclusive_group()
    cell_names.add_argument('--cell-names', action='store_true',
                            help='Use cell name headers in data file'
                            ' (csv, tsv, fcs)',
                            dest='cell_names', default=True)
    cell_names.add_argument('--no-cell-names', action='store_false',
                            help='Do not use cell names'
                            ' (csv, tsv, fcs, mtx)',
                            dest='cell_names', default=True)
    cell_names.add_argument('--cell-name-file', type=str,
                            help='Use cell name headers in FILE'
                            ' (csv, tsv, fcs, mtx)', metavar='FILE',
                            dest='cell_names', default=True)
    io_group.add_argument('--cell-axis', type=str, choices=['row', 'column'],
                          default='row',
                          help='States whether cells are on rows or columns '
                          '(csv, tsv, mtx)')
    io_group.add_argument('--gene-labels', type=str, default='both',
                          choices=['symbol', 'id', 'both'],
                          help='Choice of gene labels for 10X data'
                          ' (dir, zip, hdf5)')
    io_group.add_argument('--genome', type=str, default=None,
                          help='Genome name for 10X HDF5 data (hdf5)')
    io_group.add_argument('--metadata-channels', type=str, nargs='+',
                          default=['Time', 'Event_length', 'DNA1', 'DNA2',
                                   'Cisplatin', 'beadDist', 'bead1'],
                          help='Names of channels to remove from fcs data (fcs)', metavar='CHANNEL')

    preprocess_group = parser.add_argument_group('Preprocessing')
    cell_filter = preprocess_group.add_mutually_exclusive_group()
    cell_filter.add_argument('--min-library-size', type=int, default=2000,
                             help='Filter cells with less than COUNTS counts',
                             dest='min_library_size', metavar='COUNTS')
    cell_filter.add_argument('--no-cell-filter', action='store_false',
                             default=2000, dest='min_library_size',
                             help='Do not filter cells')
    gene_filter = preprocess_group.add_mutually_exclusive_group()
    gene_filter.add_argument('--min-cells-per-gene', type=int, default=10,
                             help='Filter genes with less than CELLS non-zero cells', dest='min_cells_per_gene',  metavar='CELLS')
    gene_filter.add_argument('--no-gene-filter', action='store_false',
                             default=2000, dest='min_cells_per_gene',
                             help='Do not filter genes')
    libnorm = preprocess_group.add_mutually_exclusive_group()
    libnorm.add_argument('--normalize', action='store_true',
                         default=True, dest='library_size_normalize',
                         help='Normalize cells by total UMI count (library size)')
    libnorm.add_argument('--no-normalize', action='store_false',
                         default=True, dest='library_size_normalize',
                         help='Do not normalize cells')
    transform = preprocess_group.add_mutually_exclusive_group()
    transform.add_argument('--transform', type=str, default='sqrt',
                           choices=['sqrt', 'log', 'arcsinh'],
                           help='Sublinear data transformation function')
    transform.add_argument('--no-transform', action='store_false',
                           default='sqrt',
                           dest='transform', help='Do not transform data')
    preprocess_group.add_argument('--pseudocount', type=float, default=1,
                                  help='Pseudocount to add to genes prior '
                                  'to log transform', metavar='PCOUNT')
    preprocess_group.add_argument('--cofactor', type=float, default=5,
                                  help='Factor by which to divide genes prior '
                                  'to arcsinh transform')

    kernel_group = parser.add_argument_group('Kernel Computation')
    kernel_group.add_argument('-k', '--knn', type=int, default=5, dest='knn',
                              help='Number of nearest neighbors on which to '
                              'build kernel')
    decay = kernel_group.add_mutually_exclusive_group()
    decay.add_argument('-a', '--decay', type=int, default=15, dest='decay',
                       help='Sets decay rate of kernel tails')
    decay.add_argument('--no-decay', action='store_false', default=15,
                       dest='decay', help='Do not use alpha decay')
    pca = kernel_group.add_mutually_exclusive_group()
    pca.add_argument('--pca', type=int, default=100, dest='n_pca',
                     help='Number of principal components to use for '
                     'neighborhoods')
    pca.add_argument('--no-pca', action='store_false', default=100,
                     dest='n_pca', help='Do not use PCA')
    kernel_group.add_argument('--knn-dist', type=str, default='euclidean',
                              help='Distance metric to use for calculating '
                              'neighborhoods. Recommended values are '
                              '"euclidean" and "cosine"', metavar='DISTANCE')
    kernel_group.add_argument('-t', '--threads', type=int, default=1,
                              help='Use THREADS threads. If -1 all CPUs are used',
                              metavar='THREADS', dest='random_state')
    kernel_group.add_argument('--seed', type=int, default=None,
                              help='Integer random seed', metavar='SEED',
                              dest='random_state')
    verbose = kernel_group.add_mutually_exclusive_group()
    verbose.add_argument('-v', '--verbose', action='store_true', default=True,
                         help='Print verbose output')
    verbose.add_argument('-q', '--quiet', action='store_false', default=True,
                         help='Do not print verbose output', dest='verbose')
    verbose.add_argument('-vv', '--debug', action='store_true', default=False,
                         help='Print debugging output', dest='debug')

    phate_group = parser.add_argument_group('PHATE')
    phate_group.add_argument('-n', '--n-components', type=int, default=2,
                             help='Number of dimensions in which the data will'
                             ' be embedded for PHATE', metavar='N')
    phate_group.add_argument('--gamma', type=float, default=1,
                             help='Informational distance constant between -1 '
                             'and 1. `gamma=1` gives the PHATE log potential, '
                             '`gamma=0` gives a square root potential')
    phate_group.add_argument('--t-phate', type=str, default='auto',
                             help='Level of diffusion for PHATE',
                             metavar='T')
    phate_group.add_argument('--output', type=str, default='phate.csv',
                             help='Output CSV file to save low-dimensional '
                             'embedding', metavar='FILE')

    args = parser.parse_args()

    if args.validate:
        tasklogger.set_level(2)
        tasklogger.log_info("Running MAGIC validation.")
        args.filename = "https://github.com/KrishnaswamyLab/scprep/raw/master/data/test_data/test_small.csv"
        args.sparse = True
        args.gene_names = True
        args.cell_names = True
        args.cell_axis = "row"
        args.gene_labels = "both"
        args.gene_labels = None
        args.metadata_channels = None
        args.min_library_size = 1
        args.min_cells_per_gene = 1
        args.library_size_normalize = True
        args.transform = 'sqrt'
        args.pseudocount = None
        args.cofactor = None
        args.knn = 3
        args.decay = 20
        args.n_pca = None
        args.knn_dist = "euclidean"
        args.n_jobs = 1
        args.random_state = 42
        args.verbose = True
        args.debug = True
        args.n_components = 2
        args.gamma = 1
        args.t_phate = "auto"
        args.output = "phate-validate.csv"

    # fix t argument
    if args.t_phate != 'auto':
        try:
            args.t_phate = int(args.t_phate)
        except TypeError:
            parser.error(
                "argument --t-phate: invalid int value: '{}'".format(args.t_phate))

    # fix debug argument
    if args.debug:
        args.verbose = 2
    del args.debug

    # store None values where appropriate
    if args.decay is False:
        args.decay = None
    if args.n_pca is False:
        args.n_pca = None
    if args.min_library_size is False:
        args.min_library_size = None
    if args.min_cells_per_gene is False:
        args.min_cells_per_gene = None

    # check for inappropriately set defaults
    try:
        filetype = check_filetype(args.filename)
    except RuntimeError as e:
        parser.error(str(e))
    if filetype not in ['csv', 'tsv', 'csv.gz', 'tsv.gz', 'fcs']:
        if '--gene-names' not in sys.argv:
            args.gene_names = None
        else:
            parser.error(
                "Cannot handle --gene-names with {} file".format(filetype[1:]))
        if '--cell-names' not in sys.argv:
            args.cell_names = None
        else:
            parser.error(
                "Cannot handle --cell-names with {} file".format(filetype[1:]))
    if filetype not in ['csv', 'tsv', 'csv.gz', 'tsv.gz', 'mtx']:
        if '--cell-axis' not in sys.argv:
            args.cell_axis = None
        else:
            parser.error(
                "Cannot handle --cell-axis with {} file".format(filetype[1:]))
    if filetype not in ['dir', 'zip', 'hdf5', 'h5']:
        if '--gene-labels' not in sys.argv:
            args.gene_labels = None
        else:
            parser.error(
                "Cannot handle --gene-labels with {} file".format(filetype[1:]))
    if filetype not in ['hdf5', 'h5']:
        if '--genome' not in sys.argv:
            args.genome = None
        else:
            parser.error(
                "Cannot handle --genome with {} file".format(filetype[1:]))
    if filetype not in ['fcs']:
        if '--metadata-channels' not in sys.argv:
            args.metadata_channels = None
        else:
            parser.error(
                "Cannot handle --metadata-channels with {} file".format(filetype[1:]))

    # check for inappropriately set parameters
    if not args.transform == 'log':
        if '--pseudocount' in sys.argv:
            parser.error(
                "Cannot handle --pseudocount with --transform {}".format(args.transform))
        else:
            args.pseudocount = None
    if not args.transform == 'arcsinh':
        if '--cofactor' in sys.argv:
            parser.error(
                "Cannot handle --cofactor with --transform {}".format(args.transform))
        else:
            args.cofactor = None

    return args


if __name__ == "__main__":
    args = parse_args()
    tasklogger.set_level(args.verbose)
    tasklogger.log_debug(
        "Running PHATE with arguments {}".format(args.__dict__))
    run_phate_from_file(**(args.__dict__))
