import os
import scprep


def check_filetype(filename):
    if os.path.isdir(filename):
        filetype = 'dir'
    elif os.path.isfile(filename) or filename.startswith("http"):
        # get fle extension
        filetype = filename.split('.')[-1]
        if filetype == 'gz':
            # check for extension before .gz
            filetype = ".".join(filename.split('.')[-2:])
    else:
        raise RuntimeError("file {} not found".format(filename))
    return filetype


def check_load_args(filetype,
                    sparse=True,
                    gene_names=None,
                    cell_names=None,
                    cell_axis=None,
                    gene_labels=None,
                    allow_duplicates=None,
                    genome=None,
                    metadata_channels=None):
    # store relevant keyword arguments
    if filetype == 'zip':
        load_fn = scprep.io.load_10X_zip
        load_kws = {'sparse': sparse,
                    'gene_labels': gene_labels}
    elif filetype == 'dir':
        load_fn = scprep.io.load_10X
        load_kws = {'sparse': sparse,
                    'gene_labels': gene_labels}
    elif filetype in ['hdf5', 'h5']:
        load_fn = scprep.io.load_10X_HDF5
        load_kws = {'sparse': sparse,
                    'gene_labels': gene_labels,
                    'genome': genome}
    elif filetype in ['tsv', 'tsv.gz']:
        load_fn = scprep.io.load_tsv
        load_kws = {'sparse': sparse,
                    'gene_names': gene_names,
                    'cell_names': cell_names,
                    'cell_axis': cell_axis}
    elif filetype in ['csv', 'csv.gz']:
        load_fn = scprep.io.load_csv
        load_kws = {'sparse': sparse,
                    'gene_names': gene_names,
                    'cell_names': cell_names,
                    'cell_axis': cell_axis}
    elif filetype == 'mtx':
        load_fn = scprep.io.load_mtx
        load_kws = {'sparse': sparse,
                    'gene_names': gene_names,
                    'cell_names': cell_names,
                    'cell_axis': cell_axis}
    elif filetype == 'fcs':
        load_fn = scprep.io.load_fcs
        load_kws = {'sparse': sparse,
                    'gene_names': gene_names,
                    'cell_names': cell_names,
                    'metadata_channels': metadata_channels}
    else:
        raise RuntimeError("filetype {} not recognized. Expected 'csv', "
                           "'tsv', 'csv.gz', 'tsv.gz', 'mtx', 'zip', 'hdf5', "
                           "'h5', 'fcs' or a "
                           "directory".format(filetype))

    # check arguments are missing where appropriate
    load_args = ['gene_names', 'cell_names', 'cell_axis',
                 'sparse', 'gene_labels', 'allow_duplicates',
                 'metadata_channels']
    for arg in load_args:
        if arg == 'sparse':
            # allow None
            pass
        elif arg in load_kws:
            assert eval(arg) is not None, \
                "Expected {} not None for filetype {}".format(arg, filetype)
        else:
            assert eval(arg) is None, \
                "Expected {} to be None for filetype {}. Got {}".format(
                    arg, filetype, eval(arg))
    return load_fn, load_kws


def check_transform_args(transform='sqrt',
                         pseudocount=None,
                         cofactor=None):
    # store relevant keyword arguments
    if transform == 'sqrt':
        transform_fn = scprep.transform.sqrt
        transform_kws = {}
    elif transform == 'log':
        transform_fn = scprep.transform.log
        transform_kws = {'cofactor': cofactor}
    elif transform == 'arcsinh':
        transform_fn = scprep.transform.arcsinh
        transform_kws = {'pseudocount': pseudocount}
    elif transform is None:
        transform_kws = {}
    else:
        raise RuntimeError("transformation {} not recognized. "
                           "Choose from ['sqrt', 'log', 'arcsinh', "
                           "None]".format(transform))

    # check arguments are missing where appropriate
    transform_args = ['pseudocount', 'cofactor']
    for arg in transform_args:
        if arg in transform_kws:
            assert eval(arg) is not None, \
                "Expected {} not None for {} transformation".format(
                    arg, transform)
        else:
            assert eval(arg) is None, \
                "Expected {} to be None for {} transformation. Got {}".format(
                    arg, transform, eval(arg))
    return transform_fn, transform_kws
