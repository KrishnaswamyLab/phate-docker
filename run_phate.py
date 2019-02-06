import pandas as pd
import scprep
import phate
import numpy as np
import tasklogger

from utils import check_filetype, check_load_args, check_transform_args


def run_phate_from_file(
        filename,
        # data loading params
        sparse=True,
        gene_names=None,
        cell_names=None,
        cell_axis=None,
        gene_labels=None,
        allow_duplicates=None,
        genome=None,
        metadata_channels=None,
        # filtering params
        min_library_size=2000,
        min_cells_per_gene=10,
        # normalization params
        library_size_normalize=True,
        transform='sqrt',
        pseudocount=None,
        cofactor=None,
        # kernel params
        knn=5,
        decay=15,
        n_pca=100,
        knn_dist='euclidean',
        n_jobs=1,
        random_state=42,
        verbose=1,
        # phate params
        n_components=2,
        t_phate='auto',
        gamma=1,
        mds_dist='euclidean',
        mds='metric',
        # output params
        output='phate.csv',
        validate=False):
    """Run PHATE and MAGIC on a file

    Parameters
    ----------
    filename : str
        Allowed types: csv, tsv, mtx, hdf5/h5 (10X format),
        directory/zip (10X format)
    sparse : bool (recommended: True for scRNAseq, False for CyTOF)
        Force data sparsity. If `None`, sparsity is determined by data type.
    gene_names : str, list or bool
        Allowed values:
        - if filetype is csv or fcs, `True` says gene names are data
        headers, `str` gives a path to a separate csv or tsv file containing
        gene names, list gives an array of gene names, `False` means
        no gene names are given
        - if filetype is mtx, `str` gives a path to a separate csv or tsv file
        containing gene names, list gives an array of gene names, or `False`
        means no gene names are given
        - if filetype is hdf5, h5, directory or zip, must be `None`.
    cell_names : str, list or bool
        Allowed values:
        - if filetype is csv or fcs, `True` says cell names are data
        headers, `str` gives a path to a separate csv or tsv file containing
        cell names, list gives an array of cell names, `False` means
        no cell names are given
        - if filetype is mtx, `str` gives a path to a separate csv or tsv file
        containing cell names, list gives an array of cell names, or `False`
        means no gene names are given
        - if filetype is hdf5, h5, directory or zip, must be `None`.
    cell_axis : {'row', 'column'}
        States whether cells are on rows or columns. If cell_axis=='row',
        data is of shape [n_cells, n_genes]. If cell_axis=='column', data is of
        shape [n_genes, n_cells]. Only valid for filetype mtx and csv
    gene_labels : {'symbol', 'id', 'both'}
        Choice of gene labels for 10X data. Recommended: 'both'
        Only valid for directory, zip, hdf5, h5
    allow_duplicates : bool
        Allow duplicate gene names in 10X data. Recommended: True
        Only valid for directory, zip, hdf5, h5
    genome : str
        Genome name. Only valid for hdf5, h5
    metadata_channels : list of str (recommended: ['Time', 'Event_length', 'DNA1', 'DNA2', 'Cisplatin', 'beadDist', 'bead1'])
        Names of channels in fcs data which are not real measurements.
        Only valid if datatype is fcs.
    min_library_size : int or `None`, optional (default: 2000)
        Cutoff for library size normalization. If `None`,
        library size filtering is not used
    min_cells_per_gene : int or `None`, optional (default: 10)
        Minimum non-zero cells for a gene to be used. If `None`,
        genes are not removed
    library_size_normalize : `bool`, optional (default: True)
        Use library size normalization
    transform : {'sqrt', 'log', 'arcsinh', None}
        How to transform the data. If `None`, no transformation is done
    pseudocount : float (recommended: 1)
        Number of pseudocounts to add to genes prior to log transformation
    cofactor : float (recommended: 5)
        Factor by which to divide genes prior to arcsinh transformation
    knn : int, optional, default: 10
        number of nearest neighbors on which to build kernel
    decay : int, optional, default: 15
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used
    n_pca : int, optional, default: 100
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        roughly log(n_samples) time.
    knn_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean', 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph.
    n_jobs : integer, optional, default: 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging.
        For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
        n_jobs = -2, all CPUs but one are used
    random_state : integer or numpy.RandomState, optional, default: None
        The generator used to initialize random PCA
        If an integer is given, it fixes the seed
        Defaults to the global `numpy` random number generator
    verbose : `int` or `boolean`, optional (default: 1)
        If `True` or `> 0`, print status messages
    n_components : int, optional, default: 2
        number of dimensions in which the data will be embedded for PHATE
    mds_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean' and 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for MDS
    mds : string, optional, default: 'metric'
        choose from ['classic', 'metric', 'nonmetric'].
        Selects which MDS algorithm is used for dimensionality reduction
    gamma : float, optional, default: 1
        Informational distance constant between -1 and 1.
        `gamma=1` gives the PHATE log potential, `gamma=0` gives
        a square root potential.
    t_phate : int, optional, default: 'auto'
        power to which the diffusion operator is powered for PHATE.
        This sets the level of diffusion. If 'auto', t is selected
        according to the knee point in the Von Neumann Entropy of
        the diffusion operator
    output : str, optional (default: 'phate.csv')
        Output CSV file to save low-dimensional embedding
    """
    # check arguments
    filetype = check_filetype(filename)
    load_fn, load_kws = check_load_args(filetype,
                                        sparse=sparse,
                                        gene_names=gene_names,
                                        cell_names=cell_names,
                                        cell_axis=cell_axis,
                                        gene_labels=gene_labels,
                                        allow_duplicates=allow_duplicates,
                                        genome=genome,
                                        metadata_channels=metadata_channels)
    transform_fn, transform_kws = check_transform_args(transform=transform,
                                                       pseudocount=pseudocount,
                                                       cofactor=cofactor)

    # load data
    # example: scprep.io.load_csv("data.csv")
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.io
    data = load_fn(filename, **load_kws)
    data = scprep.sanitize.check_numeric(data, copy=True)

    # filter data
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.filter
    if min_library_size is not None:
        data = scprep.filter.filter_library_size(data,
                                                 cutoff=min_library_size)
    if min_cells_per_gene is not None:
        data = scprep.filter.remove_rare_genes(data,
                                               min_cells=min_cells_per_gene)

    # normalize data
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.normalize
    if library_size_normalize:
        data = scprep.normalize.library_size_normalize(data)

    # transform data
    # example: data = scprep.transform.sqrt(data)
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.transform
    if transform is not None:
        data = transform_fn(data, **transform_kws)

    # run PHATE
    # https://phate.readthedocs.io/
    phate_op = phate.PHATE(knn=knn, decay=decay, t=t_phate, n_pca=n_pca,
                           knn_dist=knn_dist,
                           n_jobs=n_jobs, random_state=random_state,
                           verbose=verbose, n_components=n_components,
                           gamma=gamma, mds_dist=mds_dist, mds=mds)
    phate_data = phate_op.fit_transform(data)

    # save as csv
    phate_data = pd.DataFrame(
        phate_data,
        columns=["PHATE{}".format(i + 1) for i in range(n_components)],
        index=data.index if hasattr(data, 'index')
        else np.arange(phate_data.shape[0]))
    if cell_axis in ['col', 'column']:
        phate_data = phate_data.T
    phate_data.to_csv(output)
    if validate:
        correct_phate_data = scprep.io.load_csv(
            'https://raw.githubusercontent.com/KrishnaswamyLab/phate-docker/'
            'master/phate-validate.csv', sparse=False)
        try:
            np.testing.assert_equal(scprep.utils.toarray(phate_data),
                                    scprep.utils.toarray(correct_phate_data))
            tasklogger.log_debug(
                "Validation complete, output is equal to expected")
        except AssertionError:
            np.testing.assert_allclose(
                scprep.utils.toarray(phate_data),
                scprep.utils.toarray(correct_phate_data),
                atol=1e-14)
            tasklogger.log_debug(
                "Validation complete, output is numerically equivalent to expected")
