import logging
import os

logger = logging.getLogger(__name__)

ASSEMBLERS_TO_RUN = [
    'masurca',
    'novoplasty',
    'shovill',
    'skesa',
    'spades',
    'unicycler'
]
ASSEMBLERS_REQUIRING_APC = [
    'masurca',
    'shovill',
    'skesa',
    'spades',
    'unicycler'
]


def importFile(toil, source_path, scheme="file"):
    """
    Import the source path using the specified shecme into the file
    store.

    Parameters
    ----------
    toil : toil.common.Toil
        instance of a Toil context manager
    source_path : str
        path of the source to import
    scheme : str
        scheme used for the source URL

    Returns
    -------
    str
        id of the imported source in the file store
    """
    try:
        file_path = scheme + "://" + source_path
        logger.info("Importing file {0}".format(file_path))
        file_id = toil.importFile(file_path)

    except Exception as exc:
        logger.info("Importing file {0} failed: {1}".format(
            file_path, exc))
        file_id = None

    return file_id


def importReadFiles(toil, data_path, plate_spec, well_specs, scheme="file"):
    """
    Import the read one and two sources from the data path containing
    the specified plate and well FASTQ source (follows Addgene's
    folder naming convention, mostly).

    Parameters
    ----------
    toil : toil.common.Toil
        instance of a Toil context manager
    data_path : str
        path containing plate and well FASTQ source
    plate_spec : str
        specification of the plate
    well_spec : list
        specification of the wells
    scheme : str
        scheme used for the source URL

    Returns
    -------
    tuple of lists
        ids of the imported read one and two files in the file store
    """
    read_one_file_ids = []
    read_two_file_ids = []

    for well_spec in well_specs:

        read_one_file_ids.append(importFile(
            toil, os.path.join(
                data_path, "{0}_FASTQ".format(plate_spec),
                "{0}_{1}_R1_001.fastq.gz".format(plate_spec, well_spec)),
            scheme))

        read_two_file_ids.append(importFile(
            toil, os.path.join(
                data_path, "{0}_FASTQ".format(plate_spec),
                "{0}_{1}_R2_001.fastq.gz".format(plate_spec, well_spec)),
            scheme))

    return read_one_file_ids, read_two_file_ids


def importContigsFile(toil, data_path, scheme="file"):
    """
    Import the contigs source from the data path containing the FASTA
    source.

    Parameters
    ----------
    data_path : str
        path containing the FASTA source
    scheme : str
        scheme used for the source URL

    Returns
    -------
    str
        id of the imported contigs file in the file store
    """
    contigs_file_id = importFile(
        toil, os.path.join(data_path, "contigs.fasta"), scheme)
    return contigs_file_id


def readGlobalFile(fileStore, file_id, *cmps):
    """
    Read the file corresponding to the specified id from the file
    store into the local temporary directory using the specified path
    components

    Parameters
    ----------
    fileStore : toil.fileStore.FileStore
        instance of a Toil file store
    file_id : str
        id of the file in the file store
    cmps : list
        path components of the file read into the local temporary
        directory

    Returns
    -------
    str
        absolute path to the file in the local temporary directory
    """
    try:
        file_path = os.path.join(fileStore.localTempDir, *cmps)
        logger.info("Reading global file {0}".format(file_path))
        fileStore.readGlobalFile(file_id, userPath=file_path)

    except Exception as exc:
        logger.info("Reading global file {0} failed: {1}".format(
            file_path, exc))
        file_path = None

    return file_path


def writeGlobalFile(fileStore, *cmps):
    """
    Write the file with the specified path components from the local
    temporary directory into the file store and return the id

    Parameters
    ----------
    fileStore : toil.fileStore.FileStore
        instance of a Toil file store
    cmps : list
        path components of the file in the local temporary directory

    Returns
    -------
    str
        id of the file written into the file store
    """
    try:
        file_path = os.path.join(fileStore.localTempDir, *cmps)
        logger.info("Writing global file {0}".format(file_path))
        file_id = fileStore.writeGlobalFile(file_path)

    except Exception as exc:
        logger.info("Writing global file {0} failed {1}".format(
            file_path, exc))
        file_id = None

    return file_id


def exportFiles(toil, output_directory, job_rv):
    """
    Export files corresponding to the specified id from the file store
    using the specified name

    Parameters
    ----------
    toil : toil.common.Toil
        instance of a Toil context manager
    output_directory : string
        name of output directory
    job_rv : dict
        name, id, and file name of files to export
    """
    for name, spec in job_rv.items():
        if spec['id'] is not None:
            try:
                logger.info("Exporting file {0}".format(spec['name']))
                toil.exportFile(spec['id'], "file://" + os.path.abspath(
                    os.path.join(output_directory, spec['name'])))

            except Exception as exc:
                logger.info("Exporting file {0} failed: {1}".format(
                    spec['name'], exc))


def exportWellAssemblyFiles(toil, assembler, output_directory, well_specs,
                            well_assembly_rvs):
    """
    Export the assembler output files, and the apc output files, if
    apc required by the assembler.

    Parameters
    ----------
    toil : toil.common.Toil
        instance of a Toil context manager
    assembler : str
        name of assembler run, from ['spades', 'shovill', 'novoplasty']
    output_directory : str
        name of directory for export destination
    well_specs : list
        string specifications of the wells
    well_assembly_rvs : list
        dictionaries of assembly job return values
    """
    nW = len(well_specs)
    for iW in range(nW):
        well_output_directory = os.path.join(output_directory, well_specs[iW])
        if not os.path.exists(well_output_directory):
            os.mkdir(well_output_directory)

        # Export all assembler output files from the file store
        exportFiles(toil, well_output_directory,
                    well_assembly_rvs[iW][assembler + "_rv"])

        if assembler in ASSEMBLERS_REQUIRING_APC:
            # Export all apc output files from the file store
            exportFiles(toil, well_output_directory,
                        well_assembly_rvs[iW]['apc_rv'])