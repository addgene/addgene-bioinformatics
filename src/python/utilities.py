import os


def importFile(toil, file_name):
    """
    Import the named file into the file store.

    Parameters
    ----------
    toil : toil.common.Toil
        instance of a Toil context manager
    file_name : str
        name of the file to import

    Returns
    -------
    str
        id of the imported file in the file store
    """
    file_id = toil.importFile("file://" + os.path.abspath(file_name))

    return file_id


def importReadFiles(toil, data_directory, plate_spec, well_specs):
    """
    Import the read one and two files in the data directory for the
    specified plate and well.

    Parameters
    ----------
    toil : toil.common.Toil
        instance of a Toil context manager
    data_directory : string
        name of directory containing plate FASTQ directory
    plate_spec : str
        specification of the plate
    well_spec : list
        specification of the wells
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
                data_directory, "{0}_FASTQ".format(plate_spec),
                "{0}_{1}_R1_001.fastq.gz".format(plate_spec, well_spec))))

        read_two_file_ids.append(importFile(
            toil, os.path.join(
                data_directory, "{0}_FASTQ".format(plate_spec),
                "{0}_{1}_R2_001.fastq.gz".format(plate_spec, well_spec))))

    return read_one_file_ids, read_two_file_ids


def importContigsFile(toil, data_directory):
    """
    Import the contigs FASTA file from the specified data directory.

    Parameters
    ----------
    data_directory : string
        name of directory containing contigs FASTA file
    Returns
    -------
    str
        id of the imported contigs file in the file store
    """
    contigs_file_id = importFile(
        toil, os.path.join(data_directory, "contigs.fasta"))
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
    file_path = os.path.join(fileStore.localTempDir, *cmps)
    fileStore.readGlobalFile(file_id, userPath=file_path)

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
        file_id = fileStore.writeGlobalFile(file_path)

    except Exception as exc:
        file_id = None

    return file_id


def exportFiles(toil, job_rv, output_directory="."):
    """
    Export files corresponding to the specified id from the file store
    using the specified name

    Parameters
    ----------
    toil : toil.common.Toil
        instance of a Toil context manager
    job_rv : dict
        name, id, and file name of files to export
    """
    for name, spec in job_rv.iteritems():
        if spec['id'] is not None:
            toil.exportFile(spec['id'], "file://" + os.path.abspath(
                os.path.join(output_directory, spec['name'])))
