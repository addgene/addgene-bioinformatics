import os


def readGlobalFile(fileStore, file_id, *cmps):
    """
    Read the file corresponding to the specified id from the file
    store into the local temporary directory using the specified name

    Parameters
    ----------
    fileStore : toil.fileStore.FileStore
        instance of a Toil file store
    file_id : str
        id of the file in the file store
    *cmps : list
        path components of the file in the local temporary directory

    Returns
    -------
    file_path : str
        absolute path to the file in the local temporary directory
    """
    file_path = os.path.join(fileStore.localTempDir, *cmps)
    fileStore.readGlobalFile(file_id, userPath=file_path)

    return file_path


def writeGlobalFile(fileStore, *cmps):
    """
    Write the named file from the local temporary directory into the
    file store and return the id

    Parameters
    ----------
    fileStore : toil.fileStore.FileStore
        instance of a Toil file store
    *cmps : list
        path components of the file in the local temporary directory
    Returns
    -------
    file_id : str
        id of the file in the file store
    """
    try:
        file_path = os.path.join(fileStore.localTempDir, *cmps)
        file_id = fileStore.writeGlobalFile(file_path)

    except Exception as exc:
        file_id = None

    return file_id
