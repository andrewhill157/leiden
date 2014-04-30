import subprocess


def bsub(command, job_name, outfile, queue='hour', memory_usage='4', io_usage='5', wall_time='4:00', project_name='test'):
    """
    Submits a job on an LSF-based distributed computing system. See LSF documentation for maximum allowed values
    for input parameters.

    Args:
        command (str): command to be executed during job execution
        job_name (str): name for the job (can be anything)
        outfile (str): output file to contain outcome of job (errors, statistics, etc.)
        queue (str): name of the queue to submit job to
        memory_usage (str): number of GB memory required
        io_usage (str): number of IO units required
        wall_time (str): max-allowed time for job (4:00 for 4 hours for example)
        project_name (str): name of project job is associated with (can be anything)

    """
    memory_usage = 'rusage[mem=' + str(memory_usage) + ']'
    io_usage = 'rusage[argon_io=' + str(io_usage) + ']'

    subprocess.call(['bsub',
                     '-o', outfile,
                     '-P', project_name,
                     '-q', queue,
                     '-W', wall_time,
                     '-R', memory_usage,
                     '-R', io_usage,
                     '-J', job_name,
                     command])

