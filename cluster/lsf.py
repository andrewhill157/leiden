import subprocess


def bsub(command, job_name, outfile, queue='hour', memory_usage='4', io_usage='5', wall_time='4:00', project_name='test'):
    """
    Submits a job on an LSF-based distributed computing system. See LSF documentation for maximum allowed values
    for input parameters.

    @param command: command to be executed during job execution
    @param job_name: name for the job (can be anything)
    @param outfile: output file to contain outcome of job (errors, statistics, etc.)
    @param queue: name of the queue to submit job to
    @param memory_usage: number of GB memory required
    @param io_usage: number of IO units required
    @param wall_time: max-allowed time for job (4:00 for 4 hours for example)
    @param project_name: name of project job is associated with (can be anything)
    @return:
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

