
import pysam
import os
from threading import Thread

class AlignmentFile(pysam.AlignmentFile):
    """
    Extension of pysam AlignmentFile that enables creating AlignmentFile from SAM string.

    Attributes:
        from_text (function): Create AlignmentFile from SAM string.
        close (function): Close pipe
    """
    @classmethod
    def from_text(cls, sam_string: str):
        """
        Create AlignmentFile from SAM string.

        AlignmentFile input data must have file descriptor. io.StringIO
        and Python string do not. Write SAM string to pipe in separate
        thread, then use pipe read fd as input fd to pysam.

        Arguments:
            sam_string (str): String containing SAM data in plaintext. Must include header.
        """
        # Open a pipe with a write and read end
        pipe_fds = os.pipe()
        read_fd, write_fd = pipe_fds
        

        def stream_segs():
            # Writes SAM data into read end of pipe.
            # Runs in separate thread to avoid blocking.
            try:
                with open(write_fd, "w") as pipe:
                    pipe.write(sam_string)
            except OSError:
                pass
            except:
                raise
        try: 
            # Create writer thread
            writer_thread = Thread(target=stream_segs, daemon=True)
            writer_thread.start()

            # AlignmentFile created from pipe 
            file = AlignmentFile(read_fd)

            # Store thread and fds so they can be closed later
            file.writer_thread = writer_thread
            file.pipe_fds = pipe_fds
            return file
        except: 
            cls._close_fds(pipe_fds)
            raise
        
    def close(self):
        """
        Close pipe
        """
        try: self.close_fds(self.pipe_fds)
        except: pass

    def _close_fds(fds: list[int]):
        """
        Try to close each file descriptor (ignores exceptions)
        """
        for fd in fds:
            try: os.close(fd)
            except: pass

    def __del__(self):
        """
        Close pipe.
        """
        self.close()
