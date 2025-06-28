import logging
import os
import sys
from rich.logging import RichHandler

def setup_logger(output_dir, analysis_name):
    """
    Sets up a logger with rich formatting for clear console output and
    simultaneously saves a plain text log to a file.
    """
    log_file_path = os.path.join(output_dir, f'{analysis_name}.log')
    
    # Configure logging to write to both a file and the console via RichHandler
    logging.basicConfig(
        level="INFO",
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="[%Y-%m-%d %H:%M:%S]",
        handlers=[
            # This handler prints to the console with rich formatting
            RichHandler(rich_tracebacks=True, show_path=False, markup=True),
            # This handler prints plain text to the log file
            logging.FileHandler(log_file_path)
        ]
    )
    # Get the logger instance
    logger = logging.getLogger("rich")
    return logger
