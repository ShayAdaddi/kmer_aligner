import os
from pathlib import Path
import logging
import sys


def set_curr_in_dir(base_directory, dir_to_set):
    curr_in_dir = Path(base_directory, dir_to_set)
    validate_path_existence(curr_in_dir)
    return curr_in_dir


def set_curr_out_dir(base_directory, dir_to_set):
    curr_out_dir = Path(base_directory, dir_to_set)
    if not os.path.exists(curr_out_dir):
        curr_out_dir.parent.mkdir(parents=True, exist_ok=True)
        curr_out_dir.mkdir(parents=True, exist_ok=True)
        logging.info('directory: {} created in order to store the iceberg results.'.format(curr_out_dir))
    else:
        logging.warning('directory: {} already exists, if you repeats on analysis which you all ready done the files will be overwritten.'.format(curr_out_dir))
    return curr_out_dir


def validate_path_existence(dir_path):
    try:
        if not os.path.exists(dir_path):
            raise FileNotFoundError(f'{dir_path}')

    except FileNotFoundError as e:
        logger = logging.getLogger()
        logger.error('files directory not exist: {}'.format(dir_path), exc_info=False)
         # sys.exit()
        raise e
