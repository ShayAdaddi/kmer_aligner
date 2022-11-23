import logging
from pathlib import Path
def set_logging_config(step_name, log_file_path, writing_mode='w'):
    try:
        if log_file_path is None:
            logging.basicConfig(format=f'%(asctime)s - analyzer - %(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S')

        elif step_name == 'analyzer':
            logging.basicConfig(format=f'%(asctime)s - analyzer - %(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S',
                                level=logging.DEBUG, force=True,
                                handlers=[logging.FileHandler(log_file_path, mode=writing_mode), logging.StreamHandler()])
        else:
            logging.basicConfig(format=f'%(asctime)s - analyzer - {step_name} - %(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S',
                                level=logging.DEBUG, force=True,
                                handlers=[logging.FileHandler(log_file_path, mode=writing_mode), logging.StreamHandler()])
    except Exception as e:
        logging.error('couldnt config logging basic configuration.')
        raise e