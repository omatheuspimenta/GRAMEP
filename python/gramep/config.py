from configparser import ConfigParser

from gramep.classify_utils import classify, predict
from gramep.mutations import get_mutations, get_variants_intersection


def execute_configparser(objective: str, config_file: str):
    """
    Execute a ConfigParser to read configuration settings.

    This function uses the ConfigParser library to read configuration settings from
    the specified config file.

    Args:
        objective (str): A description of the objective of the configuration settings.
        config_file (str): The path to the configuration file to be read.

    """

    config_parser = ConfigParser()
    parameters_file = config_file
    config_parser.read(parameters_file)

    match objective:
        case 'get-mutations':
            reference_path = str(
                config_parser.get('mutations', 'reference_path')
            )
            sequence_path = str(
                config_parser.get('mutations', 'sequence_path')
            )
            annotation_path = str(
                config_parser.get('mutations', 'annotation_path')
            )
            mode = str(config_parser.get('mutations', 'mode'))
            save_path = str(config_parser.get('mutations', 'save_path'))
            word = config_parser.getint('mutations', 'word')
            step = config_parser.getint('mutations', 'step')
            snps_max = config_parser.getint('mutations', 'snps_max')
            dictonary = str(config_parser.get('mutations', 'dictonary'))
            create_report = config_parser.getboolean(
                'mutations', 'create_report'
            )
            chunk_size = config_parser.getint('mutations', 'chunk_size')

            get_mutations(
                reference_path=reference_path,
                sequence_path=sequence_path,
                annotation_path=annotation_path,
                mode=mode,
                save_path=save_path,
                word=word,
                step=step,
                snps_max=snps_max,
                dictonary=dictonary,
                create_report=create_report,
                chunk_size=chunk_size,
            )
        case 'get-intersection':
            save_path = str(config_parser.get('intersection', 'save_path'))
            intersection_seletion = str(
                config_parser.get('intersection', 'intersection_seletion')
            )

            get_variants_intersection(
                save_path=save_path,
                intersection_seletion=intersection_seletion,
            )
        case 'classify':
            word = config_parser.getint('classify', 'word')
            step = config_parser.getint('classify', 'step')
            save_path = str(config_parser.get('classify', 'save_path'))
            dir_path = str(config_parser.get('classify', 'dir_path'))
            dictonary = str(config_parser.get('classify', 'dictonary'))
            get_kmers = config_parser.getboolean('classify', 'get_kmers')
            chunk_size = config_parser.getint('classify', 'chunk_size')

            classify(
                word=word,
                step=step,
                save_path=save_path,
                dir_path=dir_path,
                get_kmers=get_kmers,
                dictonary=dictonary,
                chunk_size=chunk_size,
            )
        case 'predict':
            word = config_parser.getint('predict', 'word')
            step = config_parser.getint('predict', 'step')
            save_path = str(config_parser.get('predict', 'save_path'))
            predict_seq_path = str(
                config_parser.get('predict', 'predict_seq_path')
            )
            dir_path = str(config_parser.get('predict', 'dir_path'))
            dictonary = str(config_parser.get('predict', 'dictonary'))
            load_ranges_path = str(
                config_parser.get('predict', 'load_ranges_path')
            )
            load_model_path = str(
                config_parser.get('predict', 'load_model_path')
            )
            chunk_size = config_parser.getint('predict', 'chunk_size')

            predict(
                word=word,
                step=step,
                save_path=save_path,
                predict_seq_path=predict_seq_path,
                dir_path=dir_path,
                dictonary=dictonary,
                load_ranges_path=load_ranges_path,
                load_model_path=load_model_path,
                chunk_size=chunk_size,
            )
