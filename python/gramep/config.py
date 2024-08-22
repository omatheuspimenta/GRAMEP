from configparser import ConfigParser

from gramep.classify_utils import classify, predict
from gramep.mutations import get_mutations, get_variants_intersection
from gramep.phylogenics import get_phylogenic


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
            save_path = str(config_parser.get('mutations', 'save_path'))
            word = config_parser.getint('mutations', 'word')
            step = config_parser.getint('mutations', 'step')
            snps_max = config_parser.getint('mutations', 'snps_max')
            dictonary = str(config_parser.get('mutations', 'dictonary'))
            create_report = config_parser.getboolean(
                'mutations', 'create_report'
            )
            save_kmers = config_parser.getboolean('mutations', 'save_kmers')
            load_exclusive_kmers = config_parser.getboolean(
                'mutations', 'load_exclusive_kmers'
            )
            path_exclusive_kmers = str(
                config_parser.get('mutations', 'path_exclusive_kmers')
            )
            chunk_size = config_parser.getint('mutations', 'chunk_size')

            get_mutations(
                reference_path=reference_path,
                sequence_path=sequence_path,
                annotation_path=annotation_path,
                save_path=save_path,
                word=word,
                step=step,
                snps_max=snps_max,
                dictonary=dictonary,
                create_report=create_report,
                save_kmers=save_kmers,
                load_exclusive_kmers=load_exclusive_kmers,
                path_exclusive_kmers=path_exclusive_kmers,
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
            should_save_data = config_parser.getboolean(
                'classify', 'should_save_data'
            )
            should_save_model = config_parser.getboolean(
                'classify', 'should_save_model'
            )
            should_save_confusion_matrix = config_parser.getboolean(
                'classify', 'should_save_confusion_matrix'
            )
            chunk_size = config_parser.getint('classify', 'chunk_size')

            classify(
                word=word,
                step=step,
                save_path=save_path,
                dir_path=dir_path,
                dictonary=dictonary,
                should_save_data=should_save_data,
                should_save_model=should_save_model,
                should_save_confusion_matrix=should_save_confusion_matrix,
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
        case 'phylogenic':
            save_path = str(config_parser.get('phylogenic', 'save_path'))
            save_heatmap = config_parser.getboolean(
                'phylogenic', 'save_heatmap'
            )

            get_phylogenic(save_path=save_path, save_heatmap=save_heatmap)
