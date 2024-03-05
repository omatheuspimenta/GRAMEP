"""Message class.

This class contains the messages used in the gramep package.
"""

import logging

from rich.logging import RichHandler

# Define RESULT as a logging level
RESULT = 25
logging.addLevelName(RESULT, 'RESULT')


def result_logging(self, message, *args, **kws):
    """Define logging level."""
    self.log(RESULT, message, *args, **kws)


# Logging configuration
FORMAT = '%(message)s'
# format="[%(levelname)s]: "
logging.basicConfig(
    format=FORMAT,
    level='INFO',
    handlers=[RichHandler(show_time=False, show_path=False, markup=True)],
)

logging.Logger.result = result_logging
log_text = logging.getLogger('rich')
log_text.setLevel(20)


class Messages:
    """Related about messages."""

    def __init__(self):
        """Initialize the class."""
        pass

    def info_classifying(self):
        """Messages.

        Classifying the sequences using Random Forest
        with default parameters
        """
        log_text.info(
            ':evergreen_tree::deciduous_tree: Classifying the sequences using Random Forest with default parameters'
        )

    def info_model_saved(self, string):
        """Messages.

        Model saved in folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Model saved in %s folder', string
        )

    def info_ranges_saved(self, string):
        """Messages.

        Ranges saved in folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Ranges saved in %s folder', string
        )

    def info_kmers_saved(self, string):
        """Messages.

        Exclusive kmers saved in %s folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Exclusive kmers saved in %s folder',
            string,
        )

    def info_newick_saved(self, string):
        """Messages.

        Newick file saved in %s folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Newick file saved in %s folder',
            string,
        )

    def info_phylogenic_tree_saved(self, string):
        """Messages.

        Phylogenic tree saved in %s folder
        """
        log_text.info(
            ':white_check_mark::evergreen_tree: Phylogenic tree saved in %s folder',
            string,
        )

    def info_report_saved(self, string):
        """Messages.

        Report saved in %s folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Report saved in %s folder', string
        )

    def info_freq_saved(self, string):
        """Messages.

        Frequencies saved in %s folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Frequencies saved in %s folder',
            string,
        )

    def info_intersections_saved(self, string):
        """Messages.

        Intersections saved in folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Intersections saved in %s folder',
            string,
        )

    def info_graphic_saved(self, string):
        """Messages.

        Graphic saved in folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Graphic saved in %s folder',
            string,
        )

    def info_dataframe_saved(self, string):
        """Messages.

        Dataframe saved in folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Dataframe saved in %s folder',
            string,
        )

    def info_metrics_saved(self, string):
        """Messages.

        Metrics saved in folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Metrics saved in %s folder',
            string,
        )

    # def info_loading(self, string):
    #     """Messages.

    #     Loading
    #     """
    #     log_text.info(':hourglass: Loading... %s', string)

    def info_loading_kmers(self):
        """Messages.

        Loading
        """
        log_text.info(':hourglass: Loading kmers from sequences')

    def info_start(self):
        """Messages.

        Algorithm start
        """
        log_text.info(':checkered_flag: Starting...')

    def info_start_objetive(self, string):
        """Messages.

        Algorithm start
        """
        log_text.info(':checkered_flag: Starting %s', string)

    def info_kmers_load(self):
        """Messages.

        kmers file loaded
        """
        log_text.info(':white_check_mark::open_file_folder: kmers file loaded')

    # def info_intersections_load(self):
    #     """Messages.

    #     Intersections file loaded
    #     """
    #     log_text.info(
    #         ':white_check_mark::open_file_folder: Intersections file loaded'
    #     )

    def info_get_kmers(self):
        """Messages.

        Getting information about exclusive kmers...
        """
        log_text.info(
            ':mag_right: Getting information about exclusive kmers...'
        )

    def info_done(self):
        """Messages.

        Done
        """
        log_text.info('Done! :boom::tada::sparkles:')

    def info_removed_seq(self, string):
        """Messages.

        Removed sequence
        """
        log_text.warning(':exclamation: Removed sequence %s', string)

    def info_start_prediction(self):
        """Messages.

        Prediction start
        """
        log_text.info(':speech_balloon: Prediction start...')

    def info_prediction(self):
        """Messages.

        Predicting the sequences using Random Forest with loaded model
        """
        log_text.info(
            ':speech_balloon: Predicting the sequences using Random Forest with loaded model'
        )

    def info_predictions_saved(self, string):
        """Messages.

        Predictions saved in folder
        """
        log_text.info(
            ':white_check_mark::file_folder: Predictions saved in %s folder',
            string,
        )

    def info_entropy(self):
        """Messages.

        Entropy analysis
        """
        log_text.info(':robot::bar_chart: Entropy analysis')

    def info_processing_dataframe(self):
        """Messages.

        Processing dataframe
        """
        log_text.info('Processing dataframe')

    def result_mean_kfold(self, string):
        """Messages.

        Mean accuracy - 10 repeated 10 fold Cross Validate
        """
        log_text.result(
            ':dart: Mean accuracy - 10 repeated 10 fold Cross Validate:'
            + string
        )

    def result_confusion_matrix(self, string):
        """Messages.

        Confusion matrix saved in folder
        """
        log_text.info(
            ':red_square::blue_square:Confusion matrix saved in %s folder',
            string,
        )

    def info_creating_heatmap(self):
        """Messages.

        Creating heatmap
        """
        log_text.info('[cyan]Creating heatmap ...')

    def result_heatmap(self, string):
        """Messages.

        Heatmap saved in folder
        """
        log_text.info(
            ':fire::world_map: Heatmap saved in %s folder',
            string,
        )

    def info_creating_distance_linkage(self):
        """Messages.

        Creating distance and linkage
        """
        log_text.info('[cyan] Creating distance matrix and linkage ...')

    def info_selected_parameters(self, word, step):
        """Messages.

        Selected parameters
        """
        log_text.info(
            ':pushpin: Selected parameters: Word: %s and Step: %s', word, step
        )

    def info_grid_running(self, word, step):
        """Messages.

        Running grid search
        """
        log_text.info(
            ':crystal_ball: Running grid search with word: %s and step: %s',
            word,
            step,
        )

    def result_metrics(self, string):
        """Messages.

        Metrics
        """
        log_text.result(':bar_chart: Metrics\n' + string)

    def result_accuracy(self, string):
        """Messages.

        Accuracy
        """
        log_text.result(':dart: Accuracy:' + string)

    def warning_annotation_file(self):
        """Messages.

        Annotation file not found.
        """
        log_text.warning(':x: Annotation file not found.')

    def error_null_values(self):
        """
        Messages.

        Missing values. Please check.
        """
        log_text.error(':x: Missing values. Please check.')

    def error_dict(self):
        """
        Messages.

        Sequence dictonary error. Please check.
        """
        log_text.error(':x: Sequence dictonary error. Please check.')

    def error_no_exclusive_kmers(self):
        """
        Messages.

        No exclusive kmers were identified. Check input parameters.
        """
        log_text.error(
            ':x: No exclusive kmers were identified. Check input parameters.'
        )

    def error_different_kmers_len(self):
        """
        Messages.

        The length of the k-mers must be the same. Check input parameters.
        """
        log_text.error(
            ':x: The length of the k-mers must be the same. Check input parameters.'
        )

    def error_color_pallete(self):
        log_text.error(
            ':x: Please select one of palettes for colorPalette option'
        )

    def warning_removed_sequences(self, string):
        log_text.warning(':warning: Removed %s sequences.', string)
