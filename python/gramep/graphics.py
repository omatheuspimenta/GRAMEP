"""
graphics
========

This module provides a forked version of the 'snipit' method from \
the '[snipit](https://github.com/aineniamh/snipit)' repository.
The 'snipit' method is used for generating graphical visualizations.

Contents:
    * plot_graphic: Plot graphical visualization of variations and k-mer frequencies.
"""
from collections import defaultdict

from gramep.helpers import get_colours, next_colour
from gramep.messages import Messages
from gramep.utilrs import ref_length
from matplotlib import patches as patches
from matplotlib import pyplot as plt

message = Messages()
"""
Set the Message class for logging.
"""


def plot_graphic(
    variations: list[str],
    reference_path: str,
    freq_kmers: defaultdict[int],
    sequence_name: str,
    save_path: str,
    colour_palette: str = 'classic',
) -> None:
    """
    Plot graphical visualization of variations and k-mer frequencies.

    This function generates a graphical visualization of variations and \
    k-mer frequencies
    based on the provided inputs. The plot is saved to the specified save path.

    Args:
        variations (list): A list of variations to be plotted.
        ref_sequence (str): The reference sequence for comparison.
        freq_kmers (defaultdict[int]): A defaultdict containing k-mer frequencies.
        sequence_name (str): The name of the sequence being analyzed.
        save_path (str): The path to save the generated plot.
        colour_palette (str, optional): The colour palette to use for the plot. \
        Default is 'classic'.

    Returns:
        Message class: A message confirming the plot was saved.
    """
    sorted_variations = []
    y_level = 0
    snp_dict = defaultdict(list)
    ref_vars = {}
    length = ref_length(reference_path)
    colour_dict = get_colours(colour_palette)

    for var in variations:
        pos, variation = var.split(':')
        freq = freq_kmers[var]
        sorted_variations.append((int(pos), variation, freq))
    sorted_variations.sort()

    for variation in sorted_variations:
        y_level += 1
        pos, var, freq = variation
        x_position = pos
        ref = var[0]
        base = var[1]
        ref_vars[x_position] = ref
        snp_dict[x_position].append(
            (
                str('Mutation' + str(y_level) + '\n' + str(freq)),
                ref,
                base,
                y_level,
            )
        )

    spacing = length / (len(snp_dict) + 1)
    y_inc = (spacing * 0.8 * y_level) / length

    # if len(snp_dict) <10:
    #     width = 10
    # else:
    #     width = 0.25* len(snp_dict)

    if y_level < 5:
        height = 5
    else:
        height = y_inc * 3 + 0.5 * y_level + y_inc * 2

    width = height + (height * 0.5)

    # width and height of the figure
    fig, ax = plt.subplots(1, 1, figsize=(width, height), dpi=300)

    y_level = 0
    for snp in snp_dict:
        for sequence in snp_dict[snp]:
            mutation_name, _, _, _ = sequence
            y_level += y_inc

            # either grey or white
            col = next_colour()

            # for each record (sequence) draw a rectangle the length of the whole \
            # genome (either grey or white)
            rect = patches.Rectangle(
                (0, y_level - (0.5 * y_inc)),
                length,
                y_inc,
                alpha=0.3,
                fill=True,
                edgecolor='none',
                facecolor=col,
            )
            ax.add_patch(rect)

            # for each record add the name to the left hand side
            ax.text(
                -50, y_level, mutation_name, size=9, ha='right', va='center'
            )

    position = 0
    for snp in sorted(snp_dict):
        position += spacing

        # write text adjacent to the SNPs shown with the numeric position
        # the text alignment is toggled right/left (top/bottom considering \
        # 90-deg rotation) if the plot is flipped
        ax.text(
            position,
            y_level + (0.55 * y_inc),
            snp,
            size=9,
            ha='center',
            va='bottom',
            rotation=90,
        )

        # snp position labels
        left_of_box = position - (0.4 * spacing)
        right_of_box = position + (0.4 * spacing)

        top_polygon = y_inc * -0.7
        bottom_polygon = y_inc * -1.7

        for sequence in snp_dict[snp]:

            name, ref, var, y_pos = sequence
            bottom_of_box = (y_pos * y_inc) - (0.5 * y_inc)
            # draw box for snp
            # if recombi_out:
            #     rect = patches.Rectangle((left_of_box,bottom_of_box),\
            # spacing*0.8,  y_inc,alpha=0.5, fill=True, \
            # edgecolor='none',facecolor=colour_dict[recombi_out])
            if var in colour_dict:
                rect = patches.Rectangle(
                    (left_of_box, bottom_of_box),
                    spacing * 0.8,
                    y_inc,
                    alpha=0.5,
                    fill=True,
                    edgecolor='none',
                    facecolor=colour_dict[var.upper()],
                )
            else:
                rect = patches.Rectangle(
                    (left_of_box, bottom_of_box),
                    spacing * 0.8,
                    y_inc,
                    alpha=0.5,
                    fill=True,
                    edgecolor='none',
                    facecolor='dimgrey',
                )

            ax.add_patch(rect)

            # sequence variant text
            ax.text(
                position, y_pos * y_inc, var, size=9, ha='center', va='center'
            )

        # reference variant text

        ax.text(position, y_inc * -0.2, ref, size=9, ha='center', va='center')

        # polygon showing mapping from genome to spaced out snps
        x = [snp - 0.5, snp + 0.5, right_of_box, left_of_box, snp - 0.5]

        y = [
            bottom_polygon,
            bottom_polygon,
            top_polygon,
            top_polygon,
            bottom_polygon,
        ]
        coords = list(zip(x, y))

        # draw polygon
        poly = patches.Polygon(
            coords, alpha=0.2, fill=True, edgecolor='none', facecolor='dimgrey'
        )
        ax.add_patch(poly)

        rect = patches.Rectangle(
            (left_of_box, top_polygon),
            spacing * 0.8,
            y_inc,
            alpha=0.1,
            fill=True,
            edgecolor='none',
            facecolor='dimgrey',
        )
        ax.add_patch(rect)

    if len(snp_dict) == 0:
        # snp position labels
        left_of_box = position - (0.4 * position)
        right_of_box = position + (0.4 * position)

        top_polygon = y_inc * -0.7
        bottom_polygon = y_inc * -1.7

    # reference variant rectangle
    rect = patches.Rectangle(
        (0, (top_polygon)),
        length,
        y_inc,
        alpha=0.2,
        fill=True,
        edgecolor='none',
        facecolor='dimgrey',
    )
    ax.add_patch(rect)

    ax.text(-20, y_inc * -0.2, 'Reference', size=9, ha='right', va='center')

    ref_genome_position = y_inc * -2.7

    # reference genome rectangle
    rect = patches.Rectangle(
        (0, ref_genome_position),
        length,
        y_inc,
        alpha=0.2,
        fill=True,
        edgecolor='none',
        facecolor='dimgrey',
    )
    ax.add_patch(rect)

    for var in ref_vars:
        ax.plot(
            [var, var],
            [ref_genome_position, ref_genome_position + (y_inc * 0.98)],
            color='#cbaca4',
        )

    ax.spines['top'].set_visible(False)   ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.yticks([])

    ax.set_xlim(0, length)
    ax.set_ylim(ref_genome_position, y_level + (y_inc * 1.05))

    seq_name = sequence_name.split('/')[-1].split('.')[0]

    ax.tick_params(axis='x', labelsize=8)
    plt.xlabel('Genome position (base)', fontsize=9)
    plt.title(str(str(seq_name) + ' Mutations'), loc='left')
    plt.tight_layout()

    outputDir = save_path + '/' + seq_name + '/results.pdf'
    plt.savefig(outputDir)

    plt.close()

    return message.info_graphic_saved(str(save_path + '/' + seq_name))
