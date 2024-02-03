import matplotlib.pyplot as plt
import numpy as np

def plot_occurrences():
    """
    Plot occurrences data for major taxonomic groups.

    This function creates a bar plot to visualize the total occurrences and ENA/Mgnify occurrences
    for different major taxonomic groups. The y-axis is set to a logarithmic scale for better visualization.

    Data:
    - taxa: List of major taxonomic groups.
    - total_occurrences: Total occurrences for each taxonomic group.
    - ena_mgnify: Occurrences from ENA/Mgnify for each taxonomic group.
    """
    # Occurrence data
    taxa = ["Animals", "Plants", "Fungi", "Bacteria", "Protists", "Archaea", "Viruses", "Incertae sedis"]
    total_occurrences = [2097448406, 442531533, 38914204, 22722639, 15895213, 442031, 910025, 8014894]
    ena_mgnify = [33065, 376547, 955943, 18355383, 3098014, 335722, 0, 630700]

    # Set up figure and axis
    fig, ax1 = plt.subplots(figsize=(12, 6), dpi=1200)

    # Plot Total Occurrences and set y-axis to logarithmic scale
    color = 'blue'
    ax1.bar(taxa, total_occurrences, color=color, label='Total Occurrences')
    ax1.set_xlabel('Major Taxonomic Group')
    ax1.set_ylabel('Total Occurrences (log scale)', color=color)
    ax1.tick_params('y', colors=color)
    ax1.set_yscale('log')  # Set y-axis to logarithmic scale

    # Create a second y-axis to plot ENA/Mgnify occurrences and set y-axis to logarithmic scale
    ax2 = ax1.twinx()
    color = 'green'
    ax2.bar(taxa, ena_mgnify, color=color, alpha=0.5, label='ENA/Mgnify Occurrences')
    ax2.set_ylabel('ENA/Mgnify Occurrences (log scale)', color=color)
    ax2.tick_params('y', colors=color)
    ax2.set_yscale('log')  # Set y-axis to logarithmic scale

    # Title and legend
    plt.title('Occurrences by Major Taxonomic Group (log scale)')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.85, 0.9))

    plt.show()

if __name__ == "__main__":
    plot_occurrences()
