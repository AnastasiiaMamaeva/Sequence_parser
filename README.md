# Introduction

When reading scientific articles that analyze protein or DNA sequences, researchers frequently encounter accession numbers instead of the sequences themselves. While this practice is standard, it poses a significant inconvenience for readers aiming to replicate alignments or other in silico experiments described in the articles. Manually retrieving sequences from databases like NCBI can be time-consuming, tedious, and prone to errors.

This project aims to solve this issue by developing a powerful and user-friendly tool that automates sequence retrieval and organization, allowing researchers to focus on their analyses rather than data gathering.

# Usage

The tool is designed to be intuitive and flexible, catering to a range of user needs. Users can input raw text containing accession numbers, even if the text includes extraneous "junk" information such as row numbers, comments, or additional annotations. By specifying a few parameters, such as the database to query or whether to split sequences into multiple files or save them in a single file, the user will acquire:

1. A Clean Table: An organized table summarizing details about the sequences, including sequence type, source organism, and database information.
2. Sequence File(s): Ready-to-use sequence files in standard formats like FASTA.

This tool streamlines sequence retrieval, significantly enhancing the efficiency and accuracy of replicating published experiments. It will be a valuable resource for the bioinformatics and computational biology community, saving time and reducing errors in essential data collection steps.

--
this project is building as a part of the course: https://github.com/szabgab/wis-python-course-2024-11
