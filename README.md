# Introduction

When reading scientific articles that analyze protein or DNA sequences, researchers frequently encounter accession numbers instead of the sequences themselves. While this practice is standard, it poses a significant inconvenience for readers aiming to replicate alignments or other in silico experiments described in the articles. Manually retrieving sequences from databases like NCBI can be time-consuming, tedious, and prone to errors.

This project aims to solve this issue by developing a powerful and user-friendly tool that automates sequence retrieval and organization, allowing researchers to focus on their analyses rather than data gathering.

# Usage

The tool is designed to be intuitive and flexible, catering to a range of user needs.

To install the program, download seq_parser.py and make sure you have all requirements in your environment (see the file requirements.txt).

To run the program:

Open command terminal.
Redirect it to the directory where seq_parser.py is placed (cd FolderPath for Bash).
Run the command python seq_parser.py.
The dialog window will open.
Choose:

A file that contains sequences to parse or input it as text in the corresponding window.
Select the NCBI database from which you want to parse your data.
Select whether you want to have output sequences in one merged file or as a folder with multiple files—one file for each sequence.
Select a place you want to store your sequences.
Press "Submit" once you are done.
If you don't input a file/text or output folder, the program will kindly ask you to do so.

You will get a folder named output (or output_1, output_2, etc., if it already exists) filled with files of sequences named with accession numbers or as search_results if the user chose to save as one file.

Fetching takes some time, so the user should wait patiently until the program displays another window saying "All done!" or gives an error in case something went wrong.

You can find an example input file in the project folder.

Have a good time using the app and analyzing sequences!

--
this project is building as a part of the course: https://github.com/szabgab/wis-python-course-2024-11

