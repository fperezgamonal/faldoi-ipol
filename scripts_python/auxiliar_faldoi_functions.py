"""
    Call a binary executable
"""
import sys
import os
import shlex
import subprocess
def execute_shell_program(command_line, pr_out):
    """
        Execute a shell program.
    """
    args = shlex.split(command_line)
    #print args
    with open(pr_out, 'w') as ouput_file:
        proc = subprocess.Popen(args, stdout = ouput_file)
        proc.wait()


def cut_matching_list(input_file):
    """
      Recieve a file to reorder the colums of sift matches in the correct order
    """
    dest = input_file.split('.txt')[-2] + '_cut.txt'
    with open(input_file) as input_file, open(dest, 'w+') as dest_w:
        for line in input_file:
            new_l = line.split()
            reorder_line = '%s %s %s %s\n'%(new_l[1],
                                            new_l[0], new_l[5], new_l[4])
            dest_w.write(reorder_line)
    return dest

def cut_deep_list(input_file):
    """
      Recieve a file with matches and rescore to see what is an outlier.
    """
    dest = input_file[:-4] + '_cut.txt'

    with open(input_file) as input_file, open(dest, 'w+') as dest_w:
        for line in input_file:
            new_l = line.split()
            reorder_line = '%s %s %s %s\n'%(new_l[0],
                                            new_l[1], new_l[2], new_l[3])
            dest_w.write(reorder_line)
    return dest


def delete_outliers(input_file, thres):
    """
      Delete from the list the matches that are consider outliers based upon
      a threshold
    """
    th = float(thres)
    
    dest = input_file[:-4] + '_out.txt'

    with open(input_file) as input_file, open(dest, 'w+') as dest_w:
        for line in input_file:
            new_l = line.split()
            val = float(new_l[4])
            if val > th:
                reorder_line = '%s %s %s %s %s\n'%(new_l[0],new_l[1],
                                                new_l[2], new_l[3], new_l[4])
                dest_w.write(reorder_line)
    return dest


def joint_matches(input_file1, input_file2):
    """
      Concatenate two matches files into a single one
    """
    dest = input_file1.split('.')[0] + '_final.txt'
    with open(input_file1) as input_file, open(dest, 'w+') as dest_w:
        for line in input_file:
            new_l = line.split()
            element = '%s %s %s %s\n'%(new_l[0],new_l[1],
                                                new_l[2], new_l[3])
            dest_w.write(element)

    with open(input_file2) as input_file, open(dest, 'a') as dest_w:
        for line in input_file:
            new_l = line.split()
            element = '%s %s %s %s\n'%(new_l[0],new_l[1],
                                                new_l[2], new_l[3])
            dest_w.write(element)            
    return dest

